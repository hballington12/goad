use crate::diff::n2f_mapping_go;
use crate::result::MuellerMatrix;
use crate::{
    beam::{Beam, BeamPropagation, BeamVariant, DefaultBeamVariant},
    bins::{generate_bins, SolidAngleBin},
    diff::Mapping,
    field::Field,
    geom::{Face, Geom},
    orientation, output,
    result::{Ampl, AmplMatrix, GOComponent, Mueller, Results, ScattResult, ScatteringBin},
    settings::{load_config, Settings},
};

use nalgebra::{Complex, Matrix2, Point3, Vector3};
use pyo3::prelude::*;
use rayon::prelude::*;

#[cfg(test)]
mod tests {

    use super::*;
    use nalgebra::Complex;

    #[test]
    fn cube_inside_ico() {
        let mut geom = Geom::from_file("./examples/data/cube_inside_ico.obj").unwrap();
        geom.shapes[0].refr_index = Complex {
            // modify the refractive index of the outer shape
            re: 2.0,
            im: 0.1,
        };
        geom.shapes[1].parent_id = Some(0); // set the parent of the second shape to be the first
        assert_ne!(geom.shapes[0].refr_index, geom.shapes[1].refr_index);
        assert_eq!(
            geom.shapes[0],
            geom.shapes[geom.shapes[1].parent_id.unwrap()]
        );

        // Use default config to avoid loading local.toml which may have test-breaking settings
        let default_settings =
            crate::settings::load_default_config().expect("Failed to load default config");
        let mut problem = Problem::new(Some(geom), Some(default_settings));

        problem.propagate_next();
    }
}

/// A solvable physics problem.
#[pyclass]
#[derive(Debug, Clone)] // Added Default derive
pub struct Problem {
    pub base_geom: Geom,                // original geometry
    pub geom: Geom,                     // geometry to trace beams in
    pub beam_queue: Vec<Beam>,          // beams awaiting near-field propagation
    pub out_beam_queue: Vec<Beam>,      // beams awaiting diffraction
    pub ext_diff_beam_queue: Vec<Beam>, // beams awaiting external diffraction
    pub settings: Settings,             // runtime settings
    pub result: Results,                // results of the problem
}

#[pymethods]
impl Problem {
    #[new]
    #[pyo3(signature = (settings = None, geom = None))]
    fn py_new(settings: Option<Settings>, geom: Option<Geom>) -> Self {
        Problem::new(geom, settings)
    }

    /// Setter function for the problem settings
    #[setter]
    pub fn set_settings(&mut self, settings: Settings) {
        self.settings = settings;
    }

    /// Getter function for the problem settings
    #[getter]
    pub fn get_settings(&self) -> Settings {
        self.settings.clone()
    }

    /// Getter function for the geometry
    #[getter]
    pub fn get_geom(&self) -> Geom {
        self.geom.clone()
    }

    pub fn py_solve(&mut self) -> PyResult<()> {
        self.reset();

        // println!("{:#?}", self.settings);

        let euler = match self.settings.orientation.scheme {
            orientation::Scheme::Discrete { ref eulers } => &eulers[0].clone(),
            _ => {
                panic!("Python solve is only supperted for discrete orientation scheme")
            }
        };

        self.run(Some(euler));
        Ok(())
    }

    pub fn py_print_stats(&self) -> PyResult<()> {
        println!("{}", self.result.powers);
        Ok(())
    }

    /// Get the results object
    #[getter]
    pub fn get_results(&self) -> Results {
        self.result.clone()
    }
}

impl Problem {
    /// Creates a new `Problem` from optional `Geom` and `Settings`.
    /// If settings not provided, loads from config file.
    /// If geom not provided, loads from file using settings.geom_name.
    pub fn new(geom: Option<Geom>, settings: Option<Settings>) -> Self {
        let settings = settings.unwrap_or_else(|| load_config().expect("Failed to load config"));
        let mut geom = geom.unwrap_or_else(|| {
            Geom::from_file(&settings.geom_name).expect("Failed to load geometry")
        });
        init_geom(&settings, &mut geom);

        let bins = generate_bins(&settings.binning.scheme);
        let solution = Results::new_empty(&bins);

        let problem = Self {
            base_geom: geom.clone(),
            geom,
            beam_queue: vec![],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            settings,
            result: solution,
        };

        problem
    }

    /// Resets the problem.
    pub fn reset(&mut self) {
        self.beam_queue.clear();
        self.out_beam_queue.clear();
        self.ext_diff_beam_queue.clear();
        self.result = Results::new_empty(&self.result.bins());
        self.geom.clone_from(&self.base_geom);
    }

    /// Initialises the geometry and scales it.
    pub fn init(&mut self) {
        // Apply geometry scaling if set
        if let Some(scale) = &self.settings.geom_scale {
            self.geom.vector_scale(scale);
        }
        // Apply distortion if set
        if let Some(distortion) = self.settings.distortion {
            self.geom.distort(distortion, self.settings.seed);
        }
        self.geom.recentre();
        self.settings.scale = self.geom.rescale();
    }

    /// Illuminates the problem with a basic initial beam.
    pub fn illuminate(&mut self) {
        let scaled_wavelength = self.settings.wavelength * self.settings.scale;

        let beam = basic_initial_beam(
            &self.geom,
            scaled_wavelength,
            self.settings.medium_refr_index,
        );

        self.beam_queue.push(beam);
    }

    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new_with_field(geom: Geom, beam: Beam) -> Self {
        let settings = load_config().expect("Failed to load config");

        let bins = generate_bins(&settings.binning.scheme);
        let solution = Results::new_empty(&bins);

        Self {
            base_geom: geom.clone(),
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            settings,
            result: solution,
        }
    }

    fn diffract_outbeams(
        queue: &[Beam],
        bins: &[SolidAngleBin],
        fov_factor: Option<f32>,
    ) -> Vec<Ampl> {
        // Calculate far-field amplitudes by diffracting all outbeams in parallel
        queue
            .par_iter()
            .map(|outbeam| outbeam.diffract(&bins, fov_factor))
            .reduce(
                || vec![Matrix2::<Complex<f32>>::zeros(); bins.len()],
                |mut acc, local| {
                    acc.iter_mut().zip(local).for_each(|(a, l)| *a += l);
                    acc
                },
            )
    }

    /// Combines the external diffraction and outbeams to get the far-field solution.
    fn combine_far(&mut self) {
        for result in self.result.field_2d.iter_mut() {
            // Combine beam and external diffraction amplitudes
            let ampl_total = match (result.ampl_beam, result.ampl_ext) {
                (Some(beam), Some(ext)) => beam + ext,
                (Some(beam), None) => beam,
                (None, Some(ext)) => ext,
                (None, None) => continue, // Skip if both are None
            };

            // Store the combined amplitude
            result.ampl_total = Some(ampl_total);
        }
    }

    /// Helper to add amplitudes to results, handling None cases
    fn add_amplitudes_to_results<B: ScatteringBin>(
        results: &mut [ScattResult<B>],
        amplitudes: Vec<Ampl>,
        component: GOComponent,
    ) {
        for (result, ampl) in results.iter_mut().zip(amplitudes) {
            // Skip invalid amplitudes
            if !ampl.iter().all(|c| c.re.is_finite() && c.im.is_finite()) {
                continue;
            }

            let target_ampl = match component {
                GOComponent::Total => &mut result.ampl_total,
                GOComponent::Beam => &mut result.ampl_beam,
                GOComponent::ExtDiff => &mut result.ampl_ext,
            };

            match target_ampl.as_mut() {
                Some(existing) => *existing += ampl,
                None => *target_ampl = Some(ampl),
            }
        }
    }

    pub fn solve_far_ext_diff(&mut self) {
        let fov_factor = None; // don't truncate by field of view for external diffraction
        let bins = self.result.bins();
        let ampls = Self::diffract_outbeams(&self.ext_diff_beam_queue, &bins, fov_factor);

        Self::add_amplitudes_to_results(&mut self.result.field_2d, ampls, GOComponent::ExtDiff);
    }

    pub fn solve_far_outbeams(&mut self) {
        let ampls = match self.settings.mapping {
            Mapping::GeometricOptics => n2f_mapping_go(
                &mut self.out_beam_queue,
                &self.settings.binning,
                &self.result.bins(),
            ),
            Mapping::ApertureDiffraction => {
                let fov_factor = self.settings.fov_factor; // truncate by field of view for outbeams
                let bins = self.result.bins();
                Self::diffract_outbeams(&self.out_beam_queue, &bins, fov_factor)
            }
        };

        Self::add_amplitudes_to_results(&mut self.result.field_2d, ampls, GOComponent::Beam);
    }
    pub fn solve_far(&mut self) {
        self.solve_far_ext_diff();
        self.solve_far_outbeams();
        self.combine_far();
    }

    fn compute_mueller(&mut self) {
        for result in self.result.field_2d.iter_mut() {
            // Convert amplitude matrices to Mueller matrices
            result.mueller_total = result.ampl_total.map(|a| a.to_mueller());
            result.mueller_beam = result.ampl_beam.map(|a| a.to_mueller());
            result.mueller_ext = result.ampl_ext.map(|a| a.to_mueller());
        }
    }

    pub fn solve(&mut self) {
        self.solve_near();
        self.solve_far();
        self.compute_mueller();
        self.try_mueller_to_1d();
    }

    pub fn try_mueller_to_1d(&mut self) {
        self.result.try_mueller_to_1d(&self.settings.binning.scheme);
    }

    pub fn try_params(&mut self) {
        match self.result.compute_params(self.settings.wavelength) {
            Ok(()) => {
                // println!("Params computed successfully");
            }
            Err(e) => {
                println!("Failed to compute params: {}", e);
            }
        }
    }

    pub fn run(&mut self, euler: Option<&orientation::Euler>) {
        self.init();
        match euler {
            Some(euler) => {
                self.orient(euler);
            }
            None => {
                // No rotation
            }
        }
        self.illuminate();
        self.solve();
        self.try_mueller_to_1d();
        self.try_params();
    }

    /// Trace beams to solve the near-field problem.
    pub fn solve_near(&mut self) {
        loop {
            if self.beam_queue.len() == 0 {
                break;
            }

            let input_power = self.result.powers.input;
            let output_power = self.result.powers.output;

            if output_power / input_power > self.settings.cutoff {
                // add remaining power in beam queue to missing power due to cutoff
                self.result.powers.trnc_cop += self
                    .beam_queue
                    .iter()
                    .map(|beam| beam.power() / self.settings.scale.powi(2))
                    .sum::<f32>();
                break;
            }

            self.propagate_next();
        }
    }

    pub fn writeup(&self) {
        // Collect Mueller matrices by component type
        let mueller_total: Vec<Mueller> = self
            .result
            .field_2d
            .iter()
            .map(|field| field.mueller_total.unwrap_or_else(Mueller::zeros))
            .collect();

        let mueller_beam: Vec<Mueller> = self
            .result
            .field_2d
            .iter()
            .map(|field| field.mueller_beam.unwrap_or_else(Mueller::zeros))
            .collect();

        let mueller_ext: Vec<Mueller> = self
            .result
            .field_2d
            .iter()
            .map(|field| field.mueller_ext.unwrap_or_else(Mueller::zeros))
            .collect();

        let _ = output::write_mueller(
            &self.result.bins(),
            &mueller_total,
            "",
            &self.settings.directory,
        );
        let _ = output::write_mueller(
            &self.result.bins(),
            &mueller_beam,
            "_beam",
            &self.settings.directory,
        );
        let _ = output::write_mueller(
            &self.result.bins(),
            &mueller_ext,
            "_ext",
            &self.settings.directory,
        );
        let _ = output::write_result(&self.result, &self.settings.directory);
    }

    /// Propagates the next beam in the queue.
    pub fn propagate_next(&mut self) -> Option<BeamPropagation> {
        // Try to pop the next beam from the queue
        let Some(mut beam) = self.beam_queue.pop() else {
            return None;
        };

        // Compute the outputs by propagating the beam
        let outputs = match &mut beam.variant {
            BeamVariant::Default(..) => self.propagate_default(&mut beam),
            BeamVariant::Initial => self.propagate_initial(&mut beam),
            _ => {
                println!("Unknown beam type, returning empty outputs.");
                Vec::new()
            }
        };

        self.result.powers.absorbed += beam.absorbed_power / self.settings.scale.powi(2);
        self.result.powers.trnc_clip +=
            (beam.clipping_area - beam.csa()) * beam.power() / self.settings.scale.powi(2);

        // Process each output beam
        for output in outputs.iter() {
            let output_power = output.power() / self.settings.scale.powi(2);
            match (&beam.variant, &output.variant) {
                (BeamVariant::Default(..), BeamVariant::Default(..)) => {
                    self.insert_beam(output.clone())
                }
                (BeamVariant::Default(..), BeamVariant::OutGoing) => {
                    self.result.powers.output += output_power;
                    self.insert_outbeam(output.clone());
                }
                (BeamVariant::Initial, BeamVariant::Default(..)) => {
                    self.result.powers.input += output_power;
                    self.insert_beam(output.clone());
                }
                (BeamVariant::Initial, BeamVariant::ExternalDiff) => {
                    self.result.powers.ext_diff += output_power;
                    self.ext_diff_beam_queue.push(output.clone());
                }
                _ => {}
            }
        }
        Some(BeamPropagation::new(beam, outputs))
    }

    fn propagate_initial(&mut self, beam: &mut Beam) -> Vec<Beam> {
        match beam.propagate(
            &mut self.geom,
            self.settings.medium_refr_index,
            self.settings.beam_area_threshold(),
        ) {
            Ok((outputs, ..)) => outputs,

            Err(_) => Vec::new(),
        }
    }

    /// Propagates a beam with the default settings.
    /// Cycles through checks to decide whether to propagate the beam or not.
    fn propagate_default(&mut self, beam: &mut Beam) -> Vec<Beam> {
        // beam power is below threshold
        if beam.power() < self.settings.beam_power_threshold * self.settings.scale.powi(2) {
            self.result.powers.trnc_energy += beam.power() / self.settings.scale.powi(2);
            return Vec::new();
        }

        // beam area is below threshold
        if beam.face.data().area.unwrap() < self.settings.beam_area_threshold() {
            self.result.powers.trnc_area += beam.power() / self.settings.scale.powi(2);
            return Vec::new();
        }

        // total internal reflection considerations
        if let BeamVariant::Default(DefaultBeamVariant::Tir) = beam.variant {
            if beam.tir_count >= self.settings.max_tir {
                self.result.powers.trnc_ref += beam.power() / self.settings.scale.powi(2);
                return Vec::new();
            } else {
                return self.propagate(beam);
            }
        }

        // beam recursion over the maximum
        if beam.rec_count > self.settings.max_rec {
            self.result.powers.trnc_rec += beam.power() / self.settings.scale.powi(2);
            return Vec::new();
        }

        // else, propagate the beam
        self.propagate(beam)
    }

    fn propagate(&mut self, beam: &mut Beam) -> Vec<Beam> {
        match beam.propagate(
            &mut self.geom,
            self.settings.medium_refr_index,
            self.settings.beam_area_threshold(),
        ) {
            Ok((outputs, area_power_loss)) => {
                self.result.powers.trnc_area += area_power_loss / self.settings.scale.powi(2);
                outputs
            }
            Err(_) => {
                self.result.powers.clip_err += beam.power() / self.settings.scale.powi(2);
                Vec::new()
            }
        }
    }

    /// Inserts a beam into the beam queue such that beams with greatest power
    /// are prioritised for dequeueing. Order is ascending because beams are
    /// process by popping.
    pub fn insert_beam(&mut self, beam: Beam) {
        let pos = get_position_by_power(beam.power(), &self.beam_queue, true);
        self.beam_queue.insert(pos, beam);
    }

    /// Inserts a beam into the outbeam queue such that beams with greatest power
    /// are prioritised for dequeueing. Order is descending because outbeams are
    /// process sequentially.
    pub fn insert_outbeam(&mut self, beam: Beam) {
        let pos = get_position_by_power(beam.power(), &self.out_beam_queue, false);
        self.out_beam_queue.insert(pos, beam);
    }

    pub fn orient(&mut self, euler: &orientation::Euler) {
        if let Err(error) = self
            .geom
            .euler_rotate(euler, self.settings.orientation.euler_convention)
        {
            panic!("Error rotating geometry: {}", error);
        }
    }
}

/// Collects a 2d array as a list of lists.
/// There is probably already a function for this in ndarray.
pub fn collect_mueller(muellers: &[Mueller]) -> Vec<Vec<f32>> {
    let mut mueller_list = Vec::new();
    for mueller in muellers.iter() {
        mueller_list.push(mueller.to_vec());
    }
    mueller_list
}

/// Find the position to insert the beam using binary search.
fn get_position_by_power(value: f32, queue: &Vec<Beam>, ascending: bool) -> usize {
    queue
        .binary_search_by(|x| {
            let cmp = x
                .power()
                .partial_cmp(&value)
                .unwrap_or(std::cmp::Ordering::Equal);

            if ascending {
                cmp
            } else {
                cmp.reverse()
            }
        })
        .unwrap_or_else(|e| e)
}

/// Creates a basic initial beam for full illumination of the geometry along the z-axis.
fn basic_initial_beam(geom: &Geom, wavelength: f32, medium_refractive_index: Complex<f32>) -> Beam {
    const FAC: f32 = 1.1; // scale factor to stretch beam to cover geometry
    let bounds = geom.bounds();
    let (min, max) = (bounds.0.map(|v| v * FAC), bounds.1.map(|v| v * FAC));

    let clip_vertices = vec![
        Point3::new(max[0], max[1], max[2]),
        Point3::new(max[0], min[1], max[2]),
        Point3::new(min[0], min[1], max[2]),
        Point3::new(min[0], max[1], max[2]),
    ];

    let mut clip = Face::new_simple(clip_vertices, None, None).unwrap();
    clip.data_mut().area = Some((max[0] - min[0]) * (max[1] - min[1]));
    let mut field = Field::new_identity(Vector3::x(), -Vector3::z()).unwrap();

    // propagate field backwards so its as if the beam comes from z=0
    let dist = bounds.1[2];
    let wavenumber = 2.0 * std::f32::consts::PI / wavelength;
    let arg = -dist * wavenumber * medium_refractive_index.re;
    field.ampl *= Complex::new(arg.cos(), arg.sin());

    let beam = Beam::new_from_field(
        clip,
        -Vector3::z(),
        medium_refractive_index,
        field,
        wavelength,
    );
    beam
}

/// Initialises the geometry with the refractive indices from the settings.
/// In the future, this function will be extended to provide additional checks
/// to ensure the geometry is well-defined.
pub fn init_geom(settings: &Settings, geom: &mut Geom) {
    for shape in geom.shapes.iter_mut() {
        shape.refr_index = settings.particle_refr_index[0]; // default refr index is first value
    }
    for (i, refr_index) in settings.particle_refr_index.iter().enumerate() {
        if i >= geom.shapes.len() {
            break;
        }
        geom.shapes[i].refr_index = *refr_index;
    }
    geom.recentre();
}
