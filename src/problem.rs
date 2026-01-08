use std::f32::consts::PI;

use crate::diff::n2f_go;
use crate::field::{Ampl, AmplMatrix};
use crate::geom::load_geom;
use crate::multiproblem::{init_result, load_settings_or_default};
use crate::settings::{default_e_perp, default_prop};
use crate::{
    beam::{Beam, BeamPropagation, BeamVariant, DefaultBeamVariant},
    diff::Mapping,
    field::Field,
    geom::{Face, Geom},
    orientation, output,
    result::{GOComponent, Mueller, Results},
    settings::Settings,
    zones::ZoneType,
};

use anyhow::Result;
use log::debug;
use nalgebra::{Complex, Point3};
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rayon::prelude::*;

#[cfg(test)]
mod tests {

    use super::*;
    use nalgebra::Complex;

    #[test]
    fn backscatter_params_computed() {
        // Use default config
        let settings =
            crate::settings::load_default_config().expect("Failed to load default config");
        let geoms = Geom::load("./examples/data/hex.obj").expect("load geom");
        let mut geom = geoms[0].clone();
        init_geom(&settings, &mut geom);

        let mut problem = Problem::new(Some(geom), Some(settings)).unwrap();
        let euler = crate::orientation::Euler::new(30.0, 30.0, 0.0);
        problem.run(Some(&euler)).expect("run");

        let result = &problem.result;

        // Check that backward zone is populated
        let backward_zone = result.zones.backward_zone();
        assert!(
            backward_zone.is_some(),
            "backward zone should exist after solve"
        );

        let bs = backward_zone.unwrap().field_2d.first().unwrap();
        // S11 should be positive for any scattering
        assert!(
            bs.mueller_total[(0, 0)] > 0.0,
            "Total S11 should be positive"
        );

        // Check that backscatter params are computed
        assert!(
            result
                .params
                .backscatter_cross(&GOComponent::Total)
                .is_some(),
            "BackscatterCross should be computed"
        );
        assert!(
            result.params.lidar_ratio(&GOComponent::Total).is_some(),
            "LidarRatio should be computed"
        );
        assert!(
            result
                .params
                .depolarization_ratio(&GOComponent::Total)
                .is_some(),
            "DepolarizationRatio should be computed"
        );

        // Sanity checks on values
        let bs_cross = result
            .params
            .backscatter_cross(&GOComponent::Total)
            .unwrap();
        let lidar = result.params.lidar_ratio(&GOComponent::Total).unwrap();
        let depol = result
            .params
            .depolarization_ratio(&GOComponent::Total)
            .unwrap();

        assert!(bs_cross > 0.0, "BackscatterCross should be positive");
        assert!(lidar > 0.0, "LidarRatio should be positive");
        assert!(
            depol >= 0.0 && depol <= 1.0,
            "DepolarizationRatio should be in [0, 1], got: {}",
            depol
        );
    }

    #[test]
    fn cube_inside_ico() {
        let geoms = Geom::load("./examples/data/cube_inside_ico.obj").unwrap();
        let mut geom = geoms[0].clone();
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
        let mut problem = Problem::new(Some(geom), Some(default_settings)).unwrap();

        problem.propagate_next();
    }
}

/// A solvable physics problem.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
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

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl Problem {
    #[new]
    #[pyo3(signature = (settings = None, geom = None))]
    fn py_new(settings: Option<Settings>, geom: Option<Geom>) -> PyResult<Self> {
        Problem::new(geom, settings).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Failed to create Problem: {}", e))
        })
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

        if let Err(err) = self.run(Some(&euler)) {
            log::error!("Error running problem (will skip this solve): {}", err);
        }

        Ok(())
    }

    pub fn py_print_stats(&self) -> PyResult<()> {
        log::info!("{}", self.result.powers);
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
    pub fn new(geom: Option<Geom>, settings: Option<Settings>) -> Result<Self> {
        let settings = load_settings_or_default(settings);
        let mut geom = match geom {
            Some(g) => g,
            None => load_geom(&settings.geom_name).map_err(|e| {
                anyhow::anyhow!(
                    "Failed to load geometry file '{}': {}",
                    settings.geom_name,
                    e
                )
            })?,
        };
        init_geom(&settings, &mut geom);
        let result = init_result(&settings);

        Ok(Self {
            base_geom: geom.clone(),
            geom,
            beam_queue: vec![],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            settings,
            result,
        })
    }

    /// Resets the problem.
    pub fn reset(&mut self) {
        self.beam_queue.clear();
        self.out_beam_queue.clear();
        self.ext_diff_beam_queue.clear();
        self.result = init_result(&self.settings);
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
    pub fn illuminate(&mut self) -> Result<()> {
        let scaled_wavelength = self.settings.wavelength * self.settings.scale;

        let beam = basic_initial_beam(
            &self.geom,
            scaled_wavelength,
            self.settings.medium_refr_index,
        )?;

        self.beam_queue.push(beam);
        Ok(())
    }

    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new_with_field(geom: Geom, beam: Beam) -> Self {
        let settings = load_settings_or_default(None);
        let result = init_result(&settings);

        Self {
            base_geom: geom.clone(),
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            settings,
            result,
        }
    }

    /// Solve far field for a single zone.
    ///
    /// Determines coherence based on zone type:
    /// - Forward zones: always coherent (optical theorem)
    /// - Other zones: respect global coherence setting
    fn solve_far_zone(&mut self, component: GOComponent, zone_idx: usize) {
        let (queue, mapping, fov_factor) = match component {
            GOComponent::Beam => (
                self.out_beam_queue.clone(),
                self.settings.mapping,
                self.settings.fov_factor,
            ),
            GOComponent::ExtDiff => (
                self.ext_diff_beam_queue.clone(),
                Mapping::ApertureDiffraction,
                None,
            ),
            GOComponent::Total => {
                panic!("No such beam queue exists for GOComponent: {:?}", component)
            }
        };

        let zone = &self.result.zones.all()[zone_idx];
        let zone_type = zone.zone_type;
        let zone_scheme = zone.scheme.clone();
        let bins = zone.bins.clone();

        // Forward zones are always coherent (optical theorem)
        let use_coherence = match zone_type {
            ZoneType::Forward => true,
            _ => self.settings.coherence,
        };

        // Map beams to this zone's bins
        let map_beam_to_zone = |beam: &Beam| -> Vec<(usize, Ampl)> {
            match mapping {
                Mapping::GeometricOptics => n2f_go(&zone_scheme, &bins, beam),
                Mapping::ApertureDiffraction => beam.diffract(&bins, fov_factor),
            }
        };

        if use_coherence {
            // Coherent: accumulate amplitudes, convert to Mueller at end
            let zero_ampls: Vec<(usize, Ampl)> =
                bins.iter().map(|_| Ampl::zeros()).enumerate().collect();

            let ampls: Vec<Ampl> = queue
                .par_iter()
                .map(|beam| map_beam_to_zone(beam))
                .reduce(
                    || zero_ampls.clone(),
                    |mut acc, val| {
                        for (i, ampl) in val.into_iter() {
                            acc[i].1 += ampl;
                        }
                        acc
                    },
                )
                .into_iter()
                .map(|x| x.1)
                .collect();

            // Assign to zone's field_2d
            let zone = &mut self.result.zones.all_mut()[zone_idx];
            for (field, ampl) in zone.field_2d.iter_mut().zip(ampls) {
                match component {
                    GOComponent::Total => {
                        field.ampl_total += ampl;
                        field.mueller_total = field.ampl_total.to_mueller();
                    }
                    GOComponent::Beam => {
                        field.ampl_beam += ampl;
                        field.mueller_beam = field.ampl_beam.to_mueller();
                    }
                    GOComponent::ExtDiff => {
                        field.ampl_ext += ampl;
                        field.mueller_ext = field.ampl_ext.to_mueller();
                    }
                }
            }
        } else {
            // Incoherent: convert each beam to Mueller, then sum
            let zero_muellers: Vec<(usize, Mueller)> =
                bins.iter().map(|_| Mueller::zeros()).enumerate().collect();

            let muellers: Vec<Mueller> = queue
                .par_iter()
                .map(|beam| {
                    let ampls = map_beam_to_zone(beam);
                    ampls
                        .into_iter()
                        .map(|(i, a)| (i, a.to_mueller()))
                        .collect()
                })
                .reduce(
                    || zero_muellers.clone(),
                    |mut acc, val: Vec<(usize, Mueller)>| {
                        for (i, mueller) in val.into_iter() {
                            acc[i].1 += mueller;
                        }
                        acc
                    },
                )
                .into_iter()
                .map(|x| x.1)
                .collect();

            // Assign to zone's field_2d
            let zone = &mut self.result.zones.all_mut()[zone_idx];
            for (field, mueller) in zone.field_2d.iter_mut().zip(muellers) {
                match component {
                    GOComponent::Total => field.mueller_total += mueller,
                    GOComponent::Beam => field.mueller_beam += mueller,
                    GOComponent::ExtDiff => field.mueller_ext += mueller,
                }
            }
        }
    }

    /// Combine beam and ext_diff components for a zone.
    fn combine_far_zone(&mut self, zone_idx: usize) {
        let zone = &mut self.result.zones.all_mut()[zone_idx];
        let use_coherence = match zone.zone_type {
            ZoneType::Forward => true,
            _ => self.settings.coherence,
        };

        for field in zone.field_2d.iter_mut() {
            if use_coherence {
                field.ampl_total = field.ampl_beam + field.ampl_ext;
                field.mueller_total = field.ampl_total.to_mueller();
            } else {
                field.mueller_total = field.mueller_beam + field.mueller_ext;
            }
        }
    }

    /// Solves the far field problem by mapping the near field either by geometric optics or aperture diffraction. Optionally, choose to consider coherence between beams.
    pub fn solve_far(&mut self) {
        // Process each zone
        let num_zones = self.result.zones.len();
        for zone_idx in 0..num_zones {
            self.solve_far_zone(GOComponent::ExtDiff, zone_idx);
            self.solve_far_zone(GOComponent::Beam, zone_idx);
            self.combine_far_zone(zone_idx);
        }
    }
    /// Solve an entire problem by tracing beams in the near field, then mapping to the far field, and finally converting to 1D mueller matrices
    pub fn solve(&mut self) {
        debug!("solving near-field problem");
        self.solve_near();
        debug!("solving far-field problem");
        self.solve_far();
        debug!("computing 1d-mueller matrices");
        self.mueller_to_1d();
        debug!("computing parameters");
        self.compute_params();
    }

    pub fn compute_params(&mut self) {
        let _ = self.result.compute_params(self.settings.wavelength);
    }

    pub fn mueller_to_1d(&mut self) {
        self.result.mueller_to_1d();
    }

    pub fn run(&mut self, euler: Option<&orientation::Euler>) -> Result<()> {
        self.init();
        match euler {
            Some(euler) => {
                self.orient(euler)?;
            }
            None => {
                // No rotation
            }
        }
        self.illuminate()?;
        self.solve();
        Ok(())
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
        let output_manager = output::OutputManager::new(&self.settings, &self.result);
        let _ = output_manager.write_all();
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
                log::warn!("Unknown beam type, returning empty outputs.");
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

    pub fn orient(&mut self, euler: &orientation::Euler) -> Result<()> {
        self.geom
            .euler_rotate(euler, self.settings.orientation.euler_convention)
    }
}

// /// Collects a 2d array as a list of lists.
// /// There is probably already a function for this in ndarray.
// pub fn collect_mueller(muellers: &[Mueller]) -> Vec<Vec<f32>> {
//     let mut mueller_list = Vec::new();
//     for mueller in muellers.iter() {
//         mueller_list.push(mueller.to_vec());
//     }
//     mueller_list
// }

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
fn basic_initial_beam(
    geom: &Geom,
    wavelength: f32,
    medium_refractive_index: Complex<f32>,
) -> Result<Beam> {
    const FAC: f32 = 1.1; // scale factor to stretch beam to cover geometry
    let bounds = geom.bounds();
    let (min, max) = (bounds.0.map(|v| v * FAC), bounds.1.map(|v| v * FAC));

    let clip_vertices = vec![
        Point3::new(max[0], max[1], max[2]),
        Point3::new(max[0], min[1], max[2]),
        Point3::new(min[0], min[1], max[2]),
        Point3::new(min[0], max[1], max[2]),
    ];

    let mut clip = Face::new_simple(clip_vertices, None, None)?;
    clip.data_mut().area = Some((max[0] - min[0]) * (max[1] - min[1]));
    let mut field = Field::new_identity(default_e_perp(), default_prop())?;

    // propagate field backwards so its as if the beam comes from z=0
    let dist = bounds.1[2] * FAC;
    let wavenumber = 2.0 * PI / wavelength;
    let arg = -dist * wavenumber * medium_refractive_index.re;
    field.wind(arg);

    let beam = Beam::new_from_field(clip, medium_refractive_index, field, wavelength);
    Ok(beam)
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
