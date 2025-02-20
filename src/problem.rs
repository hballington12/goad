use crate::{
    beam::{Beam, BeamPropagation, BeamType, BeamVariant},
    bins::generate_bins,
    field::Field,
    geom::{Face, Geom},
    helpers::draw_face,
    orientation::Orientations,
    output,
    settings::{self, Settings},
};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};
use ndarray::Array2;
use rayon::prelude::*;
use std::{fmt, ops::DivAssign, ops::Add};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use pyo3::prelude::*;

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

        let mut problem = Problem::new(geom, None);

        problem.propagate_next();
    }
}

#[derive(Debug, Copy,  Clone, PartialEq)]
pub struct Powers {
    pub input: f32,       // near-field input power
    pub output: f32,      // near-field output power
    pub absorbed: f32,    // near-field absorbed power
    pub trnc_ref: f32,    // truncated power due to max reflections
    pub trnc_rec: f32,    // truncated power due to max recursions
    pub trnc_clip: f32,   // truncated power due to clipping
    pub trnc_energy: f32, // truncated power due to threshold beam power
    pub clip_err: f32, // truncated power due to clipping error
    pub trnc_area: f32,   // truncated power due to area threshold
    pub trnc_cop: f32,    // truncated power due to cutoff power
    pub ext_diff: f32,    // external diffraction power
}

impl DivAssign<f32> for Powers {
    fn div_assign(&mut self, rhs: f32) {
        self.input /= rhs;
        self.output /= rhs;
        self.absorbed /= rhs;
        self.trnc_ref /= rhs;
        self.trnc_rec /= rhs;
        self.trnc_clip /= rhs;
        self.trnc_energy /= rhs;
        self.clip_err /= rhs;
        self.trnc_area /= rhs;
        self.trnc_cop /= rhs;
        self.ext_diff /= rhs;
    }
}

impl Add for Powers {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            input: self.input + other.input,
            output: self.output + other.output,
            absorbed: self.absorbed + other.absorbed,
            trnc_ref: self.trnc_ref + other.trnc_ref,
            trnc_rec: self.trnc_rec + other.trnc_rec,
            trnc_clip: self.trnc_clip + other.trnc_clip,
            trnc_energy: self.trnc_energy + other.trnc_energy,
            clip_err: self.clip_err + other.clip_err,
            trnc_area: self.trnc_area + other.trnc_area,
            trnc_cop: self.trnc_cop + other.trnc_cop,
            ext_diff: self.ext_diff + other.ext_diff,
        }
    }
}

use std::ops::AddAssign;
use std::time::Instant;

impl AddAssign for Powers {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            input: self.input + other.input,
            output: self.output + other.output,
            absorbed: self.absorbed + other.absorbed,
            trnc_ref: self.trnc_ref + other.trnc_ref,
            trnc_rec: self.trnc_rec + other.trnc_rec,
            trnc_clip: self.trnc_clip + other.trnc_clip,
            trnc_energy: self.trnc_energy + other.trnc_energy,
            clip_err: self.clip_err + other.clip_err,
            trnc_area: self.trnc_area + other.trnc_area,
            trnc_cop: self.trnc_cop + other.trnc_cop,
            ext_diff: self.ext_diff + other.ext_diff,
        };
    }
}

impl Powers {
    pub fn new() -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            absorbed: 0.0,
            trnc_ref: 0.0,
            trnc_rec: 0.0,
            trnc_clip: 0.0,
            trnc_energy: 0.0,
            clip_err: 0.0,
            trnc_area: 0.0,
            trnc_cop: 0.0,
            ext_diff: 0.0,
        }
    }

    /// Returns the power unaccounted for.
    pub fn missing(&self) -> f32 {
        self.input
            - (self.output
                + self.absorbed
                + self.trnc_ref
                + self.trnc_rec
                // + self.trnc_clip
                + self.trnc_area    
                + self.clip_err
                + self.trnc_cop
                + self.trnc_energy)
    }
}

impl fmt::Display for Powers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Powers:")?;
        writeln!(f, "  Input:            {:.6}", self.input)?;
        writeln!(f, "  Output:           {:.6}", self.output)?;
        writeln!(f, "  Absorbed:         {:.6}", self.absorbed)?;
        writeln!(f, "  Trunc. Refl:      {:.6}", self.trnc_ref)?;
        writeln!(f, "  Trunc. Rec:       {:.6}", self.trnc_rec)?;
        // writeln!(f, "  Trunc. Clip:   {:.6}", self.trnc_clip)?;
        writeln!(f, "  Clip Err:         {:.6}", self.clip_err)?;
        writeln!(f, "  Trunc. Energy:    {:.6}", self.trnc_energy)?;
        writeln!(f, "  Trunc. Area:      {:.6}", self.trnc_area)?;
        writeln!(f, "  Trunc. Cop:       {:.6}", self.trnc_cop)?;
        writeln!(f, "  Other:            {:.6}", self.missing())?;
        writeln!(f, "  External Diff:    {:.6}", self.ext_diff)
    }
}

/// A solvable physics problem.
#[pyclass]
#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct Problem {
    pub geom: Geom,                       // geometry to trace beams in
    pub beam_queue: Vec<Beam>,            // beams awaiting near-field propagation
    pub out_beam_queue: Vec<Beam>,        // beams awaiting diffraction
    pub ext_diff_beam_queue: Vec<Beam>,   // beams awaiting external diffraction
    pub powers: Powers,                   // different power contributions
    pub bins: Vec<(f32, f32)>,            // bins for far-field diffraction
    pub ampl: Vec<Matrix2<Complex<f32>>>, // total amplitude in far-field
    pub settings: Settings,               // runtime settings
    scale_factor: f32, // scaling factor for geometry
}

#[pymethods]
impl Problem {
    #[new]
    fn py_new(geom: Geom, settings: Settings) -> Self {
        // println!("Geometry: {:#?}", geom);
        Problem::new(geom, Some(settings))
    }

    pub fn py_solve(&mut self) -> PyResult<()>{
        self.init();
        self.solve_near();
        self.solve_far();
        Ok(())
    }

    pub fn py_print_stats(&self) -> PyResult<()> {
        println!("{}", self.powers);
        Ok(())
    }
}

impl Problem {
    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new(mut geom: Geom, settings: Option<Settings>) -> Self {

        let mut settings = settings.unwrap_or_else(settings::load_config);

        let bins = generate_bins(&settings.binning.scheme);
        let total_ampl_far_field =
            vec![Matrix2::<Complex<f32>>::zeros(); bins.len()];
        
        // rescale geometry so the max dimension is 1
        geom.recentre();
        let scale_factor = geom.rescale();
        settings.wavelength = settings.wavelength * scale_factor;
        settings.beam_power_threshold *= scale_factor.powi(2); // power scales with area
        settings.beam_area_threshold_fac *= scale_factor.powi(2); // area scales with length^2

        let problem = Self {
            geom,
            beam_queue: vec![],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            powers: Powers::new(),
            bins: bins,
            ampl: total_ampl_far_field,
            settings,
            scale_factor
        };

        problem
    }


    pub fn init(&mut self) {
        let beam = basic_initial_beam(&self.geom, self.settings.wavelength, self.settings.medium_refr_index);
        self.beam_queue.push(beam);
    }

    /// Resets the problem and reilluminates it with a basic initial beam.
    pub fn reset(&mut self) {
        self.beam_queue.clear();
        self.out_beam_queue.clear();    
        self.ext_diff_beam_queue.clear();
        self.powers = Powers::new();
        self.ampl.iter_mut().for_each(|a| a.fill(Complex::ZERO));

        let beam = basic_initial_beam(&self.geom, self.settings.wavelength, self.settings.medium_refr_index);

        self.beam_queue.push(beam);
    }

    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new_with_field(geom: Geom, beam: Beam) -> Self {
        let  settings = settings::load_config();

        let bins = generate_bins(&settings.binning.scheme);
        let total_ampl_far_field =
            vec![Matrix2::<Complex<f32>>::zeros(); bins.len()];

        Self {
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            powers: Powers::new(),
            bins: bins,
            ampl: total_ampl_far_field,
            settings,
            scale_factor: 1.0,
        }
    }

    fn diffract_outbeams(
        queue: &mut Vec<Beam>,
        bins: &[(f32, f32)],
        total_ampl_far_field: &mut [Matrix2<Complex<f32>>],
    ) {
        let ampl_far_field = queue
            .par_iter()
            .map(|outbeam| {
                outbeam.diffract(bins)
            })
            .reduce(
                || vec![Matrix2::<Complex<f32>>::zeros(); bins.len()],
                |mut acc, local| {
                    for (a, l) in acc.iter_mut().zip(local) {
                        *a += l;
                    }
                    acc
                },
            );

        // TODO move this outside the function
        for (i, ampl) in ampl_far_field.iter().enumerate() {
            total_ampl_far_field[i] += ampl;
        }
    }

    pub fn solve_far_ext_diff(&mut self) {
        Self::diffract_outbeams(
            &mut self.ext_diff_beam_queue,
            &self.bins,
            &mut self.ampl,
        );
    }

    pub fn solve_far_outbeams(&mut self) {
        Self::diffract_outbeams(
            &mut self.out_beam_queue,
            &self.bins,
            &mut self.ampl,
        );
    }
    pub fn solve_far(&mut self) {
        self.solve_far_ext_diff();
        self.solve_far_outbeams();
    }

    pub fn solve(&mut self) {
        self.solve_near();
        self.solve_far();
    }

    /// Trace beams to solve the near-field problem.
    pub fn solve_near(&mut self) {
        loop {
            if self.beam_queue.len() == 0 {
                break;
            }

            if self.powers.output / self.powers.input > self.settings.total_power_cutoff {
                // add remaining power in beam queue to missing power due to cutoff
                self.powers.trnc_cop += self.beam_queue.iter().map(|beam| beam.power() / self.scale_factor.powi(2)).sum::<f32>();
                break;
            }

            self.propagate_next();
        }
    }

    pub fn writeup(&self) {
        let mueller = output::ampl_to_mueller(&self.bins, &self.ampl);
        let _ = output::writeup(&self.bins, &mueller);
    }

    /// Propagates the next beam in the queue.
    pub fn propagate_next(&mut self) -> Option<BeamPropagation> {
        // Try to pop the next beam from the queue
        let Some(mut beam) = self.beam_queue.pop() else {
            println!("No beams left to pop!");
            return None;
        };

        // Compute the outputs by propagating the beam
        let outputs = match &mut beam.type_ {
            BeamType::Default => {
                // truncation conditions
                if beam.power() < self.settings.beam_power_threshold {
                    self.powers.trnc_energy += beam.power() / self.scale_factor.powi(2);
                    Vec::new()
                } else if beam.face.data().area.unwrap() < self.settings.beam_area_threshold() {
                    self.powers.trnc_area += beam.power() / self.scale_factor.powi(2);
                    Vec::new()
                } else if beam.variant == Some(BeamVariant::Tir) {
                    if beam.tir_count > self.settings.max_tir {
                        self.powers.trnc_ref += beam.power() / self.scale_factor.powi(2);
                        Vec::new()
                    } else {
                        match beam.propagate(
                            &mut self.geom,
                            self.settings.medium_refr_index,
                            self.settings.beam_area_threshold()
                        ) {
                            Ok((outputs,area_power_loss)) => {
                                self.powers.trnc_area += area_power_loss / self.scale_factor.powi(2);
                                outputs},
                            Err(_) => {
                                self.powers.clip_err += beam.power() / self.scale_factor.powi(2);
                                Vec::new()
                            }
                        }
                    }
                } else if beam.rec_count > self.settings.max_rec {
                    self.powers.trnc_rec += beam.power() / self.scale_factor.powi(2);
                    Vec::new()
                } else {
                    match beam.propagate(
                        &mut self.geom,
                        self.settings.medium_refr_index,
                        self.settings.beam_area_threshold()
                    ) {
                        Ok((outputs,area_power_loss)) => {
                                self.powers.trnc_area += area_power_loss / self.scale_factor.powi(2);
                                outputs},
                        Err(_) => {
                            self.powers.clip_err += beam.power() / self.scale_factor.powi(2);
                            Vec::new()
                        }
                    }
                }
            }
            BeamType::Initial => match beam.propagate(
                &mut self.geom,
                self.settings.medium_refr_index,
                self.settings.beam_area_threshold()
            ) {
                Ok((outputs,..)) => {
                outputs},

                Err(_) => {
                    Vec::new()
                }
            },
            _ => {
                println!("Unknown beam type, returning empty outputs.");
                Vec::new()
            }
        };

        self.powers.absorbed += beam.absorbed_power / self.scale_factor.powi(2);
        self.powers.trnc_clip += (beam.clipping_area - beam.csa()).abs() * beam.power() / self.scale_factor.powi(2);

        // Process each output beam
        for output in outputs.iter() {
            let output_power = output.power() / self.scale_factor.powi(2);
            match (&beam.type_, &output.type_) {
                (BeamType::Default, BeamType::Default) => self.insert_beam(output.clone()),
                (BeamType::Default, BeamType::OutGoing) => {
                    self.powers.output += output_power;
                    self.insert_outbeam(output.clone());
                }
                (BeamType::Initial, BeamType::Default) => {
                    self.powers.input += output_power;
                    self.insert_beam(output.clone());
                }
                (BeamType::Initial, BeamType::ExternalDiff) => {
                    self.powers.ext_diff += output_power;
                    self.ext_diff_beam_queue.push(output.clone());
                }
                _ => {}
            }
        }
        Some(BeamPropagation::new(beam, outputs))
    }

    /// Draws a `BeamPropagation` on top of a `Geom`.
    pub fn draw_propagation(&self, propagation: &BeamPropagation) {
        // draw the geometry
        for shape in &self.geom.shapes {
            for face in &shape.faces {
                draw_face(face, GREEN, 4.0);
            }
        }
        propagation.draw();
    }

    /// Inserts a beam into the beam queue such that beams with greatest power
    /// are prioritised for dequeueing.
    pub fn insert_beam(&mut self, beam: Beam) {
        let value = beam.power() ;

        // Find the position to insert the beam using binary search
        let pos = self
            .beam_queue
            .binary_search_by(|x| {
                x.power()
                    .partial_cmp(&value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(|e| e);

        // Insert the beam at the determined position
        self.beam_queue.insert(pos, beam);

        // Or just push
        // self.beam_queue.push(beam);
    }

    /// Inserts a beam into the outbeam queue such that beams with greatest power
    /// are prioritised for dequeueing.
    pub fn insert_outbeam(&mut self, beam: Beam) {
        let value = beam.power();

        // Find the position to insert the beam using binary search
        let pos = self
            .out_beam_queue
            .binary_search_by(|x| {
                value
                    .partial_cmp(&x.power())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(|e| e);

        // Insert the beam at the determined position
        self.out_beam_queue.insert(pos, beam);

        // Or just push
        // self.out_beam_queue.push(beam);
    }
}

/// Creates a basic initial beam for full illumination of the geometry along the z-axis.
fn basic_initial_beam(geom: &Geom, wavelength: f32, medium_refractive_index: Complex<f32>) -> Beam {
    const FAC: f32 = 1.1;
    let bounds = geom.bounds();
    let (min, max) = (bounds.0.map(|v| v * FAC), bounds.1.map(|v| v * FAC));    

    let clip_vertices = vec![
        Point3::new(max[0], max[1], max[2]),
        Point3::new(max[0], min[1], max[2]),
        Point3::new(min[0], min[1], max[2]),
        Point3::new(min[0], max[1], max[2]),
    ];

    let mut clip = Face::new_simple(clip_vertices, None).unwrap();
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

/// A problem for a single geometry with multiple orientations.
#[derive(Debug, PartialEq)] // Added Default derive
pub struct MultiProblem {
    pub geom: Geom,
    pub problems: Vec<Problem>,
    pub orientations: Orientations,
    pub bins: Vec<(f32, f32)>,            // bins for far-field diffraction
    pub mueller: Array2::<f32>, // total mueller in far-field
    pub settings: Settings,               // runtime settings
    pub powers: Powers,                   // different power contributions
}

impl MultiProblem {
    /// Creates a new `MultiOrientProblem` from a `settings: Settings` configuration.
    pub fn new(settings: Settings) -> Self {
        let mut geom = Geom::from_file(&settings.geom_name).unwrap();

        init_geom(&settings, &mut geom);

        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let problems = Vec::new();
        let bins = generate_bins(&settings.binning.scheme);
        let  mueller = Array2::<f32>::zeros((bins.len(), 16));
        let powers = Powers::new();

        Self { geom, problems, orientations, bins, mueller, settings, powers }
    }

    /// Solves a `MultiOrientProblem` by averaging over the problems.
    pub fn solve(&mut self) {

        let start = Instant::now();
        println!("Solving problem...");

        // init a base problem that can be reset
        let problem_base = Problem::new(self.geom.clone(), Some(self.settings.clone()));
        // let mut problem = problem_base.clone();

        let m = MultiProgress::new();
        let n = self.orientations.num_orientations;
        let pb = m.add(ProgressBar::new(n as u64));
        pb.set_style(
            ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] {bar:40.green/blue} {pos:>5}/{len:5} {msg} ETA: {eta_precise}",
            )
            .unwrap()
            .progress_chars("█▇▆▅▄▃▂▁")
        );
        pb.set_message("orientation".to_string());

        (self.mueller,self.powers) = self.orientations.eulers.par_iter().map(|(a,b,g)| {
            let mut problem = problem_base.clone();
            if let Err(error) = problem.geom.euler_rotate(*a, *b, *g) {
                panic!("Error rotating geometry: {}", error);
            }

            problem.init();
            problem.solve();

            let mueller = output::ampl_to_mueller(&self.bins, &problem.ampl);

            pb.inc(1);
            (mueller, problem.powers)
        }).reduce(|| (Array2::<f32>::zeros((self.bins.len(), 16)), Powers::new()), |(mut acc_mueller, mut acc_powers), (local_mueller, local_powers)| {
            for (a, l) in acc_mueller.iter_mut().zip(local_mueller) {
                *a += l;
            }
            acc_powers += local_powers;
            (acc_mueller, acc_powers)
        });

        // normalise
        for mut row in self.mueller.outer_iter_mut() {
            for val in row.iter_mut() {
                *val /= self.orientations.num_orientations as f32;
            }
        }
        self.powers /= self.orientations.num_orientations as f32;
        let end = Instant::now();
        let duration = end.duration_since(start);
        let time_per_orientation = duration / self.orientations.num_orientations as u32;

        println!(
            "Time taken: {:.2?}, Time per orientation: {:.2?}",
            duration, time_per_orientation
        );

        pb.finish_with_message(format!("(done)"));
        println!("Average {}", self.powers);
    }

    pub fn writeup(&self) {
        let _ = output::writeup(&self.bins, &self.mueller);
    }


}

/// Initialises the geometry with the refractive indices from the settings.
/// In the future, this function will be extended to provide additional checks
/// to ensure the geometry is well-defined.
fn init_geom(settings: &Settings, geom: &mut Geom) {
    for shape in geom.shapes.iter_mut() {
        shape.refr_index = settings.particle_refr_index[0]; // default refr index is first value
    }
    for (i,refr_index) in settings.particle_refr_index.iter().enumerate() {
        if i >= geom.shapes.len() {
            break;
        }
        geom.shapes[i].refr_index = *refr_index;
    }
    geom.recentre();
}