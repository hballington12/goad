use crate::{
    beam::{Beam, BeamPropagation, BeamType, BeamVariant},
    bins,
    field::Field,
    geom::{Face, Geom},
    helpers::draw_face,
    output,
    settings::{self, Settings},
};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};
use rayon::prelude::*;
use std::{fmt, sync::Arc};

use clap::Parser;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

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

        let mut problem = Problem::new(geom);

        problem.propagate_next();
    }
}

#[derive(Debug, PartialEq)]
pub struct Powers {
    pub input: f32,       // near-field input power
    pub output: f32,      // near-field output power
    pub absorbed: f32,    // near-field absorbed power
    pub trnc_ref: f32,    // truncated power due to max reflections
    pub trnc_rec: f32,    // truncated power due to max recursions
    pub trnc_clip: f32,   // truncated power due to clipping
    pub trnc_energy: f32, // truncated power due to threshold beam power
    pub trnc_area: f32,   // truncated power due to area threshold
    pub ext_diff: f32,    // external diffraction power
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
            trnc_area: 0.0,
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
                + self.trnc_energy)
    }
}

impl fmt::Display for Powers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Powers:")?;
        writeln!(f, "  Input:            {:.6}", self.input)?;
        writeln!(f, "  Output:           {:.6}", self.output)?;
        writeln!(f, "  Absorbed:         {:.6}", self.absorbed)?;
        writeln!(f, "  Truncated Ref:    {:.6}", self.trnc_ref)?;
        writeln!(f, "  Truncated Rec:    {:.6}", self.trnc_rec)?;
        // writeln!(f, "  Truncated Clip:   {:.6}", self.trnc_clip)?;
        writeln!(f, "  Truncated Energy: {:.6}", self.trnc_energy)?;
        writeln!(f, "  Truncated Area:   {:.6}", self.trnc_area)?;
        writeln!(f, "  Other:            {:.6}", self.missing())?;
        writeln!(f, "  External Diff:    {:.6}", self.ext_diff)
    }
}

/// A solvable physics problem.
#[derive(Debug, PartialEq)] // Added Default derive
pub struct Problem {
    pub geom: Geom,                       // geometry to trace beams in
    pub beam_queue: Vec<Beam>,            // beams awaiting near-field propagation
    pub out_beam_queue: Vec<Beam>,        // beams awaiting diffraction
    pub ext_diff_beam_queue: Vec<Beam>,   // beams awaiting external diffraction
    pub powers: Powers,                   // different power contributions
    pub bins: Vec<(f32, f32)>,            // bins for far-field diffraction
    pub ampl: Vec<Matrix2<Complex<f32>>>, // total amplitude in far-field
    pub settings: Settings,               // runtime settings
}

impl Problem {
    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new(geom: Geom) -> Self {
        let mut settings = settings::load_config();
        // let args = settings::CliArgs::parse();

        println!("Settings: {:#?}", settings);

        let theta_phi_combinations = bins::generate_theta_phi_combinations(
            settings.far_field_resolution,
            settings.far_field_resolution,
        );
        let total_ampl_far_field =
            vec![Matrix2::<Complex<f32>>::zeros(); theta_phi_combinations.len()];

        // rescale geometry so the max dimension is 1
        // let mut geom = geom.clone();
        // let scaling_factor = geom.rescale();
        // settings.wavelength = settings.wavelength / scaling_factor;

        let beam = basic_initial_beam(&geom, settings.wavelength, settings.medium_refr_index);

        let problem = Self {
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            powers: Powers::new(),
            bins: theta_phi_combinations,
            ampl: total_ampl_far_field,
            settings,
        };

        problem
    }

    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new_with_field(geom: Geom, beam: Beam) -> Self {
        let mut settings = settings::load_config();
        // let args = settings::CliArgs::parse();

        let theta_phi_combinations = bins::generate_theta_phi_combinations(
            settings.far_field_resolution,
            settings.far_field_resolution,
        );
        let total_ampl_far_field =
            vec![Matrix2::<Complex<f32>>::zeros(); theta_phi_combinations.len()];

        Self {
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            powers: Powers::new(),
            bins: theta_phi_combinations,
            ampl: total_ampl_far_field,
            settings,
        }
    }

    fn diffract_outbeams(
        queue: &mut Vec<Beam>,
        theta_phi_combinations: &[(f32, f32)],
        total_ampl_far_field: &mut [Matrix2<Complex<f32>>],
        description: &str,
    ) {
        let m = MultiProgress::new();

        let n = queue.len();
        let pb = m.add(ProgressBar::new(n as u64));
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] {bar:40.green/blue} {pos:>5}/{len:5} {msg}",
            )
            .unwrap()
            .progress_chars("➤➤➤➤➤➤➤➤"),
        );
        pb.set_message(description.to_string());
        let pb2 = m.add(ProgressBar::new(1000 as u64));
        pb2.set_style(
            ProgressStyle::with_template(
                "{spinner:.red} [{elapsed_precise}] {bar:40.yellow/blue} {percent:>6.1.white}%     {msg}",
            )
            .unwrap()
            .progress_chars("➤➤➤➤➤➤➤➤"),
        );
        pb2.set_message("power diffracted".to_string());

        let total_power = queue.iter().map(|beam| beam.power()).sum::<f32>();

        let ampl_far_field = queue
            .par_iter()
            .map(|outbeam| {
                let outbeam_power = outbeam.power();
                pb.inc(1);
                pb2.inc((outbeam_power / total_power * 1000.0) as u64);
                outbeam.diffract(theta_phi_combinations)
            })
            .reduce(
                || vec![Matrix2::<Complex<f32>>::zeros(); theta_phi_combinations.len()],
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

        queue.clear();
        let _ = output::writeup(&theta_phi_combinations, &total_ampl_far_field);
        pb.finish_with_message(format!("{} (done)", description));
        pb2.finish_with_message(format!("power diffracted (%)"));
        // total write
    }

    pub fn solve_far_ext_diff(&mut self) {
        Self::diffract_outbeams(
            &mut self.ext_diff_beam_queue,
            &self.bins,
            &mut self.ampl,
            "external diffraction",
        );
    }

    pub fn solve_far_outbeams(&mut self) {
        Self::diffract_outbeams(
            &mut self.out_beam_queue,
            &self.bins,
            &mut self.ampl,
            "outgoing beams",
        );
    }
    pub fn solve_far(&mut self) {
        self.solve_far_ext_diff();
        self.solve_far_outbeams();
    }

    /// Trace beams to solve the near-field problem.
    pub fn solve_near(&mut self) {
        loop {
            if self.beam_queue.len() == 0 {
                println!("all beams traced...");
                break;
            }

            if self.powers.output / self.powers.input > self.settings.total_power_cutoff {
                println!("cut off power out reached...");
                break;
            }

            self.propagate_next();
        }

        println!("{}", self.powers);
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
                    self.powers.trnc_energy += beam.power();
                    Vec::new()
                } else if beam.face.data().area.unwrap() < self.settings.beam_area_threshold() {
                    self.powers.trnc_area += beam.power();
                    Vec::new()
                } else if beam.variant == Some(BeamVariant::Tir) {
                    if beam.tir_count > self.settings.max_tir {
                        self.powers.trnc_ref += beam.power();
                        Vec::new()
                    } else {
                        beam.propagate(
                            &mut self.geom,
                            self.settings.medium_refr_index,
                        )
                    }
                } else if beam.rec_count > self.settings.max_rec {
                    self.powers.trnc_rec += beam.power();
                    Vec::new()
                } else {
                    beam.propagate(
                        &mut self.geom,
                        self.settings.medium_refr_index,
                    )
                }
            }
            BeamType::Initial => beam.propagate(
                &mut self.geom,
                self.settings.medium_refr_index,
            ),
            _ => {
                println!("Unknown beam type, returning empty outputs.");
                Vec::new()
            }
        };

        self.powers.absorbed += beam.absorbed_power;
        self.powers.trnc_clip += (beam.clipping_area - beam.csa()).abs() * beam.power();

        // Process each output beam
        for output in outputs.iter() {
            let output_power = output.power();
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
        let value = beam.power();

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
