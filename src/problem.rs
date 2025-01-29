use crate::{
    beam::{Beam, BeamPropagation, BeamType, BeamVariant},
    bins, config, diff,
    geom::{Face, Geom},
    helpers::draw_face,
    output,
};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2};
use std::fmt;

#[cfg(test)]
mod tests {

    use super::*;
    use nalgebra::{Complex, Vector3};

    #[test]
    fn cube_inside_ico() {
        use nalgebra::Point3;

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

        let projection = Vector3::new(0.0, 0.0, -1.0);
        let e_perp = Vector3::x(); // choose e_perp along x-axis for now

        let lower_left = vec![-10.0, -3.0];
        let upper_right = vec![10.0, 3.0];
        let clip_vertices = vec![
            Point3::new(lower_left[0], upper_right[1], 10.0),
            Point3::new(lower_left[0], lower_left[1], 10.0),
            Point3::new(upper_right[0], lower_left[1], 10.0),
            Point3::new(upper_right[0], upper_right[1], 10.0),
        ];
        let mut clip = Face::new_simple(clip_vertices, None).unwrap();
        clip.data_mut().area =
            Some((upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]));

        let mut problem = Problem::new(
            geom,
            Beam::new_initial(clip, projection, Complex::new(1.31, 0.1), e_perp).unwrap(),
        );

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
        writeln!(f, "  Unaccounted:      {:.6}", self.missing())?;
        writeln!(f, "  External Diff:    {:.6}", self.ext_diff)
    }
}

/// A solvable physics problem.
#[derive(Debug, PartialEq)] // Added Default derive
pub struct Problem {
    pub geom: Geom,                     // geometry to trace beams in
    pub beam_queue: Vec<Beam>,          // beams awaiting near-field propagation
    pub out_beam_queue: Vec<Beam>,      // beams awaiting diffraction
    pub ext_diff_beam_queue: Vec<Beam>, // beams awaiting external diffraction
    pub powers: Powers,                 // different power contributions
}

impl Problem {
    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new(geom: Geom, beam: Beam) -> Self {
        Self {
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            powers: Powers::new(),
        }
    }

    pub fn solve_far(&mut self) {
        println!("solving far-field problem...");
        // set up bins
        let theta_phi_combinations = bins::generate_theta_phi_combinations();
        let mut total_ampl_far_field =
            vec![Matrix2::<Complex<f32>>::zeros(); theta_phi_combinations.len()];

        loop {
            if self.out_beam_queue.len() == 0 {
                println!("no beams remaining to solve far-field problem...");
                break;
            }

            let outbeam = self.out_beam_queue.pop().unwrap();
            match &outbeam.face {
                Face::Simple(face) => {
                    // println!("simple face, ready for diffraction...");
                    // println!("Press any key to start...");
                    // let _ = std::io::stdin().read_line(&mut String::new());

                    let verts = &face.exterior;
                    let ampl = outbeam.field.ampl;
                    let prop = outbeam.prop;
                    let vk7 = outbeam.field.e_perp;
                    // println!("Vertices: {:?}", verts);
                    // println!("Amplitude: {:?}", ampl);
                    // println!("Propagation: {:?}", prop);
                    // println!("E_perp: {:?}", vk7);
                    let ampl_far_field =
                        diff::diffraction(verts, ampl, prop, vk7, &theta_phi_combinations);

                    for (i, ampl) in ampl_far_field.iter().enumerate() {
                        total_ampl_far_field[i] += ampl;
                    }

                    let _ = output::writeup(&theta_phi_combinations, &total_ampl_far_field);
                    // let _ = output::writeup(&theta_phi_combinations, &ampl_far_field);
                    // println!("done.");
                }
                Face::Complex { .. } => {
                    println!("complex face not supported yet...");
                }
            }
        }
    }

    /// Trace beams to solve the near-field problem.
    pub fn solve_near(&mut self) {
        loop {
            if self.beam_queue.len() == 0 {
                println!("all beams traced...");
                break;
            }

            if self.powers.output / self.powers.input > config::TOTAL_POWER_CUTOFF {
                println!("cut off power out reached...");
                break;
            }

            self.propagate_next();
        }

        println!("done.");
        println!("{}", self.powers);
        println!(
            "number of beams awaiting diffraction: {}",
            self.out_beam_queue.len()
        );
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
                if beam.power() < config::BEAM_POWER_THRESHOLD {
                    self.powers.trnc_energy += beam.power();
                    Vec::new()
                } else if beam.variant == Some(BeamVariant::Tir) {
                    if beam.tir_count > config::MAX_TIR {
                        self.powers.trnc_ref += beam.power();
                        Vec::new()
                    } else {
                        beam.propagate(&mut self.geom)
                    }
                } else if beam.rec_count > config::MAX_REC {
                    self.powers.trnc_rec += beam.power();
                    Vec::new()
                } else {
                    beam.propagate(&mut self.geom)
                }
            }
            BeamType::Initial => beam.propagate(&mut self.geom),
            _ => Vec::new(),
        };

        self.powers.absorbed += beam.absorbed_power;
        self.powers.trnc_clip += (beam.clipping_area - beam.csa()).abs() * beam.power();

        // Process each output beam
        for output in &outputs {
            match (&beam.type_, &output.type_) {
                // Handle Default beams
                (BeamType::Default, BeamType::Default) => {
                    self.insert_beam(output.clone());
                }
                (BeamType::Default, BeamType::OutGoing) => {
                    self.powers.output += output.power();
                    self.insert_outbeam(output.clone());
                }

                // Handle Initial beams
                (BeamType::Initial, BeamType::Default) => {
                    self.powers.input += output.power();
                    self.insert_beam(output.clone());
                }
                (BeamType::Initial, BeamType::ExternalDiff) => {
                    self.powers.ext_diff += output.power();
                    self.ext_diff_beam_queue.push(output.clone());
                }

                _ => {} // Ignore other cases
            }
        }

        // Create and return the propagation result
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
                x.power()
                    .partial_cmp(&value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(|e| e);

        // Insert the beam at the determined position
        self.out_beam_queue.insert(pos, beam);

        // Or just push
        // self.beam_queue.push(beam);
    }
}
