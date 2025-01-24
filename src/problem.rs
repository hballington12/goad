use crate::{
    beam::{Beam, BeamPropagation},
    clip::Clipping,
    geom::{Face, Geom},
    helpers::draw_face,
};
use macroquad::prelude::*;

#[cfg(test)]
mod tests {
    use std::cell::Ref;

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

/// A solvable physics problem.
#[derive(Debug, PartialEq)] // Added Default derive
pub struct Problem {
    pub geom: Geom,                // geometry to trace beams in
    pub beam_queue: Vec<Beam>,     // beams awaiting near-field propagation
    pub out_beam_queue: Vec<Beam>, // beams awaiting diffraction
    pub power_in: f32,             // power in, excluding remainders from any initial beams
    pub power_out: f32,            // power out in the near-field
}

impl Problem {
    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new(geom: Geom, beam: Beam) -> Self {
        println!("geom struct: {:?}", geom.containment_graph);
        Self {
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            power_in: 0.0,
            power_out: 0.0,
        }
    }

    /// Propagates the next beam in the queue.
    pub fn propagate_next(&mut self) -> Option<BeamPropagation> {
        // try to pop the next beam
        if let Some(mut beam) = self.beam_queue.pop() {
            let outputs = beam.propagate(&mut self.geom);

            // Process each output beam
            for output in outputs.iter() {
                match (&beam, output) {
                    // Handle Default beams
                    (Beam::Default { .. }, Beam::Default { .. }) => {
                        self.insert_beam(output.clone());
                    }
                    (Beam::Default { .. }, Beam::OutGoing(..)) => {
                        self.power_out += output.data().power();
                        self.out_beam_queue.push(output.clone());
                    }

                    // Handle Initial beams
                    (Beam::Initial(..), Beam::Default { .. }) => {
                        self.power_in += output.data().power();
                        self.insert_beam(output.clone());
                    }

                    _ => {}
                }
            }

            let propagation = BeamPropagation::new(beam, outputs);
            // println!(
            //     "this power in: {}, out: {}, conservation: {}",
            //     propagation.input_power(),
            //     propagation.output_power(),
            //     propagation.output_power() / propagation.input_power()
            // );
            // println!(
            //     "number of output beams in this propagation: {}",
            //     propagation.outputs.len()
            // );
            // println!("{:?}", propagation.input);
            println!(
                "total power in: {}, out: {}, conservation: {}",
                self.power_in,
                self.power_out,
                self.power_out / self.power_in
            );
            Some(propagation)
        } else {
            println!("no beams left to pop!");
            None
        }
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
        let value = beam.data().power();

        // Find the position to insert the beam using binary search
        let pos = self
            .beam_queue
            .binary_search_by(|x| {
                x.data()
                    .power()
                    .partial_cmp(&value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(|e| e);

        // Insert the beam at the determined position
        self.beam_queue.insert(pos, beam);

        // Or just push
        // self.beam_queue.push(beam);
    }
}
