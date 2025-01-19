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
    use crate::geom::RefrIndex;
    use nalgebra::Vector3;

    #[test]
    fn cube_inside_ico() {
        use nalgebra::Point3;

        let mut geom = Geom::from_file("./examples/data/cube_inside_ico.obj").unwrap();
        geom.shapes[0].refr_index = RefrIndex {
            // modify the refractive index of the outer shape
            real: 2.0,
            imag: 0.1,
        };
        geom.shapes[1].parent_id = Some(0); // set the parent of the second shape to be the first
        assert_ne!(geom.shapes[0].refr_index, geom.shapes[1].refr_index);
        assert_eq!(
            geom.shapes[0],
            geom.shapes[geom.shapes[1].parent_id.unwrap()]
        );

        let projection = Vector3::new(0.0, 0.0, -1.0);

        let lower_left = vec![-10.0, -3.0];
        let upper_right = vec![10.0, 3.0];
        let clip_vertices = vec![
            Point3::new(lower_left[0], upper_right[1], 10.0),
            Point3::new(lower_left[0], lower_left[1], 10.0),
            Point3::new(upper_right[0], lower_left[1], 10.0),
            Point3::new(upper_right[0], upper_right[1], 10.0),
        ];
        let clip = Face::new_simple(clip_vertices, None);

        let mut problem = Problem::new(
            geom,
            Beam::new_initial(clip, projection, RefrIndex::new(1.31, 0.1)),
        );

        problem.propagate_next();
        assert!(false)
    }
}

/// A solvable physics problem.
#[derive(Debug, PartialEq)] // Added Default derive
pub struct Problem {
    pub geom: Geom,            // geometry to trace beams in
    pub beam_queue: Vec<Beam>, // beams which need to be propagated
}

impl Problem {
    /// Creates a new `Problem` from a `Geom` and an initial `Beam`.
    pub fn new(geom: Geom, beam: Beam) -> Self {
        Self {
            geom,
            beam_queue: vec![beam],
        }
    }

    /// Propagates the next beam in the queue.
    pub fn propagate_next(&mut self) -> Option<BeamPropagation> {
        // try to pop the next beam
        if let Some(beam) = self.beam_queue.pop() {
            let propagation = self.propagate(beam);
            Some(propagation)
        } else {
            panic!("Tried to pop() beam but there were no beams to pop.");
            None
        }
    }

    /// Propagates a beam.
    fn propagate(&mut self, beam: Beam) -> BeamPropagation {
        // do beam propagation

        println!("-------------------");
        println!("propagating beam...");

        // make clipping object from the beam
        match beam {
            // if default beam, normal
            Beam::Default { data, .. } => {
                let mut propagation = BeamPropagation::new(data);
                propagation.propagate(&mut self.geom);
                self.beam_queue.extend(propagation.outputs.clone());
                propagation
            }
            // if initial beam, as normal but discard remainders
            Beam::Initial(data) => {
                // create a beam propagation
                let mut propagation = BeamPropagation::new(data);
                propagation.propagate(&mut self.geom); // propagate the beam

                // clone below not strictly needed, but probably acceptable so the propagation can be returned for debugging
                self.beam_queue.extend(propagation.outputs.clone()); // add outputs to the beam queue
                propagation // return
            }
            // if outgoing, panic
            Beam::OutGoing(_) => panic!("Outgoing beams cannot be propagated."),
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
}
