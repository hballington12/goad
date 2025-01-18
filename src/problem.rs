use crate::{
    beam::{Beam, BeamPropagation},
    clip::Clipping,
    geom::{Face, Geom},
};

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
    pub fn propagate_next(&mut self) {
        // try to pop the next beam
        if let Some(beam) = self.beam_queue.pop() {
            self.propagate(beam);
        } else {
            panic!("Tried to pop() beam but there were no beams to pop.");
        }
    }

    fn propagate(&mut self, mut beam: Beam) {
        // do beam propagation

        // make clipping object from the beam
        match beam {
            // if default beam, normal
            Beam::Default { mut data, .. } => {
                let mut clipping = Clipping::new(&mut self.geom, &mut data.face, &data.proj);
                clipping.clip(); // do the clipping -> contains the intersections
                println!("{}", clipping.stats.unwrap());
            }
            // if initial beam, as normal but discard remainders
            Beam::Initial(mut beam_data) => {
                // create a beam propagation
                let mut propagation = BeamPropagation::new(beam_data);
                propagation.propagate(&mut self.geom); // propagate the beam

                // // deal with intersections
                // for (i, intsn) in clipping.intersections.into_iter().enumerate() {
                //     // for each intersection:
                //     // filter some conditions to determine if and how new beams should be made
                //     // ie. area, max recursions, max reflections, energy, total internal reflection
                //     let area = intsn
                //         .data()
                //         .area
                //         .expect("Face area of intersection should always be Some().");

                //     println!("Intersection face areas: {:?}", area);

                //     if area < THRESHOLD_AREA {
                //         continue;
                //     }

                //     let map = clipping.source_mapping[i]; // determined which shape, face in the geometry the intersection was at
                //                                           // let refr_index = self.geom.shapes[map.0].refr_index; // use the shape in the geometry to retrieve the refractive index

                //     // do some stuff to determine new propagation directions,
                //     // electric field, amplitude matrix, etc
                //     // use refractive index tree to determine external refr index when
                //     // a beam intersects with the same shape from which it originates

                //     // always create reflected beam (sort propagation later)
                //     // let proj = beam.data().proj.clone(); // use the same propagation vector for now

                //     // let reflected_beam = Beam::new_default(intsn, proj, refr_index);
                // }
                // deal with remainders
            }
            // if outgoing, panic
            Beam::OutGoing(_) => panic!("Outgoing beams cannot be propagated."),
        }
    }
}
