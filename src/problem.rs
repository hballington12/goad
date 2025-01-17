use crate::{
    beam::Beam,
    clip::Clipping,
    geom::{Face, Geom},
};

#[cfg(test)]
mod tests {}

const THRESHOLD_AREA: f32 = 0.1; // minimum area for new beam to propagate

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
            Beam::Default(mut data) => {
                let mut clipping = Clipping::new(&mut self.geom, &mut data.face, &data.proj);
                clipping.clip(); // do the clipping -> contains the intersections
                println!("{}", clipping.stats.unwrap());
            }
            //if initial beam, as normal but discard remainders
            Beam::Initial(mut data) => {
                let mut clipping = Clipping::new(&mut self.geom, &mut data.face, &data.proj);
                clipping.clip(); // do the clipping -> contains the intersections
                assert!(
                    clipping.stats.clone().unwrap().total_consvtn > 0.999,
                    "Error, poor total energy conservation in clipping: {:?}",
                    clipping.stats.clone().unwrap().total_consvtn
                );
                println!("{}", clipping.stats.unwrap());
                // deal with intersections
                for (i, intsn) in clipping.intersections.into_iter().enumerate() {
                    // for each intersection:
                    // filter some conditions to determine if and how new beams should be made
                    // ie. area, max recursions, max reflections, energy, total internal reflection
                    let area = intsn
                        .data()
                        .area
                        .expect("Face area of intersection should always be Some().");

                    println!("Intersection face areas: {:?}", area);

                    if area < THRESHOLD_AREA {
                        continue;
                    }

                    let map = clipping.source_mapping[i]; // determined which shape, face in the geometry the intersection was at
                    let refr_index = self.geom.shapes[map.0].refr_index; // use the shape in the geometry to retrieve the refractive index

                    // do some stuff to determine new propagation directions,
                    // electric field, amplitude matrix, etc
                    // use refractive index tree to determine external refr index when
                    // a beam intersects with the same shape from which it originates

                    // always create reflected beam (sort propagation later)
                    let proj = beam.data().proj.clone(); // use the same propagation vector for now

                    let reflected_beam = Beam::new_default(intsn, proj, refr_index);
                }
                // deal with remainders
            }
            // if outgoing, panic
            Beam::OutGoing(_) => panic!("Outgoing beams cannot be propagated."),
        }
    }
}
