use crate::{beam::Beam, clip::Clipping, geom::Geom};

#[cfg(test)]
mod tests {}

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

    fn propagate(&mut self, beam: Beam) {
        // do beam propagation

        // make clipping object from the beam
        match beam {
            // if default beam, normal
            Beam::Default(mut data) => {
                let mut clipping = Clipping::new(&mut self.geom, &mut data.face, &data.proj);
                clipping.clip(); // do the clipping -> contains the intersections
                println!("{}", clipping.stats.unwrap());
                // deal with intersections
                for intsn in clipping.intersections {
                    // for each intersection:
                    // filter some conditions to determine if and how new beams should be made
                    // ie. area, max recursions, max reflections, energy, total internal reflection
                }
                // deal with remainders
            }
            //if initial beam, as normal but discard remainders
            Beam::Initial(mut data) => {
                let mut clipping = Clipping::new(&mut self.geom, &mut data.face, &data.proj);
                clipping.clip(); // do the clipping -> contains the intersections
                println!("{}", clipping.stats.unwrap());
            }
            // if outgoing, panic
            Beam::OutGoing(_) => panic!("Outgoing beams cannot be propagated."),
        }
    }
}
