use nalgebra::Vector3;

use crate::geom::Face;
use crate::geom::RefrIndex;

#[derive(Debug, PartialEq)]
pub enum Beam {
    Initial(BeamData),  // an initial beam to be traced in the near-field
    Default(BeamData),  // a beam to be traced in the near-field
    OutGoing(BeamData), // a beam to be mapped to the far-field
}

impl Beam {
    pub fn new_initial(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Beam::Default(BeamData::new(face, proj, refr_index))
    }
    pub fn new_default(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Beam::Default(BeamData::new(face, proj, refr_index))
    }
    pub fn to_outgoing(beam: Beam) -> Beam {
        match beam {
            Beam::Default(data) => Beam::OutGoing(data),
            Beam::OutGoing(_) => panic!(
                "You probably don't want to convert an outgoing beam to another outgoing beam."
            ),
            Beam::Initial(_) => {
                panic!("The initial beam should not be converted to an outgoing beam.")
            }
        }
    }
}

#[derive(Debug, PartialEq)] // Added Default derive
pub struct BeamData {
    pub face: Face,
    pub proj: Vector3<f32>,
    pub refr_index: RefrIndex,
}

/// Creates a new beam
impl BeamData {
    pub fn new(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Self {
            face,
            proj,
            refr_index,
        }
    }
}
