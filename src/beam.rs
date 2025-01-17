use nalgebra::Vector3;

use crate::clip::Clipping;
use crate::config;
use crate::geom::Face;
use crate::geom::Geom;
use crate::geom::RefrIndex;

#[derive(Debug, PartialEq)]
pub struct BeamPropagation {
    pub source: BeamData,
    pub refr_index: RefrIndex,
    pub sinks: Vec<Beam>,
}

impl BeamPropagation {
    /// Makes a new `BeamPropagation` struct.
    pub fn new(source: BeamData) -> Self {
        let refr_index = source.refr_index.clone();
        Self {
            source,
            refr_index,
            sinks: Vec::new(),
        }
    }

    pub fn propagate(&mut self, geom: &mut Geom) {
        // create a clipping to clip the beam against the geometry
        let mut clipping = Clipping::new(geom, &mut self.source.face, &self.source.proj);
        clipping.clip(); // do the clipping -> contains the intersections

        assert!(
            clipping.stats.clone().unwrap().total_consvtn > 0.999,
            "Error, poor total energy conservation in clipping: {:?}",
            clipping.stats.clone().unwrap().total_consvtn
        );
        println!("{}", clipping.stats.unwrap());

        // remove intersections with very small area
        let intersections: Vec<_> = clipping
            .intersections
            .into_iter()
            .filter(|x| x.data().area.unwrap() > config::BEAM_AREA_THRESHOLD)
            .collect();

        // create transmitted beams
        let transmitted_beams: Vec<_> = intersections
            .iter()
            .filter_map(|x| {
                // apply filtering here -> return None if fail
                // total internal reflection -> None
                // max recursions reached -> None
                // determine refractive index
                // if x.

                // determine transmitted propagation direction
                let proj = self.source.proj;
                // determine other BeamData values here later...
                Some(Beam::new_default(x.clone(), proj, self.refr_index))
            })
            .collect();

        // create reflected beams
        let reflected_beams: Vec<_> = intersections
            .into_iter() // can consume the intersections here to avoid cloning
            .filter_map(|x| {
                // apply filtering here -> return None if fail
                // total internal reflection && max tir reached -> None
                // max recursions reached -> None
                // refractive index of reflected beam is always the same as the source
                // determine reflected propagation direction
                let proj = self.source.proj;
                // determine other BeamData values here later...
                Some(Beam::new_default(x, proj, self.refr_index))
            })
            .collect();

        // println!("number of reflected beams: {:?}", reflected_beams.len());
        // println!("number of transmitted beams: {:?}", reflected_beams.len());
        // println!("reflected beams: {:?}", reflected_beams);

        // add the reflected and transmitted beams to the sinks
        self.sinks.extend(reflected_beams);
        self.sinks.extend(transmitted_beams);
    }
}

#[derive(Debug, PartialEq)]
pub enum Beam {
    Initial(BeamData),  // an initial beam to be traced in the near-field
    Default(BeamData),  // a beam to be traced in the near-field
    OutGoing(BeamData), // a beam to be mapped to the far-field
}

impl Beam {
    pub fn new_initial(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Beam::Initial(BeamData::new(face, proj, refr_index))
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
                panic!("An initial beam should not be converted to an outgoing beam.")
            }
        }
    }
    pub fn data(&self) -> &BeamData {
        match self {
            Beam::Initial(data) => data,
            Beam::Default(data) => data,
            Beam::OutGoing(data) => data,
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
