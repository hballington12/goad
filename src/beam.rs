use macroquad::prelude::*;
use nalgebra::Vector3;

use crate::clip::Clipping;
use crate::config;
use crate::geom::Face;
use crate::geom::Geom;
use crate::geom::RefrIndex;
use crate::helpers::draw_face;

#[derive(Debug, Clone, PartialEq)]
pub struct BeamPropagation {
    pub input: BeamData,
    pub refr_index: RefrIndex,
    pub outputs: Vec<Beam>,
}

impl BeamPropagation {
    /// Makes a new `BeamPropagation` struct.
    pub fn new(source: BeamData) -> Self {
        let refr_index = source.refr_index.clone();
        Self {
            input: source,
            refr_index,
            outputs: Vec::new(),
        }
    }

    /// Computes the propagation of a `BeamPropagation`, yielding the output
    /// beams which can then be dealt with as needed.
    pub fn propagate(&mut self, geom: &mut Geom) {
        // create a clipping to clip the beam against the geometry
        let mut clipping = Clipping::new(geom, &mut self.input.face, &self.input.proj);
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
            .filter(|x| {
                println!("{:?}", x.face.data().area.unwrap());
                x.face.data().area.unwrap() > config::BEAM_AREA_THRESHOLD
            })
            .collect();

        println!("number of 'good' intersections: {}", intersections.len());

        // create transmitted beams
        let transmitted_beams: Vec<_> = intersections
            .iter()
            .filter_map(|x| {
                // apply filtering here -> return None if fail
                // total internal reflection -> None
                // max recursions reached -> None
                // determine refractive index
                let sink_shape_id = x.shape_id;
                // sinks have no parent id at the moment, perhaps we need to make some changes
                let source_shape_id = self.input.face.data().parent_id.unwrap_or_else(|| 999);
                // println!(
                //     "sink shape id: {}, source shape id: {}",
                //     sink_shape_id, source_shape_id
                // );

                // determine transmitted propagation direction
                let proj = self.input.proj;
                // determine other BeamData values here later...
                Some(Beam::new_default(
                    x.face.clone(),
                    proj,
                    self.refr_index,
                    x.shape_id,
                ))
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
                let proj = self.input.proj;
                // determine other BeamData values here later...
                Some(Beam::new_default(x.face, proj, self.refr_index, x.shape_id))
            })
            .collect();

        // println!("number of reflected beams: {:?}", reflected_beams.len());
        // println!("number of transmitted beams: {:?}", reflected_beams.len());
        // println!("reflected beams: {:?}", reflected_beams);

        // add the reflected and transmitted beams to the sinks
        println!("note: reflected beams are currently omitted");
        // self.outputs.extend(reflected_beams);
        println!("adding {} beams to the outputs", transmitted_beams.len());
        self.outputs.extend(transmitted_beams);
    }

    /// Draws a `Beam Propagation`
    pub fn draw(&self) {
        // draw the input
        draw_face(&self.input.face, YELLOW, 4.0);
        // draw the outputs
        for beam in &self.outputs {
            draw_face(&beam.data().face, BLUE, 4.0);
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Beam {
    Initial(BeamData), // an initial beam to be traced in the near-field
    Default {
        data: BeamData,  // a beam to be traced in the near-field
        shape_id: usize, // each default beam should also be associated with a shape in a geometry
    },
    OutGoing(BeamData), // a beam to be mapped to the far-field
}

impl Beam {
    pub fn new_initial(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Beam::Initial(BeamData::new(face, proj, refr_index))
    }
    pub fn new_default(
        face: Face,
        proj: Vector3<f32>,
        refr_index: RefrIndex,
        shape_id: usize,
    ) -> Self {
        Beam::Default {
            data: BeamData::new(face, proj, refr_index),
            shape_id,
        }
    }
    pub fn to_outgoing(beam: Beam) -> Beam {
        match beam {
            Beam::Default { data, .. } => Beam::OutGoing(data),
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
            Beam::Default { data, .. } => data,
            Beam::OutGoing(data) => data,
        }
    }
}

#[derive(Debug, Clone, PartialEq)] // Added Default derive
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
