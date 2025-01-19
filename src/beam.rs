use geo::Coord;
use macroquad::prelude::*;
use macroquad::ui::Id;
use nalgebra::Vector3;

use crate::clip::Clipping;
use crate::config;
use crate::geom::Face;
use crate::geom::Geom;
use crate::geom::RefrIndex;
use crate::helpers::draw_face;
use crate::helpers::lines_to_screen;

#[derive(Debug, Clone, PartialEq)]
pub struct BeamPropagation {
    pub input: Beam,
    pub refr_index: RefrIndex,
    pub outputs: Vec<Beam>,
}

impl BeamPropagation {
    /// Makes a new `BeamPropagation` struct, which represents a beam propagation.
    pub fn new(input: Beam, outputs: Vec<Beam>) -> Self {
        let refr_index = input.data().refr_index.clone();
        Self {
            input,
            refr_index,
            outputs,
        }
    }

    /// Draws a `Beam Propagation`
    pub fn draw(&self) {
        // draw the input
        draw_face(&self.input.data().face, YELLOW, 4.0);
        // draw the outputs
        for beam in &self.outputs {
            draw_face(&beam.data().face, BLUE, 4.0);
        }
        // draw lines from the outputs to the input
        let line_strings: Vec<_> = self
            .outputs
            .iter()
            .map(|x| {
                let output_mid = x.data().face.midpoint();
                let input_mid = self.input.data().face.midpoint();
                let vec = input_mid - output_mid;
                let input_normal = self.input.data().face.data().normal;
                let norm_dist_to_plane = vec.dot(&input_normal);
                let dist_to_plane =
                    norm_dist_to_plane / (input_normal.dot(&self.input.data().proj));
                // ray cast along propagation direction
                let intsn = output_mid + dist_to_plane * self.input.data().proj;
                vec![
                    Coord {
                        x: output_mid.coords.x,
                        y: output_mid.coords.y,
                    },
                    Coord {
                        x: intsn.coords.x,
                        y: intsn.coords.y,
                    },
                ]
            })
            .collect();

        lines_to_screen(line_strings, RED, 2.0);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Beam {
    Initial(BeamData), // an initial beam to be traced in the near-field
    Default {
        data: BeamData, // a beam to be traced in the near-field
    },
    OutGoing(BeamData), // a beam to be mapped to the far-field
}

impl Beam {
    pub fn new_initial(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Beam::Initial(BeamData::new(face, proj, refr_index))
    }
    pub fn new_default(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        Beam::Default {
            data: BeamData::new(face, proj, refr_index),
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

    /// Computes the propagation of a `Beam`, yielding the output
    /// beams which can then be dealt with as needed.
    pub fn propagate(&mut self, geom: &mut Geom) -> Vec<Beam> {
        println!("===================");
        println!("propagating beam...");
        let mut outputs = Vec::new();
        match self {
            Beam::Initial(data) => {
                let outputs_beams = Self::process_beam(geom, data);
                println!("adding {} beams to the outputs", outputs_beams.len());
                outputs.extend(outputs_beams);
            }
            Beam::Default { data } => {
                let output_beams = Self::process_beam(geom, data);

                println!("adding {} beams to the outputs", output_beams.len());
                outputs.extend(output_beams);
            }
            Beam::OutGoing(..) => {
                panic!("tried to propagate an outgoing beam, which is not yet supported.");
            }
        }
        outputs
    }

    fn process_beam(geom: &mut Geom, data: &mut BeamData) -> Vec<Beam> {
        let mut output_beams = Vec::new();
        let refr_index_in = data.refr_index;
        match data.face.data().parent_id {
            Some(id) => {
                println!("beam in originated from shape #{}", id)
            }
            None => {
                println!("beam did not originate from a shape")
            }
        }

        let mut clipping = Clipping::new(geom, &mut data.face, &data.proj);
        clipping.clip(); // do the clipping -> contains the intersections

        assert!(
            clipping.stats.clone().unwrap().total_consvtn > 0.99,
            "Error, poor total energy conservation in clipping: {:?}",
            clipping.stats.clone().unwrap().total_consvtn
        );
        println!("{}", clipping.stats.unwrap());

        // remove intersections with below threshold area
        let intersections: Vec<_> = clipping
            .intersections
            .into_iter()
            .filter(|x| {
                println!("intersection... {:?}", x.data().area.unwrap());
                x.data().area.unwrap() > config::BEAM_AREA_THRESHOLD
            })
            .collect();

        // remove remainders with below threshold area
        let remainder_beams: Vec<_> = clipping
            .remaining
            .into_iter()
            .filter_map(|x| {
                println!("remainder... {:?}", x.data().area.unwrap());
                if x.data().area.unwrap() > config::BEAM_AREA_THRESHOLD {
                    None
                } else {
                    Some(Beam::OutGoing(BeamData {
                        face: x,
                        proj: data.proj,
                        refr_index: data.refr_index,
                    }))
                }
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
                let sink_shape_id = x
                    .data()
                    .parent_id
                    .expect("Sink face should always have a parent shape id.");

                println!("intersection was with shape #{}", sink_shape_id);
                // sinks have no parent id at the moment, perhaps we need to make some changes
                let source_shape_id = data.face.data().parent_id.unwrap_or_else(|| 999);
                // println!(
                //     "sink shape id: {}, source shape id: {}",
                //     sink_shape_id, source_shape_id
                // );

                // determine transmitted propagation direction
                let proj = data.proj;
                // determine other BeamData values here later...
                Some(Beam::new_default(x.clone(), proj, data.refr_index))
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
                let proj = data.proj;
                // determine other BeamData values here later...
                Some(Beam::new_default(x, proj, data.refr_index))
            })
            .collect();

        output_beams.extend(reflected_beams);
        output_beams.extend(transmitted_beams);
        output_beams.extend(remainder_beams);

        output_beams
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
