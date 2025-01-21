use std::cell::Ref;
use std::result;

use geo::Coord;
use macroquad::conf;
use macroquad::prelude::*;
use macroquad::ui::Id;
use nalgebra::Vector3;

use crate::beam;
use crate::clip::Clipping;
use crate::config;
use crate::geom::Face;
use crate::geom::Geom;
use crate::geom::RefrIndex;
use crate::helpers::draw_face;
use crate::helpers::lines_to_screen;
use crate::snell::get_sin_theta_t;

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
                    norm_dist_to_plane / (input_normal.dot(&self.input.data().prop));
                // ray cast along propagation direction
                let intsn = output_mid + dist_to_plane * self.input.data().prop;
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
                println!("beam was outgoing. skipping...");
                // panic!("tried to propagate an outgoing beam, which is not yet supported.");
            }
        }
        outputs
    }

    fn process_beam(geom: &mut Geom, beam_data: &mut BeamData) -> Vec<Beam> {
        let mut output_beams = Vec::new();
        let n1 = beam_data.refr_index;
        let input_shape_id = beam_data.face.data().shape_id; // Option(shape id)

        match input_shape_id {
            Some(id) => {
                println!("beam in originated from shape #{}", id)
            }
            None => {
                println!("beam did not originate from a shape")
            }
        }

        let mut clipping = Clipping::new(geom, &mut beam_data.face, &beam_data.prop);
        clipping.clip(); // do the clipping -> contains the intersections

        // assert!(
        //     clipping.stats.clone().unwrap().total_consvtn > 0.99,
        //     "Error, poor total energy conservation in clipping: {:?}",
        //     clipping.stats.clone().unwrap().total_consvtn
        // );
        println!("{}", clipping.stats.unwrap());

        let remainder_beams: Vec<_> = clipping
            .remaining
            .into_iter()
            .filter_map(|x| {
                if x.data().area.unwrap() > config::BEAM_AREA_THRESHOLD {
                    None // remove below threshold
                } else {
                    Some(Beam::OutGoing(BeamData {
                        face: x,
                        prop: beam_data.prop,
                        refr_index: beam_data.refr_index,
                    }))
                }
            })
            .collect();

        // remove intersections with below threshold area
        let intersections: Vec<_> = clipping
            .intersections
            .into_iter()
            .filter(|x| {
                println!("intersection... {:?}", x.data().area.unwrap());
                x.data().area.unwrap() > config::BEAM_AREA_THRESHOLD
            })
            .collect();

        println!("number of 'good' intersections: {}", intersections.len());

        // create  beams
        let beams: Vec<_> = intersections
            .iter()
            .filter_map(|x| {
                let normal = x.data().normal;
                let theta_i = normal.dot(&beam_data.prop).abs().acos();
                let n2 = Self::get_n2(geom, x.data().shape_id.unwrap(), input_shape_id);

                let transmitted = if theta_i > (n2.real / n1.real).asin() {
                    None
                } else {
                    let stt = get_sin_theta_t(theta_i, n1, n2); // sin(theta_t)
                    let prop = Self::get_refraction_vector(&normal, &beam_data.prop, stt);
                    assert!(beam_data.prop.dot(&prop) > 0.0);
                    Some(Beam::new_default(x.clone(), prop, n2))
                };

                let reflected = {
                    let prop = Self::get_reflection_vector(&normal, &beam_data.prop);
                    Some(Beam::new_default(x.clone(), prop, n1))
                };

                // max recursions reached -> None

                Some((reflected, transmitted))

                // determine other BeamData values here later...
            })
            .into_iter()
            .flat_map(|(refl, trans)| refl.into_iter().chain(trans))
            .collect();

        output_beams.extend(beams);
        output_beams.extend(remainder_beams);

        output_beams
    }

    // Helper function to get the refractive index n2 based on the conditions
    fn get_n2(geom: &Geom, output_shape_id: usize, input_shape_id: Option<usize>) -> RefrIndex {
        let ni = geom.shapes[output_shape_id].refr_index; // refr index inside x

        let no = geom
            .containment_graph
            .get_parent(output_shape_id)
            .map_or(config::MEDIUM_REFR_INDEX, |parent_id| {
                geom.shapes[parent_id].refr_index
            });

        let going_into = match input_shape_id {
            None => true, // If no input shape, it's going into the second medium
            Some(id) if id == output_shape_id => false, // If input and output shape are the same, it's going out
            Some(id) => geom.containment_graph.get_parent(output_shape_id) == Some(id), // Check for parent-child relationship
        };

        if going_into {
            ni
        } else {
            no
        }
    }

    /// Returns a transmitted propagation vector, where `stt` is the sine of the angle of transmission.
    fn get_refraction_vector(norm: &Vector3<f32>, prop: &Vector3<f32>, stt: f32) -> Vector3<f32> {
        if stt < 0.0001 {
            return *prop;
        }
        // upward facing normal
        let n = if norm.dot(&prop) > 0.0 {
            *norm
        } else {
            *norm * -1.0
        };

        let w = n.dot(&prop); // cos theta_i
        let wvz = w / w.abs();
        let wf = stt;
        let ctt = (1.0 - stt.powi(2)).sqrt();
        let mut result = wf * (prop - w * n) + wvz * ctt * n;

        result.normalize_mut();

        result
    }

    fn get_reflection_vector(norm: &Vector3<f32>, prop: &Vector3<f32>) -> Vector3<f32> {
        // upward facing normal
        let n = if norm.dot(&prop) > 0.0 {
            *norm
        } else {
            *norm * -1.0
        };
        let w = n.dot(&prop); // cos theta_i
        let mut result = prop - 2.0 * w * n;
        result.normalize_mut();
        assert!((result.dot(&n) - w) < 0.01);
        result
    }
}

#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct BeamData {
    pub face: Face,
    pub prop: Vector3<f32>,
    pub refr_index: RefrIndex,
}

/// Creates a new beam
impl BeamData {
    pub fn new(face: Face, proj: Vector3<f32>, refr_index: RefrIndex) -> Self {
        let prop = proj.normalize();
        Self {
            face,
            prop: prop,
            refr_index,
        }
    }
}
