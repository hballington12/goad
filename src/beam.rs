use std::f32::consts::PI;

use geo::Coord;
use macroquad::prelude::*;
use nalgebra::Complex;
use nalgebra::ComplexField;
use nalgebra::Matrix2;
use nalgebra::Vector3;

use crate::beam;
use crate::clip::Clipping;
use crate::config;
use crate::field;
use crate::field::Field;
use crate::fresnel;
use crate::geom::Face;
use crate::geom::Geom;
use crate::helpers::draw_face;
use crate::helpers::lines_to_screen;
use crate::snell::get_theta_t;

#[derive(Debug, Clone, PartialEq)]
pub struct BeamPropagation {
    pub input: Beam,
    pub refr_index: Complex<f32>,
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
        let input_mid = self.input.data().face.data().midpoint;
        // draw lines from the outputs to the input
        let mut line_strings: Vec<_> = self
            .outputs
            .iter()
            .map(|x| {
                let output_mid = x.data().face.midpoint();
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

        // draw a small line in the direction of propagation
        let length = 3.0;
        let propagation_line = vec![vec![
            Coord {
                x: input_mid.coords.x,
                y: input_mid.coords.y,
            },
            Coord {
                x: input_mid.coords.x + self.input.data().prop.x,
                y: input_mid.coords.y + self.input.data().prop.y,
            },
        ]];

        lines_to_screen(line_strings, RED, 2.0);
        lines_to_screen(propagation_line, MAGENTA, 5.0);
    }

    pub fn input_power(&self) -> f32 {
        self.input.data().power()
    }

    pub fn output_power(&self) -> f32 {
        let total = self
            .outputs
            .iter()
            .fold(0.0, |acc, x| acc + x.data().power());

        total
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Beam {
    Initial(BeamData), // an initial beam to be traced in the near-field
    Default {
        data: BeamData,       // a beam to be traced in the near-field
        variant: BeamVariant, // whether the beam was a refl, refraction, or tir,
    },
    OutGoing(BeamData), // a beam to be mapped to the far-field
}

impl Beam {
    /// Creates a new initial field. The amplitude matrix is the identity matrix
    /// with the specified perpendicular field vector.
    pub fn new_initial(
        face: Face,
        prop: Vector3<f32>,
        refr_index: Complex<f32>,
        e_perp: Vector3<f32>,
    ) -> Self {
        Beam::Initial(BeamData::new(
            face,
            prop,
            refr_index,
            0,
            0,
            Field::new_identity(e_perp, prop),
        ))
    }
    pub fn new_default(
        face: Face,
        proj: Vector3<f32>,
        refr_index: Complex<f32>,
        rec_count: i32,
        tir_count: i32,
        variant: BeamVariant,
        field: Field,
    ) -> Self {
        Beam::Default {
            data: BeamData::new(face, proj, refr_index, rec_count, tir_count, field),
            variant,
        }
    }
    pub fn new_outgoing(beam_data: &BeamData) -> Beam {
        Beam::OutGoing(beam_data.clone())
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
            Beam::Default { data, variant } => {
                if data.rec_count < config::MAX_REC
                    || (*variant == BeamVariant::Tir && data.tir_count < config::MAX_TIR)
                {
                    let output_beams = Self::process_beam(geom, data);
                    println!("adding {} beams to the outputs", output_beams.len());
                    outputs.extend(output_beams);
                } else {
                    println!(
                        "beam, variant: {:?}, trunacted with rec: {}, tir: {}",
                        variant, data.rec_count, data.tir_count
                    );
                }
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

        let mut clipping = Clipping::new(geom, &mut beam_data.face, &beam_data.prop);
        clipping.clip(); // do the clipping -> contains the intersections

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
                        rec_count: beam_data.rec_count,
                        tir_count: beam_data.tir_count,
                        field: beam_data.field.clone(),
                        power_in: 0.0,
                        power_out: 0.0,
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

        // create  beams
        let beams: Vec<_> = intersections
            .iter()
            .filter_map(|x| {
                let normal = x.data().normal;
                let theta_i = normal.dot(&beam_data.prop).abs().acos();
                let n2 = Self::get_n2(geom, x.data().shape_id.unwrap(), input_shape_id);
                let e_perp = if normal.dot(&beam_data.prop) > 0.001 {
                    normal.cross(&beam_data.prop).normalize() // new e_perp
                } else {
                    -beam_data.field.e_perp
                };
                let rot = Field::rotation_matrix(beam_data.field.e_perp, e_perp, beam_data.prop)
                    .map(|x| nalgebra::Complex::new(x, 0.0)); // rotation matrix
                let mut ampl = rot * beam_data.field.ampl.clone();
                // let mut ampl = beam_data.field.ampl.clone();
                let dist = (x.midpoint() - beam_data.face.data().midpoint).dot(&beam_data.prop); // z-distance
                let arg = dist * config::WAVENO * n1.re; // optical path length
                ampl *= Complex::new(arg.cos(), arg.sin()); //  apply distance phase factor
                let arg = -2.0 * config::WAVENO * n1.im * dist.sqrt(); // absorption
                ampl *= Complex::new(arg.cos(), arg.sin()); //  apply absorption factor

                let refracted = if theta_i > (n2.re / n1.re).asin() {
                    // if total internal reflection
                    None
                } else {
                    let theta_t = get_theta_t(theta_i, n1, n2); // sin(theta_t)
                    let prop =
                        Self::get_refraction_vector(&normal, &beam_data.prop, theta_i, theta_t);
                    let fresnel = fresnel::refr(n1, n2, theta_i, theta_t);
                    let refr_ampl = fresnel * ampl.clone();

                    debug_assert!(beam_data.prop.dot(&prop) > 0.0);
                    debug_assert!((prop.dot(&normal).abs() - theta_t.cos()).abs() < 0.01);

                    Some(Beam::new_default(
                        x.clone(),
                        prop,
                        n2,
                        beam_data.rec_count + 1,
                        beam_data.tir_count,
                        BeamVariant::Refr,
                        Field::new(e_perp, prop, refr_ampl),
                    ))
                };

                let reflected = {
                    let prop = Self::get_reflection_vector(&normal, &beam_data.prop);

                    debug_assert!((prop.dot(&normal) - theta_i.cos()) < 0.01);

                    if theta_i > (n2.re / n1.re).asin() {
                        // if total internal reflection
                        let fresnel = -Matrix2::identity().map(|x| nalgebra::Complex::new(x, 0.0));
                        let refl_ampl = fresnel * ampl;

                        Some(Beam::new_default(
                            x.clone(),
                            prop,
                            n1,
                            beam_data.rec_count + 1,
                            beam_data.tir_count + 1,
                            BeamVariant::Tir,
                            Field::new(e_perp, prop, refl_ampl),
                        ))
                    } else {
                        let theta_t = get_theta_t(theta_i, n1, n2); // sin(theta_t)
                        let fresnel = fresnel::refl(n1, n2, theta_i, theta_t);
                        let refl_ampl = fresnel * ampl;

                        Some(Beam::new_default(
                            x.clone(),
                            prop,
                            n1,
                            beam_data.rec_count + 1,
                            beam_data.tir_count,
                            BeamVariant::Refl,
                            Field::new(e_perp, prop, refl_ampl),
                        ))
                    }
                };

                Some((reflected, refracted))

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
    fn get_n2(geom: &Geom, output_shape_id: usize, input_shape_id: Option<usize>) -> Complex<f32> {
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
    fn get_refraction_vector(
        norm: &Vector3<f32>,
        prop: &Vector3<f32>,
        theta_i: f32,
        theta_t: f32,
    ) -> Vector3<f32> {
        if theta_t.sin() < 0.0001 {
            return *prop;
        }
        // upward facing normal
        let n = if norm.dot(&prop) > 0.0 {
            *norm
        } else {
            *norm * -1.0
        };

        let alpha = PI - theta_t;
        let a = (theta_t - theta_i).sin() / theta_i.sin();
        let b = alpha.sin() / theta_i.sin();

        let mut result = b * prop - a * n;

        result.normalize_mut();

        debug_assert!((theta_t.cos() - result.dot(&norm).abs()).abs() < 0.01);

        result
    }

    fn get_reflection_vector(norm: &Vector3<f32>, prop: &Vector3<f32>) -> Vector3<f32> {
        // upward facing normal
        let n = if norm.dot(&prop) > 0.0 {
            *norm
        } else {
            *norm * -1.0
        };
        let cti = n.dot(&prop); // cos theta_i
        let mut result = prop - 2.0 * cti * n;
        result.normalize_mut();
        assert!((result.dot(&n) - cti) < 0.01);
        result
    }
}

/// Contains information about a beam.
#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct BeamData {
    pub face: Face,
    pub prop: Vector3<f32>,
    pub refr_index: Complex<f32>,
    pub rec_count: i32,
    pub tir_count: i32,
    pub field: Field,
    pub power_in: f32,
    pub power_out: f32,
}

/// Creates a new beam
impl BeamData {
    pub fn new(
        face: Face,
        prop: Vector3<f32>,
        refr_index: Complex<f32>,
        rec_count: i32,
        tir_count: i32,
        field: Field,
    ) -> Self {
        let prop = prop.normalize();
        Self {
            face,
            prop,
            refr_index,
            rec_count,
            tir_count,
            field,
            power_in: 0.0,
            power_out: 0.0,
        }
    }

    /// Returns the cross sectional area of the beam.
    pub fn csa(&self) -> f32 {
        let area = self.face.data().area.unwrap();
        let norm = self.face.data().normal;
        let cosine = self.prop.dot(&norm).abs();

        area * cosine
    }

    /// Returns the power of a beam.
    pub fn power(&self) -> f32 {
        self.field.intensity() * self.refr_index.re * self.csa()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BeamVariant {
    Refl, // refraction
    Refr, // reflection
    Tir,  // total internal reflection
}
