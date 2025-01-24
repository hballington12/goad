use anyhow::Result;
use std::f32::consts::PI;
use std::fmt::Error;

use geo::Coord;
use macroquad::prelude::*;
use nalgebra::Complex;
use nalgebra::Matrix2;
use nalgebra::Point3;
use nalgebra::Vector3;

use crate::clip::Clipping;
use crate::config;
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

        // // draw lines from the outputs midpoints to the input
        // let line_strings: Vec<_> = self
        //     .outputs
        //     .iter()
        //     .map(|x| Self::get_line(&x.data().face.midpoint(), &self.input))
        //     .collect();

        // draw lines from all vertices of outputs to the input
        let mut line_strings = Vec::new();
        for output in &self.outputs {
            for vertex in &output.data().face.data().exterior {
                line_strings.push(Self::get_line(&vertex, &self.input));
            }
        }

        // draw a small line in the direction of propagation
        let length = 1.0;
        let propagation_line = vec![vec![
            Coord {
                x: input_mid.coords.x,
                y: input_mid.coords.y,
            },
            Coord {
                x: input_mid.coords.x + self.input.data().prop.x * length,
                y: input_mid.coords.y + self.input.data().prop.y * length,
            },
        ]];

        // draw a small line in the direction of normal
        let length = 1.5;
        let normal_line = vec![vec![
            Coord {
                x: input_mid.coords.x,
                y: input_mid.coords.y,
            },
            Coord {
                x: input_mid.coords.x + self.input.data().face.data().normal.x * length,
                y: input_mid.coords.y + self.input.data().face.data().normal.y * length,
            },
        ]];

        lines_to_screen(line_strings, RED, 2.0);
        lines_to_screen(propagation_line, MAGENTA, 5.0);
        lines_to_screen(normal_line, WHITE, 2.5);
    }

    fn get_line(point: &Point3<f32>, input: &Beam) -> Vec<Coord<f32>> {
        let output_mid = point;
        let input_mid = input.data().face.data().midpoint;
        let vec = input_mid - output_mid;
        let input_normal = input.data().face.data().normal;
        let norm_dist_to_plane = vec.dot(&input_normal);
        let dist_to_plane = norm_dist_to_plane / (input_normal.dot(&input.data().prop));
        // ray cast along propagation direction
        let intsn = output_mid + dist_to_plane * input.data().prop;
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
    ) -> Result<Self> {
        let field = Field::new_identity(e_perp, prop)?;
        Ok(Beam::Initial(BeamData::new(
            face, prop, refr_index, 0, 0, field,
        )))
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
                if data.power() > config::BEAM_POWER_THRESHOLD
                    && (data.rec_count < config::MAX_REC
                        || (*variant == BeamVariant::Tir && data.tir_count < config::MAX_TIR))
                {
                    let output_beams = Self::process_beam(geom, data);
                    println!("adding {} beams to the outputs", output_beams.len());
                    outputs.extend(output_beams);
                } else {
                    println!(
                        "beam, variant: {:?}, trunacted with rec: {}, tir: {}, power: {}",
                        variant,
                        data.rec_count,
                        data.tir_count,
                        data.power()
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

    /// Processes a beam. The beam is propagated, the remainders, reflected,
    /// and refracted beams are computed and output.
    fn process_beam(geom: &mut Geom, beam_data: &mut BeamData) -> Vec<Beam> {
        let mut clipping = Clipping::new(geom, &mut beam_data.face, &beam_data.prop);
        let _ = clipping.clip();

        let (intersections, remainders) = (
            filter_faces(clipping.intersections),
            filter_faces(clipping.remaining),
        );

        let remainder_beams = remainders_to_beams(beam_data, remainders);
        let beams = create_beams(geom, beam_data, intersections);

        let mut output_beams = Vec::new();
        output_beams.extend(beams);
        output_beams.extend(remainder_beams);
        output_beams
    }
}

/// Returns a transmitted propagation vector, where `stt` is the sine of the angle of transmission.
pub fn get_refraction_vector(
    norm: &Vector3<f32>,
    prop: &Vector3<f32>,
    theta_i: f32,
    theta_t: f32,
) -> Vector3<f32> {
    if theta_t.sin() < config::COLINEAR_THRESHOLD {
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

    debug_assert!((theta_t.cos() - result.dot(&norm).abs()).abs() < config::COLINEAR_THRESHOLD);

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
    assert!((result.dot(&n) - cti) < config::COLINEAR_THRESHOLD);
    result
}

fn create_beams(geom: &mut Geom, beam_data: &mut BeamData, intersections: Vec<Face>) -> Vec<Beam> {
    let n1 = beam_data.refr_index;

    // create  beams
    intersections
        .iter()
        .filter_map(|face| {
            let normal = face.data().normal;
            let theta_i = normal.dot(&beam_data.prop).abs().acos();
            let n2 = get_n2(geom, beam_data, face, normal);
            let e_perp = get_e_perp(normal, &beam_data);
            let rot = get_rotation_matrix(&beam_data, e_perp);
            let ampl = match get_ampl(&beam_data, rot, face, n1) {
                Ok(ampl) => ampl,
                Err(_) => return None, // skip this intersection if get_ampl() returns error
            };

            let refracted =
                match create_refracted(face, ampl, e_perp, normal, beam_data, theta_i, n1, n2) {
                    Ok(beam) => beam,
                    Err(_) => {
                        // count skipped beams
                        None
                    }
                };

            let reflected =
                match create_reflected(face, ampl, e_perp, normal, beam_data, theta_i, n1, n2) {
                    Ok(beam) => beam,
                    Err(_) => {
                        // count skipped beams
                        None
                    }
                };

            Some((reflected, refracted))

            // determine other BeamData values here later...
        })
        .into_iter()
        .flat_map(|(refl, trans)| refl.into_iter().chain(trans))
        .collect()
}

/// Takes an amplitude matrix from the input beam data, rotates it into the new
/// scattering plane using the rotation matrix `rot`, computes the distance to
/// the intersection `face`, and applies the corresponding phase and absorption
/// factors.
fn get_ampl(
    beam_data: &BeamData,
    rot: Matrix2<Complex<f32>>,
    face: &Face,
    n1: Complex<f32>,
) -> Result<Matrix2<Complex<f32>>> {
    let mut ampl = rot * beam_data.field.ampl.clone();
    let dist = (face.midpoint() - beam_data.face.data().midpoint).dot(&beam_data.prop); // z-distance

    if dist < 0.0 {
        return Err(anyhow::anyhow!("distance less than 0: {}", dist));
    }

    let arg = dist * config::WAVENO * n1.re; // optical path length
    ampl *= Complex::new(arg.cos(), arg.sin()); //  apply distance phase factor
    let exp_absorption = (-2.0 * config::WAVENO * n1.im * dist.sqrt()).exp(); // absorption
    ampl *= Complex::new(exp_absorption, 0.0); //  apply absorption factor
    Ok(ampl)
}

/// Returns a rotation matrix for rotating from the plane perpendicular to e_perp
/// in `beam_data` to the plane perpendicular to `e_perp`.
fn get_rotation_matrix(beam_data: &BeamData, e_perp: Vector3<f32>) -> Matrix2<Complex<f32>> {
    Field::rotation_matrix(beam_data.field.e_perp, e_perp, beam_data.prop)
        .map(|x| nalgebra::Complex::new(x, 0.0))
}

/// Determines the new `e_perp` vector for an intersection at a `face``.
fn get_e_perp(normal: Vector3<f32>, beam_data: &BeamData) -> Vector3<f32> {
    if normal.dot(&beam_data.prop).abs() > 1.0 - config::COLINEAR_THRESHOLD {
        -beam_data.field.e_perp
    } else {
        normal.cross(&beam_data.prop).normalize() // new e_perp
    }
}

/// Determines the refractive index of the second medium when a beam intersects
/// with a face.
fn get_n2(
    geom: &mut Geom,
    beam_data: &mut BeamData,
    face: &Face,
    normal: Vector3<f32>,
) -> Complex<f32> {
    let id = face.data().shape_id.unwrap();
    if normal.dot(&beam_data.prop) < 0.0 {
        geom.shapes[id].refr_index
    } else {
        geom.n_out(id)
    }
}

/// Creates a new reflected beam
fn create_reflected(
    face: &Face,
    ampl: Matrix2<Complex<f32>>,
    e_perp: Vector3<f32>,
    normal: Vector3<f32>,
    beam_data: &BeamData,
    theta_i: f32,
    n1: Complex<f32>,
    n2: Complex<f32>,
) -> Result<Option<Beam>> {
    let prop = get_reflection_vector(&normal, &beam_data.prop);

    debug_assert!((prop.dot(&normal) - theta_i.cos()) < config::COLINEAR_THRESHOLD);
    debug_assert!(!Field::ampl_intensity(&ampl).is_nan());

    if theta_i > (n2.re / n1.re).asin() {
        // if total internal reflection
        let fresnel = -Matrix2::identity().map(|x| nalgebra::Complex::new(x, 0.0));
        let refl_ampl = fresnel * ampl;
        debug_assert!(!Field::ampl_intensity(&refl_ampl).is_nan());

        Ok(Some(Beam::new_default(
            face.clone(),
            prop,
            n1,
            beam_data.rec_count + 1,
            beam_data.tir_count + 1,
            BeamVariant::Tir,
            Field::new(e_perp, prop, refl_ampl).unwrap(),
        )))
    } else {
        let theta_t = get_theta_t(theta_i, n1, n2); // sin(theta_t)
        let fresnel = fresnel::refl(n1, n2, theta_i, theta_t);
        let refl_ampl = fresnel * ampl;

        Ok(Some(Beam::new_default(
            face.clone(),
            prop,
            n1,
            beam_data.rec_count + 1,
            beam_data.tir_count,
            BeamVariant::Refl,
            Field::new(e_perp, prop, refl_ampl)?,
        )))
    }
}

/// Creates a new refracted beam.
fn create_refracted(
    face: &Face,
    ampl: Matrix2<Complex<f32>>,
    e_perp: Vector3<f32>,
    normal: Vector3<f32>,
    beam_data: &BeamData,
    theta_i: f32,
    n1: Complex<f32>,
    n2: Complex<f32>,
) -> Result<Option<Beam>> {
    if theta_i > (n2.re / n1.re).asin() {
        // if total internal reflection
        Ok(None)
    } else {
        let theta_t = get_theta_t(theta_i, n1, n2); // sin(theta_t)
        let prop = get_refraction_vector(&normal, &beam_data.prop, theta_i, theta_t);
        let fresnel = fresnel::refr(n1, n2, theta_i, theta_t);
        let refr_ampl = fresnel * ampl.clone();

        debug_assert!(beam_data.prop.dot(&prop) > 0.0);
        debug_assert!((prop.dot(&normal).abs() - theta_t.cos()).abs() < config::COLINEAR_THRESHOLD);

        Ok(Some(Beam::new_default(
            face.clone(),
            prop,
            n2,
            beam_data.rec_count + 1,
            beam_data.tir_count,
            BeamVariant::Refr,
            Field::new(e_perp, prop, refr_ampl)?,
        )))
    }
}

/// Filter faces by threshold area.
fn filter_faces(faces: Vec<Face>) -> Vec<Face> {
    // remove intersections with below threshold area
    faces
        .into_iter()
        .filter(|x| x.data().area.unwrap() > config::BEAM_AREA_THRESHOLD)
        .collect()
}

/// Converts the remainder faces from a clipping into beams with the same field
/// properties as the original beam.
fn remainders_to_beams(beam_data: &mut BeamData, remainders: Vec<Face>) -> Vec<Beam> {
    let remainder_beams: Vec<_> = remainders
        .into_iter()
        .filter_map(|remainder| {
            Some(Beam::OutGoing(BeamData {
                face: remainder,
                prop: beam_data.prop,
                refr_index: beam_data.refr_index,
                rec_count: beam_data.rec_count,
                tir_count: beam_data.tir_count,
                field: beam_data.field.clone(),
            }))
        })
        .collect();
    remainder_beams
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
