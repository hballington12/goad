use anyhow::Result;
use std::f32::consts::PI;

use geo::Coord;
#[cfg(feature = "macroquad")]
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};

use crate::{
    clip::Clipping,
    diff,
    field::Field,
    fresnel,
    geom::{Face, Geom},
    helpers,
    settings,
    snell::get_theta_t,
};

#[derive(Debug, Clone, PartialEq)]
pub struct BeamPropagation {
    pub input: Beam,
    pub refr_index: Complex<f32>,
    pub outputs: Vec<Beam>,
}

impl BeamPropagation {
    /// Makes a new `BeamPropagation` struct, which represents a beam propagation.
    pub fn new(input: Beam, outputs: Vec<Beam>) -> Self {
        let refr_index = input.refr_index.clone();
        Self {
            input,
            refr_index,
            outputs,
        }
    }

    /// Draws a `Beam Propagation`
    #[cfg(feature = "macroquad")]
    pub fn draw(&self) {
        // draw the input
        #[cfg(feature = "macroquad")]
        helpers::draw_face(&self.input.face, YELLOW, 4.0);
        // draw the outputs
        for beam in &self.outputs {
            if beam.type_ == BeamType::Default {
                #[cfg(feature = "macroquad")]
                helpers::draw_face(&beam.face, BLUE, 4.0);
            } else if beam.type_ == BeamType::OutGoing {
                #[cfg(feature = "macroquad")]
                helpers::draw_face(&beam.face, RED, 3.0);
            }
        }
        let input_mid = self.input.face.data().midpoint;

        // // draw lines from the outputs midpoints to the input
        // let line_strings: Vec<_> = self
        //     .outputs
        //     .iter()
        //     .map(|x| Self::get_line(&x.data().face.midpoint(), &self.input))
        //     .collect();

        // draw lines from all vertices of outputs to the input
        let mut line_strings = Vec::new();
        for output in &self.outputs {
            for vertex in &output.face.data().exterior {
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
                x: input_mid.coords.x + self.input.prop.x * length,
                y: input_mid.coords.y + self.input.prop.y * length,
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
                x: input_mid.coords.x + self.input.face.data().normal.x * length,
                y: input_mid.coords.y + self.input.face.data().normal.y * length,
            },
        ]];

        #[cfg(feature = "macroquad")]
        helpers::lines_to_screen(line_strings, RED, 2.0);
        #[cfg(feature = "macroquad")]
        helpers::lines_to_screen(propagation_line, MAGENTA, 5.0);
        #[cfg(feature = "macroquad")]
        helpers::lines_to_screen(normal_line, WHITE, 2.5);
    }

    fn get_line(point: &Point3<f32>, input: &Beam) -> Vec<Coord<f32>> {
        let output_mid = point;
        let input_mid = input.face.data().midpoint;
        let vec = input_mid - output_mid;
        let input_normal = input.face.data().normal;
        let norm_dist_to_plane = vec.dot(&input_normal);
        let dist_to_plane = norm_dist_to_plane / (input_normal.dot(&input.prop));
        // ray cast along propagation direction
        let intsn = output_mid + dist_to_plane * input.prop;
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
        self.input.power()
    }

    pub fn output_power(&self) -> f32 {
        let total = self.outputs.iter().fold(0.0, |acc, x| acc + x.power());

        total
    }
}

impl Beam {
    /// Creates a new initial field. The amplitude matrix is the identity matrix
    /// with the specified perpendicular field vector.
    pub fn new_initial(
        face: Face,
        prop: Vector3<f32>,
        refr_index: Complex<f32>,
        e_perp: Vector3<f32>,
        wavelength: f32,
    ) -> Result<Self> {
        let field = Field::new_identity(e_perp, prop)?;
        Ok(Beam::new(
            face,
            prop,
            refr_index,
            0,
            0,
            field,
            None,
            BeamType::Initial,
            wavelength,
        ))
    }

    pub fn new_from_field(
        face: Face,
        prop: Vector3<f32>,
        refr_index: Complex<f32>,
        field: Field,
        wavelength: f32,
    ) -> Self {
        Beam::new(
            face,
            prop,
            refr_index,
            0,
            0,
            field,
            None,
            BeamType::Initial,
            wavelength,
        )
    }

    /// Processes data from a beam. The beam is propagated, the remainders, reflected,
    /// and refracted beams are computed and output.
    pub fn propagate(
        &mut self,
        geom: &mut Geom,
        medium_refr_index: Complex<f32>,
        area_threshold: f32,
    ) -> Result<(Vec<Beam>, f32)> {
        let mut clipping = Clipping::new(geom, &mut self.face, &self.prop);
        clipping.clip(area_threshold)?;

        self.clipping_area = match clipping.stats {
            Some(stats) => stats.intersection_area + stats.remaining_area,
            _ => 0.0,
        };

        let (intersections, remainders) = (
            clipping.intersections.into_iter().collect(),
            clipping.remaining.into_iter().collect(),
        );

        let remainder_beams = self.remainders_to_beams(remainders, medium_refr_index);
        let beams = self.create_beams(geom, intersections, medium_refr_index);

        let mut output_beams = Vec::new();
        output_beams.extend(beams);
        output_beams.extend(remainder_beams);
        let output_power = output_beams.iter().fold(0.0, |acc, x| acc + x.power());
        let power_loss = self.power() - self.absorbed_power - output_power;

        Ok((output_beams, power_loss))
    }

    fn create_beams(
        &mut self,
        geom: &mut Geom,
        intersections: Vec<Face>,
        medium_refr_index: Complex<f32>,
    ) -> Vec<Beam> {
        let n1 = self.refr_index;

        let mut outputs = Vec::new();
        for face in &intersections {
            let normal = face.data().normal;
            let theta_i = normal.dot(&self.prop).abs().acos();
            let n2 = get_n2(geom, self, face, normal, medium_refr_index);
            let e_perp = get_e_perp(normal, &self);
            let rot = get_rotation_matrix(&self, e_perp);
            let (ampl, absorbed_intensity) = get_ampl(&self, rot, face, n1);

            self.absorbed_power +=
                absorbed_intensity * face.data().area.unwrap() * theta_i.cos() * n1.re;

            if self.type_ == BeamType::Initial {
                let external_diff = Beam::new(
                    face.clone(),
                    self.prop,
                    n1,
                    self.rec_count + 1,
                    self.tir_count,
                    Field::new(e_perp, self.prop, ampl).unwrap(),
                    None,
                    BeamType::ExternalDiff,
                    self.wavelength,
                );
                outputs.push(external_diff);
            }

            // untracked energy leaks can occur here if the amplitude matrix contains NaN values
            let refracted =
                create_refracted(face, ampl, e_perp, normal, self, theta_i, n1, n2).unwrap_or(None);
            let reflected =
                create_reflected(face, ampl, e_perp, normal, self, theta_i, n1, n2).unwrap_or(None);

            if refracted.is_some() {
                outputs.push(refracted.unwrap().clone());
            }
            if reflected.is_some() {
                outputs.push(reflected.unwrap().clone());
            }
        }

        outputs
    }

    /// Uses the earcut function from the geom crate to convert a beam with
    /// a complex face into beams with simple faces. The medium refractive index
    /// is required to map the phase.
    fn earcut(beam: &Beam, medium_refr_index: Complex<f32>) -> Vec<Beam> {
        let mut outputs = Vec::new();
        let midpoint = beam.face.data().midpoint;
        match &beam.face {
            Face::Simple(_) => outputs.push(beam.clone()),
            Face::Complex { .. } => {
                let faces = Face::earcut(&beam.face);
                for face in faces {
                    let dist = (face.data().midpoint - midpoint).dot(&beam.prop);
                    let arg = dist * beam.wavenumber() * medium_refr_index.re;
                    let ampl = beam.field.ampl.clone() * Complex::new(arg.cos(), arg.sin());

                    let new_beam = Beam::new(
                        face,
                        beam.prop,
                        beam.refr_index,
                        beam.rec_count,
                        beam.tir_count,
                        Field::new(beam.field.e_perp, beam.prop, ampl).unwrap(),
                        beam.variant.clone(),
                        beam.type_.clone(),
                        beam.wavelength,
                    );

                    outputs.push(new_beam);
                }
            }
        }
        outputs
    }
}

/// Returns a transmitted propagation vector, where `stt` is the sine of the angle of transmission.
pub fn get_refraction_vector(
    norm: &Vector3<f32>,
    prop: &Vector3<f32>,
    theta_i: f32,
    theta_t: f32,
) -> Vector3<f32> {
    if theta_t.sin() < settings::COLINEAR_THRESHOLD {
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

    debug_assert!((theta_t.cos() - result.dot(&norm).abs()).abs() < settings::COLINEAR_THRESHOLD);

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
    assert!((result.dot(&n) - cti) < settings::COLINEAR_THRESHOLD);
    result
}

/// Takes an amplitude matrix from the input beam data, rotates it into the new
/// scattering plane using the rotation matrix `rot`, computes the distance to
/// the intersection `face`, and applies the corresponding phase and absorption
/// factors.
fn get_ampl(
    beam: &Beam,
    rot: Matrix2<Complex<f32>>,
    face: &Face,
    n1: Complex<f32>,
) -> (Matrix2<Complex<f32>>, f32) {
    let mut ampl = rot * beam.field.ampl.clone();

    let dist = (face.midpoint() - beam.face.data().midpoint).dot(&beam.prop); // z-distance
    let wavenumber = beam.wavenumber();

    let arg = dist * wavenumber * n1.re; // optical path length
    ampl *= Complex::new(arg.cos(), arg.sin()); //  apply distance phase factor

    let dist_sqrt = dist.signum() * dist.abs().sqrt(); // TODO: improve this

    let absorbed_intensity = Field::ampl_intensity(&ampl)
        * (1.0 - (-2.0 * wavenumber * n1.im * dist_sqrt).exp().powi(2));

    let exp_absorption = (-2.0 * wavenumber * n1.im * dist_sqrt).exp(); // absorption

    ampl *= Complex::new(exp_absorption, 0.0); //  apply absorption factor

    (ampl, absorbed_intensity)
}

/// Returns a rotation matrix for rotating from the plane perpendicular to e_perp
/// in `beam` to the plane perpendicular to `e_perp`.
fn get_rotation_matrix(beam: &Beam, e_perp: Vector3<f32>) -> Matrix2<Complex<f32>> {
    Field::rotation_matrix(beam.field.e_perp, e_perp, beam.prop)
        .map(|x| nalgebra::Complex::new(x, 0.0))
}

/// Determines the new `e_perp` vector for an intersection at a `face``.
fn get_e_perp(normal: Vector3<f32>, beam: &Beam) -> Vector3<f32> {
    if normal.dot(&beam.prop).abs() > 1.0 - settings::COLINEAR_THRESHOLD {
        -beam.field.e_perp
    } else {
        normal.cross(&beam.prop).normalize() // new e_perp
    }
}

/// Determines the refractive index of the second medium when a beam intersects
/// with a face.
fn get_n2(
    geom: &mut Geom,
    beam: &mut Beam,
    face: &Face,
    normal: Vector3<f32>,
    medium_refr_index: Complex<f32>,
) -> Complex<f32> {
    let id = face.data().shape_id.unwrap();
    if normal.dot(&beam.prop) < 0.0 {
        geom.shapes[id].refr_index
    } else {
        geom.n_out(id, medium_refr_index)
    }
}

/// Creates a new reflected beam
fn create_reflected(
    face: &Face,
    ampl: Matrix2<Complex<f32>>,
    e_perp: Vector3<f32>,
    normal: Vector3<f32>,
    beam: &Beam,
    theta_i: f32,
    n1: Complex<f32>,
    n2: Complex<f32>,
) -> Result<Option<Beam>> {
    let prop = get_reflection_vector(&normal, &beam.prop);

    debug_assert!((prop.dot(&normal) - theta_i.cos()) < settings::COLINEAR_THRESHOLD);
    debug_assert!(!Field::ampl_intensity(&ampl).is_nan());

    if theta_i > (n2.re / n1.re).asin() {
        // if total internal reflection
        let fresnel = -Matrix2::identity().map(|x| nalgebra::Complex::new(x, 0.0));
        let refl_ampl = fresnel * ampl;
        debug_assert!(!Field::ampl_intensity(&refl_ampl).is_nan());

        Ok(Some(Beam::new(
            face.clone(),
            prop,
            n1,
            beam.rec_count + 1,
            beam.tir_count + 1,
            Field::new(e_perp, prop, refl_ampl)?,
            Some(BeamVariant::Tir),
            BeamType::Default,
            beam.wavelength,
        )))
    } else {
        let theta_t = get_theta_t(theta_i, n1, n2)?; // sin(theta_t)
        let fresnel = fresnel::refl(n1, n2, theta_i, theta_t);
        let refl_ampl = fresnel * ampl;

        Ok(Some(Beam::new(
            face.clone(),
            prop,
            n1,
            beam.rec_count + 1,
            beam.tir_count,
            Field::new(e_perp, prop, refl_ampl)?,
            Some(BeamVariant::Refl),
            BeamType::Default,
            beam.wavelength,
        )))
    }
}

/// Creates a new refracted beam.
fn create_refracted(
    face: &Face,
    ampl: Matrix2<Complex<f32>>,
    e_perp: Vector3<f32>,
    normal: Vector3<f32>,
    beam: &Beam,
    theta_i: f32,
    n1: Complex<f32>,
    n2: Complex<f32>,
) -> Result<Option<Beam>> {
    if theta_i >= (n2.re / n1.re).asin() {
        // if total internal reflection
        Ok(None)
    } else {
        let theta_t = get_theta_t(theta_i, n1, n2)?; // sin(theta_t)
        let prop = get_refraction_vector(&normal, &beam.prop, theta_i, theta_t);
        let fresnel = fresnel::refr(n1, n2, theta_i, theta_t);
        let refr_ampl = fresnel * ampl.clone();

        debug_assert!(beam.prop.dot(&prop) > 0.0);
        debug_assert!(
            (prop.dot(&normal).abs() - theta_t.cos()).abs() < settings::COLINEAR_THRESHOLD
        );

        Ok(Some(Beam::new(
            face.clone(),
            prop,
            n2,
            beam.rec_count + 1,
            beam.tir_count,
            Field::new(e_perp, prop, refr_ampl)?,
            Some(BeamVariant::Refr),
            BeamType::Default,
            beam.wavelength,
        )))
    }
}

/// Converts the remainder faces from a clipping into beams with the same field
/// properties as the original beam.
impl Beam {
    fn remainders_to_beams(
        &mut self,
        remainders: Vec<Face>,
        medium_refr_index: Complex<f32>,
    ) -> Vec<Beam> {
        // need to account for distance along propagation direction from
        // midpoint of remainder to midpoint of original face. Propagate
        // the field back or forward by this distance.
        let self_midpoint = self.face.data().midpoint;
        let remainder_beams: Vec<_> = remainders
            .into_iter()
            .filter_map(|remainder| {
                let dist = (remainder.data().midpoint - self_midpoint).dot(&self.prop);
                let arg = dist * self.wavenumber() * medium_refr_index.re;
                // let arg: f32 = 0.0;
                let ampl = self.field.ampl.clone() * Complex::new(arg.cos(), arg.sin());
                Some(Beam::new(
                    remainder,
                    self.prop,
                    self.refr_index,
                    self.rec_count,
                    self.tir_count,
                    Field::new(self.field.e_perp, self.prop, ampl).unwrap(),
                    None,
                    BeamType::OutGoing,
                    self.wavelength,
                ))
            })
            .collect();

        // Also convert any complex faces into simple faces
        let mut output_beams = Vec::new();
        for beam in remainder_beams {
            output_beams.extend(Beam::earcut(&beam, medium_refr_index));
        }
        output_beams
    }
}

/// Contains information about a beam.
#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct Beam {
    pub face: Face,
    pub prop: Vector3<f32>,
    pub refr_index: Complex<f32>,
    pub rec_count: i32,
    pub tir_count: i32,
    pub field: Field,
    pub absorbed_power: f32,          // power absorbed by the medium
    pub clipping_area: f32,           // total area accounted for by intersections and remainders
    pub variant: Option<BeamVariant>, // variant of beam, e.g. reflection, refraction, total internal reflection
    pub type_: BeamType, // type of beam, e.g. initial, default, outgoing, external diff
    pub wavelength: f32,
}

/// Creates a new beam
impl Beam {
    pub fn new(
        face: Face,
        prop: Vector3<f32>,
        refr_index: Complex<f32>,
        rec_count: i32,
        tir_count: i32,
        field: Field,
        variant: Option<BeamVariant>,
        type_: BeamType,
        wavelength: f32,
    ) -> Self {
        let prop = prop.normalize();
        Self {
            face,
            prop,
            refr_index,
            rec_count,
            tir_count,
            field,
            absorbed_power: 0.0,
            clipping_area: 0.0,
            variant,
            type_,
            wavelength,
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

    pub fn wavenumber(&self) -> f32 {
        2.0 * PI / self.wavelength
    }

    pub fn diffract(
        &self,
        theta_phi_combinations: &[(f32, f32)],
        fov_factor: Option<f32>,
    ) -> Vec<Matrix2<Complex<f32>>> {
        match &self.face {
            Face::Simple(face) => {
                let verts = &face.exterior;
                let ampl = self.field.ampl;
                let prop = self.prop;
                let vk7 = self.field.e_perp;
                diff::diffraction(
                    verts,
                    ampl,
                    prop,
                    vk7,
                    &theta_phi_combinations,
                    self.wavenumber(),
                    fov_factor,
                )
            }
            Face::Complex { .. } => {
                println!("complex face not supported yet...");
                vec![Matrix2::zeros(); theta_phi_combinations.len()]
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BeamVariant {
    Refl, // refraction
    Refr, // reflection
    Tir,  // total internal reflection
}

#[derive(Debug, Clone, PartialEq)]
pub enum BeamType {
    Initial,
    Default,
    OutGoing,
    ExternalDiff,
}
