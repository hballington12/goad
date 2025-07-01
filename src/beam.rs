//! Electromagnetic beam representation and propagation for geometric optics.
//!
//! This module implements the core electromagnetic field segments (beams) that
//! propagate through particle geometries in the GOAD method. Each beam carries
//! field amplitude information, geometric cross-section data, and propagation
//! state for tracking reflection, refraction, and absorption processes.
//!
//! The beam system provides:
//! - Polarized electromagnetic field representation
//! - Surface interaction through Fresnel equations
//! - Geometric clipping and intersection handling
//! - Power conservation tracking
//! - Far-field diffraction calculation
//!
//! # Key Components
//!
//! - [`Beam`]: Individual electromagnetic field segments
//! - [`BeamPropagation`]: Results from beam-surface interactions
//! - [`BeamType`] and [`BeamVariant`]: Classification and origin tracking
//! - Surface interaction functions (reflection, refraction, TIR)
//! - Triangulation for complex apertures

use anyhow::Result;
use std::f32::consts::PI;

use geo::Coord;
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};

use crate::{
    clip::Clipping,
    diff,
    field::Field,
    fresnel,
    geom::{Face, Geom},
    helpers::{draw_face, lines_to_screen},
    settings,
    snell::get_theta_t,
};

/// Result of a single beam propagation step through geometry.
/// 
/// **Context**: Beam propagation involves complex interactions with particle surfaces
/// that can split a single input [`Beam`] into multiple output beams through reflection,
/// refraction, and diffraction. This structure captures the complete result of one
/// propagation step for analysis and visualization.
/// 
/// **How it Works**: Stores the input [`Beam`] that was propagated, the medium refractive
/// index used for calculations, and all output beams produced by surface interactions.
/// Used for debugging, visualization, and power conservation verification.
#[derive(Debug, Clone, PartialEq)]
pub struct BeamPropagation {
    pub input: Beam,              // The [`Beam`] that was propagated
    pub refr_index: Complex<f32>, // Medium refractive index used
    pub outputs: Vec<Beam>,       // Resulting [`Beam`] objects from interactions
}

impl BeamPropagation {
    /// Creates a propagation result from input beam and output beams.
    /// 
    /// **Context**: After beam propagation completes, the results need to be
    /// packaged for analysis, visualization, or debugging. This constructor
    /// captures the complete propagation step.
    /// 
    /// **How it Works**: Copies the refractive index from the input beam
    /// and stores both input and output beams for later analysis.
    /// 
    /// # Example
    /// ```rust
    /// Some(BeamPropagation::new(beam, outputs))
    /// ```
    pub fn new(input: Beam, outputs: Vec<Beam>) -> Self {
        let refr_index = input.refr_index.clone();
        Self {
            input,
            refr_index,
            outputs,
        }
    }

    /// Renders the propagation for visualization.
    /// 
    /// **Context**: Interactive debugging and educational visualization require
    /// showing both the input beam and resulting output beams from surface
    /// interactions, along with propagation directions.
    /// 
    /// **How it Works**: Draws the input beam face in yellow, output beam faces
    /// in different colors based on type, and lines connecting output vertices
    /// to the input beam to show the interaction relationships.
    pub fn draw(&self) {
        // draw the input
        draw_face(&self.input.face, YELLOW, 4.0);
        // draw the outputs
        for beam in &self.outputs {
            if beam.type_ == BeamType::Default {
                draw_face(&beam.face, BLUE, 4.0);
            } else if beam.type_ == BeamType::OutGoing {
                draw_face(&beam.face, RED, 3.0);
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

        lines_to_screen(line_strings, RED, 2.0);
        lines_to_screen(propagation_line, MAGENTA, 5.0);
        lines_to_screen(normal_line, WHITE, 2.5);
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

    /// Returns the power of the input beam.
    /// 
    /// **Context**: Power conservation analysis requires tracking power flow
    /// through each propagation step to verify energy conservation and
    /// identify power losses.
    /// 
    /// **How it Works**: Calls the power method on the input beam.
    pub fn input_power(&self) -> f32 {
        self.input.power()
    }

    /// Returns the total power of all output beams.
    /// 
    /// **Context**: Power conservation verification requires summing the power
    /// of all output beams to compare with the input power and identify
    /// absorption or numerical losses.
    /// 
    /// **How it Works**: Sums the power of all beams in the outputs vector.
    pub fn output_power(&self) -> f32 {
        let total = self.outputs.iter().fold(0.0, |acc, x| acc + x.power());

        total
    }
}

impl Beam {
    /// Creates an initial incident beam with identity amplitude matrix.
    /// 
    /// **Context**: Scattering simulations begin with an incident electromagnetic
    /// wave that must be represented as a beam with appropriate field properties.
    /// Initial beams use identity amplitude matrices to represent unscattered
    /// incident radiation.
    /// 
    /// **How it Works**: Constructs a [`crate::field::Field`] with identity amplitude matrix using
    /// the specified perpendicular polarization vector, then creates a beam
    /// with [`BeamType::Initial`] and zero interaction counters.
    /// 
    /// # Example
    /// ```rust
    /// let beam = Beam::new_initial(
    ///     clip,
    ///     projection,                    // propagation direction
    ///     Complex::new(1.00, 0.0),      // refractive index
    ///     e_perp,                       // perpendicular polarization vector
    ///     wavelength,                   // wavelength
    /// ).unwrap();
    /// ```
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

    /// Creates a beam with a pre-constructed electromagnetic field.
    /// 
    /// **Context**: Some simulation setups require specific [`crate::field::Field`] configurations
    /// that differ from the standard identity matrix initialization. This
    /// constructor allows direct specification of field properties.
    /// 
    /// **How it Works**: Uses the provided [`crate::field::Field`] directly without modification,
    /// creating an initial-type beam with the specified field properties.
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

    /// Propagates the beam through geometry to compute surface interactions.
    /// 
    /// **Context**: Beam propagation is the core of near-field electromagnetic simulation.
    /// Each beam must be tested for intersections with particle surfaces, and
    /// appropriate reflection and refraction beams must be generated according
    /// to electromagnetic boundary conditions.
    /// 
    /// **How it Works**: Uses geometric [`crate::clip::Clipping`] to find surface intersections,
    /// creates new beams for reflected and refracted components based on [`crate::fresnel`]
    /// equations, handles remainder beams that don't interact with surfaces,
    /// and tracks power absorption. Returns output beams and power loss information.
    /// 
    /// # Example
    /// ```rust
    /// match beam.propagate(
    ///     &mut self.geom,
    ///     self.settings.medium_refr_index,
    ///     self.settings.beam_area_threshold(),
    /// ) {
    ///     Ok((outputs, area_power_loss)) => {
    ///         self.result.powers.trnc_area += area_power_loss / self.settings.scale.powi(2);
    ///         outputs
    ///     }
    ///     Err(_) => {
    ///         self.result.powers.clip_err += beam.power() / self.settings.scale.powi(2);
    ///         Vec::new()
    ///     }
    /// }
    /// ```
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

    /// Converts a beam with complex face geometry into triangular beams.
    /// 
    /// **Context**: Some geometric faces contain holes or complex polygonal shapes
    /// that cannot be processed directly by diffraction algorithms. These faces
    /// must be triangulated into simpler elements while preserving electromagnetic
    /// field continuity.
    /// 
    /// **How it Works**: Uses ear clipping triangulation to split complex faces
    /// into triangular faces, creates new beams for each triangle with appropriate
    /// phase corrections based on the distance from the original face midpoint.
    /// Simple faces are passed through unchanged.
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

/// A segment of electromagnetic radiation with defined cross-section and field properties.
/// 
/// **Context**: Electromagnetic field propagation through particles requires discretizing
/// the field into manageable segments that can interact with surfaces independently.
/// Each beam carries electromagnetic [`crate::field::Field`] information, geometric cross-section data,
/// and propagation state for tracking reflection, refraction, and absorption processes.
/// 
/// **How it Works**: A [`Beam`] combines a geometric [`crate::geom::Face`] (defining the cross-sectional area),
/// a propagation direction, electromagnetic [`crate::field::Field`] properties, and interaction history.
/// The field contains amplitude and polarization information, while counters track
/// interaction depth for convergence control. The beam can be classified by [`BeamType`]
/// (initial, internal, outgoing) and [`BeamVariant`] (reflected, refracted, total internal reflection).
#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct Beam {
    pub face: Face,                   // [`crate::geom::Face`] defining the cross-sectional area
    pub prop: Vector3<f32>,           // propagation direction vector
    pub refr_index: Complex<f32>,     // refractive index of current medium
    pub rec_count: i32,               // recursion depth counter for convergence
    pub tir_count: i32,               // total internal reflection event counter
    pub field: Field,                 // electromagnetic [`crate::field::Field`] properties
    pub absorbed_power: f32,          // power absorbed by the medium
    pub clipping_area: f32,           // total area accounted for by intersections and remainders
    pub variant: Option<BeamVariant>, // variant of beam, e.g. reflection, refraction, total internal reflection
    pub type_: BeamType,              // [`BeamType`] classification (initial, default, outgoing, external diff)
    pub wavelength: f32,              // wavelength in current medium
}

impl Beam {
    /// Creates a new beam with specified properties.
    /// 
    /// **Context**: Beam creation occurs throughout the simulation as initial beams
    /// are established and surface interactions generate new reflected and refracted
    /// beams. Each beam requires electromagnetic [`crate::field::Field`] properties, geometric information,
    /// and interaction history for proper propagation.
    /// 
    /// **How it Works**: Normalizes the propagation direction and initializes all
    /// beam properties including [`crate::geom::Face`] geometry, [`crate::field::Field`] information, interaction counters,
    /// and [`BeamType`] classification. Sets absorbed power and clipping area to zero for new beams.
    /// 
    /// # Example
    /// ```rust
    /// let new_beam = Beam::new(
    ///     face,
    ///     beam.prop,
    ///     beam.refr_index,
    ///     beam.rec_count,
    ///     beam.tir_count,
    ///     Field::new(beam.field.e_perp, beam.prop, ampl).unwrap(),
    ///     beam.variant.clone(),
    ///     beam.type_.clone(),
    ///     beam.wavelength,
    /// );
    /// ```
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

    /// Returns the effective cross-sectional area for power calculations.
    /// 
    /// **Context**: Beam power calculations require the geometric cross-sectional
    /// area projected along the propagation direction. This accounts for oblique
    /// incidence where the effective area differs from the face area.
    /// 
    /// **How it Works**: Multiplies the face area by the cosine of the angle
    /// between the face normal and propagation direction to get the projected area.
    pub fn csa(&self) -> f32 {
        let area = self.face.data().area.unwrap();
        let norm = self.face.data().normal;
        let cosine = self.prop.dot(&norm).abs();

        area * cosine
    }

    /// Calculates the electromagnetic power carried by the beam.
    /// 
    /// **Context**: Power tracking is essential for energy conservation verification,
    /// threshold-based beam termination, and scattering cross-section calculations.
    /// Power determines the relative importance of beams in the simulation.
    /// 
    /// **How it Works**: Multiplies field intensity by the real part of the
    /// refractive index and the effective cross-sectional area to get total power.
    pub fn power(&self) -> f32 {
        self.field.intensity() * self.refr_index.re * self.csa()
    }

    /// Returns the wavenumber for the beam's wavelength.
    /// 
    /// **Context**: Phase calculations for beam propagation and interference
    /// require the wavenumber, which relates wavelength to spatial frequency.
    /// 
    /// **How it Works**: Calculates 2π/wavelength to get the wavenumber
    /// in the vacuum (refractive index effects are handled elsewhere).
    pub fn wavenumber(&self) -> f32 {
        2.0 * PI / self.wavelength
    }

    /// Computes far-field diffraction amplitudes at specified angles.
    /// 
    /// **Context**: Converting near-field beams to far-field scattering requires
    /// computing diffraction integrals over the beam cross-section. This transformation
    /// connects the geometric optics beam representation to scattering observables.
    /// 
    /// **How it Works**: Calls the diffraction module to compute Fourier transforms
    /// of the field distribution over the beam face vertices at each requested
    /// scattering angle. Returns amplitude matrices for polarization analysis.
    /// 
    /// # Example
    /// ```rust
    /// let ampl_far_field = queue
    ///     .par_iter()
    ///     .map(|outbeam| outbeam.diffract(bins, fov_factor))
    ///     .reduce(/* combine results */);
    /// ```
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

/// Classification of beam based on surface interaction mechanism.
/// 
/// **Context**: When beams interact with particle surfaces, they can undergo
/// different physical processes depending on the interface properties and
/// incident angles. Each mechanism requires different field calculations.
/// 
/// **How it Works**: Variants track the physical origin of the beam to enable
/// appropriate handling in power tracking and result analysis.
#[derive(Debug, Clone, PartialEq)]
pub enum BeamVariant {
    Refl, // reflection
    Refr, // refraction
    Tir,  // total internal reflection
}

/// Classification of beam based on its role in the simulation pipeline.
/// 
/// **Context**: Different stages of the simulation require different handling
/// of beams. Initial beams need special processing, internal beams require
/// threshold checks, and outgoing beams undergo far-field diffraction.
/// 
/// **How it Works**: Types determine which processing path the beam follows
/// and how its power contributions are categorized in the results.
#[derive(Debug, Clone, PartialEq)]
pub enum BeamType {
    Initial,      // Incident beam from illumination
    Default,      // Internal beam undergoing propagation
    OutGoing,     // Beam exiting the particle
    ExternalDiff, // Beam contributing to external diffraction
}
