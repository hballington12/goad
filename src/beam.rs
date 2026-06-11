use anyhow::Result;
use std::f32::consts::PI;
use std::sync::atomic::{AtomicU64, Ordering};

use nalgebra::{Complex, Matrix4, Vector3};

use crate::{
    clip::Clipping,
    field::Field,
    fresnel,
    geom::{Face, Geom},
    settings,
    snell::get_theta_t,
};

/// Process-unique identifier for a beam. Assigned at construction in
/// `Beam::new`; preserved across `Clone` (a clone is the same logical beam).
pub type BeamId = u64;

static BEAM_ID_COUNTER: AtomicU64 = AtomicU64::new(0);

fn next_beam_id() -> BeamId {
    BEAM_ID_COUNTER.fetch_add(1, Ordering::Relaxed)
}

/// Where an output beam goes after a propagation event. Shared by the
/// production batch path (via `categorise_outputs`) and the recording path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputKind {
    /// Re-propagated as a new internal beam.
    NearField,
    /// Exits geometry — feeds the far-field diffraction stage.
    OutGoing,
    /// External (silhouette) diffraction beam from the initial pass.
    ExternalDiff,
}

/// Pure classification of an (input, output) variant pair. Returns None for
/// combinations the solver drops on the floor (e.g. anything emitted by an
/// OutGoing beam — which never happens but the match must be total).
pub fn classify_output(input: &BeamVariant, output: &BeamVariant) -> Option<OutputKind> {
    match (input, output) {
        (BeamVariant::Default(..), BeamVariant::Default(..))
        | (BeamVariant::Initial, BeamVariant::Default(..)) => Some(OutputKind::NearField),
        (BeamVariant::Default(..), BeamVariant::OutGoing) => Some(OutputKind::OutGoing),
        (BeamVariant::Initial, BeamVariant::ExternalDiff) => Some(OutputKind::ExternalDiff),
        _ => None,
    }
}

impl Beam {
    /// Determines the refractive index of the second medium when a beam intersects with a face.
    fn get_n2(
        &self,
        geom: &Geom,
        face: &Face,
        normal: Vector3<f32>,
        medium_refr_index: Complex<f32>,
    ) -> Complex<f32> {
        let id = face.data().shape_id.unwrap();
        if normal.dot(&self.field.prop()) < 0.0 {
            geom.shapes[id].refr_index
        } else {
            geom.n_out(id, medium_refr_index)
        }
    }

    /// Determines the new `e_perp` vector for an intersection at a `face`.
    fn get_e_perp(&self, normal: &Vector3<f32>) -> Vector3<f32> {
        let dot = normal.dot(&self.field.prop());
        let e_perp = if dot.abs() > 1.0 - settings::COLINEAR_THRESHOLD {
            -self.field.e_perp()
        } else {
            normal.cross(&self.field.prop()).normalize() // new e_perp
        };
        if dot > 0.0 {
            -e_perp
        } else {
            e_perp
        }
    }

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
        let rec = 0;
        let tir = 0;
        Ok(Beam::new(
            face,
            refr_index,
            rec,
            tir,
            field,
            BeamVariant::Initial,
            wavelength,
        ))
    }

    pub fn new_from_field(
        face: Face,
        refr_index: Complex<f32>,
        field: Field,
        wavelength: f32,
    ) -> Self {
        let rec = 0;
        let tir = 0;
        Beam::new(
            face,
            refr_index,
            rec,
            tir,
            field,
            BeamVariant::Initial,
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
        let prop = self.field.prop();
        let mut clipping = Clipping::new(geom, &mut self.face, &prop);
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
        let intersection_beams =
            self.intersections_to_beams(geom, intersections, medium_refr_index);

        let mut output_beams = Vec::new();
        output_beams.extend(intersection_beams);
        output_beams.extend(remainder_beams);
        let output_power = output_beams.iter().fold(0.0, |acc, x| acc + x.power());
        let power_loss = self.power() - self.absorbed_power - output_power;

        Ok((output_beams, power_loss))
    }

    fn intersections_to_beams(
        &mut self,
        geom: &mut Geom,
        intersections: Vec<Face>,
        medium_refr_index: Complex<f32>,
    ) -> Vec<Beam> {
        let n1 = self.refr_index;
        let mut outputs = Vec::new();
        for face in &intersections {
            let normal = face.data().normal;
            let theta_i = normal.dot(&self.field.prop()).abs().clamp(-1.0, 1.0).acos();
            let n2 = self.get_n2(geom, face, normal, medium_refr_index);
            let e_perp = self.get_e_perp(&normal);
            let mut field = self.field.new_from_e_perp(&e_perp);

            let dist = (face.midpoint() - self.face.data().midpoint)
                .dot(&self.field.prop())
                .abs(); // z-distance
            let wavenumber = self.wavenumber();
            // Compute absorbed_intensity before applying the attenuation —
            // intensity is invariant under phase winding so it can be read
            // after `propagate` increments phase but before it attenuates.
            let absorbed_intensity =
                field.intensity() * (1.0 - (-wavenumber * n1.im * dist).exp().powi(2));
            field.propagate(dist, wavenumber, n1);
            self.absorbed_power +=
                absorbed_intensity * face.data().area.unwrap() * theta_i.cos() * n1.re;

            if self.variant == BeamVariant::Initial {
                if let Ok(flipped_face) = face.flipped() {
                    let external_diff = Beam::new(
                        flipped_face,
                        n1,
                        self.rec_count + 1,
                        self.tir_count,
                        field.clone(),
                        BeamVariant::ExternalDiff,
                        self.wavelength,
                    );
                    outputs.push(external_diff);
                }
            }

            // untracked energy leaks can occur here if the amplitude matrix contains NaN values
            let refracted = self
                .create_refracted(face, theta_i, n1, n2, &field)
                .unwrap_or(None);
            let reflected = self
                .create_reflected(face, theta_i, n1, n2, &field)
                .unwrap_or(None);

            if refracted.is_some() {
                outputs.push(refracted.unwrap());
            }
            if reflected.is_some() {
                outputs.push(reflected.unwrap());
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
                    let dist = (face.data().midpoint - midpoint).dot(&beam.field.prop());
                    let arg = dist * beam.wavenumber() * medium_refr_index.re;
                    let mut field = beam.field.clone();
                    field.wind(arg);

                    let new_beam = Beam::new(
                        face,
                        beam.refr_index,
                        beam.rec_count,
                        beam.tir_count,
                        field,
                        beam.variant.clone(),
                        beam.wavelength,
                    );

                    outputs.push(new_beam);
                }
            }
        }
        outputs
    }

    /// Computes the polar scattering angle of a beam in degrees
    fn get_polar_angle(&self) -> f32 {
        ((-self.field.prop()[2]).acos().to_degrees()).abs() // compute outgoing theta
    }

    /// Computes the azimuthal scattering angle of a beam in degrees. Returns a value in the range [0, 360)
    fn get_azimuthal_angle(&self) -> f32 {
        let kx = self.field.prop()[0];
        let ky = self.field.prop()[1];
        let mut phi = ky.atan2(kx).to_degrees();
        if phi < 0.0 {
            phi += 360.0
        }
        phi
    }

    /// Returns the polar and azimuthal scattering angles of a beam in degrees.
    pub fn get_scattering_angles(&self) -> (f32, f32) {
        let theta = self.get_polar_angle();
        let phi = self.get_azimuthal_angle();
        (theta, phi)
    }
}

/// Converts the remainder faces from a clipping into beams with the same field
/// properties as the original beam.
impl Beam {
    /// Returns a transmitted propagation vector, where `stt` is the sine of the angle of transmission.
    fn get_refraction_vector(
        &self,
        norm: &Vector3<f32>,
        theta_i: f32,
        theta_t: f32,
    ) -> Vector3<f32> {
        let prop = self.field.prop();
        if theta_t.sin() < settings::COLINEAR_THRESHOLD {
            return prop;
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

        debug_assert!(
            (theta_t.cos() - result.dot(&norm).abs()).abs() < settings::COLINEAR_THRESHOLD
        );

        result
    }

    fn get_reflection_vector(&self, norm: &Vector3<f32>) -> Vector3<f32> {
        let prop = self.field.prop();
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
    /// Creates a new reflected beam
    fn create_reflected(
        &self,
        face: &Face,
        theta_i: f32,
        n1: Complex<f32>,
        n2: Complex<f32>,
        field_in: &Field,
    ) -> Result<Option<Beam>> {
        let normal = face.data().normal;
        let prop = self.get_reflection_vector(&normal);
        let mut field = field_in.clone();
        field.set_prop(prop);

        debug_assert!((field.prop().dot(&normal) - theta_i.cos()) < settings::COLINEAR_THRESHOLD);
        debug_assert!(!Field::ampl_intensity(&field.ampl()).is_nan());

        let fresnel = fresnel::refl(n1, n2, theta_i);
        field.matmul(&fresnel);
        debug_assert!(!Field::ampl_intensity(&field.ampl()).is_nan());

        debug_assert!(!Field::ampl_intensity(&field.ampl()).is_nan());

        if theta_i > (n2.re / n1.re).asin() {
            Ok(Some(Beam::new(
                face.clone(),
                n1,
                self.rec_count,
                self.tir_count + 1,
                field,
                BeamVariant::Default(DefaultBeamVariant::Tir),
                self.wavelength,
            )))
        } else {
            Ok(Some(Beam::new(
                face.clone(),
                n1,
                self.rec_count + 1,
                self.tir_count,
                field,
                BeamVariant::Default(DefaultBeamVariant::Refl),
                self.wavelength,
            )))
        }
    }

    /// Creates a new refracted beam.
    fn create_refracted(
        &self,
        face: &Face,
        theta_i: f32,
        n1: Complex<f32>,
        n2: Complex<f32>,
        field_in: &Field,
    ) -> Result<Option<Beam>> {
        let mut field = field_in.clone();
        let normal = face.data().normal;
        if theta_i >= (n2.re / n1.re).asin() {
            // if total internal reflection
            Ok(None)
        } else {
            let theta_t = get_theta_t(theta_i, n1, n2)?; // sin(theta_t)
            let prop = self.get_refraction_vector(&normal, theta_i, theta_t);
            let fresnel = fresnel::refr(n1, n2, theta_i);

            field.set_prop(prop);
            field.matmul(&fresnel);

            debug_assert!(field.prop().dot(&prop) > 0.0);
            debug_assert!(
                (field.prop().dot(&normal).abs() - theta_t.cos()).abs()
                    < settings::COLINEAR_THRESHOLD
            );

            Ok(Some(Beam::new(
                face.clone(),
                n2,
                self.rec_count + 1,
                self.tir_count,
                field,
                BeamVariant::Default(DefaultBeamVariant::Refr),
                self.wavelength,
            )))
        }
    }

    fn remainders_to_beams(
        &self,
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
                let dist = (remainder.data().midpoint - self_midpoint).dot(&self.field.prop());
                let arg = dist * self.wavenumber() * medium_refr_index.re;
                let mut field = self.field.clone();
                field.wind(arg);

                Some(Beam::new(
                    remainder,
                    self.refr_index,
                    self.rec_count,
                    self.tir_count,
                    field,
                    BeamVariant::OutGoing,
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
#[derive(Debug, Clone)]
pub struct Beam {
    /// Process-unique identity. Set in `new`; preserved by `Clone` so cross-
    /// thread / cross-batch clones are recognised as the same logical beam.
    pub id: BeamId,
    pub face: Face,
    pub refr_index: Complex<f32>,
    pub rec_count: i32,
    pub tir_count: i32,
    pub field: Field,
    pub absorbed_power: f32,  // power absorbed by the medium
    pub clipping_area: f32,   // total area accounted for by intersections and remainders
    pub variant: BeamVariant, // type of beam, e.g. initial, default, outgoing, external diff
    pub wavelength: f32,
}

// Equality is over physical content; identity (`id`) is excluded so a clone
// (which intentionally shares the parent's id) still compares equal to its
// source.
impl PartialEq for Beam {
    fn eq(&self, other: &Self) -> bool {
        self.face == other.face
            && self.refr_index == other.refr_index
            && self.rec_count == other.rec_count
            && self.tir_count == other.tir_count
            && self.field == other.field
            && self.absorbed_power == other.absorbed_power
            && self.clipping_area == other.clipping_area
            && self.variant == other.variant
            && self.wavelength == other.wavelength
    }
}

/// Creates a new beam
impl Beam {
    pub fn new(
        face: Face,
        refr_index: Complex<f32>,
        rec_count: i32,
        tir_count: i32,
        field: Field,
        variant: BeamVariant,
        wavelength: f32,
    ) -> Self {
        Self {
            id: next_beam_id(),
            face,
            refr_index,
            rec_count,
            tir_count,
            field,
            absorbed_power: 0.0,
            clipping_area: 0.0,
            variant,
            wavelength,
        }
    }

    /// Returns the cross sectional area of the beam.
    pub fn csa(&self) -> f32 {
        let area = self.face.data().area.unwrap();
        let norm = self.face.data().normal;
        let cosine = self.field.prop().dot(&norm).abs();

        area * cosine
    }

    /// Returns the power of a beam.
    pub fn power(&self) -> f32 {
        self.field.intensity() * self.refr_index.re * self.csa()
    }

    pub fn wavenumber(&self) -> f32 {
        2.0 * PI / self.wavelength
    }

    /// Returns a new Beam with the given 4x4 transformation matrix applied.
    /// The transformation is applied to the face (vertices, normal, midpoint)
    /// and the field (prop and e_perp vectors are rotated by the upper-left 3x3).
    pub fn transformed(&self, transform: &Matrix4<f32>) -> Result<Self> {
        // Clone and transform the face
        let mut new_face = self.face.clone();
        new_face.transform(transform)?;

        // Extract the 3x3 rotation part from the 4x4 matrix
        let rot3 = transform.fixed_view::<3, 3>(0, 0).into_owned();

        // Rotate the field
        let new_field = self.field.rotated(&rot3);

        Ok(Self {
            id: self.id,
            face: new_face,
            refr_index: self.refr_index,
            rec_count: self.rec_count,
            tir_count: self.tir_count,
            field: new_field,
            absorbed_power: self.absorbed_power,
            clipping_area: self.clipping_area,
            variant: self.variant.clone(),
            wavelength: self.wavelength,
        })
    }

}

#[derive(Debug, Clone, PartialEq)]
pub enum DefaultBeamVariant {
    Refl, // refraction
    Refr, // reflection
    Tir,  // total internal reflection
}

#[derive(Debug, Clone, PartialEq)]
pub enum BeamVariant {
    Initial,
    Default(DefaultBeamVariant),
    OutGoing,
    ExternalDiff,
}
