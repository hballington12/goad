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

/// Regression switch (gate G2 of the lossy-media rework): when the
/// environment variable `GOAD_LEGACY_INTERFACE` is set, interfaces use
/// the pre-rework path (Macke `get_theta_t`, angle-test TIR, and
/// `cos_theta_t`-form Fresnel), which must reproduce pre-rework output
/// byte-identically. Removed at stage 6; see
/// docs/theory/lossy_media_implementation.md.
fn legacy_interface() -> bool {
    static LEGACY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *LEGACY.get_or_init(|| std::env::var("GOAD_LEGACY_INTERFACE").is_ok_and(|v| v != "0"))
}

/// Exact integral of `exp(g . (r - c))` over a planar polygon face with
/// centroid `c`: the ramp-weighted area `D_eff / cos(theta)` of the
/// theory note, reducing to the plain area as `g -> 0`. Evaluated by
/// the real-argument analogue of the stable Fraunhofer edge sum: with
/// the in-plane axis `u` chosen along `g`, Green's theorem gives
/// `(1/|g|) * sum_j dy_j exp(|g| x_j) phi(|g| dx_j)` with
/// `phi(t) = (e^t - 1)/t`, singularity-free (PGOM computes the same
/// quantity as its `PYEXP` contour form).
fn ramp_weighted_area(face: &Face, g: Vector3<f32>) -> f32 {
    let data = face.data();
    let area = data.area.unwrap();
    let normal = data.normal;
    // in-plane component of g
    let g_t = g - g.dot(&normal) * normal;
    let gm = g_t.norm();
    let verts = &data.exterior;
    let r_max = verts
        .iter()
        .map(|p| (p - data.midpoint).norm())
        .fold(0.0_f32, f32::max);
    if gm * r_max < 1e-4 {
        // relative error O((g R)^2) below the f32 budget
        return area;
    }
    let u = g_t / gm;
    let v = normal.cross(&u);

    let n = verts.len();
    let mut sum = 0.0_f32;
    for j in 0..n {
        let a = verts[j] - data.midpoint;
        let b = verts[(j + 1) % n] - data.midpoint;
        let x0 = a.dot(&u);
        let dy = (b - a).dot(&v);
        let t = gm * (b.dot(&u) - x0);
        let phi = if t.abs() < 1e-3 {
            1.0 + t * (0.5 + t / 6.0)
        } else {
            (t.exp() - 1.0) / t
        };
        sum += dy * (gm * x0).exp() * phi;
    }
    // Winding relative to the normal decides the sign of the contour
    // sum; the integral of a positive function is positive.
    (sum / gm).abs()
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

    /// The perpendicular basis bivector `h1 = n x w / sqrt((n x w).(n x w))`
    /// for an intersection at a face: the generalisation of
    /// [`Self::get_e_perp`] to inhomogeneous waves (Dupertuis eq. 13).
    /// Complex only when the wave's tangential attenuation direction
    /// differs from its tangential propagation direction; reduces to
    /// the real path exactly when the absorption vector is zero,
    /// keeping lossless media bit-stable. Sign conventions mirror
    /// [`Self::get_e_perp`].
    fn get_e_perp_bivector(&self, normal: &Vector3<f32>) -> (Vector3<f32>, Vector3<f32>) {
        if self.field.absorption() == Vector3::zeros()
            && self.field.e_perp_im() == Vector3::zeros()
        {
            return (self.get_e_perp(normal), Vector3::zeros());
        }
        let dot = normal.dot(&self.field.prop());
        if dot.abs() > 1.0 - settings::COLINEAR_THRESHOLD {
            return (-self.field.e_perp(), -self.field.e_perp_im());
        }
        let c_re = normal.cross(&(self.field.n_app() * self.field.prop()));
        let c_im = normal.cross(&self.field.absorption());
        // bilinear normalisation: h = c / sqrt(c . c)
        let cc = Complex::new(
            c_re.dot(&c_re) - c_im.dot(&c_im),
            2.0 * c_re.dot(&c_im),
        );
        let inv = cc.sqrt().inv();
        let h_re = inv.re * c_re - inv.im * c_im;
        let h_im = inv.re * c_im + inv.im * c_re;
        if dot > 0.0 {
            (-h_re, -h_im)
        } else {
            (h_re, h_im)
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
        let mut field = Field::new_identity(e_perp, prop)?;
        field.sync_wave_homogeneous(refr_index);
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
        mut field: Field,
        wavelength: f32,
    ) -> Self {
        field.sync_wave_homogeneous(refr_index);
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

        self.clipping_area = clipping
            .stats
            .as_ref()
            .map(|stats| stats.intersection_area + stats.remaining_area);

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
        let power_loss = self.power() - self.absorbed_power - self.interface_mixed - output_power;

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
            let mut field = if legacy_interface() {
                let e_perp = self.get_e_perp(&normal);
                self.field.new_from_e_perp(&e_perp)
            } else {
                let (e_re, e_im) = self.get_e_perp_bivector(&normal);
                self.field.new_from_bivector_basis(&e_re, &e_im)
            };

            let wavenumber = self.wavenumber();
            // Compute absorbed_intensity before applying the attenuation:
            // intensity is invariant under phase winding, and the
            // amplitude factor is known ahead of the propagation.
            let intensity_before = field.intensity();
            let absorbed_intensity;
            if legacy_interface() {
                let dist = (face.midpoint() - self.face.data().midpoint)
                    .dot(&self.field.prop())
                    .abs(); // z-distance
                absorbed_intensity =
                    intensity_before * (1.0 - (-wavenumber * n1.im * dist).exp().powi(2));
                field.propagate(dist, wavenumber, n1);
            } else {
                let d = face.midpoint() - self.face.data().midpoint;
                field.propagate_displacement(&d, wavenumber);
                // The absorbed power is booked from the exact tube
                // ledger below, not from a centroid intensity change,
                // which would conflate dissipation with transverse
                // ramp re-referencing.
                absorbed_intensity = 0.0;
            }
            if legacy_interface() {
                self.absorbed_power +=
                    absorbed_intensity * face.data().area.unwrap() * theta_i.cos() * n1.re;
            } else if field.absorption() != Vector3::zeros() {
                // Exact tube ledger (theory note, power section): the
                // absorbed power on this leg is the difference between
                // the power entering the flux tube at the parent face
                // and the power arriving at this face. Both are ramp
                // integrals over this face: the arriving intensity
                // carries the wave's full ramp `-2kA`, while the
                // no-dissipation transport of the parent's own face
                // ramp has `-2k(A - (A.e) grad_l)` with
                // `grad_l = n_p / (e . n_p)` the path-length gradient,
                // and the centroid prefactors differ by the transverse
                // (non-dissipative) part of the re-referencing only,
                // so nothing here overflows for strong absorption.
                let a = field.absorption();
                let e = self.field.prop();
                let n_p = self.face.data().normal;
                let d = face.midpoint() - self.face.data().midpoint;
                let a_e = a.dot(&e);
                let e_np = e.dot(&n_p);
                let l_c = d.dot(&n_p) / e_np; // centroid path length
                let grad_l = n_p / e_np;
                let cos_c = theta_i.cos();
                let n_app = field.n_app();

                let p_arr = field.intensity()
                    * n_app
                    * cos_c
                    * ramp_weighted_area(face, -2.0 * wavenumber * a);
                // transverse-only re-referencing factor: exp(-2k (A.d - (A.e) l_c))
                let trans = (-2.0 * wavenumber * (a.dot(&d) - a_e * l_c)).exp();
                let p_in = intensity_before
                    * trans
                    * n_app
                    * cos_c
                    * ramp_weighted_area(face, -2.0 * wavenumber * (a - a_e * grad_l));
                self.absorbed_power += p_in - p_arr;
            }

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
            let (refracted, reflected) = if legacy_interface() {
                (
                    self.create_refracted(face, theta_i, n1, n2, &field)
                        .unwrap_or(None),
                    self.create_reflected(face, theta_i, n1, n2, &field)
                        .unwrap_or(None),
                )
            } else {
                // Orient the normal from the incident medium into medium 2.
                let normal_t = if normal.dot(&self.field.prop()) > 0.0 {
                    normal
                } else {
                    -normal
                };
                let tr = crate::snell::transfer(
                    field.n_app() * field.prop(),
                    field.absorption(),
                    normal_t,
                    n2,
                );
                (
                    self.create_refracted_general(face, &tr, n1, n2, &field),
                    self.create_reflected_general(face, &tr, normal_t, n1, n2, &field),
                )
            };

            // Book the interface mixed term: for inhomogeneous incidence
            // R + T = 1 + M (Dupertuis II eq. 60), so the children's
            // powers need not sum to the arriving power; the imbalance
            // is interference flux at the interface, tracked in its own
            // ledger entry rather than polluting truncation counters.
            // Skipped for the initial beam (homogeneous incidence from
            // the lossless host has M = 0 exactly, and its outputs
            // include the Babinet external-diffraction duplicates).
            if !legacy_interface() && self.variant != BeamVariant::Initial {
                let deff =
                    ramp_weighted_area(face, -2.0 * wavenumber * field.absorption());
                let p_in = field.intensity() * field.n_app() * theta_i.cos() * deff;
                let p_out = refracted.as_ref().map_or(0.0, |b| b.power())
                    + reflected.as_ref().map_or(0.0, |b| b.power());
                self.interface_mixed += p_in - p_out;
            }

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
                    let d = face.data().midpoint - midpoint;
                    let mut field = beam.field.clone();
                    if legacy_interface() {
                        let arg = d.dot(&beam.field.prop()) * beam.wavenumber() * medium_refr_index.re;
                        field.wind(arg);
                    } else {
                        // Re-reference phase and (for ramped beams)
                        // amplitude to the sub-face centroid.
                        field.propagate_displacement(&d, beam.wavenumber());
                    }

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

    /// Creates the reflected beam from the reduced-wave-vector transfer.
    ///
    /// The propagation mirror is unchanged from the legacy path; what
    /// generalises is the wave state (the absorption vector mirrors
    /// with the frame) and the Fresnel coefficients, which take the
    /// complex normal components. TIR is no longer an angle test:
    /// a bounce counts as TIR when the transmitted normal flux
    /// `Re(kappa_t)` is negligible, which reduces exactly to the
    /// classical condition for lossless media and remains continuous
    /// as `n' -> 0`.
    fn create_reflected_general(
        &self,
        face: &Face,
        tr: &crate::snell::Transfer,
        normal_t: Vector3<f32>,
        n1: Complex<f32>,
        n2: Complex<f32>,
        field_in: &Field,
    ) -> Option<Beam> {
        let normal = face.data().normal;
        let prop = self.get_reflection_vector(&normal);
        let mut field = field_in.clone();
        field.set_prop(prop);

        // Mirror the absorption vector about the interface plane.
        let a = field.absorption();
        let a_n = a.dot(&normal_t);
        field.set_wave(field.n_app(), a - 2.0 * a_n * normal_t);

        let fresnel = fresnel::refl_general(tr.kappa_i, tr.kappa_t, n1, n2);
        field.matmul(&fresnel);
        debug_assert!(!Field::ampl_intensity(&field.ampl()).is_nan());

        let is_tir = tr.kappa_t.re <= settings::TIR_FLUX_EPS * tr.kappa_t.norm();
        let (rec, tir, variant) = if is_tir {
            (
                self.rec_count,
                self.tir_count + 1,
                BeamVariant::Default(DefaultBeamVariant::Tir),
            )
        } else {
            (
                self.rec_count + 1,
                self.tir_count,
                BeamVariant::Default(DefaultBeamVariant::Refl),
            )
        };
        Some(Beam::new(
            face.clone(),
            n1,
            rec,
            tir,
            field,
            variant,
            self.wavelength,
        ))
    }

    /// Creates the transmitted beam from the reduced-wave-vector
    /// transfer, or `None` when the transmitted wave carries no normal
    /// flux (the TIR limit). Direction, apparent index, and absorption
    /// vector all come from the transmitted `w2`; the Fresnel
    /// coefficients take the complex normal components.
    fn create_refracted_general(
        &self,
        face: &Face,
        tr: &crate::snell::Transfer,
        n1: Complex<f32>,
        n2: Complex<f32>,
        field_in: &Field,
    ) -> Option<Beam> {
        if tr.kappa_t.re <= settings::TIR_FLUX_EPS * tr.kappa_t.norm() {
            return None;
        }
        let n_app = tr.w_re.norm();
        let prop = tr.w_re / n_app;

        let mut field = field_in.clone();
        field.set_prop(prop);
        field.set_wave(n_app, tr.w_im);

        let fresnel = fresnel::refr_general(tr.kappa_i, tr.kappa_t, n1, n2);
        field.matmul(&fresnel);
        debug_assert!(!Field::ampl_intensity(&field.ampl()).is_nan());

        Some(Beam::new(
            face.clone(),
            n2,
            self.rec_count + 1,
            self.tir_count,
            field,
            BeamVariant::Default(DefaultBeamVariant::Refr),
            self.wavelength,
        ))
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
                let d = remainder.data().midpoint - self_midpoint;
                let mut field = self.field.clone();
                if legacy_interface() {
                    let arg = d.dot(&self.field.prop()) * self.wavenumber() * medium_refr_index.re;
                    field.wind(arg);
                } else {
                    field.propagate_displacement(&d, self.wavenumber());
                }

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
    /// Interface energy imbalance accumulated by this beam's Fresnel
    /// splittings: per interface, R + T = 1 + M for inhomogeneous
    /// incidence (Dupertuis II eq. 60). Booked to the powers ledger
    /// as `interface_mixed`. Zero for homogeneous incidence.
    pub interface_mixed: f32,
    /// Total area accounted for by intersections and remainders.
    /// `None` until a clip has successfully completed, so a failed or
    /// never-run clip cannot masquerade as a real zero-area result.
    pub clipping_area: Option<f32>,
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
            && self.interface_mixed == other.interface_mixed
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
        // The field's wave state is the caller's responsibility: fresh
        // fields are synced homogeneous at the true creation sites
        // (new_initial, new_from_field), transferred fields carry the
        // state set by create_{refracted,reflected}_general, and clones
        // (earcut, remainders, external diffraction) inherit it.
        debug_assert!(
            legacy_interface() || field.dispersion_residual(refr_index) < 1e-3,
            "beam field wave state inconsistent with its medium"
        );
        Self {
            id: next_beam_id(),
            face,
            refr_index,
            rec_count,
            tir_count,
            field,
            absorbed_power: 0.0,
            interface_mixed: 0.0,
            clipping_area: None,
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

    /// Returns the power of a beam: the normal Poynting flux integrated
    /// over the beam cross-section. For an inhomogeneous wave the flux
    /// factor is the apparent index `N` and the cross-section carries
    /// the amplitude ramp, `D_eff = cos(theta) * integral of
    /// exp(-2 k A_t . (r - c)) dA` (theory note, power section); both
    /// reduce exactly to the legacy `n.re * area * cos(theta)` when the
    /// absorption vector vanishes.
    pub fn power(&self) -> f32 {
        if legacy_interface() {
            return self.field.intensity() * self.refr_index.re * self.csa();
        }
        let a = self.field.absorption();
        if a == Vector3::zeros() {
            return self.field.intensity() * self.field.n_app() * self.csa();
        }
        let g = -2.0 * self.wavenumber() * a;
        let cosine = self.field.prop().dot(&self.face.data().normal).abs();
        self.field.intensity() * self.field.n_app() * cosine * ramp_weighted_area(&self.face, g)
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
            interface_mixed: self.interface_mixed,
            clipping_area: self.clipping_area,
            variant: self.variant.clone(),
            wavelength: self.wavelength,
        })
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    /// The ramp-weighted area must match the closed form on a unit
    /// square, be invariant under rigid rotation, ignore the
    /// out-of-plane component of `g`, and reduce to the area at `g = 0`.
    #[test]
    fn ramp_weighted_area_matches_analytic() {
        let verts = vec![
            Point3::new(0.0_f32, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let face = Face::new_simple(verts, None, None).unwrap();

        assert!((ramp_weighted_area(&face, Vector3::zeros()) - 1.0).abs() < 1e-6);

        // g along x, centroid at x = 1/2:
        // integral = 2 sinh(g/2) / g
        let g = 1.7_f32;
        let expected = 2.0 * (g / 2.0).sinh() / g;
        let got = ramp_weighted_area(&face, Vector3::new(g, 0.0, 0.0));
        assert!((got - expected).abs() < 1e-5, "{got} vs {expected}");

        // out-of-plane component must not contribute
        let got_z = ramp_weighted_area(&face, Vector3::new(g, 0.0, 3.0));
        assert!((got_z - expected).abs() < 1e-5);

        // rigid rotation invariance
        let rot = nalgebra::Rotation3::from_euler_angles(0.4_f32, -0.9, 1.3);
        let rverts: Vec<Point3<f32>> = face
            .data()
            .exterior
            .iter()
            .map(|p| rot * p)
            .collect();
        let rface = Face::new_simple(rverts, None, None).unwrap();
        let rg = rot * Vector3::new(g, 0.0, 0.0);
        let got_r = ramp_weighted_area(&rface, rg);
        assert!((got_r - expected).abs() < 1e-4, "{got_r} vs {expected}");
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
