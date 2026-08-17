use crate::result::{ApproxEq, Mueller};
#[cfg(debug_assertions)]
use crate::settings;

use anyhow::Result;
use std::fmt::Debug;

use nalgebra::{Complex, Matrix2, Matrix3, RealField, Vector3};

#[cfg(test)]
mod tests {

    use super::*;
    use crate::settings;

    #[test]
    fn identity_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let field = Field::new_identity(e_perp, prop).unwrap();
        assert!((field.intensity() - 1.0).abs() < settings::COLINEAR_THRESHOLD);
    }

    /// The homogeneous wave state satisfies the dispersion relation for
    /// its medium, transforms correctly under rotation, and is zero
    /// absorption in a lossless medium.
    #[test]
    fn wave_state_homogeneous_dispersion() {
        let prop = Vector3::new(0.6_f32, 0.0, 0.8);
        let e_perp = Vector3::new(0.0_f32, 1.0, 0.0);
        let mut field = Field::new_identity(e_perp, prop).unwrap();

        let lossless = Complex::new(1.31, 0.0);
        field.sync_wave_homogeneous(lossless);
        assert_eq!(field.absorption(), Vector3::zeros());
        assert_eq!(field.n_app(), 1.31);
        assert!(field.dispersion_residual(lossless) < 1e-5);

        let lossy = Complex::new(1.31, 0.1);
        field.sync_wave_homogeneous(lossy);
        assert!((field.absorption() - 0.1 * prop).norm() < 1e-7);
        assert!(field.dispersion_residual(lossy) < 1e-5);

        // Rotation carries the absorption vector with the frame.
        let rot = nalgebra::Rotation3::from_euler_angles(0.3_f32, 1.1, -0.4);
        let rotated = field.rotated(rot.matrix());
        assert!(rotated.dispersion_residual(lossy) < 1e-4);
        assert!((rotated.absorption() - 0.1 * rotated.prop()).norm() < 1e-6);
    }

    /// For a homogeneous wave (absorption parallel to prop) and real
    /// bases, the bivector rotation must agree with the real one: the
    /// complex factors cancel identically.
    #[test]
    fn bivector_rotation_reduces_to_real() {
        let prop = Vector3::new(0.6_f32, 0.0, 0.8);
        let e_perp = Vector3::y();
        let mut field = Field::new_identity(e_perp, prop).unwrap();
        field.sync_wave_homogeneous(Complex::new(1.31, 0.1));
        field.set_ampl(Ampl::new(
            Complex::new(0.8, 0.1),
            Complex::new(-0.2, 0.05),
            Complex::new(0.1, -0.3),
            Complex::new(0.7, 0.2),
        ));

        // A unit vector perpendicular to prop, distinct from e_perp.
        let target = prop.cross(&Vector3::y()).normalize();

        let real = field.new_from_e_perp(&target);
        let biv = field.new_from_bivector_basis(&target, &Vector3::zeros());
        for (a, b) in real.ampl().iter().zip(biv.ampl().iter()) {
            assert!((a - b).norm() < 1e-6, "{a} vs {b}");
        }
    }

    /// Vector-displacement propagation must reduce to the legacy
    /// scalar form for homogeneous waves, and must implement the ramp
    /// semantics for inhomogeneous ones: transverse displacement
    /// changes amplitude but not phase when `A` is transverse.
    #[test]
    fn propagate_displacement_semantics() {
        let prop = Vector3::new(0.6_f32, 0.0, 0.8);
        let e_perp = Vector3::y();
        let k = 11.8;
        let n = Complex::new(1.31, 0.01);

        // Homogeneous reduction, displacement along prop.
        let mut a = Field::new_identity(e_perp, prop).unwrap();
        a.sync_wave_homogeneous(n);
        let mut b = a.clone();
        let dist = 3.7_f32;
        a.propagate(dist, k, n);
        b.propagate_displacement(&(dist * prop), k);
        assert!((a.phase() - b.phase()).abs() < 1e-3);
        assert!((a.intensity() - b.intensity()).abs() / a.intensity() < 1e-5);

        // Transverse ramp: A perpendicular to prop (an exiting beam).
        let mut c = Field::new_identity(e_perp, prop).unwrap();
        let a_vec = 0.02_f32 * Vector3::y();
        c.set_wave(1.0, a_vec);
        let i0 = c.intensity();
        // along prop: pure phase, no attenuation
        let f = c.propagate_displacement(&(2.0 * prop), k);
        assert_eq!(f, 1.0);
        assert!((c.intensity() - i0).abs() < 1e-7);
        let p0 = c.phase();
        // transverse: pure amplitude re-reference, no phase
        let f = c.propagate_displacement(&Vector3::new(0.0, 1.5, 0.0), k);
        let expected = (-k * 0.02 * 1.5_f32).exp();
        assert!((f - expected).abs() < 1e-6);
        assert!((c.intensity() / i0 - expected * expected).abs() < 1e-4);
        assert!((c.phase() - p0).abs() < 1e-6);
    }

    /// Rotating into a complex bivector basis and straight back must be
    /// the identity: the rotation matrices are mutually inverse under
    /// the bilinear algebra.
    #[test]
    fn bivector_rotation_round_trip() {
        // An inhomogeneous wave: oblique propagation, absorption along
        // the entry normal (the first-entry configuration).
        let th = 40.0_f32.to_radians();
        let prop = Vector3::new(th.sin(), 0.0, th.cos());
        let w = (1.31_f32 * prop, 0.1_f32 * Vector3::<f32>::z());

        // Bivector bases from two facet normals: h = n x w, bilinearly
        // normalised (h . w = 0 exactly by construction).
        let make_basis = |n: Vector3<f32>| {
            let c_re: Vector3<f32> = n.cross(&w.0);
            let c_im: Vector3<f32> = n.cross(&w.1);
            let cc = Complex::new(c_re.dot(&c_re) - c_im.dot(&c_im), 2.0 * c_re.dot(&c_im));
            let inv: Complex<f32> = cc.sqrt().inv();
            (inv.re * c_re - inv.im * c_im, inv.re * c_im + inv.im * c_re)
        };
        let e1 = make_basis(Vector3::new(0.3, 0.2, 0.93).normalize());
        let e2 = make_basis(Vector3::new(-0.1, 0.5, 0.86).normalize());

        let m12 = Field::rotation_matrix_bivector(e1, e2, w);
        let m21 = Field::rotation_matrix_bivector(e2, e1, w);
        let prod = m12 * m21;
        let id = Matrix2::<Complex<f32>>::identity();
        for (a, b) in prod.iter().zip(id.iter()) {
            assert!((a - b).norm() < 1e-5, "round trip not identity: {prod}");
        }
        // The bases themselves must be genuinely complex here, or the
        // test proves nothing.
        assert!(e1.1.norm() > 1e-3 && e2.1.norm() > 1e-3);
    }

    #[test]
    fn field_partial_eq_with_tolerance() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();

        // Create two nearly identical fields
        let field1 = Field::new_identity(e_perp, prop).unwrap();

        // Create a slightly different field
        let e_perp2 = Vector3::new(1.0 + 1e-6, 0.0, 0.0).normalize();
        let mut field2 = Field::new_identity(e_perp2, prop).unwrap();
        field2.phase = 1e-6;

        // Should be equal within tolerance (1e-5)
        assert_eq!(field1, field2);

        // Create a field with larger differences
        let e_perp3 = Vector3::new(1.0 + 1e-4, 0.0, 0.0).normalize();
        let mut field3 = Field::new_identity(e_perp3, prop).unwrap();
        field3.phase = 1e-4;

        // Should not be equal (difference exceeds tolerance)
        assert_ne!(field1, field3);
    }
}

pub type Ampl = Matrix2<Complex<f32>>;

pub trait AmplMatrix {
    fn s11(&self) -> f32;
    fn s12(&self) -> f32;
    fn s13(&self) -> f32;
    fn s14(&self) -> f32;
    fn s21(&self) -> f32;
    fn s22(&self) -> f32;
    fn s23(&self) -> f32;
    fn s24(&self) -> f32;
    fn s31(&self) -> f32;
    fn s32(&self) -> f32;
    fn s33(&self) -> f32;
    fn s34(&self) -> f32;
    fn s41(&self) -> f32;
    fn s42(&self) -> f32;
    fn s43(&self) -> f32;
    fn s44(&self) -> f32;
    fn to_mueller(&self) -> Mueller;
    fn zeros() -> Self;
    fn identity() -> Self;
    fn is_valid(&self) -> bool;
}

impl AmplMatrix for Ampl {
    fn s11(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj()
            + self[(0, 1)] * self[(0, 1)].conj()
            + self[(1, 0)] * self[(1, 0)].conj()
            + self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s12(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj() - self[(0, 1)] * self[(0, 1)].conj()
            + self[(1, 0)] * self[(1, 0)].conj()
            - self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s13(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() + self[(1, 1)] * self[(1, 0)].conj()).re
    }
    fn s14(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() - self[(1, 1)] * self[(1, 0)].conj()).im
    }
    fn s21(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj() + self[(0, 1)] * self[(0, 1)].conj()
            - self[(1, 0)] * self[(1, 0)].conj()
            - self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s22(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj()
            - self[(0, 1)] * self[(0, 1)].conj()
            - self[(1, 0)] * self[(1, 0)].conj()
            + self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s23(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() - self[(1, 1)] * self[(1, 0)].conj()).re
    }
    fn s24(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() + self[(1, 1)] * self[(1, 0)].conj()).im
    }
    fn s31(&self) -> f32 {
        (self[(0, 0)] * self[(1, 0)].conj() + self[(1, 1)] * self[(0, 1)].conj()).re
    }
    fn s32(&self) -> f32 {
        (self[(0, 0)] * self[(1, 0)].conj() - self[(1, 1)] * self[(0, 1)].conj()).re
    }
    fn s33(&self) -> f32 {
        (self[(0, 0)] * self[(1, 1)].conj() + self[(0, 1)] * self[(1, 0)].conj()).re
    }
    fn s34(&self) -> f32 {
        (self[(0, 0)] * self[(1, 1)].conj() + self[(1, 0)] * self[(0, 1)].conj()).im
    }
    fn s41(&self) -> f32 {
        (self[(1, 0)] * self[(0, 0)].conj() + self[(1, 1)] * self[(0, 1)].conj()).im
    }
    fn s42(&self) -> f32 {
        (self[(1, 0)] * self[(0, 0)].conj() - self[(1, 1)] * self[(0, 1)].conj()).im
    }
    fn s43(&self) -> f32 {
        (self[(1, 1)] * self[(0, 0)].conj() - self[(0, 1)] * self[(1, 0)].conj()).im
    }
    fn s44(&self) -> f32 {
        (self[(1, 1)] * self[(0, 0)].conj() - self[(0, 1)] * self[(1, 0)].conj()).re
    }
    fn to_mueller(&self) -> Mueller {
        Mueller::new(
            self.s11(),
            self.s12(),
            self.s13(),
            self.s14(),
            self.s21(),
            self.s22(),
            self.s23(),
            self.s24(),
            self.s31(),
            self.s32(),
            self.s33(),
            self.s34(),
            self.s41(),
            self.s42(),
            self.s43(),
            self.s44(),
        )
    }
    fn zeros() -> Self {
        Self::zeros()
    }
    fn identity() -> Self {
        Self::identity()
    }
    fn is_valid(&self) -> bool {
        self.iter().all(|c| c.re.is_finite() && c.im.is_finite())
    }
}

/// Essentially represents a plane wave field with phase, amplitude matrix, and corresponding electric field vectors.
#[derive(Debug, Clone)]
pub struct Field {
    ampl: Ampl,
    e_perp: Vector3<f32>,
    phase: f32,
    prop: Vector3<f32>,
    /// Apparent refractive index `N` of the wave: the magnitude of the
    /// real part of the reduced wave vector `w = N e + iA`. For a
    /// homogeneous wave (the only kind the tracer produces until the
    /// lossy-media transfer lands) this equals `refr_index.re`.
    /// See docs/theory/lossy_media_beam_tracing.typ.
    n_app: f32,
    /// Absorption vector `A`: the imaginary part of the reduced wave
    /// vector, i.e. the gradient of the log-amplitude in units of the
    /// wavenumber. Amplitude varies as `exp(-k A . r)`. Homogeneous
    /// waves have `A = refr_index.im * prop`; a beam in a lossless
    /// medium has `A = 0`.
    absorption: Vector3<f32>,
    /// Imaginary part of the perpendicular polarisation basis vector:
    /// the full basis is the bivector `e_perp + i e_perp_im`
    /// (Dupertuis eq. 13, `h1 = n x w / sqrt((n x w).(n x w))`),
    /// which is genuinely complex when the wave's tangential phase and
    /// attenuation directions differ at the facet that created it.
    /// Zero for every homogeneous wave, which is what keeps the real
    /// fast paths exact for lossless media.
    e_perp_im: Vector3<f32>,
}

impl PartialEq for Field {
    /// Compare two Field structs with a tolerance of 1e-5
    fn eq(&self, other: &Self) -> bool {
        const TOLERANCE: f32 = 1e-5;

        // Check amplitude matrices using ApproxEq trait
        self.ampl.approx_eq(&other.ampl, TOLERANCE)
            // Check vectors component-wise with tolerance
            && (self.e_perp.x - other.e_perp.x).abs() < TOLERANCE
            && (self.e_perp.y - other.e_perp.y).abs() < TOLERANCE
            && (self.e_perp.z - other.e_perp.z).abs() < TOLERANCE
            && (self.prop.x - other.prop.x).abs() < TOLERANCE
            && (self.prop.y - other.prop.y).abs() < TOLERANCE
            && (self.prop.z - other.prop.z).abs() < TOLERANCE
            // Check phase with tolerance
            && (self.phase - other.phase).abs() < TOLERANCE
    }
}

impl Field {
    /// Multiply the field amplitude by a complex-valued rotation matrix
    pub fn matmul(&mut self, rot: &Matrix2<Complex<f32>>) {
        self.ampl = rot * self.ampl;
    }

    /// Multiply the field amplitude by a real-valued factor
    pub fn mul(&mut self, fac: f32) {
        self.ampl *= Complex::new(fac, 0.0);
    }

    /// Increment the phase by a real value
    pub fn wind(&mut self, arg: f32) {
        self.phase += arg;
    }

    /// Propagate this field forward by `distance` through a medium of
    /// (possibly complex) refractive index `n` at wavenumber `k`. Phase
    /// is wound by `distance · k · n.re`; amplitude is attenuated by
    /// `exp(-k · n.im · distance)`.
    pub fn propagate(&mut self, distance: f32, wavenumber: f32, n: Complex<f32>) {
        self.wind(distance * wavenumber * n.re);
        self.mul((-wavenumber * n.im * distance).exp());
    }

    /// Propagate this field by the vector displacement `d` using the
    /// stored wave state: phase is wound by `k N (prop . d)` and the
    /// amplitude is scaled by `exp(-k A . d)`. This is the general form
    /// of [`Self::propagate`]: for a homogeneous wave (`A = n' prop`)
    /// and `d` along `prop` the two agree, while for an inhomogeneous
    /// wave the attenuation depends on the full displacement, not just
    /// its along-`prop` component (`A . d != (A . prop)(prop . d)`).
    /// Also serves as the amplitude re-referencing when a beam's
    /// reference point moves within its cross-section (earcut,
    /// remainders), where a transverse ramp changes the amplitude at
    /// the new centroid. Returns the amplitude factor applied, which
    /// callers use for absorbed-power tallies.
    pub fn propagate_displacement(&mut self, d: &Vector3<f32>, wavenumber: f32) -> f32 {
        let dist = self.prop.dot(d);
        self.wind(dist * wavenumber * self.n_app);
        let a_d = self.absorption.dot(d);
        let factor = (-wavenumber * a_d).exp();
        if factor != 1.0 {
            self.mul(factor);
        }
        factor
    }

    // Getters

    /// Returns the propagation vector of the field
    pub fn prop(&self) -> Vector3<f32> {
        self.prop
    }

    /// Returns the amplitude of the field
    pub fn ampl(&self) -> Ampl {
        let phase = self.phase;
        self.ampl_wo_phase() * Complex::new(phase.cos(), phase.sin())
    }

    /// Returns the amplitude of the field without phase due to path length
    pub fn ampl_wo_phase(&self) -> Ampl {
        self.ampl
    }

    /// Returns the phase of the field
    pub fn phase(&self) -> f32 {
        self.phase
    }

    /// Returns the apparent refractive index `N` of the wave (the
    /// magnitude of the real part of the reduced wave vector).
    pub fn n_app(&self) -> f32 {
        self.n_app
    }

    /// Returns the absorption vector `A` (imaginary part of the reduced
    /// wave vector): amplitude varies across space as `exp(-k A . r)`.
    pub fn absorption(&self) -> Vector3<f32> {
        self.absorption
    }

    /// Sets the wave state explicitly (used by the interface transfer,
    /// which produces inhomogeneous waves in absorbing media).
    pub fn set_wave(&mut self, n_app: f32, absorption: Vector3<f32>) {
        self.n_app = n_app;
        self.absorption = absorption;
    }

    /// Sets the wave state to the homogeneous wave of the given medium:
    /// `N = n.re`, `A = n.im * prop`. This is exact for every beam the
    /// tracer currently produces; the lossy-media interface transfer
    /// will replace these values with the inhomogeneous ones.
    pub fn sync_wave_homogeneous(&mut self, refr_index: Complex<f32>) {
        self.n_app = refr_index.re;
        self.absorption = refr_index.im * self.prop;
        debug_assert!(self.dispersion_residual(refr_index) < 1e-4);
    }

    /// Residual of the dispersion relation `w . w = n_tilde^2` for the
    /// stored wave state, with `w = N prop + i A`. Zero (to rounding)
    /// for any physical wave in the medium of index `refr_index`.
    pub fn dispersion_residual(&self, refr_index: Complex<f32>) -> f32 {
        let re = self.n_app * self.prop;
        let a = self.absorption;
        let ww_re = re.dot(&re) - a.dot(&a);
        let ww_im = 2.0 * re.dot(&a);
        let n2 = refr_index * refr_index;
        ((ww_re - n2.re).powi(2) + (ww_im - n2.im).powi(2)).sqrt()
    }

    /// Returns the perpendicular component of the electric field
    pub fn e_perp(&self) -> Vector3<f32> {
        self.e_perp
    }

    /// Returns the parallel component of the electric field
    pub fn e_par(&self) -> Vector3<f32> {
        self.e_perp.cross(&self.prop).normalize()
    }

    // Setters

    pub fn set_ampl(&mut self, ampl0: Ampl) {
        self.ampl = ampl0;
    }

    pub fn set_phase(&mut self, phase: f32) {
        self.phase = phase;
    }

    pub fn set_e_perp(&mut self, e_perp: &Vector3<f32>) {
        self.e_perp = *e_perp;
        // A plain (real) basis assignment discards any bivector part.
        self.e_perp_im = Vector3::zeros();
    }

    /// Returns the imaginary part of the perpendicular basis bivector.
    pub fn e_perp_im(&self) -> Vector3<f32> {
        self.e_perp_im
    }

    pub fn set_prop(&mut self, prop: Vector3<f32>) {
        self.prop = prop;
    }

    /// Returns a rotation matrix for rotating from the current plane perpendicular to the plane perpendicular to `e_perp`.
    pub fn get_rotation_matrix(&self, e_perp: &Vector3<f32>) -> Matrix2<Complex<f32>> {
        Field::rotation_matrix(self.e_perp, *e_perp, self.prop)
            .map(|x| nalgebra::Complex::new(x, 0.0))
    }

    pub fn new_from_e_perp(&self, e_perp: &Vector3<f32>) -> Self {
        let rot = self.get_rotation_matrix(e_perp);
        let mut field = self.clone();
        field.set_e_perp(e_perp);
        field.matmul(&rot);
        field
    }

    /// Rotates the amplitude into a new (generally complex) bivector
    /// basis `e_out = e_out_re + i e_out_im`, the generalisation of
    /// [`Self::new_from_e_perp`] to inhomogeneous waves. When neither
    /// the current basis, the target basis, nor the wave carries an
    /// imaginary part, this reduces exactly to the real path (same
    /// code, so lossless media stay bit-stable).
    pub fn new_from_bivector_basis(
        &self,
        e_out_re: &Vector3<f32>,
        e_out_im: &Vector3<f32>,
    ) -> Self {
        if self.e_perp_im == Vector3::zeros()
            && *e_out_im == Vector3::zeros()
            && self.absorption == Vector3::zeros()
        {
            return self.new_from_e_perp(e_out_re);
        }
        let rot = Self::rotation_matrix_bivector(
            (self.e_perp, self.e_perp_im),
            (*e_out_re, *e_out_im),
            (self.n_app * self.prop, self.absorption),
        );
        let mut field = self.clone();
        field.e_perp = *e_out_re;
        field.e_perp_im = *e_out_im;
        field.matmul(&rot);
        field
    }

    /// Bilinear (unconjugated) dot product of two bivectors.
    fn cdot(
        a: &(Vector3<f32>, Vector3<f32>),
        b: &(Vector3<f32>, Vector3<f32>),
    ) -> Complex<f32> {
        Complex::new(
            a.0.dot(&b.0) - a.1.dot(&b.1),
            a.0.dot(&b.1) + a.1.dot(&b.0),
        )
    }

    /// Cross product of two bivectors.
    fn ccross(
        a: &(Vector3<f32>, Vector3<f32>),
        b: &(Vector3<f32>, Vector3<f32>),
    ) -> (Vector3<f32>, Vector3<f32>) {
        (
            a.0.cross(&b.0) - a.1.cross(&b.1),
            a.0.cross(&b.1) + a.1.cross(&b.0),
        )
    }

    /// Scales a bivector by a complex scalar.
    fn cscale(s: Complex<f32>, v: &(Vector3<f32>, Vector3<f32>)) -> (Vector3<f32>, Vector3<f32>) {
        (s.re * v.0 - s.im * v.1, s.re * v.1 + s.im * v.0)
    }

    /// The 2x2 basis-rotation matrix between two bivector bases sharing
    /// the reduced wave vector `w`, the complex generalisation of
    /// [`Self::rotation_matrix`]. All products are bilinear (no
    /// conjugation), per the bivector algebra of Dupertuis paper I.
    /// The second basis vector is `w x h1`, bilinearly normalised
    /// (which absorbs the `1/n_tilde` of Dupertuis eq. 13 since
    /// `(w x h1).(w x h1) = n_tilde^2` for a bilinear-unit `h1`).
    pub fn rotation_matrix_bivector(
        e_in: (Vector3<f32>, Vector3<f32>),
        e_out: (Vector3<f32>, Vector3<f32>),
        w: (Vector3<f32>, Vector3<f32>),
    ) -> Matrix2<Complex<f32>> {
        let evo2_raw = Self::ccross(&w, &e_in);
        let evo2_norm = Self::cdot(&evo2_raw, &evo2_raw).sqrt();
        let evo2 = Self::cscale(evo2_norm.inv(), &evo2_raw);

        let dot1 = Self::cdot(&e_out, &e_in);
        let dot2 = Self::cdot(&e_out, &evo2);

        let result = Matrix2::new(dot1, dot2, -dot2, dot1);
        let det = dot1 * dot1 + dot2 * dot2;
        result / det.sqrt()
    }

    /// Returns a new Field with prop and e_perp rotated by the given 3x3 rotation matrix.
    /// The amplitude and phase remain unchanged since they represent the relationship
    /// between polarization components, which is preserved under coordinate rotation.
    pub fn rotated(&self, rot: &Matrix3<f32>) -> Self {
        let new_prop = (rot * self.prop).normalize();
        let new_e_perp = (rot * self.e_perp).normalize();

        Self {
            prop: new_prop,
            e_perp: new_e_perp,
            ampl: self.ampl,
            phase: self.phase,
            n_app: self.n_app,
            absorption: rot * self.absorption,
            e_perp_im: rot * self.e_perp_im,
        }
    }

    /// Creates a new unit electric field with the given input perpendicular
    /// and propagation vectors.
    pub fn new_identity(e_perp: Vector3<f32>, prop: Vector3<f32>) -> Result<Self> {
        #[cfg(debug_assertions)]
        {
            let norm_e_perp_diff = e_perp.norm() - 1.0;
            if norm_e_perp_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!("e-perp is not normalised: {:?}", e_perp));
            }

            let norm_prop_diff = prop.norm() - 1.0;
            if norm_prop_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                    "propagation vector is not normalised: {:?}",
                    prop
                ));
            }

            let dot_product = e_perp.dot(&prop);
            if dot_product.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                "e-perp and propagation vector are not perpendicular, e_perp is: {:?}, prop is: {:?}, dot product is: {:?}",
                e_perp,
                prop,
                dot_product
            ));
            }
        }

        let field = Self {
            ampl: Matrix2::identity(),
            e_perp,
            phase: 0.0,
            prop,
            // Homogeneous vacuum defaults; Beam::new re-syncs the wave
            // state to its medium via sync_wave_homogeneous.
            n_app: 1.0,
            absorption: Vector3::zeros(),
            e_perp_im: Vector3::zeros(),
        };

        Ok(field)
    }

    /// Creates an electric field with the given input perpendicular field
    /// vector, propagation vector, and amplitude matrix.
    pub fn new(e_perp: Vector3<f32>, prop: Vector3<f32>, ampl0: Ampl, phase: f32) -> Result<Self> {
        #[cfg(debug_assertions)]
        {
            let norm_e_perp_diff = e_perp.norm() - 1.0;
            if norm_e_perp_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!("e-perp is not normalised: {:?}", e_perp));
            }

            let norm_prop_diff = prop.norm() - 1.0;
            if norm_prop_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                    "propagation vector is not normalised: {:?}",
                    prop
                ));
            }

            let dot_product = e_perp.dot(&prop);
            if dot_product.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                "e-perp and propagation vector are not perpendicular, e_perp is: {:?}, prop is: {:?}, dot product is: {:?}",
                e_perp,
                prop,
                dot_product
            ));
            }
        }
        let field = Self {
            ampl: ampl0,
            e_perp,
            phase,
            prop,
            n_app: 1.0,
            absorption: Vector3::zeros(),
            e_perp_im: Vector3::zeros(),
        };

        Ok(field)
    }

    /// Returns the 2x2 rotation matrix for rotating an amplitude matrix
    /// about the propagation vector `prop` from
    /// the plane perpendicular to `e_perp_in` to the plane perpendicular to
    /// `e_perp_out`.
    pub fn rotation_matrix<T: RealField + std::marker::Copy>(
        e_perp_in: Vector3<T>,
        e_perp_out: Vector3<T>,
        prop: Vector3<T>,
    ) -> Matrix2<T> {
        let dot1 = e_perp_out.dot(&e_perp_in);
        let evo2 = prop.cross(&e_perp_in).normalize();
        let dot2 = e_perp_out.dot(&evo2);

        let result = Matrix2::new(dot1, dot2.clone(), -dot2, dot1.clone());
        let det = result.determinant();

        result / det.abs().sqrt()
    }

    /// Returns the field intensity.
    pub fn intensity(&self) -> f32 {
        0.5 * self.ampl.norm_squared()
    }

    pub fn ampl_intensity(ampl: &Matrix2<Complex<f32>>) -> f32 {
        0.5 * ampl.norm_squared()
    }
}
