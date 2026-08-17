use anyhow::Result;
use nalgebra::{Complex, Vector3};

/// Result of the reduced-wave-vector interface transfer.
///
/// The transmitted wave is `w2 = w_t + kappa_t * normal` with `w_t` the
/// (complex) tangential part of the incident reduced wave vector,
/// conserved across the interface. `kappa_i` and `kappa_t` are the
/// complex normal components that feed the generalised Fresnel
/// coefficients. See docs/theory/lossy_media_beam_tracing.typ.
pub struct Transfer {
    /// Real part of the transmitted reduced wave vector, `N2 * e2`.
    pub w_re: Vector3<f32>,
    /// Imaginary part of the transmitted reduced wave vector, the
    /// absorption vector `A2`.
    pub w_im: Vector3<f32>,
    /// Complex normal component of the transmitted wave (`gamma`).
    pub kappa_t: Complex<f32>,
    /// Complex normal component of the incident wave.
    pub kappa_i: Complex<f32>,
}

/// Interface transfer of the reduced wave vector `w = N e + iA`:
/// tangential continuity plus one complex square root,
/// `gamma = sqrt(n2^2 - w_t . w_t)`, taking the continuation-rule
/// branch (`re + im > 0`), which is the analytic continuation from the
/// beam's real-angle spectrum and is continuous with the lossless
/// limit on both sides of the critical angle (spectrally validated in
/// examples/lossy_media_validation/spectral_branch_check.py).
///
/// `normal` must point from the incident medium into medium 2 (along
/// transmission). This subsumes Snell's law: the transmitted direction
/// is `w_re` normalised, and the Macke effective-index construction in
/// [`get_theta_t`] is its homogeneous-incidence special case (verified
/// to 1e-15 in the validation harness).
pub fn transfer(
    w_re: Vector3<f32>,
    w_im: Vector3<f32>,
    normal: Vector3<f32>,
    n2: Complex<f32>,
) -> Transfer {
    let ki_re = w_re.dot(&normal);
    let ki_im = w_im.dot(&normal);
    let wt_re = w_re - ki_re * normal;
    let wt_im = w_im - ki_im * normal;

    let wt_sq = Complex::new(
        wt_re.dot(&wt_re) - wt_im.dot(&wt_im),
        2.0 * wt_re.dot(&wt_im),
    );
    let g = (n2 * n2 - wt_sq).sqrt();
    let gamma = if g.re + g.im < 0.0 { -g } else { g };

    Transfer {
        w_re: wt_re + gamma.re * normal,
        w_im: wt_im + gamma.im * normal,
        kappa_t: gamma,
        kappa_i: Complex::new(ki_re, ki_im),
    }
}

#[cfg(test)]
mod tests {

    use nalgebra::Complex;

    use super::*;
    use std::f32::consts::PI;

    fn homogeneous_w(n: Complex<f32>, prop: Vector3<f32>) -> (Vector3<f32>, Vector3<f32>) {
        (n.re * prop, n.im * prop)
    }

    /// The transfer must agree with the Macke effective-index Snell
    /// construction for homogeneous incidence (its special case;
    /// verified to 1e-15 in f64 by the validation harness).
    #[test]
    fn transfer_matches_macke_for_homogeneous_incidence() {
        let normal = Vector3::z();
        for (n1, n2) in [
            (Complex::new(1.0, 0.0), Complex::new(1.31, 0.1)),
            (Complex::new(1.0, 0.0), Complex::new(1.5, 0.5)),
            (Complex::new(1.31, 0.01), Complex::new(1.5, 0.2)),
        ] {
            for theta_deg in [10.0_f32, 30.0, 50.0, 70.0] {
                let th = theta_deg.to_radians();
                let prop = Vector3::new(th.sin(), 0.0, th.cos());
                let (w_re, w_im) = homogeneous_w(n1, prop);
                let tr = transfer(w_re, w_im, normal, n2);
                let sin_t = (tr.w_re - tr.w_re.dot(&normal) * normal).norm() / tr.w_re.norm();
                let sin_t_macke = get_theta_t(th, n1, n2).unwrap().sin();
                assert!(
                    (sin_t - sin_t_macke).abs() < 2e-6,
                    "{n1} -> {n2} at {theta_deg}: {sin_t} vs {sin_t_macke}"
                );
            }
        }
    }

    /// Invariants of the transfer: tangential continuity, the
    /// dispersion relation in medium 2, and the structure of the
    /// transmitted absorption vector (normal at first entry from a
    /// lossless host, transverse on exit into one).
    #[test]
    fn transfer_invariants() {
        let normal = Vector3::z();
        let th = 40.0_f32.to_radians();
        let prop = Vector3::new(th.sin(), 0.0, th.cos());

        // first entry: vacuum -> absorbing
        let n2 = Complex::new(1.31, 0.1);
        let (w_re, w_im) = homogeneous_w(Complex::new(1.0, 0.0), prop);
        let tr = transfer(w_re, w_im, normal, n2);
        let tang = tr.w_im - tr.w_im.dot(&normal) * normal;
        assert!(tang.norm() < 1e-7, "first-entry A must be normal to the facet");
        let ww_re = tr.w_re.dot(&tr.w_re) - tr.w_im.dot(&tr.w_im);
        let ww_im = 2.0 * tr.w_re.dot(&tr.w_im);
        let n2sq = n2 * n2;
        assert!((ww_re - n2sq.re).abs() < 1e-4 && (ww_im - n2sq.im).abs() < 1e-4);

        // exit: the interior inhomogeneous wave into vacuum
        let tr2 = transfer(tr.w_re, tr.w_im, normal, Complex::new(1.0, 0.0));
        let cos_ae = tr2.w_im.dot(&tr2.w_re);
        assert!(
            cos_ae.abs() < 1e-5,
            "exit-ramp absorption must be transverse, got A.we = {cos_ae}"
        );
    }

    /// Golden value from the validation harness: homogeneous ice-like
    /// wave (1.31 + 0.01i) at 60 degrees onto vacuum. The continuation
    /// branch gives gamma = -0.018329 + 0.536037i (the principal root
    /// is its negation, the wrong branch above critical).
    #[test]
    fn transfer_continuation_branch_above_critical() {
        let normal = Vector3::z();
        let th = 60.0_f32.to_radians();
        let prop = Vector3::new(th.sin(), 0.0, th.cos());
        let (w_re, w_im) = homogeneous_w(Complex::new(1.31, 0.01), prop);
        let tr = transfer(w_re, w_im, normal, Complex::new(1.0, 0.0));
        assert!((tr.kappa_t.re - -0.018329).abs() < 1e-4, "{}", tr.kappa_t);
        assert!((tr.kappa_t.im - 0.536037).abs() < 1e-4, "{}", tr.kappa_t);
        assert!((tr.kappa_i.re - 0.655).abs() < 1e-5);
        assert!((tr.kappa_i.im - 0.005).abs() < 1e-6);
    }

    #[test]
    fn normal_incidence_same_media() {
        let theta_i = 0.0;
        let m1 = Complex::new(1.0, 0.0);
        let m2 = m1;
        let theta_t = get_theta_t(theta_i, m1, m2).unwrap();
        assert!(theta_i - theta_t < 0.01)
    }

    #[test]
    fn normal_incidence() {
        let theta_i = 0.0;
        let m1 = Complex::new(1.0, 0.0);
        let m2 = Complex::new(1.31, 0.0);
        let theta_t = get_theta_t(theta_i, m1, m2).unwrap();
        let abs_difference = (theta_i - theta_t).abs();
        assert!(abs_difference < f32::EPSILON)
    }

    #[test]
    fn angle30_incidence() {
        let theta_i = 30.0 * PI / 180.0;
        let m1 = Complex::new(1.0, 0.0);
        let m2 = Complex::new(1.31, 0.0);
        let theta_t = get_theta_t(theta_i, m1, m2).unwrap();
        let abs_difference = (theta_t - 0.3916126).abs();
        assert!(abs_difference < 0.001)
    }

    #[test]
    fn absorbing_test() {
        let theta_i = 1.17773;
        let m1 = Complex::new(1.0, 0.0);
        let m2 = Complex::new(1.5, 0.1);
        let theta_t = get_theta_t(theta_i, m1, m2).unwrap();
        let abs_difference = (theta_t - 0.662387).abs();
        assert!(abs_difference < 0.001)
    }
}

/// Returns the sine of the transmitted angle according to Snell's Law.
/// Port from Fortran code rt_c.f90, Macke 1996.
/// All angles are in radians.
pub fn get_theta_t(theta_i: f32, m1: Complex<f32>, m2: Complex<f32>) -> Result<f32> {
    if m1 == m2 {
        return Ok(theta_i);
    }

    let k1 = m1.im / m1.re; // imag(inc) / real(inc)
    let k2 = m2.im / m2.re; // imag(trans) / real(trans)
    let krel = (k2 - k1) / (1.0 + k1 * k2);
    let nrel = m2.re / m1.re * (1.0 + k1 * k2) / (1.0 + k1 * k1);

    let ref1 = nrel * nrel;
    let ref2 = krel * krel;
    let ref3 = (1.0 + ref2) * (1.0 + ref2);
    let ref6 = ref1 * ref3 / ((1.0 + krel * k2) * (1.0 + krel * k2));

    let sintiq = (theta_i).sin().powi(2);
    let ref4 = 1.0 - (1.0 - ref2) / ref1 / ref3 * sintiq;
    let ref5 = 2.0 * krel / ref1 / ref3 * sintiq;

    let q4 = ref4 * ref4 + ref5 * ref5;
    let q2 = q4.sqrt();

    let test1 = (ref4 / q2).acos() / 2.0;

    let g = test1;

    let ref7 = (g.cos() - k2 * g.sin()) * (g.cos() - k2 * g.sin());
    let rnstar = (sintiq + ref6 * q2 * ref7).sqrt();

    let theta_t = (theta_i.sin() / rnstar).asin();

    if theta_t.is_nan() {
        Err(anyhow::anyhow!("theta_t is NaN"))
    } else {
        Ok(theta_t)
    }
}
