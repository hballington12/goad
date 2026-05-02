use nalgebra::{Complex, Matrix2, Vector2};

/// Returns the diagonal 2x2 Fresnel reflection matrix for an interface between
/// media of refractive index `n1` and `n2` at incident angle `theta_i`.
///
/// `cos(theta_t)` is computed internally as a complex number from Snell's law,
/// so the same formula is correct both below and above the critical angle.
/// Above critical, the complex sqrt produces a non-zero imaginary part for
/// `cos(theta_t)`, which gives the s and p reflection coefficients different
/// complex phases — the s/p phase difference is what produces non-zero
/// `S34`/`S43` Mueller elements from a TIR event.
///
/// Diagonal entries:
/// - `[0, 0]` = `r_p` (parallel polarization)
/// - `[1, 1]` = `r_s` (perpendicular polarization)
pub fn refl(n1: Complex<f32>, n2: Complex<f32>, theta_i: f32) -> Matrix2<Complex<f32>> {
    let cti = Complex::from(theta_i.cos());
    let ctt = cos_theta_t(n1, n2, theta_i);
    let f11 = (n2 * cti - n1 * ctt) / (n1 * ctt + n2 * cti);
    let f22 = (n1 * cti - n2 * ctt) / (n1 * cti + n2 * ctt);
    Matrix2::from_diagonal(&Vector2::new(f11, f22))
}

/// Returns the diagonal 2x2 Fresnel transmission matrix. Only meaningful below
/// the critical angle (TIR has no transmitted ray); the formula uses a complex
/// `cos(theta_t)` so above-critical inputs do not produce NaNs.
pub fn refr(n1: Complex<f32>, n2: Complex<f32>, theta_i: f32) -> Matrix2<Complex<f32>> {
    let cti = Complex::from(theta_i.cos());
    let ctt = cos_theta_t(n1, n2, theta_i);
    let f11 = (Complex::from(2.0) * n1 * cti) / (n1 * ctt + n2 * cti);
    let f22 = (Complex::from(2.0) * n1 * cti) / (n1 * cti + n2 * ctt);
    Matrix2::from_diagonal(&Vector2::new(f11, f22))
}

/// `cos(theta_t)` from Snell's law, as a complex value.
fn cos_theta_t(n1: Complex<f32>, n2: Complex<f32>, theta_i: f32) -> Complex<f32> {
    let sti = theta_i.sin();
    let ratio = n1 / n2;
    let sin_sq_tt = ratio * ratio * Complex::from(sti * sti);
    let cos_sq_tt = Complex::from(1.0) - sin_sq_tt;
    cos_sq_tt.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: Complex<f32>, b: Complex<f32>, tol: f32) -> bool {
        (a.re - b.re).abs() < tol && (a.im - b.im).abs() < tol
    }

    #[test]
    fn sub_critical_matches_real_formula() {
        let n1 = Complex::new(1.0, 0.0);
        let n2 = Complex::new(1.31, 0.0);
        let theta_i = 30f32.to_radians();
        let m = refl(n1, n2, theta_i);
        assert!(m[(0, 0)].im.abs() < 1e-6);
        assert!(m[(1, 1)].im.abs() < 1e-6);
    }

    #[test]
    fn tir_unit_modulus() {
        // Ice to air at 60 deg internal — well above critical (~49.8 deg).
        // Energy conservation forces |r_s| = |r_p| = 1 under TIR.
        let n1 = Complex::new(1.31, 0.0);
        let n2 = Complex::new(1.0, 0.0);
        let theta_i = 60f32.to_radians();
        let m = refl(n1, n2, theta_i);
        assert!((m[(0, 0)].norm() - 1.0).abs() < 1e-5);
        assert!((m[(1, 1)].norm() - 1.0).abs() < 1e-5);
    }

    #[test]
    fn tir_phase_difference_nonzero() {
        // Above critical, r_p and r_s should have different complex phases.
        let n1 = Complex::new(1.31, 0.0);
        let n2 = Complex::new(1.0, 0.0);
        let theta_i = 60f32.to_radians();
        let m = refl(n1, n2, theta_i);
        let delta = m[(0, 0)].arg() - m[(1, 1)].arg();
        assert!(
            delta.abs() > 0.1,
            "s/p phase difference too small: {}",
            delta
        );
    }

    #[test]
    fn at_critical_angle_real() {
        // At the critical angle, cos(theta_t) = 0 so r_s = r_p = 1.
        let n1 = Complex::new(1.31, 0.0);
        let n2 = Complex::new(1.0, 0.0);
        let critical = (1.0_f32 / 1.31).asin();
        let m = refl(n1, n2, critical);
        assert!(approx(m[(0, 0)], Complex::new(1.0, 0.0), 1e-3));
        assert!(approx(m[(1, 1)], Complex::new(1.0, 0.0), 1e-3));
    }
}
