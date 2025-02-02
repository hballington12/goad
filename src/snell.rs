use anyhow::Result;
use nalgebra::Complex;

#[cfg(test)]
mod tests {

    use nalgebra::Complex;

    use super::*;
    use std::f32::consts::PI;

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
