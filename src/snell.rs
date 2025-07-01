//! Generalized Snell's law for complex refractive indices.
//!
//! This module implements the generalized form of Snell's law that handles
//! complex refractive indices for absorbing and dispersive materials. This
//! extension is essential for electromagnetic scattering simulations involving
//! realistic materials with absorption and dispersion effects.
//!
//! The implementation provides:
//! - Generalized Snell's law for complex media
//! - Proper handling of absorption and dispersion
//! - Numerical stability for edge cases
//! - Validation against limiting cases
//! - Support for highly absorbing materials
//!
//! # Mathematical Foundation
//!
//! Based on the generalized Snell's law formulation from Macke (1996):
//! - Decomposition of complex refractive indices
//! - Proper treatment of phase and amplitude relationships
//! - Numerical algorithms for stable computation
//! - Validation of physical constraints

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

/// Computes transmitted angle using generalized Snell's law for complex refractive indices.
/// 
/// **Context**: Classical Snell's law applies to real refractive indices, but electromagnetic
/// scattering simulations often involve absorbing materials with complex refractive indices.
/// The generalized form accounts for both the real and imaginary parts of the refractive
/// index, ensuring proper treatment of absorption and dispersion effects.
/// 
/// **How it Works**: Implements the algorithm from Macke (1996) that handles complex
/// refractive indices by decomposing them into real and imaginary components, computing
/// relative parameters, and solving the generalized form of Snell's law. Includes
/// numerical safeguards against invalid angle calculations.
/// 
/// # Example
/// ```rust,no_run
/// // let theta_t = get_theta_t(theta_i, n1, n2)?;
/// // let prop = get_refraction_vector(&normal, &beam.prop, theta_i, theta_t);
/// ```
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
