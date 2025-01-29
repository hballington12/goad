use nalgebra::{Complex, Matrix2, Vector2};

/// Returns the matrix representation of the Fresnel equations for reflection.
///
/// # Parameters
/// - `n1`: Refractive index of the first medium (complex).
/// - `n2`: Refractive index of the second medium (complex).
/// - `theta_i`: Incident angle (in radians).
/// - `theta_t`: Transmission angle (in radians).
///
/// # Returns
/// A 2x2 diagonal matrix representing the Fresnel reflection coefficients.
pub fn refl(
    n1: Complex<f32>,
    n2: Complex<f32>,
    theta_i: f32,
    theta_t: f32,
) -> Matrix2<Complex<f32>> {
    let cti = theta_i.cos();
    let ctt = theta_t.cos();
    let f11 = (n2 * cti - n1 * ctt) / (n1 * ctt + n2 * cti);
    let f22 = (n1 * cti - n2 * ctt) / (n1 * cti + n2 * ctt);
    Matrix2::from_diagonal(&Vector2::new(f11, f22))
}

/// Returns the matrix representation of the Fresnel equations for refraction.
///
/// # Parameters
/// - `n1`: Refractive index of the first medium (complex).
/// - `n2`: Refractive index of the second medium (complex).
/// - `theta_i`: Incident angle (in radians).
/// - `theta_t`: Transmission angle (in radians).
///
/// # Returns
/// A 2x2 diagonal matrix representing the Fresnel transmission coefficients.
pub fn refr(
    n1: Complex<f32>,
    n2: Complex<f32>,
    theta_i: f32,
    theta_t: f32,
) -> Matrix2<Complex<f32>> {
    let cti = theta_i.cos();
    let ctt = theta_t.cos();
    let f11 = (2.0 * n1 * cti) / (n1 * ctt + n2 * cti);
    let f22 = (2.0 * n1 * cti) / (n1 * cti + n2 * ctt);
    Matrix2::from_diagonal(&Vector2::new(f11, f22))
}
