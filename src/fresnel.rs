//! Fresnel equations for electromagnetic boundary conditions.
//!
//! This module implements the Fresnel equations that govern electromagnetic
//! wave reflection and transmission at material interfaces. These equations
//! provide the exact electromagnetic boundary conditions required for computing
//! field amplitudes in reflection and refraction processes.
//!
//! The Fresnel calculations provide:
//! - Reflection coefficients for s and p polarizations
//! - Transmission coefficients with impedance matching
//! - Complex refractive index support for absorbing materials
//! - Matrix representation for direct amplitude computation
//! - Proper treatment of polarization-dependent scattering
//!
//! # Physical Foundation
//!
//! Based on Maxwell's equations at material boundaries:
//! - Continuity of tangential electric field
//! - Continuity of tangential magnetic field
//! - Proper impedance relationships
//! - Conservation of electromagnetic energy

use nalgebra::{Complex, Matrix2, Vector2};

/// Computes Fresnel reflection coefficients for electromagnetic surface interactions.
/// 
/// **Context**: When electromagnetic waves encounter interfaces between materials
/// with different refractive indices, the reflected field amplitudes depend on
/// polarization, incident angle, and material properties. The Fresnel equations
/// provide the exact electromagnetic boundary conditions for these interactions.
/// 
/// **How it Works**: Calculates reflection coefficients separately for s-polarized
/// (perpendicular) and p-polarized (parallel) field components using the classic
/// Fresnel formulas. Returns a diagonal matrix with these coefficients for direct
/// multiplication with field amplitude matrices.
/// 
/// # Example
/// ```rust
/// let fresnel = fresnel::refl(n1, n2, theta_i, theta_t);
/// let refl_ampl = fresnel * ampl;
/// ```
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

/// Computes Fresnel transmission coefficients for electromagnetic surface interactions.
/// 
/// **Context**: Transmitted (refracted) electromagnetic fields at material interfaces
/// require different amplitude scaling than reflected fields. The transmission
/// coefficients account for impedance matching between media and ensure power
/// conservation at the interface.
/// 
/// **How it Works**: Applies Fresnel transmission formulas for both s-polarized
/// and p-polarized components. The coefficients account for both the change in
/// field amplitude and the impedance difference between the two media.
/// 
/// # Example
/// ```rust
/// let fresnel = fresnel::refr(n1, n2, theta_i, theta_t);
/// let refr_ampl = fresnel * ampl.clone();
/// ```
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
