//! Tests for aperture diffraction polarization preservation.
//!
//! These tests verify that for forward scattering (theta ≈ 0), the diffracted beam:
//! 1. Has approximately identity amplitude matrix (scaled by Fraunhofer factor)
//! 2. Preserves the physical electric field polarization direction

use nalgebra::{Complex, Matrix2, Point3, Vector3};
use std::f32::consts::PI;

use goad::bins::{AngleBin, SolidAngleBin};
use goad::diff::n2f_aperture_diffraction;

const TOL: f32 = 1e-3;

// =============================================================================
// Helper functions
// =============================================================================

/// Creates a square aperture in the xy plane centered at origin
fn square_aperture() -> Vec<Point3<f32>> {
    vec![
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 0.0),
        Point3::new(1.0, -1.0, 0.0),
    ]
}

/// Creates a forward scattering bin (theta=0, phi=0)
fn forward_bin() -> SolidAngleBin {
    // Small bin centered at theta=0 (forward direction)
    let theta_bin = AngleBin::new(0.0, 1.0);
    let phi_bin = AngleBin::new(0.0, 1.0);
    SolidAngleBin::new(theta_bin, phi_bin)
}

/// Check if amplitude matrix is proportional to identity (diagonal with equal elements)
fn is_proportional_to_identity(ampl: &Matrix2<Complex<f32>>, tol: f32) -> bool {
    // Off-diagonal elements should be near zero
    let off_diag_small = ampl[(0, 1)].norm() < tol && ampl[(1, 0)].norm() < tol;

    // Diagonal elements should be approximately equal
    let diag_ratio = if ampl[(1, 1)].norm() > 1e-10 {
        (ampl[(0, 0)] / ampl[(1, 1)] - Complex::new(1.0, 0.0)).norm()
    } else {
        0.0
    };

    off_diag_small && diag_ratio < tol
}

// =============================================================================
// Tests: Forward diffraction preserves polarization for various e_perp directions
// =============================================================================

#[test]
fn test_forward_diffraction_e_perp_along_x() {
    // e_perp along +x, propagation along +z
    let e_perp = Vector3::new(1.0, 0.0, 0.0);
    let prop = Vector3::new(0.0, 0.0, 1.0);

    let verts = square_aperture();
    let ampl = Matrix2::<Complex<f32>>::identity();
    let bins = vec![forward_bin()];
    let wavenumber = 2.0 * PI;

    let result = n2f_aperture_diffraction(&verts, ampl, prop, e_perp, &bins, wavenumber, None);

    assert!(!result.is_empty(), "Should produce output for forward bin");
    let forward_ampl = &result[0];

    // Check amplitude is proportional to identity
    assert!(
        is_proportional_to_identity(forward_ampl, TOL),
        "Forward amplitude should be proportional to identity, got: {:?}",
        forward_ampl
    );
}

#[test]
fn test_forward_diffraction_e_perp_along_y() {
    // e_perp along +y, propagation along +z
    let e_perp = Vector3::new(0.0, 1.0, 0.0);
    let prop = Vector3::new(0.0, 0.0, 1.0);

    let verts = square_aperture();
    let ampl = Matrix2::<Complex<f32>>::identity();
    let bins = vec![forward_bin()];
    let wavenumber = 2.0 * PI;

    let result = n2f_aperture_diffraction(&verts, ampl, prop, e_perp, &bins, wavenumber, None);

    assert!(!result.is_empty(), "Should produce output for forward bin");
    let forward_ampl = &result[0];

    assert!(
        is_proportional_to_identity(forward_ampl, TOL),
        "Forward amplitude should be proportional to identity, got: {:?}",
        forward_ampl
    );
}

#[test]
fn test_forward_diffraction_e_perp_diagonal_xy() {
    // e_perp along diagonal in xy plane, propagation along +z
    let e_perp = Vector3::new(1.0, 1.0, 0.0).normalize();
    let prop = Vector3::new(0.0, 0.0, 1.0);

    let verts = square_aperture();
    let ampl = Matrix2::<Complex<f32>>::identity();
    let bins = vec![forward_bin()];
    let wavenumber = 2.0 * PI;

    let result = n2f_aperture_diffraction(&verts, ampl, prop, e_perp, &bins, wavenumber, None);

    assert!(!result.is_empty(), "Should produce output for forward bin");
    let forward_ampl = &result[0];

    assert!(
        is_proportional_to_identity(forward_ampl, TOL),
        "Forward amplitude should be proportional to identity, got: {:?}",
        forward_ampl
    );
}

#[test]
fn test_forward_diffraction_e_perp_negative_x() {
    // e_perp along -x, propagation along +z
    let e_perp = Vector3::new(-1.0, 0.0, 0.0);
    let prop = Vector3::new(0.0, 0.0, 1.0);

    let verts = square_aperture();
    let ampl = Matrix2::<Complex<f32>>::identity();
    let bins = vec![forward_bin()];
    let wavenumber = 2.0 * PI;

    let result = n2f_aperture_diffraction(&verts, ampl, prop, e_perp, &bins, wavenumber, None);

    assert!(!result.is_empty(), "Should produce output for forward bin");
    let forward_ampl = &result[0];

    assert!(
        is_proportional_to_identity(forward_ampl, TOL),
        "Forward amplitude should be proportional to identity, got: {:?}",
        forward_ampl
    );
}

#[test]
fn test_forward_diffraction_prop_along_negative_z() {
    // e_perp along +x, propagation along -z (incoming from +z side)
    let e_perp = Vector3::new(1.0, 0.0, 0.0);
    let prop = Vector3::new(0.0, 0.0, -1.0);

    let verts = square_aperture();
    let ampl = Matrix2::<Complex<f32>>::identity();

    // For prop along -z, "forward" is theta=180 degrees
    let theta_bin = AngleBin::new(179.0, 180.0);
    let phi_bin = AngleBin::new(0.0, 1.0);
    let backward_bin = SolidAngleBin::new(theta_bin, phi_bin);
    let bins = vec![backward_bin];

    let wavenumber = 2.0 * PI;

    let result = n2f_aperture_diffraction(&verts, ampl, prop, e_perp, &bins, wavenumber, None);

    assert!(!result.is_empty(), "Should produce output for forward bin");
    let forward_ampl = &result[0];

    assert!(
        is_proportional_to_identity(forward_ampl, TOL),
        "Forward amplitude should be proportional to identity, got: {:?}",
        forward_ampl
    );
}

#[test]
fn test_forward_diffraction_oblique_prop() {
    // Oblique propagation: prop tilted 45 degrees from z in xz plane
    let prop = Vector3::new(1.0, 0.0, 1.0).normalize();
    // e_perp must be perpendicular to prop, choose y direction
    let e_perp = Vector3::new(0.0, 1.0, 0.0);

    let verts = square_aperture();
    let ampl = Matrix2::<Complex<f32>>::identity();

    // Forward scattering for this prop is at theta=45, phi=0
    let theta_bin = AngleBin::new(44.0, 46.0);
    let phi_bin = AngleBin::new(-1.0, 1.0);
    let forward_bin = SolidAngleBin::new(theta_bin, phi_bin);
    let bins = vec![forward_bin];

    let wavenumber = 2.0 * PI;

    let result = n2f_aperture_diffraction(&verts, ampl, prop, e_perp, &bins, wavenumber, None);

    assert!(!result.is_empty(), "Should produce output for forward bin");
    let forward_ampl = &result[0];

    assert!(
        is_proportional_to_identity(forward_ampl, TOL),
        "Forward amplitude should be proportional to identity, got: {:?}",
        forward_ampl
    );
}

// =============================================================================
// Tests: Verify different e_perp choices give consistent results
// =============================================================================

#[test]
fn test_different_e_perp_same_intensity() {
    // For unpolarized light, different e_perp choices should give same total intensity
    let prop = Vector3::new(0.0, 0.0, 1.0);
    let verts = square_aperture();
    let bins = vec![forward_bin()];
    let wavenumber = 2.0 * PI;

    let e_perp_x = Vector3::new(1.0, 0.0, 0.0);
    let e_perp_y = Vector3::new(0.0, 1.0, 0.0);
    let e_perp_diag = Vector3::new(1.0, 1.0, 0.0).normalize();

    let ampl = Matrix2::<Complex<f32>>::identity();

    let result_x = n2f_aperture_diffraction(&verts, ampl, prop, e_perp_x, &bins, wavenumber, None);
    let result_y = n2f_aperture_diffraction(&verts, ampl, prop, e_perp_y, &bins, wavenumber, None);
    let result_diag =
        n2f_aperture_diffraction(&verts, ampl, prop, e_perp_diag, &bins, wavenumber, None);

    // Compute intensities (Frobenius norm squared of amplitude matrices)
    let intensity_x = result_x[0].norm_squared();
    let intensity_y = result_y[0].norm_squared();
    let intensity_diag = result_diag[0].norm_squared();

    // All should give the same intensity
    assert!(
        (intensity_x - intensity_y).abs() / intensity_x < TOL,
        "Intensity should be same for e_perp_x and e_perp_y: {} vs {}",
        intensity_x,
        intensity_y
    );
    assert!(
        (intensity_x - intensity_diag).abs() / intensity_x < TOL,
        "Intensity should be same for e_perp_x and e_perp_diag: {} vs {}",
        intensity_x,
        intensity_diag
    );
}
