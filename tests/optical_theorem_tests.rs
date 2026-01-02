//! Tests for optical theorem validation.
//!
//! These tests verify that for a simple plane geometry with GO disabled (max_rec=0, max_tir=0):
//! 1. The optical theorem extinction is ~2x the integrated extinction (diffraction only)
//! 2. The optical theorem extinction is ~2x the input power
//! 3. Split and unsplit plane geometries give identical results
//! 4. Extinction cross section is wavelength-independent (for pure diffraction)
//!
//! Additional tests verify the forward scattering amplitude matrix properties:
//! 5. The forward amplitude matrix is diagonal (off-diagonal elements < 1% of diagonal)
//! 6. The diagonal elements have +90 degree phase (within 1 degree)

use goad::{
    bins::Scheme,
    multiproblem::MultiProblem,
    orientation::{Euler, Orientation, Scheme as OrientationScheme},
    params::Param,
    result::GOComponent,
    settings,
    zones::{ZoneConfig, ZoneType},
};
use nalgebra::{Complex, Matrix2};
use num_complex::Complex32;
use std::f32::consts::PI;

mod helpers;

/// Helper struct to hold extinction results from a simulation run
struct ExtinctionResults {
    input_power: f32,
    ext_cross_integrated: f32,
    ext_cross_optical_theorem: f32,
}

/// Run a simulation with the given geometry and wavelength, returning extinction results
fn run_plane_simulation(geom_name: &str, wavelength: f32) -> ExtinctionResults {
    let mut settings = settings::load_default_config().unwrap();

    // Set geometry and wavelength
    settings.geom_name = geom_name.to_string();
    settings.wavelength = wavelength;

    // Disable geometric optics (only external diffraction)
    settings.max_rec = 0;
    settings.max_tir = 0;

    // Set refractive indices
    settings.medium_refr_index = Complex32::new(1.0, 0.0);
    settings.particle_refr_index = vec![Complex32::new(1.31, 0.0)];

    // Set up zone scheme with fine resolution near forward direction
    settings.zones = vec![ZoneConfig::new(Scheme::Interval {
        thetas: vec![0.0, 0.001, 0.01, 0.1, 1.0, 180.0],
        theta_spacings: vec![0.0001, 0.001, 0.01, 0.1, 1.0],
        phis: vec![0.0, 360.0],
        phi_spacings: vec![2.0],
    })];

    // Single orientation
    settings.orientation = Orientation {
        scheme: OrientationScheme::Discrete {
            eulers: vec![Euler::new(0.0, 10.0, 0.0)],
        },
        euler_convention: goad::orientation::EulerConvention::ZYZ,
    };

    // Run simulation
    let mut multiproblem =
        MultiProblem::new(None, Some(settings)).expect("Failed to create MultiProblem");
    multiproblem.solve();

    // Extract results
    let input_power = multiproblem.result.powers.input;

    let full_zone = multiproblem
        .result
        .zones
        .iter()
        .find(|z| z.zone_type == ZoneType::Full)
        .expect("No full zone found");

    let forward_zone = multiproblem
        .result
        .zones
        .iter()
        .find(|z| z.zone_type == ZoneType::Forward)
        .expect("No forward zone found");

    let ext_cross_integrated = full_zone
        .params
        .get(&Param::ExtCross, &GOComponent::Total)
        .expect("No ExtCross_Total in full zone");

    let ext_cross_optical_theorem = forward_zone
        .params
        .get(&Param::ExtCrossOpticalTheorem, &GOComponent::Total)
        .expect("No ExtCrossOpticalTheorem_Total in forward zone");

    ExtinctionResults {
        input_power,
        ext_cross_integrated,
        ext_cross_optical_theorem,
    }
}

#[test]
fn test_optical_theorem_plane_extinction() {
    // Run simulations for three cases
    let results_plane = run_plane_simulation("examples/data/plane_xy_rect.obj", 0.532);
    let results_split = run_plane_simulation("examples/data/plane_xy_rect_split.obj", 0.532);
    let results_half_wl = run_plane_simulation("examples/data/plane_xy_rect_split.obj", 0.266);

    // === Check 1: Ratio of optical theorem to integrated extinction ≈ 2.0 (within 10%) ===
    let ratio_plane = results_plane.ext_cross_optical_theorem / results_plane.ext_cross_integrated;
    let ratio_split = results_split.ext_cross_optical_theorem / results_split.ext_cross_integrated;
    let ratio_half_wl =
        results_half_wl.ext_cross_optical_theorem / results_half_wl.ext_cross_integrated;

    assert!(
        (ratio_plane - 2.0).abs() / 2.0 < 0.10,
        "Plane: optical theorem / integrated ratio should be ~2.0, got {:.3}",
        ratio_plane
    );
    assert!(
        (ratio_split - 2.0).abs() / 2.0 < 0.10,
        "Split plane: optical theorem / integrated ratio should be ~2.0, got {:.3}",
        ratio_split
    );
    assert!(
        (ratio_half_wl - 2.0).abs() / 2.0 < 0.10,
        "Half wavelength: optical theorem / integrated ratio should be ~2.0, got {:.3}",
        ratio_half_wl
    );

    // === Check 2: Optical theorem extinction ≈ 2x input power (within 10%) ===
    let two_input_plane = 2.0 * results_plane.input_power;
    let two_input_split = 2.0 * results_split.input_power;
    let two_input_half_wl = 2.0 * results_half_wl.input_power;

    assert!(
        (results_plane.ext_cross_optical_theorem - two_input_plane).abs() / two_input_plane < 0.10,
        "Plane: optical theorem ext should be ~2x input power. Got {:.3}, expected {:.3}",
        results_plane.ext_cross_optical_theorem,
        two_input_plane
    );
    assert!(
        (results_split.ext_cross_optical_theorem - two_input_split).abs() / two_input_split < 0.10,
        "Split: optical theorem ext should be ~2x input power. Got {:.3}, expected {:.3}",
        results_split.ext_cross_optical_theorem,
        two_input_split
    );
    assert!(
        (results_half_wl.ext_cross_optical_theorem - two_input_half_wl).abs() / two_input_half_wl
            < 0.10,
        "Half wavelength: optical theorem ext should be ~2x input power. Got {:.3}, expected {:.3}",
        results_half_wl.ext_cross_optical_theorem,
        two_input_half_wl
    );

    // === Check 3: Split vs unsplit should match within 0.1% ===
    let integrated_diff = (results_plane.ext_cross_integrated - results_split.ext_cross_integrated)
        .abs()
        / results_plane.ext_cross_integrated;
    let optical_theorem_diff =
        (results_plane.ext_cross_optical_theorem - results_split.ext_cross_optical_theorem).abs()
            / results_plane.ext_cross_optical_theorem;

    assert!(
        integrated_diff < 0.001,
        "Split vs unsplit integrated extinction should match within 0.1%. Diff: {:.4}%",
        integrated_diff * 100.0
    );
    assert!(
        optical_theorem_diff < 0.001,
        "Split vs unsplit optical theorem extinction should match within 0.1%. Diff: {:.4}%",
        optical_theorem_diff * 100.0
    );

    // === Check 4: Half wavelength should match original within 10% ===
    let wl_integrated_diff =
        (results_split.ext_cross_integrated - results_half_wl.ext_cross_integrated).abs()
            / results_split.ext_cross_integrated;
    let wl_optical_theorem_diff =
        (results_split.ext_cross_optical_theorem - results_half_wl.ext_cross_optical_theorem).abs()
            / results_split.ext_cross_optical_theorem;

    assert!(
        wl_integrated_diff < 0.10,
        "Half wavelength integrated extinction should match original within 10%. Diff: {:.2}%",
        wl_integrated_diff * 100.0
    );
    assert!(
        wl_optical_theorem_diff < 0.10,
        "Half wavelength optical theorem extinction should match original within 10%. Diff: {:.2}%",
        wl_optical_theorem_diff * 100.0
    );
}

/// Run a simulation and return the forward scattering amplitude matrix
fn get_forward_amplitude(geom_name: &str) -> Matrix2<Complex<f32>> {
    let mut settings = settings::load_default_config().unwrap();

    // Set geometry
    settings.geom_name = geom_name.to_string();
    settings.wavelength = 0.532;

    // Disable geometric optics (only external diffraction)
    settings.max_rec = 0;
    settings.max_tir = 0;

    // Set refractive indices
    settings.medium_refr_index = Complex32::new(1.0, 0.0);
    settings.particle_refr_index = vec![Complex32::new(1.31, 0.0)];

    // Set up zone scheme with fine resolution near forward direction
    settings.zones = vec![ZoneConfig::new(Scheme::Interval {
        thetas: vec![0.0, 0.001, 0.01, 0.1, 1.0, 180.0],
        theta_spacings: vec![0.0001, 0.001, 0.01, 0.1, 1.0],
        phis: vec![0.0, 360.0],
        phi_spacings: vec![2.0],
    })];

    // Single orientation (same as extinction test)
    settings.orientation = Orientation {
        scheme: OrientationScheme::Discrete {
            eulers: vec![Euler::new(0.0, 10.0, 0.0)],
        },
        euler_convention: goad::orientation::EulerConvention::ZYZ,
    };

    // Run simulation
    let mut multiproblem =
        MultiProblem::new(None, Some(settings)).expect("Failed to create MultiProblem");
    multiproblem.solve();

    // Get forward zone amplitude
    let forward_zone = multiproblem
        .result
        .zones
        .iter()
        .find(|z| z.zone_type == ZoneType::Forward)
        .expect("No forward zone found");

    forward_zone.field_2d[0].ampl_total
}

/// Check that amplitude matrix is diagonal (off-diagonal < 1% of diagonal magnitude)
fn check_diagonal(ampl: &Matrix2<Complex<f32>>, name: &str) {
    let diag_mag = (ampl[(0, 0)].norm() + ampl[(1, 1)].norm()) / 2.0;
    let off_diag_01 = ampl[(0, 1)].norm();
    let off_diag_10 = ampl[(1, 0)].norm();

    let ratio_01 = off_diag_01 / diag_mag;
    let ratio_10 = off_diag_10 / diag_mag;

    assert!(
        ratio_01 < 0.01,
        "{}: off-diagonal (0,1) should be < 1% of diagonal. Got {:.4}%",
        name,
        ratio_01 * 100.0
    );
    assert!(
        ratio_10 < 0.01,
        "{}: off-diagonal (1,0) should be < 1% of diagonal. Got {:.4}%",
        name,
        ratio_10 * 100.0
    );
}

/// Check that diagonal elements have +90 degree phase (within 1 degree)
fn check_phase_90(ampl: &Matrix2<Complex<f32>>, name: &str) {
    let phase_00 = ampl[(0, 0)].arg() * 180.0 / PI;
    let phase_11 = ampl[(1, 1)].arg() * 180.0 / PI;

    // Phase should be +90 degrees (or -270, which is equivalent)
    // Normalize to [-180, 180] range
    let phase_diff_00 = ((phase_00 - 90.0 + 180.0) % 360.0) - 180.0;
    let phase_diff_11 = ((phase_11 - 90.0 + 180.0) % 360.0) - 180.0;

    assert!(
        phase_diff_00.abs() < 1.0,
        "{}: S2 (0,0) phase should be +90 deg. Got {:.2} deg (diff: {:.2} deg)",
        name,
        phase_00,
        phase_diff_00
    );
    assert!(
        phase_diff_11.abs() < 1.0,
        "{}: S1 (1,1) phase should be +90 deg. Got {:.2} deg (diff: {:.2} deg)",
        name,
        phase_11,
        phase_diff_11
    );
}

#[test]
fn test_forward_amplitude_matrix_properties() {
    // Test with plane geometry
    let ampl_plane = get_forward_amplitude("examples/data/plane_xy_rect.obj");
    check_diagonal(&ampl_plane, "plane_xy_rect");
    check_phase_90(&ampl_plane, "plane_xy_rect");

    // Test with hex geometry
    let ampl_hex = get_forward_amplitude("examples/data/hex.obj");
    check_diagonal(&ampl_hex, "hex");
    check_phase_90(&ampl_hex, "hex");

    // Test with concave geometries
    let ampl_concave1 = get_forward_amplitude("examples/data/concave1.obj");
    check_diagonal(&ampl_concave1, "concave1");
    check_phase_90(&ampl_concave1, "concave1");

    let ampl_concave2 = get_forward_amplitude("examples/data/concave2.obj");
    check_diagonal(&ampl_concave2, "concave2");
    check_phase_90(&ampl_concave2, "concave2");
}
