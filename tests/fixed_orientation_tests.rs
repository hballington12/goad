use goad::{
    bins::{self, BinningScheme},
    multiproblem::MultiProblem,
    orientation::Euler,
    result::MuellerMatrix,
    settings,
};
use helpers::{compare_results, load_reference_mueller};
use num_complex::Complex32;

pub mod helpers;
// Tolerance for comparing Mueller matrix elements
const FRAC_TOL: f32 = 1e-4; // fractional error
const ABS_TOL: f32 = 1e4; // absolute error

#[test]
fn fixed_hex_30_30_30() {
    let mut settings = settings::load_default_config().unwrap();
    // Reduce binning for faster testing
    settings.binning = BinningScheme {
        scheme: bins::Scheme::new_simple(19, 19),
    };
    settings.orientation = goad::orientation::Orientation {
        scheme: goad::orientation::Scheme::Discrete {
            eulers: vec![Euler::new(30.0, 30.0, 30.0)],
        },
        euler_convention: goad::orientation::EulerConvention::ZYZ,
    };

    let mut multiproblem =
        MultiProblem::new(None, Some(settings)).expect("Failed to create MultiProblem");
    multiproblem.solve();

    let result = multiproblem
        .result
        .field_2d
        .iter()
        .map(|m| m.mueller_total.to_vec())
        .collect::<Vec<Vec<f32>>>();
    let reference = load_reference_mueller("fixed_hex_30_30_30_mueller_scatgrid").unwrap();
    compare_results(result, reference, FRAC_TOL, ABS_TOL).unwrap();
}

#[test]
fn fixed_hex_30_20_20() {
    let mut settings = settings::load_default_config().unwrap();
    // Reduce binning for faster testing
    settings.binning = BinningScheme {
        scheme: bins::Scheme::new_simple(19, 19),
    };
    settings.orientation = goad::orientation::Orientation {
        scheme: goad::orientation::Scheme::Discrete {
            eulers: vec![Euler::new(30.0, 20.0, 20.0)],
        },
        euler_convention: goad::orientation::EulerConvention::ZYZ,
    };
    // Change the refractive index
    settings.particle_refr_index = vec![Complex32::new(1.3117, 0.1)];

    let mut multiproblem =
        MultiProblem::new(None, Some(settings)).expect("Failed to create MultiProblem");
    multiproblem.solve();

    let result = multiproblem
        .result
        .field_2d
        .iter()
        .map(|m| m.mueller_total.to_vec())
        .collect::<Vec<Vec<f32>>>();
    let reference = load_reference_mueller("fixed_hex_30_20_20_mueller_scatgrid").unwrap();
    compare_results(result, reference, FRAC_TOL, ABS_TOL).unwrap();
}
