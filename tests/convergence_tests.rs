use goad::{
    bins::{self, BinningScheme},
    convergence::Convergence,
    multiproblem::MultiProblem,
    orientation::{Euler, EulerConvention, Orientation, Scheme as OrientScheme},
    params::Param,
    result::MuellerMatrix,
    settings,
};
use helpers::compare_results;

pub mod helpers;

// Tolerance for comparing Mueller matrix elements between solvers
const FRAC_TOL: f32 = 1e-4;
const ABS_TOL: f32 = 1e4;

/// Compare Convergence vs MultiProblem with the same orientations.
/// Both should produce identical results when given identical inputs.
#[test]
fn convergence_vs_multiproblem_identical() {
    let mut settings = settings::load_default_config().unwrap();

    // Use a small fixed set of orientations for reproducibility
    settings.binning = BinningScheme {
        scheme: bins::Scheme::new_simple(9, 9),
    };
    settings.orientation = Orientation {
        scheme: OrientScheme::Discrete {
            eulers: vec![
                Euler::new(0.0, 0.0, 0.0),
                Euler::new(30.0, 30.0, 30.0),
                Euler::new(60.0, 45.0, 15.0),
                Euler::new(90.0, 60.0, 45.0),
            ],
        },
        euler_convention: EulerConvention::ZYZ,
    };
    settings.seed = Some(42); // Fixed seed for reproducibility
    settings.quiet = true;

    // Solve with MultiProblem
    let mut multiproblem =
        MultiProblem::new(None, Some(settings.clone())).expect("Failed to create MultiProblem");
    multiproblem.solve();

    // Solve with Convergence (set target to match orientation count)
    let mut convergence =
        Convergence::new(None, Some(settings)).expect("Failed to create Convergence");
    convergence.add_target(Param::Asymmetry, 0.001); // tight target to ensure all 4 run
    convergence.max_orientations = 4; // match the 4 discrete orientations
    convergence.solve().unwrap();

    // Compare results
    let mp_result: Vec<Vec<f32>> = multiproblem
        .result
        .field_2d
        .iter()
        .map(|m| m.mueller_total.to_vec())
        .collect();

    let conv_result: Vec<Vec<f32>> = convergence
        .mean()
        .field_2d
        .iter()
        .map(|m| m.mueller_total.to_vec())
        .collect();

    compare_results(conv_result, mp_result, FRAC_TOL, ABS_TOL)
        .expect("Convergence and MultiProblem results should match");
}
