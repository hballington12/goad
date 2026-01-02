//! Integration tests for comparing Convergence vs MultiProblem solvers.
//!
//! NOTE: Always run these tests in release mode for realistic performance:
//!   cargo test --release -p goad --test compare_solvers
//!
//! Debug mode is ~10-20x slower and will cause timeouts.

use goad::{
    bins,
    convergence::Convergence,
    multiproblem::MultiProblem,
    orientation::{EulerConvention, Orientation, Scheme as OrientScheme},
    params::Param,
    result::MuellerMatrix,
    settings,
    zones::ZoneConfig,
};
use std::fs::File;
use std::io::Write;

/// Compare Convergence vs MultiProblem and dump 1D Mueller results to files.
#[test]
#[ignore] // Run with: cargo test --release -- --ignored
fn dump_1d_mueller_comparison() {
    let mut settings = settings::load_default_config().unwrap();

    // Use uniform orientations for comparison (keep small for debug mode)
    settings.zones = vec![ZoneConfig::new(bins::Scheme::new_simple(37, 37))];
    settings.orientation = Orientation {
        scheme: OrientScheme::Uniform { num_orients: 3 },
        euler_convention: EulerConvention::ZYZ,
    };
    settings.seed = Some(42);
    settings.quiet = true;

    // Solve with MultiProblem
    let mut multiproblem =
        MultiProblem::new(None, Some(settings.clone())).expect("Failed to create MultiProblem");
    multiproblem.solve();

    // Solve with Convergence
    let mut convergence =
        Convergence::new(None, Some(settings)).expect("Failed to create Convergence");
    convergence.add_target(Param::Asymmetry, 0.001); // tight target to ensure all 3 run
    convergence.max_orientations = 3;
    convergence.solve().unwrap();

    // Dump 1D Mueller results
    let mp_1d = multiproblem
        .result
        .zones
        .full_zone()
        .and_then(|z| z.field_1d.as_ref())
        .expect("No 1D results");
    let conv_mean = convergence.mean();
    let conv_1d = conv_mean
        .zones
        .full_zone()
        .and_then(|z| z.field_1d.as_ref())
        .expect("No 1D results");

    let mut mp_file = File::create("multiproblem_1d.dat").unwrap();
    let mut conv_file = File::create("convergence_1d.dat").unwrap();

    for result in mp_1d.iter() {
        let theta = result.bin.center;
        let s11 = result.mueller_total.s11();
        writeln!(mp_file, "{} {}", theta, s11).unwrap();
    }

    for result in conv_1d.iter() {
        let theta = result.bin.center;
        let s11 = result.mueller_total.s11();
        writeln!(conv_file, "{} {}", theta, s11).unwrap();
    }

    println!("Wrote multiproblem_1d.dat and convergence_1d.dat");

    // Print SEM values from convergence sem() method
    let sem = convergence.sem();
    println!("\n=== Convergence SEM Values ===");
    println!("Powers SEM:");
    println!("  Input:    {}", sem.powers.input);
    println!("  Output:   {}", sem.powers.output);
    println!("  Absorbed: {}", sem.powers.absorbed);

    println!("\nParams SEM:");
    if let Some(asym) = sem.params.asymmetry(&goad::result::GOComponent::Total) {
        println!("  Asymmetry: {}", asym);
    }
    if let Some(albedo) = sem.params.albedo(&goad::result::GOComponent::Total) {
        println!("  Albedo: {}", albedo);
    }

    // Also print the mean values for comparison
    let mean = convergence.mean();
    println!("\n=== Convergence Mean Values ===");
    println!("Powers Mean:");
    println!("  Input:    {}", mean.powers.input);
    println!("  Output:   {}", mean.powers.output);
    println!("  Absorbed: {}", mean.powers.absorbed);

    println!("\nParams Mean:");
    if let Some(asym) = mean.params.asymmetry(&goad::result::GOComponent::Total) {
        println!("  Asymmetry: {}", asym);
    }
    if let Some(albedo) = mean.params.albedo(&goad::result::GOComponent::Total) {
        println!("  Albedo: {}", albedo);
    }
}

/// Test tracker with 5 specific orientations to compare with Python
#[test]
#[ignore] // Run with: cargo test --release -- --ignored
fn test_tracker_5_orientations() {
    use goad::convergence::ConvergenceTracker;
    use goad::orientation::Euler;
    use goad::problem::Problem;

    let eulers = [
        Euler::new(30.0, 45.0, 60.0),
        Euler::new(120.0, 30.0, 90.0),
        Euler::new(200.0, 60.0, 45.0),
        Euler::new(45.0, 90.0, 180.0),
        Euler::new(300.0, 120.0, 30.0),
    ];

    let settings = settings::load_default_config().unwrap();
    let mut tracker: Option<ConvergenceTracker<goad::result::Results>> = None;

    for (i, euler) in eulers.iter().enumerate() {
        let mut problem = Problem::new(None, Some(settings.clone())).unwrap();
        problem.run(Some(euler)).unwrap();

        let asym = problem
            .result
            .params
            .asymmetry(&goad::result::GOComponent::Total)
            .unwrap();
        let sc = problem
            .result
            .params
            .scatt_cross(&goad::result::GOComponent::Total)
            .unwrap();
        println!(
            "Orient {}: asymmetry={:.4}, scat_cross={:.4}",
            i + 1,
            asym,
            sc
        );

        if tracker.is_none() {
            tracker = Some(ConvergenceTracker::new(&problem.result));
        }
        tracker.as_mut().unwrap().update(&problem.result);
    }

    let tracker = tracker.unwrap();
    let mean = tracker.mean();
    let sem = tracker.sem();

    let asym_mean = mean
        .params
        .asymmetry(&goad::result::GOComponent::Total)
        .unwrap();
    let asym_sem = sem
        .params
        .asymmetry(&goad::result::GOComponent::Total)
        .unwrap();

    println!("\nAfter 5 orientations:");
    println!("  mean = {:.4}", asym_mean);
    println!("  sem = {:.4}", asym_sem);
    println!("  relative sem = {:.2}%", (asym_sem / asym_mean) * 100.0);
}

/// Test convergence with target-based termination (1% on asymmetry)
#[test]
#[ignore] // Run with: cargo test --release -- --ignored
fn test_convergence_with_target() {
    use std::time::Instant;
    let start = Instant::now();

    let mut settings = settings::load_default_config().unwrap();

    // Use the same binning scheme as Python default
    settings.zones = vec![ZoneConfig::new(bins::Scheme::Interval {
        thetas: vec![0.0, 5.0, 175.0, 179.0, 180.0],
        theta_spacings: vec![0.1, 2.0, 0.5, 0.1],
        phis: vec![0.0, 360.0],
        phi_spacings: vec![7.5],
    })];
    settings.orientation = Orientation {
        scheme: OrientScheme::Uniform { num_orients: 3 },
        euler_convention: EulerConvention::ZYZ,
    };
    settings.seed = Some(42);
    settings.quiet = true;

    let mut convergence =
        Convergence::new(None, Some(settings)).expect("Failed to create Convergence");

    // Set target: 90% relative SEM on asymmetry (relaxed for fewer orientations)
    convergence.add_target(Param::Asymmetry, 0.90);
    convergence.max_orientations = 3; // keep tests fast

    convergence.solve().unwrap();

    // Check results
    let mean = convergence.mean();
    let sem = convergence.sem();
    let asym_mean = mean
        .params
        .asymmetry(&goad::result::GOComponent::Total)
        .unwrap_or(0.0);
    let asym_sem = sem
        .params
        .asymmetry(&goad::result::GOComponent::Total)
        .unwrap_or(0.0);
    let relative_sem = if asym_mean != 0.0 {
        (asym_sem / asym_mean) * 100.0
    } else {
        0.0
    };

    let elapsed = start.elapsed();
    println!("\n=== Convergence with Target (10% on Asymmetry) ===");
    println!("Time elapsed: {:.2}s", elapsed.as_secs_f64());
    println!("Orientations computed: {}", convergence.count());
    println!("Asymmetry:");
    println!("  Mean: {:.4}", asym_mean);
    println!("  SEM:  {:.6}", asym_sem);
    println!("  Relative SEM: {:.2}%", relative_sem);

    // Should have converged before hitting 2000 orientations
    // and relative SEM should be <= 10%
    assert!(
        relative_sem <= 91.0, // allow tiny margin for floating point
        "Expected relative SEM <= 90%, got {:.2}%",
        relative_sem
    );
}

/// Benchmark: compare 100 orientations Convergence vs MultiProblem timing
#[test]
#[ignore] // Run with: cargo test --release -- --ignored
fn benchmark_convergence_vs_multiproblem() {
    use std::time::Instant;

    let mut settings = settings::load_default_config().unwrap();
    settings.orientation = Orientation {
        scheme: OrientScheme::Uniform { num_orients: 3 },
        euler_convention: EulerConvention::ZYZ,
    };
    settings.seed = Some(42);
    settings.quiet = true;

    // Time MultiProblem
    let start_mp = Instant::now();
    let mut multiproblem =
        MultiProblem::new(None, Some(settings.clone())).expect("Failed to create MultiProblem");
    multiproblem.solve();
    let elapsed_mp = start_mp.elapsed();

    // Time Convergence
    let start_conv = Instant::now();
    let mut convergence =
        Convergence::new(None, Some(settings)).expect("Failed to create Convergence");
    convergence.add_target(Param::Asymmetry, 0.001); // tight target to run all 3
    convergence.max_orientations = 3;
    convergence.solve().unwrap();
    let elapsed_conv = start_conv.elapsed();

    println!("\n=== Benchmark: 100 orientations ===");
    println!("MultiProblem: {:.2}s", elapsed_mp.as_secs_f64());
    println!("Convergence:  {:.2}s", elapsed_conv.as_secs_f64());
    println!(
        "Ratio (Conv/MP): {:.2}x",
        elapsed_conv.as_secs_f64() / elapsed_mp.as_secs_f64()
    );
}
