use goad::{
    bins::{self, BinningScheme},
    convergence::Convergence,
    multiproblem::MultiProblem,
    orientation::{EulerConvention, Orientation, Scheme as OrientScheme},
    result::MuellerMatrix,
    settings,
};
use std::fs::File;
use std::io::Write;

/// Compare Convergence vs MultiProblem and dump 1D Mueller results to files.
#[test]
fn dump_1d_mueller_comparison() {
    let mut settings = settings::load_default_config().unwrap();

    // Use uniform orientations for comparison
    settings.binning = BinningScheme {
        scheme: bins::Scheme::new_simple(37, 37),
    };
    settings.orientation = Orientation {
        scheme: OrientScheme::Uniform { num_orients: 100 },
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
    convergence.convergence_target = 100;
    convergence.solve();

    // Dump 1D Mueller results
    let mp_1d = multiproblem
        .result
        .field_1d
        .as_ref()
        .expect("No 1D results");
    let conv_1d = convergence.result.field_1d.as_ref().expect("No 1D results");

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

    // Print SEM values from convergence error field
    println!("\n=== Convergence SEM Values ===");
    println!("Powers SEM:");
    println!("  Input:    {}", convergence.error.powers.input);
    println!("  Output:   {}", convergence.error.powers.output);
    println!("  Absorbed: {}", convergence.error.powers.absorbed);

    println!("\nParams SEM:");
    if let Some(asym) = convergence
        .error
        .params
        .asymmetry(&goad::result::GOComponent::Total)
    {
        println!("  Asymmetry: {}", asym);
    }
    if let Some(albedo) = convergence
        .error
        .params
        .albedo(&goad::result::GOComponent::Total)
    {
        println!("  Albedo: {}", albedo);
    }

    // Also print the mean values for comparison
    println!("\n=== Convergence Mean Values ===");
    println!("Powers Mean:");
    println!("  Input:    {}", convergence.result.powers.input);
    println!("  Output:   {}", convergence.result.powers.output);
    println!("  Absorbed: {}", convergence.result.powers.absorbed);

    println!("\nParams Mean:");
    if let Some(asym) = convergence
        .result
        .params
        .asymmetry(&goad::result::GOComponent::Total)
    {
        println!("  Asymmetry: {}", asym);
    }
    if let Some(albedo) = convergence
        .result
        .params
        .albedo(&goad::result::GOComponent::Total)
    {
        println!("  Albedo: {}", albedo);
    }
}
