// --8<-- [start:convergence]
fn main() {
    use goad::convergence::Convergence;
    use goad::params::Param;
    use goad::result::GOComponent;
    use goad::settings;

    // Load default settings
    let settings = settings::load_default_config().unwrap();

    // Create a convergence solver
    let mut convergence = Convergence::new(None, Some(settings)).unwrap();

    // Set convergence target: 1% relative SEM on asymmetry parameter
    convergence.add_target(Param::Asymmetry, 0.03);
    convergence.add_target(Param::Asymmetry, 0.02);
    convergence.add_target(Param::ScatCross, 0.02);
    convergence.add_target(Param::ExtCrossOpticalTheorem, 0.02);
    convergence.add_target(Param::DepolarizationRatio, 0.1);
    convergence.add_target(Param::LidarRatio, 0.1);

    // Optionally set max orientations as safety cap (default is 100k)
    convergence.max_orientations = 500;

    // Solve - will terminate when target is reached or max_orientations hit
    convergence.solve().unwrap();

    // Print results (using mean() and sem() methods)
    let mean = convergence.mean();
    let sem = convergence.sem();
    let asym = mean.params.asymmetry(&GOComponent::Total).unwrap_or(0.0);
    let asym_sem = sem.params.asymmetry(&GOComponent::Total).unwrap_or(0.0);
    let relative_sem = (asym_sem / asym.abs()) * 100.0;

    println!("Orientations computed: {}", convergence.count());
    println!(
        "Asymmetry: {:.4} +/- {:.4} ({:.2}% relative SEM)",
        asym, asym_sem, relative_sem
    );

    convergence.writeup();
}
// --8<-- [end:convergence]
