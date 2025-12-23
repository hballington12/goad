// --8<-- [start:convergence]
fn main() {
    use goad::convergence::Convergence;
    use goad::orientation::{EulerConvention, Orientation, Scheme};
    use goad::params::Param;
    use goad::result::GOComponent;
    use goad::settings;

    // Load default settings and configure for convergence
    let mut settings = settings::load_default_config().unwrap();

    // Use uniform random orientations (enough for convergence)
    settings.orientation = Orientation {
        scheme: Scheme::Uniform { num_orients: 2000 },
        euler_convention: EulerConvention::ZYZ,
    };

    // Create a convergence solver
    let mut convergence = Convergence::new(None, Some(settings)).unwrap();

    // Set convergence target: 1% relative SEM on asymmetry parameter
    convergence.add_target(Param::Asymmetry, 0.01);

    // Optionally set max orientations as safety cap (default is 100k)
    convergence.max_orientations = 2000;

    // Solve - will terminate when target is reached or max_orientations hit
    convergence.solve();

    // Print results
    let asym = convergence
        .result
        .params
        .asymmetry(&GOComponent::Total)
        .unwrap_or(0.0);
    let asym_sem = convergence
        .error
        .params
        .asymmetry(&GOComponent::Total)
        .unwrap_or(0.0);
    let relative_sem = (asym_sem / asym.abs()) * 100.0;

    println!("Orientations computed: {}", convergence.count());
    println!(
        "Asymmetry: {:.4} +/- {:.4} ({:.2}% relative SEM)",
        asym, asym_sem, relative_sem
    );
}
// --8<-- [end:convergence]
