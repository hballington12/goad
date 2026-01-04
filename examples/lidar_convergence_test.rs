//! Test script for lidar ratio convergence with zones=[]
//! This script logs the indicatif mean value to a file as it runs.

fn main() {
    use goad::convergence::Convergence;
    use goad::orientation::{Orientation, Scheme};
    use goad::params::Param;
    use goad::result::GOComponent;
    use goad::settings;

    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Load default settings (uses default geometry and config)
    let mut settings = settings::load_default_config().unwrap();

    // Clear zones to use only auto-generated Forward/Backward
    settings.zones = vec![];

    // Use Sobol sampling for convergence (not discrete!)
    settings.orientation = Orientation {
        scheme: Scheme::Sobol { num_orients: 1 },
        euler_convention: settings.orientation.euler_convention,
    };

    // Create a convergence solver
    let mut convergence = Convergence::new(None, Some(settings)).unwrap();

    // Enable logging of mean values during convergence
    convergence.set_log_file("lidar_convergence_running.csv");

    // Set convergence target: 1% relative SEM on lidar ratio only
    convergence.add_target(Param::LidarRatio, 0.01);

    // Set max orientations as safety cap
    convergence.max_orientations = 500;

    // Solve - will terminate when target is reached or max_orientations hit
    convergence.solve().unwrap();

    // Get final results
    let mean = convergence.mean();
    let sem = convergence.sem();

    // Extract lidar ratio values
    let lidar = mean.params.lidar_ratio(&GOComponent::Total).unwrap_or(0.0);
    let lidar_sem = sem.params.lidar_ratio(&GOComponent::Total).unwrap_or(0.0);
    let relative_sem = if lidar.abs() > 1e-10 {
        (lidar_sem / lidar.abs()) * 100.0
    } else {
        0.0
    };

    println!("\n=== Final Results ===");
    println!("Orientations computed: {}", convergence.count());
    println!(
        "Lidar Ratio: {:.4} +/- {:.4} ({:.2}% relative SEM)",
        lidar, lidar_sem, relative_sem
    );

    println!("\nRunning log written to lidar_convergence_running.csv");

    // Optional: write full output
    convergence.writeup();
}
