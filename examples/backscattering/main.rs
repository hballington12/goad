//! Rust equivalent of run_simulation.py for debugging memory issues.
//!
//! Run with: cargo run --release --example backscattering

use goad::convergence::Convergence;
use goad::params::Param;
use goad::settings;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let output_dir = PathBuf::from("examples/backscattering");
    let geometry_file = output_dir
        .join("roughened_edge15.0_sigma0p1_merge1p0_hexagonal_column_l5.0_r25.0_6904c8.obj");

    // Load settings from defaults and override
    let mut settings = settings::load_default_config()?;
    settings.geom_name = geometry_file.to_string_lossy().to_string();
    settings.wavelength = 0.532;
    settings.particle_refr_index = vec![nalgebra::Complex::new(1.31, 0.0)];
    settings.zones = vec![]; // empty zones for backscatter only
    settings.max_tir = 20;
    settings.beam_area_threshold_fac = 0.01;
    settings.cutoff = 0.999;
    settings.beam_power_threshold = 0.01;
    settings.directory = output_dir.clone();
    settings.quiet = false;

    // Create convergence solver
    let mut convergence = Convergence::new(None, Some(settings))?;

    // Add targets: 5% error on lidar ratio and depolarization ratio
    convergence.add_target(Param::LidarRatio, 0.05);
    convergence.add_target(Param::DepolarizationRatio, 0.05);

    // Run the simulation
    convergence.solve()?;

    // Save results
    convergence.writeup();

    println!("Simulation complete!");
    println!("Results saved to: {:?}", output_dir);

    Ok(())
}
