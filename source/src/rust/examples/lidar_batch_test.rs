//! Rust replication of the Python run_test_batch.py script.
//! Tests plate_50um configurations to compare with Python results.

use std::path::PathBuf;
use std::time::Instant;

use goad::convergence::Convergence;
use goad::orientation::{EulerConvention, Orientation, Scheme};
use goad::params::Param;
use goad::result::GOComponent;
use goad::settings::Settings;

use nalgebra::Complex;

// Physical parameters from Python config
const WAVELENGTH_UM: f32 = 0.355; // ATLID wavelength
const REFR_INDEX_RE: f32 = 1.3249; // Ice at 355nm
const REFR_INDEX_IM: f32 = 2e-11; // Absorption (negligible)

// Convergence settings
const CONVERGENCE_THRESHOLD: f32 = 0.05; // 5% relative error
const MAX_ORIENTATIONS: usize = 100_000;

// GOAD settings for backscatter
const MAX_TIR: i32 = 20;
const BEAM_AREA_THRESHOLD_FAC: f32 = 0.01;
const BEAM_POWER_THRESHOLD: f32 = 0.01;
const CUTOFF: f32 = 0.999;

fn run_single_config(
    geom_dir: &PathBuf,
    results_dir: &PathBuf,
    config_name: &str,
) -> Result<ConfigResult, String> {
    // Find all .obj files in the geometry directory
    let obj_files: Vec<_> = std::fs::read_dir(geom_dir)
        .map_err(|e| format!("Failed to read dir: {}", e))?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "obj"))
        .collect();

    if obj_files.is_empty() {
        return Err(format!("No OBJ files in {:?}", geom_dir));
    }

    println!("\n{}", "=".repeat(60));
    println!("Running: {}", config_name);
    println!("Geometries: {} files", obj_files.len());
    println!("{}", "=".repeat(60));

    // Create settings matching Python
    let settings = Settings {
        wavelength: WAVELENGTH_UM,
        medium_refr_index: Complex::new(1.0, 0.0),
        particle_refr_index: vec![Complex::new(REFR_INDEX_RE, REFR_INDEX_IM)],
        geom_name: geom_dir.to_string_lossy().to_string(),
        zones: vec![], // Empty zones = auto Forward/Backward
        max_tir: MAX_TIR,
        max_rec: 10,
        beam_area_threshold_fac: BEAM_AREA_THRESHOLD_FAC,
        beam_power_threshold: BEAM_POWER_THRESHOLD,
        cutoff: CUTOFF,
        directory: results_dir.clone(),
        orientation: Orientation {
            scheme: Scheme::Sobol { num_orients: 1 },
            euler_convention: EulerConvention::ZYZ,
        },
        seed: None,
        scale: 1.0,
        distortion: None,
        geom_scale: None,
        fov_factor: None,
        mapping: goad::diff::Mapping::ApertureDiffraction,
        output: goad::settings::OutputConfig {
            settings_json: true,
            mueller_2d: false,
            mueller_1d: false,
            mueller_components: goad::settings::MuellerComponentConfig {
                total: true,
                beam: false,
                external: false,
            },
        },
        coherence: true,
        quiet: false,
        binning: None,
    };

    // Create convergence solver
    let mut convergence = Convergence::new(None, Some(settings)).map_err(|e| e.to_string())?;

    // Enable logging
    let log_file = results_dir.join("convergence_log.csv");
    convergence.set_log_file(&log_file);

    // Add target for LidarRatio only (faster)
    convergence.add_target(Param::LidarRatio, CONVERGENCE_THRESHOLD);
    convergence.max_orientations = MAX_ORIENTATIONS;

    // Run with timing
    let start = Instant::now();
    convergence.solve().map_err(|e| e.to_string())?;
    let elapsed = start.elapsed();

    // Save results
    convergence.writeup();

    // Extract results
    let mean = convergence.mean();
    let sem = convergence.sem();

    let lidar_ratio = mean.params.lidar_ratio(&GOComponent::Total).unwrap_or(0.0);
    let lidar_ratio_sem = sem.params.lidar_ratio(&GOComponent::Total).unwrap_or(0.0);
    let depol_ratio = mean
        .params
        .depolarization_ratio(&GOComponent::Total)
        .unwrap_or(0.0);
    let depol_ratio_sem = sem
        .params
        .depolarization_ratio(&GOComponent::Total)
        .unwrap_or(0.0);

    println!("\nResults:");
    println!(
        "  Time: {:.1}s ({:.1} min)",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() / 60.0
    );
    println!("  Orientations: {}", convergence.count());
    println!("  S = {:.2} +/- {:.2} sr", lidar_ratio, lidar_ratio_sem);
    println!(
        "  delta = {:.1} +/- {:.1} %",
        depol_ratio * 100.0,
        depol_ratio_sem * 100.0
    );

    Ok(ConfigResult {
        config: config_name.to_string(),
        time_seconds: elapsed.as_secs_f64(),
        num_orientations: convergence.count(),
        lidar_ratio,
        lidar_ratio_sem,
        depol_ratio,
        depol_ratio_sem,
    })
}

#[derive(Debug)]
struct ConfigResult {
    config: String,
    time_seconds: f64,
    num_orientations: usize,
    lidar_ratio: f32,
    lidar_ratio_sem: f32,
    depol_ratio: f32,
    depol_ratio_sem: f32,
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let base_dir = PathBuf::from(
        "/Users/ixguard/Documents/work/agents/studies/radiative-transfer/earthcare-atlid-validation/goad-lidar-study",
    );
    let geom_base = base_dir.join("geometries");
    let results_base = base_dir.join("results_rust");

    // Test configurations: plate 50um only (matching Python test batch)
    let test_configs = [
        ("plate_50um_pristine", "plate", 50.0, "pristine"),
        ("plate_50um_light", "plate", 50.0, "light"),
        ("plate_50um_moderate", "plate", 50.0, "moderate"),
        ("plate_50um_severe", "plate", 50.0, "severe"),
    ];

    println!("{}", "=".repeat(60));
    println!("GOAD Rust Test Batch: 50um Plates");
    println!("Convergence target: {:.0}%", CONVERGENCE_THRESHOLD * 100.0);
    println!("{}", "=".repeat(60));

    let mut results = Vec::new();
    let total_start = Instant::now();

    for (config_name, _habit, _size, _roughness) in &test_configs {
        let geom_dir = geom_base.join(config_name);
        let results_dir = results_base.join(config_name);
        std::fs::create_dir_all(&results_dir).ok();

        match run_single_config(&geom_dir, &results_dir, config_name) {
            Ok(result) => results.push(result),
            Err(e) => {
                println!("ERROR: {}", e);
            }
        }
    }

    let total_elapsed = total_start.elapsed();

    // Summary
    println!("\n{}", "=".repeat(60));
    println!("SUMMARY");
    println!("{}", "=".repeat(60));
    println!(
        "\n{:<25} {:>8} {:>8} {:>12} {:>12}",
        "Config", "Time", "Orient", "S (sr)", "delta (%)"
    );
    println!("{}", "-".repeat(70));

    for r in &results {
        let time_str = format!("{:.0}s", r.time_seconds);
        let s_str = format!("{:.1}+/-{:.1}", r.lidar_ratio, r.lidar_ratio_sem);
        let d_str = format!(
            "{:.1}+/-{:.1}",
            r.depol_ratio * 100.0,
            r.depol_ratio_sem * 100.0
        );
        println!(
            "{:<25} {:>8} {:>8} {:>12} {:>12}",
            r.config, time_str, r.num_orientations, s_str, d_str
        );
    }

    println!("{}", "-".repeat(70));
    println!(
        "Total time: {:.0}s ({:.1} min)",
        total_elapsed.as_secs_f64(),
        total_elapsed.as_secs_f64() / 60.0
    );

    // Save results as JSON
    let output_file = results_base.join("test_batch_results.json");
    if let Ok(json) = serde_json::to_string_pretty(&serde_json::json!({
        "results": results.iter().map(|r| serde_json::json!({
            "config": r.config,
            "time_seconds": r.time_seconds,
            "num_orientations": r.num_orientations,
            "lidar_ratio_sr": r.lidar_ratio,
            "lidar_ratio_sem": r.lidar_ratio_sem,
            "depol_ratio": r.depol_ratio,
            "depol_ratio_sem": r.depol_ratio_sem,
        })).collect::<Vec<_>>(),
        "total_time_seconds": total_elapsed.as_secs_f64(),
    })) {
        std::fs::write(&output_file, json).ok();
        println!("\nResults saved to: {:?}", output_file);
    }
}
