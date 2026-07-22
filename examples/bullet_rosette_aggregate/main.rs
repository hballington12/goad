//! Rust reproduction of run_simulation.py, to debug the clipping-routine bug
//! surfaced by the large 3-bullet-rosette aggregate (trnc_clip wildly out of range).
//!
//! Settings are built to match goad-py's `Settings()` defaults exactly (NOT
//! `settings::load_default_config()`, which reads config/default.toml and has
//! different defaults - discrete orientation, looser thresholds, etc.)
//!
//! Run with: cargo run --release --example bullet_rosette_aggregate

use goad::bins;
use goad::geom::Geom;
use goad::multiproblem::MultiProblem;
use goad::orientation::{Orientation, Scheme};
use goad::settings::constants::{
    self, DEFAULT_BEAM_AREA_THRESHOLD_FAC, DEFAULT_BEAM_POWER_THRESHOLD, DEFAULT_COHERENCE,
    DEFAULT_CUTOFF, DEFAULT_EULER_ORDER, DEFAULT_MAPPING, DEFAULT_MAX_REC, DEFAULT_MAX_TIR,
    DEFAULT_MEDIUM_REFR_INDEX,
};
use goad::settings::Settings;
use goad::zones::ZoneConfig;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let output_dir = PathBuf::from("examples/bullet_rosette_aggregate");
    let geometry_file = output_dir.join("aggregate-2p9mm_tri.obj");

    let settings = Settings {
        wavelength: 0.2,
        beam_power_threshold: DEFAULT_BEAM_POWER_THRESHOLD,
        beam_area_threshold_fac: DEFAULT_BEAM_AREA_THRESHOLD_FAC,
        cutoff: DEFAULT_CUTOFF,
        medium_refr_index: DEFAULT_MEDIUM_REFR_INDEX,
        orientation: Orientation {
            scheme: Scheme::Sobol { num_orients: 1 },
            euler_convention: DEFAULT_EULER_ORDER,
        },
        max_rec: DEFAULT_MAX_REC,
        max_tir: DEFAULT_MAX_TIR,
        zones: vec![ZoneConfig::new(bins::Scheme::Interval {
            thetas: vec![0.0, 5.0, 175.0, 179.0, 180.0],
            theta_spacings: vec![0.1, 2.0, 0.5, 0.1],
            phis: vec![0.0, 360.0],
            phi_spacings: vec![7.5],
        })],
        binning: None,
        seed: Some(42),
        scale: Some(1.0),
        distortion: None,
        geom_scale: None,
        directory: output_dir.clone(),
        fov_factor: None,
        mapping: DEFAULT_MAPPING,
        output: constants::default_output_config(),
        coherence: DEFAULT_COHERENCE,
        quiet: false,
    };

    let geoms = Geom::load(
        geometry_file.to_str().expect("non-utf8 geometry path"),
        vec![nalgebra::Complex::new(1.39, 0.0)],
    )?;

    let mut multiproblem = MultiProblem::new(geoms, Some(settings))?;
    multiproblem.solve();

    println!("{:#?}", multiproblem.get_results().powers);

    Ok(())
}
