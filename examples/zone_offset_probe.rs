//! Probe the sensitivity of the forward/backward zone results to
//! ZONE_THETA_OFFSET by evaluating zero-width custom zones at a sweep of
//! theta values in a single run, including theta exactly 0 and exactly 180.
//!
//! Run with: cargo run --release --example zone_offset_probe

use goad::bins;
use goad::diff::Mapping;
use goad::geom::Geom;
use goad::multiproblem::MultiProblem;
use goad::orientation::{Euler, EulerConvention, Orientation, Scheme as OrientScheme};
use goad::settings::constants::{
    self, DEFAULT_BEAM_AREA_THRESHOLD_FAC, DEFAULT_BEAM_POWER_THRESHOLD, DEFAULT_CUTOFF,
    DEFAULT_MAX_REC, DEFAULT_MAX_TIR, DEFAULT_MEDIUM_REFR_INDEX,
};
use goad::settings::Settings;
use goad::zones::ZoneConfig;
use std::f32::consts::PI;
use std::path::PathBuf;

fn zero_width_zone(theta: f32) -> ZoneConfig {
    ZoneConfig::with_label(
        format!("probe_{theta}"),
        bins::Scheme::Custom {
            bins: vec![[[theta, theta], [0.0, 0.0]]],
            file: None,
        },
    )
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    let forward_thetas = [0.0_f32, 1e-4, 1e-3, 0.01, 0.05, 0.1];
    let backward_thetas = [179.9_f32, 179.95, 179.99, 179.999, 179.9999, 180.0];

    let mut zones: Vec<ZoneConfig> = Vec::new();
    for &t in forward_thetas.iter().chain(backward_thetas.iter()) {
        zones.push(zero_width_zone(t));
    }

    let wavelength = 0.532_f32;
    let settings = Settings {
        wavelength,
        beam_power_threshold: DEFAULT_BEAM_POWER_THRESHOLD,
        beam_area_threshold_fac: DEFAULT_BEAM_AREA_THRESHOLD_FAC,
        cutoff: DEFAULT_CUTOFF,
        medium_refr_index: DEFAULT_MEDIUM_REFR_INDEX,
        orientation: Orientation {
            scheme: OrientScheme::Sobol { num_orients: 300 },
            euler_convention: EulerConvention::ZYZ,
        },
        max_rec: DEFAULT_MAX_REC,
        max_tir: DEFAULT_MAX_TIR,
        zones,
        binning: None,
        seed: Some(42),
        scale: Some(1.0),
        distortion: None,
        geom_scale: None,
        directory: PathBuf::from("tmp/zone_offset_probe"),
        fov_factor: None,
        mapping: Mapping::ApertureDiffraction,
        output: constants::default_output_config(),
        coherence: true,
        quiet: true,
    };

    let geoms = Geom::load(
        "examples/data/hex.obj",
        vec![nalgebra::Complex::new(1.31, 0.0)],
    )?;

    let mut multiproblem = MultiProblem::new(geoms, Some(settings))?;
    multiproblem.solve();

    let results = multiproblem.get_results();
    let k = 2.0 * PI / wavelength;

    println!("theta       S11             S22             ext_cross_OT    ampl(0,0)");
    for &t in forward_thetas.iter().chain(backward_thetas.iter()) {
        let label = format!("probe_{t}");
        let Some(zone) = results.zones.get(&label) else {
            println!("{t:<10}  <zone missing>");
            continue;
        };
        let field = &zone.field_2d[0];
        let s11 = field.mueller_total[(0, 0)];
        let s22 = field.mueller_total[(1, 1)];
        let s2 = field.ampl_total[(0, 0)];
        let ext_ot = -s2.im * 4.0 * PI / k.powi(2);
        println!(
            "{t:<10}  {s11:<14.6e}  {s22:<14.6e}  {ext_ot:<14.6e}  {re:.6e}{im:+.6e}i",
            re = s2.re,
            im = s2.im,
        );
    }

    Ok(())
}
