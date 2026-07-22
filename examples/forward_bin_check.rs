//! Reconcile the 'forward' zone output with the first bin of the full zone.
//!
//! A user reported the forward-zone beam and ext values are lower than the
//! full-zone phase function first bin by a factor ~0.15 when both sit at
//! theta = 0.01 deg. Hypothesis: the 1D phase function is integrated over
//! phi (carrying a factor of 2 pi vs a differential value), while the
//! forward zone is a differential point sample; 1/(2 pi) = 0.159.
//!
//! Run with: cargo run --release --example forward_bin_check

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
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    // Full zone: first theta bin [0, 0.02] centred on 0.01 deg, 4-deg phi bins
    let full = ZoneConfig::with_label(
        "full_fine_forward",
        bins::Scheme::Interval {
            thetas: vec![0.0, 0.02, 176.02, 180.0],
            theta_spacings: vec![0.02, 2.0, 1.99],
            phis: vec![0.0, 360.0],
            phi_spacings: vec![4.0],
        },
    );
    // Point-sample zone matching the email comparison: theta 0.01, phi 0
    let probe = ZoneConfig::with_label(
        "probe_0p01",
        bins::Scheme::Custom {
            bins: vec![[[0.01, 0.01], [0.0, 0.0]]],
            file: None,
        },
    );

    let settings = Settings {
        wavelength: 0.532,
        beam_power_threshold: DEFAULT_BEAM_POWER_THRESHOLD,
        beam_area_threshold_fac: DEFAULT_BEAM_AREA_THRESHOLD_FAC,
        cutoff: DEFAULT_CUTOFF,
        medium_refr_index: DEFAULT_MEDIUM_REFR_INDEX,
        orientation: Orientation {
            scheme: OrientScheme::Discrete {
                eulers: vec![Euler::new(30.0, 30.0, 0.0)],
            },
            euler_convention: EulerConvention::ZYZ,
        },
        max_rec: DEFAULT_MAX_REC,
        max_tir: DEFAULT_MAX_TIR,
        zones: vec![full, probe],
        binning: None,
        seed: Some(42),
        scale: Some(1.0),
        distortion: None,
        geom_scale: None,
        directory: PathBuf::from("tmp/forward_bin_check"),
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

    let mut mp = MultiProblem::new(geoms, Some(settings))?;
    mp.solve();

    let results = mp.get_results();

    let full_zone = results.zones.get("full_fine_forward").expect("full zone");
    let probe_zone = results.zones.get("probe_0p01").expect("probe zone");

    // Probe: differential point sample at (0.01, 0)
    let p = &probe_zone.field_2d[0];
    let (p_beam, p_ext, p_tot) = (
        p.mueller_beam[(0, 0)],
        p.mueller_ext[(0, 0)],
        p.mueller_total[(0, 0)],
    );

    // Full zone 2D: mean over phi of the first theta row
    let first_theta = full_zone.field_2d[0].bin.theta;
    let row: Vec<_> = full_zone
        .field_2d
        .iter()
        .filter(|f| f.bin.theta == first_theta)
        .collect();
    let n = row.len() as f32;
    let m_beam = row.iter().map(|f| f.mueller_beam[(0, 0)]).sum::<f32>() / n;
    let m_ext = row.iter().map(|f| f.mueller_ext[(0, 0)]).sum::<f32>() / n;
    let m_tot = row.iter().map(|f| f.mueller_total[(0, 0)]).sum::<f32>() / n;

    // Full zone 1D: phi-integrated first bin
    let f1d = &full_zone.field_1d.as_ref().expect("1d")[0];
    let (i_beam, i_ext, i_tot) = (
        f1d.mueller_beam[(0, 0)],
        f1d.mueller_ext[(0, 0)],
        f1d.mueller_total[(0, 0)],
    );

    println!("probe (theta=0.01, phi=0) point sample:");
    println!("  S11 beam {p_beam:.6e}  ext {p_ext:.6e}  total {p_tot:.6e}");
    println!("full zone first theta row ({} phi bins), mean over phi:", row.len());
    println!("  S11 beam {m_beam:.6e}  ext {m_ext:.6e}  total {m_tot:.6e}");
    println!("full zone 1D first bin (phi-integrated):");
    println!("  S11 beam {i_beam:.6e}  ext {i_ext:.6e}  total {i_tot:.6e}");
    println!();
    println!("ratio probe / full-2D-mean:   beam {:.4}  ext {:.4}  total {:.4}",
        p_beam / m_beam, p_ext / m_ext, p_tot / m_tot);
    println!("ratio probe / full-1D:        beam {:.4}  ext {:.4}  total {:.4}",
        p_beam / i_beam, p_ext / i_ext, p_tot / i_tot);
    println!("1 / (2 pi)                  = {:.4}", 1.0 / (2.0 * std::f32::consts::PI));
    println!("ratio full-1D / full-2D-mean: beam {:.4} (expect ~2 pi = {:.4})",
        i_beam / m_beam, 2.0 * std::f32::consts::PI);

    Ok(())
}
