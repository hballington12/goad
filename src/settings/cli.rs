use clap::{Args, Parser};
use nalgebra::Complex;
use std::path::PathBuf;

use crate::bins::{self, BinningScheme};
use crate::diff::Mapping;
use crate::orientation::{Euler, EulerConvention, Orientation, Scheme};
use crate::settings::{Settings, DEFAULT_EULER_ORDER};

#[derive(Parser, Debug)]
#[command(version, about = "GOAD - Geometric Optics with Aperture Diffraction")]
#[command(author = "Harry Ballington")]
pub struct CliArgs {
    #[command(flatten)]
    pub propagation: PropagationArgs,

    #[command(flatten)]
    pub material: MaterialArgs,

    #[command(flatten)]
    pub orientation: OrientationArgs,

    #[command(flatten)]
    pub binning: BinningArgs,

    /// Random seed for reproducibility.
    #[arg(short, long)]
    pub seed: Option<u64>,

    /// Output directory for simulation results.
    #[arg(long)]
    pub dir: Option<PathBuf>,

    /// Set the field of view truncation factor for diffraction of beams.
    #[arg(long)]
    pub fov_factor: Option<f32>,

    /// Set the near to far field mapping method.
    #[arg(long, value_enum)]
    pub mapping: Option<Mapping>,

    /// Set the coherence settings for near to far field mapping
    #[arg(long, value_enum)]
    pub coherence: Option<bool>,
}

/// Beam propagation parameters
#[derive(Args, Debug)]
pub struct PropagationArgs {
    /// Wavelength in units of the geometry.
    #[arg(short, long)]
    pub w: Option<f32>,

    /// Minimum beam power threshold for propagation.
    #[arg(long)]
    pub bp: Option<f32>,

    /// Minimum area factor threshold for beam propagation.
    #[arg(long)]
    pub baf: Option<f32>,

    /// Total power cutoff fraction (0.0-1.0).
    #[arg(long)]
    pub cop: Option<f32>,

    /// Maximum recursion depth for beam tracing.
    #[arg(long)]
    pub rec: Option<i32>,

    /// Maximum allowed total internal reflections.
    #[arg(long)]
    pub tir: Option<i32>,
}

/// Material and geometry parameters
#[derive(Args, Debug)]
pub struct MaterialArgs {
    /// Path to geometry file (.obj format).
    #[arg(short, long)]
    pub geo: Option<String>,

    /// Surrounding medium refractive index.
    #[arg(long)]
    pub ri0: Option<Complex<f32>>,

    /// Particle refractive indices, space-separated.
    #[arg(short, long, value_parser, num_args = 1.., value_delimiter = ' ')]
    pub ri: Option<Vec<Complex<f32>>>,

    /// Distortion factor for the geometry.
    #[arg(long)]
    pub distortion: Option<f32>,

    /// Geometry scale factors for each axis (x, y, z).
    #[arg(long, value_parser, num_args = 1..=3, value_delimiter = ' ')]
    pub geom_scale: Option<Vec<f32>>,
}

/// Orientation parameters
#[derive(Args, Debug)]
pub struct OrientationArgs {
    /// Use uniform random orientation scheme.
    #[arg(long, group = "orientation")]
    pub uniform: Option<usize>,

    /// Use discrete orientation scheme with specified Euler angles (degrees).
    #[arg(long, value_parser = parse_euler_angles, num_args = 1.., value_delimiter = ' ', group = "orientation")]
    pub discrete: Option<Vec<Euler>>,

    /// Specify Euler angle convention for orientation.
    #[arg(long, value_parser = parse_euler_convention)]
    pub euler: Option<EulerConvention>,
}

/// Output binning parameters
#[derive(Args, Debug)]
pub struct BinningArgs {
    /// Use simple equal-spacing binning scheme.
    #[arg(long, num_args = 2, value_delimiter = ' ', group = "binning")]
    pub simple: Option<Vec<usize>>,

    /// Enable interval binning scheme with variable spacing.
    #[arg(long, group = "binning")]
    pub interval: bool,

    /// Theta angle bins for interval binning (degrees).
    #[arg(long, requires = "interval", num_args = 3.., value_delimiter = ' ')]
    pub theta: Option<Vec<f32>>,

    /// Phi angle bins for interval binning (degrees).
    #[arg(long, requires = "interval", num_args = 3.., value_delimiter = ' ')]
    pub phi: Option<Vec<f32>>,

    /// Path to custom binning scheme file.
    #[arg(long)]
    pub custom: Option<String>,
}

fn parse_euler_angles(s: &str) -> Result<Euler, String> {
    println!("Parsing Euler angles: '{}'", s);

    let angles: Vec<&str> = s.split(',').collect();
    if angles.len() != 3 {
        return Err(format!(
            "Invalid Euler angle format: '{}'. Expected 'alpha,beta,gamma'",
            s
        ));
    }

    let alpha = angles[0]
        .parse::<f32>()
        .map_err(|_| format!("Failed to parse alpha angle: {}", angles[0]))?;
    let beta = angles[1]
        .parse::<f32>()
        .map_err(|_| format!("Failed to parse beta angle: {}", angles[1]))?;
    let gamma = angles[2]
        .parse::<f32>()
        .map_err(|_| format!("Failed to parse gamma angle: {}", angles[2]))?;

    println!("Parsed Euler angles: {}, {}, {}", alpha, beta, gamma);

    Ok(Euler::new(alpha, beta, gamma))
}

fn parse_interval_specification(values: &[f32]) -> Result<(Vec<f32>, Vec<f32>), String> {
    if values.len() < 3 {
        return Err(format!(
            "Insufficient values for interval specification: need at least 3, got {}",
            values.len()
        ));
    }

    let mut positions = Vec::new();
    let mut spacings = Vec::new();

    positions.push(values[0]);

    for i in (1..values.len() - 1).step_by(2) {
        let step = values[i];
        let pos = values[i + 1];

        if step <= 0.0 {
            return Err(format!("Step size must be positive. Got {}", step));
        }

        if pos < *positions.last().unwrap() {
            return Err(format!(
                "Positions must be monotonically increasing. Got {} after {}",
                pos,
                positions.last().unwrap()
            ));
        }

        spacings.push(step);
        positions.push(pos);
    }

    if values.len() % 2 == 0 {
        return Err("Interval specification must have an odd number of values".to_string());
    }

    Ok((positions, spacings))
}

fn parse_euler_convention(s: &str) -> Result<EulerConvention, String> {
    match s.to_uppercase().as_str() {
        "XYZ" => Ok(EulerConvention::XYZ),
        "XZY" => Ok(EulerConvention::XZY),
        "YXZ" => Ok(EulerConvention::YXZ),
        "YZX" => Ok(EulerConvention::YZX),
        "ZXY" => Ok(EulerConvention::ZXY),
        "ZYX" => Ok(EulerConvention::ZYX),
        "XYX" => Ok(EulerConvention::XYX),
        "XZX" => Ok(EulerConvention::XZX),
        "YXY" => Ok(EulerConvention::YXY),
        "YZY" => Ok(EulerConvention::YZY),
        "ZXZ" => Ok(EulerConvention::ZXZ),
        "ZYZ" => Ok(EulerConvention::ZYZ),
        _ => Err(format!("Invalid Euler convention: '{}'", s)),
    }
}

pub fn update_settings_from_cli(config: &mut Settings) {
    let args = CliArgs::parse();

    if let Some(wavelength) = args.propagation.w {
        config.wavelength = wavelength;
    }
    if let Some(medium) = args.material.ri0 {
        config.medium_refr_index = medium;
    }
    if let Some(particle) = args.material.ri {
        config.particle_refr_index = particle;
    }
    if let Some(geo) = args.material.geo {
        config.geom_name = geo;
    }
    if let Some(mp) = args.propagation.bp {
        config.beam_power_threshold = mp;
    }
    if let Some(maf) = args.propagation.baf {
        config.beam_area_threshold_fac = maf;
    }
    if let Some(cop) = args.propagation.cop {
        config.cutoff = cop;
    }
    if let Some(rec) = args.propagation.rec {
        config.max_rec = rec;
    }
    if let Some(tir) = args.propagation.tir {
        config.max_tir = tir;
    }

    let euler_convention = args.orientation.euler.unwrap_or(DEFAULT_EULER_ORDER);

    if let Some(num_orients) = args.orientation.uniform {
        config.orientation = Orientation {
            scheme: Scheme::Uniform { num_orients },
            euler_convention,
        };
    } else if let Some(eulers) = args.orientation.discrete {
        config.orientation = Orientation {
            scheme: Scheme::Discrete { eulers },
            euler_convention,
        };
    } else if let Some(convention) = args.orientation.euler {
        config.orientation.euler_convention = convention;
    }

    if let Some(custom_path) = &args.binning.custom {
        config.binning = BinningScheme {
            scheme: bins::Scheme::Custom {
                bins: vec![],
                file: Some(custom_path.clone()),
            },
        };
    } else if let Some(simple_bins) = &args.binning.simple {
        if simple_bins.len() == 2 {
            let num_theta = simple_bins[0];
            let num_phi = simple_bins[1];
            config.binning = BinningScheme {
                scheme: bins::Scheme::new_simple(num_theta, num_phi),
            };
        } else {
            eprintln!(
                "Warning: Simple binning requires exactly two values. Using default binning."
            );
        }
    } else if args.binning.interval {
        let mut valid_binning = true;

        let (thetas, theta_spacings) = if let Some(theta_values) = &args.binning.theta {
            match parse_interval_specification(theta_values) {
                Ok(result) => result,
                Err(err) => {
                    eprintln!("Error in theta specification: {}", err);
                    valid_binning = false;
                    (vec![], vec![])
                }
            }
        } else {
            eprintln!("Warning: Interval binning requires --theta parameter.");
            valid_binning = false;
            (vec![], vec![])
        };

        let (phis, phi_spacings) = if let Some(phi_values) = &args.binning.phi {
            match parse_interval_specification(phi_values) {
                Ok(result) => result,
                Err(err) => {
                    eprintln!("Error in phi specification: {}", err);
                    valid_binning = false;
                    (vec![], vec![])
                }
            }
        } else {
            eprintln!("Warning: Interval binning requires --phi parameter.");
            valid_binning = false;
            (vec![], vec![])
        };

        if valid_binning {
            config.binning = BinningScheme {
                scheme: bins::Scheme::Interval {
                    thetas,
                    theta_spacings,
                    phis,
                    phi_spacings,
                },
            };
        } else {
            eprintln!("Warning: Could not create interval binning. Using default binning.");
        }
    }

    if let Some(dir) = args.dir {
        config.directory = dir;
    }

    if let Some(distortion) = args.material.distortion {
        config.distortion = Some(distortion);
    }

    if let Some(fov_factor) = args.fov_factor {
        config.fov_factor = Some(fov_factor);
    }
    if let Some(mapping) = args.mapping {
        config.mapping = mapping;
    }

    if let Some(geom_scale) = args.material.geom_scale {
        if geom_scale.len() != 3 {
            panic!("Geometry scale must have exactly 3 values (x, y, z)");
        } else {
            config.geom_scale = Some(geom_scale);
        }
    }
}
