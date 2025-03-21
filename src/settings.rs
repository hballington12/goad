use anyhow::Result;
use clap::Args;
use clap::Parser;
use config::{Config, Environment, File};
use nalgebra::Complex;
use pyo3::prelude::*;
use serde::Deserialize;
use std::env;
use std::fmt;

use crate::bins;
use crate::orientation::Euler;
use crate::{bins::BinningScheme, orientation::*};

/// Minimum distance for vertices to be considered the same.
pub const VERTEX_MERGE_DISTANCE: f32 = 0.001;
/// Scaling factor for integer coordinates during clipping.
pub const CLIP_TOLERANCE: f32 = 1e6;
/// Minimum absolute value of the dot product of two vectors to be considered colinear.
pub const COLINEAR_THRESHOLD: f32 = 0.001;
/// Minimum vector length (in geometry units) to be considered non-degenerate.
pub const VEC_LENGTH_THRESHOLD: f32 = 0.01;
/// Minimum distance traversed by ray to intersection. Intersections closer than this are ignored.
pub const RAYCAST_MINIMUM_DISTANCE: f32 = 0.01;
/// Distance to far-field in multiples of wavelength. 1e3 - 1e5 is a good range for single precision arithmetic.
pub const RADIUS: f32 = 1e4;
/// Tolerance for diffraction computations, used to avoid divide by zero errors.
pub const DIFF_EPSILON: f32 = 1e-3;
/// Minimum dx or dy in diffraction computation.
pub const DIFF_DMIN: f32 = 1e-5;
/// Tolerance for kxx or kyy in diffraction computation.
pub const KXY_EPSILON: f32 = 1e-3;
/// Small perturbation for propagation distance to reduce errors in diffraction
pub const PROP_PERTURBATION: f32 = 1e-5;
/// Default Euler angle order for the discrete orientation scheme.
pub const DEFAULT_EULER_ORDER: EulerConvention = EulerConvention::ZYZ;

/// Runtime configuration for the application.
#[pyclass]
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Settings {
    pub wavelength: f32,
    pub beam_power_threshold: f32,
    pub beam_area_threshold_fac: f32,
    pub cutoff: f32,
    pub medium_refr_index: Complex<f32>,
    pub particle_refr_index: Vec<Complex<f32>>,
    pub orientation: Orientation,
    pub geom_name: String,
    pub max_rec: i32,
    pub max_tir: i32,
    pub binning: BinningScheme,
    pub seed: Option<u64>,
    #[serde(default = "default_scale_factor")]
    pub scale: f32,
}

fn default_scale_factor() -> f32 {
    1.0
}

#[pymethods]
impl Settings {
    #[new]
    fn py_new(
        wavelength: f32,
        beam_power_threshold: f32,
        beam_area_threshold_fac: f32,
        cutoff: f32,
        medium_refr_index_re: f32,
        medium_refr_index_im: f32,
        particle_refr_index_re: f32,
        particle_refr_index_im: f32,
        geom_name: String,
        max_rec: i32,
        max_tir: i32,
        theta_res: usize,
        phi_res: usize,
        euler: Vec<f32>,
    ) -> Self {
        let medium_refr_index = Complex::new(medium_refr_index_re, medium_refr_index_im);
        let particle_refr_index =
            vec![Complex::new(particle_refr_index_re, particle_refr_index_im)];
        let orientation: Orientation = Orientation {
            scheme: Scheme::Discrete {
                eulers: vec![Euler::new(euler[0], euler[1], euler[2])],
            },
            euler_convention: EulerConvention::XYZ,
        };
        let binning = BinningScheme {
            scheme: bins::Scheme::Simple {
                num_theta: theta_res,
                num_phi: phi_res,
            },
        };
        Settings {
            wavelength,
            beam_power_threshold,
            beam_area_threshold_fac,
            cutoff,
            medium_refr_index,
            particle_refr_index,
            orientation,
            geom_name,
            max_rec,
            max_tir,
            binning,
            seed: None,
            scale: 1.0,
        }
    }

    /// Set the euler angles
    #[setter]
    fn set_euler(&mut self, euler: Vec<f32>) {
        self.orientation = Orientation {
            scheme: Scheme::Discrete {
                eulers: vec![Euler::new(euler[0], euler[1], euler[2])],
            },
            euler_convention: EulerConvention::XYZ,
        };
    }

    /// Get the euler angle, assuming the orientation scheme is discrete
    #[getter]
    fn get_euler(&self) -> Vec<f32> {
        match &self.orientation.scheme {
            Scheme::Discrete { eulers } => vec![eulers[0].alpha, eulers[0].beta, eulers[0].gamma],
            _ => vec![0.0, 0.0, 0.0],
        }
    }
}

impl Settings {
    pub fn beam_area_threshold(&self) -> f32 {
        self.wavelength * self.wavelength * self.beam_area_threshold_fac * self.scale.powi(2)
    }
}

pub fn load_default_config() -> Result<Settings> {
    let goad_dir = retrieve_project_root();
    let default_config_file = goad_dir.join("config/default.toml");

    let settings: Config = Config::builder()
        .add_source(File::from(default_config_file).required(true))
        .build()
        .unwrap_or_else(|err| {
            eprintln!("Error loading configuration: {}", err);
            std::process::exit(1);
        });

    let config: Settings = settings.try_deserialize().unwrap_or_else(|err| {
        eprintln!("Error deserializing configuration: {}", err);
        std::process::exit(1);
    });

    validate_config(&config);

    Ok(config)
}

pub fn load_config() -> Result<Settings> {
    // Try to find the project directory in different ways
    let goad_dir = retrieve_project_root();

    let default_config_file = goad_dir.join("config/default.toml");
    let local_config = goad_dir.join("config/local.toml");

    // Check if local config exists, if not use default
    let config_file = if local_config.exists() {
        // println!("Using local configuration: {:?}", local_config);
        local_config
    } else {
        // println!("Using default configuration: {:?}", default_config_file);
        default_config_file
    };

    let settings: Config = Config::builder()
        .add_source(File::from(config_file).required(true))
        .add_source(Environment::with_prefix("goad"))
        .build()
        .unwrap_or_else(|err| {
            eprintln!("Error loading configuration: {}", err);
            std::process::exit(1);
        });

    // println!("config: {:#?}", default_settings);

    let mut config: Settings = settings.try_deserialize().unwrap_or_else(|err| {
        eprintln!("Error deserializing configuration: {}", err);
        std::process::exit(1);
    });

    // Parse command-line arguments and override values
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

    // Handle orientation schemes
    if let Some(num_orients) = args.orientation.uniform {
        config.orientation = Orientation {
            scheme: Scheme::Uniform { num_orients },
            euler_convention: DEFAULT_EULER_ORDER,
        };
    } else if let Some(eulers) = args.orientation.discrete {
        config.orientation = Orientation {
            scheme: Scheme::Discrete { eulers },
            euler_convention: DEFAULT_EULER_ORDER,
        };
    }

    // Handle binning scheme
    if let Some(simple_bins) = &args.binning.simple {
        if simple_bins.len() == 2 {
            let num_theta = simple_bins[0];
            let num_phi = simple_bins[1];
            config.binning = BinningScheme {
                scheme: bins::Scheme::Simple { num_theta, num_phi },
            };
        } else {
            eprintln!(
                "Warning: Simple binning requires exactly two values. Using default binning."
            );
        }
    } else if args.binning.interval {
        let mut valid_binning = true;

        // Parse theta intervals
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

        // Parse phi intervals
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

    validate_config(&config);

    println!("{:#?}", config);

    Ok(config)
}

/// Retrieve the project root directory.
/// This function tries to find the project root directory in different ways:
/// 1. If the CARGO_MANIFEST_DIR environment variable is set, use it.
/// 2. If the GOAD_ROOT_DIR environment variable is set, use it.
/// 3. If the "config" subdirectory is found in the exectuable directory or any of its parents, use it.
/// If none of these methods work, the function will panic.
fn retrieve_project_root() -> std::path::PathBuf {
    let goad_dir = if let Ok(manifest_dir) = env::var("CARGO_MANIFEST_DIR") {
        // When running through cargo (e.g. cargo run, cargo test)
        std::path::PathBuf::from(manifest_dir)
    } else if let Ok(path) = env::var("GOAD_ROOT_DIR") {
        // Allow explicit configuration via environment variable
        std::path::PathBuf::from(path)
    } else {
        // Fallback: try to find the nearest directory containing a "config" subdirectory
        // Start from the executable directory and walk upward
        let exe_path = env::current_exe().expect("Failed to get current executable path");
        let mut current_dir = exe_path
            .parent()
            .expect("Failed to get executable directory")
            .to_path_buf();
        let mut found = false;

        while !found && current_dir.parent().is_some() {
            if current_dir.join("config").is_dir() {
                found = true;
            } else {
                current_dir = current_dir.parent().unwrap().to_path_buf();
            }
        }

        if found {
            current_dir
        } else {
            panic!("Could not find project root directory");
        }
    };
    goad_dir
}

fn validate_config(config: &Settings) {
    assert!(
        config.beam_area_threshold_fac > 1e-5,
        "Beam area threshold factor must be greater than 1e-5"
    );
    assert!(config.wavelength > 0.0, "Wavelength must be greater than 0");
}

#[derive(Parser, Debug)]
#[command(version, about = "GOAD - Geometric Optics with Aperture Diffraction")]
#[command(author = "Harry Ballington")]
#[command(help_template = "{about}\n{author}\n\nUsage: {usage}\n\n{all-args}{after-help}")]
#[command(after_help = "\x1b[1;36mEXAMPLES:\x1b[0m
    \x1b[32m# Run with a specific wavelength and geometry file\x1b[0m
    \x1b[36mgoad -w 0.5 --geo geometry.obj\x1b[0m
    
    \x1b[32m# Run with a specific refractive index and random orientations\x1b[0m
    \x1b[36mgoad --ri 1.31+0.01i --uniform 100\x1b[0m
    
    \x1b[32m# Run over discrete orientations with an interval binning scheme\x1b[0m
    \x1b[36mgoad --discrete 0.0,0.0,0.0 90.0,0.0,0.0 --interval \\\x1b[0m
    \x1b[36m     --theta 0 1 10 2 180 --phi 0 2 180\x1b[0m

    \x1b[32m# Run inside a medium other than air\x1b[0m
    \x1b[36mgoad --ri0 1.5+0.0i\x1b[0m

    \x1b[32m# Run with multiple shapes with different refractive indices\x1b[0m
    \x1b[36mgoad --ri 1.31+0.0i 1.5+0.1i --geo geometries.obj\x1b[0m
    ")]
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
    /// Omit for a randomized seed.
    #[arg(short, long)]
    pub seed: Option<u64>,
}

/// Beam propagation parameters - control how beams are traced through the geometry
#[derive(Args, Debug)]
pub struct PropagationArgs {
    /// Wavelength in units of the geometry.
    /// Should be larger than the smallest feature in the geometry.
    #[arg(short, long)]
    pub w: Option<f32>,

    /// Minimum beam power threshold for propagation.
    /// Beams with less power than this will be truncated.
    #[arg(long)]
    pub bp: Option<f32>,

    /// Minimum area factor threshold for beam propagation.
    /// The actual area threshold is wavelength² × factor.
    /// Prevents geometric optics from modeling sub-wavelength beams.
    #[arg(long)]
    pub baf: Option<f32>,

    /// Total power cutoff fraction (0.0-1.0).
    /// Simulation stops when this fraction of input power is accounted for.
    /// Set to 1.0 to disable and trace all beams to completion.
    #[arg(long)]
    pub cop: Option<f32>,

    /// Maximum recursion depth for beam tracing.
    /// Typical values: 8-15. Higher values rarely improve results
    /// when reasonable beam power thresholds are set.
    #[arg(long)]
    pub rec: Option<i32>,

    /// Maximum allowed total internal reflections.
    /// Prevents infinite TIR loops by truncating beams
    /// after this many TIR events.
    #[arg(long)]
    pub tir: Option<i32>,
}

/// Material and geometry parameters - control the physical properties of the simulation
#[derive(Args, Debug)]
pub struct MaterialArgs {
    /// Path to geometry file (.obj format).
    /// Contains all input shapes for the simulation.
    #[arg(short, long)]
    pub geo: Option<String>,

    /// Surrounding medium refractive index.
    /// Format: "re+im" (e.g., "1.3117+0.0001i").
    #[arg(long)]
    pub ri0: Option<Complex<f32>>,

    /// Particle refractive indices, space-separated.
    /// Each shape in the geometry is assigned a refractive index.
    /// If fewer values than shapes are provided, the first value is reused.
    #[arg(short, long, value_parser, num_args = 1.., value_delimiter = ' ')]
    pub ri: Option<Vec<Complex<f32>>>,
}

/// Orientation parameters - control how the particle is oriented relative to the incident beam
#[derive(Args, Debug)]
pub struct OrientationArgs {
    /// Use uniform random orientation scheme.
    /// The value specifies the number of random orientations.
    #[arg(long, group = "orientation")]
    pub uniform: Option<usize>,

    /// Use discrete orientation scheme with specified Euler angles (degrees).
    /// Format: alpha1,beta1,gamma1 alpha2,beta2,gamma2 ...
    #[arg(long, value_parser = parse_euler_angles, num_args = 1.., value_delimiter = ' ', group = "orientation")]
    pub discrete: Option<Vec<Euler>>,
}

/// Output binning parameters - control how scattered light is binned by angle
#[derive(Args, Debug)]
pub struct BinningArgs {
    /// Use simple equal-spacing binning scheme.
    /// Format: <num_theta_bins> <num_phi_bins>
    #[arg(long, num_args = 2, value_delimiter = ' ', group = "binning")]
    pub simple: Option<Vec<usize>>,

    /// Enable interval binning scheme with variable spacing.
    /// Allows fine binning in regions of interest like forward/backward scattering.
    #[arg(long, group = "binning")]
    pub interval: bool,

    /// Theta angle bins for interval binning (degrees).
    /// Format: start step1 mid1 step2 mid2 ... stepN end
    /// Example: 0 1 10 2 180 = 0° to 10° in 1° steps, then 10° to 180° in 2° steps
    #[arg(long, requires = "interval", num_args = 3.., value_delimiter = ' ')]
    pub theta: Option<Vec<f32>>,

    /// Phi angle bins for interval binning (degrees).
    /// Format: start step1 mid1 step2 mid2 ... stepN end
    /// Example: 0 2 180 = 0° to 180° in 2° steps
    #[arg(long, requires = "interval", num_args = 3.., value_delimiter = ' ')]
    pub phi: Option<Vec<f32>>,
}

/// Parse a string of Euler angles in the format "alpha,beta,gamma"
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

/// Parse interval specification in the format:
/// start step1 mid1 step2 mid2 ... stepN end
/// This returns two vectors:
/// 1. positions: [start, mid1, mid2, ..., end]
/// 2. spacings: [step1, step2, ..., stepN]
fn parse_interval_specification(values: &[f32]) -> Result<(Vec<f32>, Vec<f32>), String> {
    if values.len() < 3 {
        return Err(format!(
            "Insufficient values for interval specification: need at least 3, got {}",
            values.len()
        ));
    }

    // For values [start, step1, mid1, step2, mid2, ..., stepN, end],
    // we need to extract positions and spacings
    let mut positions = Vec::new();
    let mut spacings = Vec::new();

    // Add the start position
    positions.push(values[0]);

    // Process pairs of (step, position) until the last value
    for i in (1..values.len() - 1).step_by(2) {
        let step = values[i];
        let pos = values[i + 1];

        // Validate step
        if step <= 0.0 {
            return Err(format!("Step size must be positive. Got {}", step));
        }

        // Validate monotonicity
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

    // Check if there are an odd number of values (required for valid format)
    if values.len() % 2 == 0 {
        return Err("Interval specification must have an odd number of values".to_string());
    }

    Ok((positions, spacings))
}

impl fmt::Display for Settings {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Settings:
  - Wavelength: {:.6}
  - Beam Power Threshold: {:.6}
  - Total Power Cutoff: {:.6}
  - Medium Refractive Index: {:.6} + {:.6}i
  - Particle Refractive Indices: {:?}
  - Max Rec: {}
  - Max TIR: {}
  ",
            self.wavelength,
            self.beam_power_threshold,
            self.cutoff,
            self.medium_refr_index.re,
            self.medium_refr_index.im,
            self.particle_refr_index,
            self.max_rec,
            self.max_tir,
        )
    }
}
