use anyhow::Result;
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
        println!("Using local configuration: {:?}", local_config);
        local_config
    } else {
        println!("Using default configuration: {:?}", default_config_file);
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

    if let Some(wavelength) = args.w {
        config.wavelength = wavelength;
    }
    if let Some(medium) = args.ri0 {
        config.medium_refr_index = medium;
    }
    if let Some(particle) = args.ri {
        config.particle_refr_index = particle;
    }
    if let Some(geo) = args.geo {
        config.geom_name = geo;
    }
    if let Some(mp) = args.bp {
        config.beam_power_threshold = mp;
    }
    if let Some(maf) = args.baf {
        config.beam_area_threshold_fac = maf;
    }
    if let Some(cop) = args.cop {
        config.cutoff = cop;
    }
    if let Some(rec) = args.rec {
        config.max_rec = rec;
    }
    if let Some(tir) = args.tir {
        config.max_tir = tir;
    }

    // Handle orientation schemes
    if let Some(num_orients) = args.uniform {
        config.orientation = Orientation {
            scheme: Scheme::Uniform { num_orients },
            euler_convention: DEFAULT_EULER_ORDER,
        };
    } else if let Some(eulers) = args.discrete {
        config.orientation = Orientation {
            scheme: Scheme::Discrete { eulers },
            euler_convention: DEFAULT_EULER_ORDER,
        };
    }

    // Handle binning scheme
    if let Some(simple_bins) = &args.simple {
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
    } else if args.interval {
        let mut valid_binning = true;

        // Parse theta intervals
        let (thetas, theta_spacings) = if let Some(theta_values) = &args.theta {
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
        let (phis, phi_spacings) = if let Some(phi_values) = &args.phi {
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
pub struct CliArgs {
    /// Wavelength in units of the geometry.
    #[arg(short, long)]
    w: Option<f32>,

    /// Minimum absolute beam power threshold for new beams to propagate.
    #[arg(long)]
    bp: Option<f32>,

    /// Minimum area factor for new beams to propagate. The actual area threshold is
    /// calculated as `wavelength^2 * factor`.
    #[arg(long)]
    baf: Option<f32>,

    /// Cutoff power. The total acceptable output power per orientation before beam propagation is terminated.
    /// Once this threshold is reached, the near-field simulation will stop.
    #[arg(long)]
    cop: Option<f32>,

    /// The maximum number of recursions before a beam is truncated.
    #[arg(long)]
    rec: Option<i32>,

    /// The maximum number of total internal reflections before a beam is truncated.
    #[arg(long)]
    tir: Option<i32>,

    /// File path to the input geometry. All input shapes should be defined in this file.
    /// Currently, only the Wavefront .obj format is supported.
    #[arg(short, long)]
    geo: Option<String>,

    /// The refractive index of the surrounding medium.
    #[arg(long)]
    ri0: Option<Complex<f32>>,

    /// The refractive index of the particle/s, separated by spaces.
    /// If multiple values are provided, each shape in the geometry will be assigned a refractive index.
    /// If fewer values are provided than the number of shapes, the first value will be used for the remaining shapes.
    #[arg(short, long, value_parser, num_args = 1.., value_delimiter = ' ')]
    ri: Option<Vec<Complex<f32>>>,

    /// Random seed for the simulation.
    #[arg(short, long)]
    seed: Option<u64>,

    /// Use uniform orientation scheme with specified number of orientations
    #[arg(long, group = "orientation")]
    uniform: Option<usize>,

    /// Use discrete orientation scheme with specified Euler angles (in degrees)
    /// Format: alpha1,beta1,gamma1 alpha2,beta2,gamma2 ...
    #[arg(long, value_parser = parse_euler_angles, num_args = 1.., value_delimiter = ' ', group = "orientation")]
    discrete: Option<Vec<Euler>>,

    /// Use a simple binning scheme with the specified numbers of theta and phi bins.
    #[arg(long, num_args = 2, value_delimiter = ' ', group = "binning")]
    simple: Option<Vec<usize>>,

    /// Use an interval binning scheme with specified intervals for theta and phi in degrees.
    #[arg(long, group = "binning")]
    interval: bool,

    /// Theta angles for interval binning (must be used with --interval)
    /// Format: start step1 mid1 step2 mid2 ... stepN end
    #[arg(long, requires = "interval", num_args = 3.., value_delimiter = ' ')]
    theta: Option<Vec<f32>>,

    /// Phi angles for interval binning (must be used with --interval)
    /// Format: start step1 mid1 step2 mid2 ... stepN end
    #[arg(long, requires = "interval", num_args = 3.., value_delimiter = ' ')]
    phi: Option<Vec<f32>>,
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
