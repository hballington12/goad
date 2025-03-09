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
    pub total_power_cutoff: f32,
    pub medium_refr_index: Complex<f32>,
    pub particle_refr_index: Vec<Complex<f32>>,
    pub orientation: OrientationScheme,
    pub geom_name: String,
    pub max_rec: i32,
    pub max_tir: i32,
    pub binning: BinningScheme,
    pub seed: Option<u64>,
    #[serde(default = "default_scale_factor")]
    pub scale_factor: f32,
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
        total_power_cutoff: f32,
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
        let orientation: OrientationScheme = OrientationScheme {
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
            total_power_cutoff,
            medium_refr_index,
            particle_refr_index,
            orientation,
            geom_name,
            max_rec,
            max_tir,
            binning,
            seed: None,
            scale_factor: 1.0,
        }
    }

    /// Set the euler angles
    #[setter]
    fn set_euler(&mut self, euler: Vec<f32>) {
        self.orientation = OrientationScheme {
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
        self.wavelength * self.wavelength * self.beam_area_threshold_fac * self.scale_factor.powi(2)
    }
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

    let default_settings: Config = Config::builder()
        .add_source(File::from(config_file).required(true))
        .add_source(Environment::with_prefix("goad"))
        .build()
        .unwrap_or_else(|err| {
            eprintln!("Error loading configuration: {}", err);
            std::process::exit(1);
        });

    // println!("config: {:#?}", default_settings);

    let mut config: Settings = default_settings.try_deserialize().unwrap_or_else(|err| {
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
        config.total_power_cutoff = cop;
    }
    if let Some(rec) = args.rec {
        config.max_rec = rec;
    }
    if let Some(tir) = args.tir {
        config.max_tir = tir;
    }
    if let Some(orient) = args.orient {
        config.orientation = OrientationScheme {
            scheme: orient,
            euler_convention: EulerConvention::ZYZ,
        };
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

    /// Orientation scheme for the simulation.
    #[command(subcommand)]
    orient: Option<Scheme>,

    /// Random seed for the simulation.
    #[arg(short, long)]
    seed: Option<u64>,
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
            self.total_power_cutoff,
            self.medium_refr_index.re,
            self.medium_refr_index.im,
            self.particle_refr_index,
            self.max_rec,
            self.max_tir,
        )
    }
}
