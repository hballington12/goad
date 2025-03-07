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
/// Distance to far-field. 1e3 - 1e5 is a good range for single precision arithmetic.
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
    ) -> Self {
        let medium_refr_index = Complex::new(medium_refr_index_re, medium_refr_index_im);
        let particle_refr_index =
            vec![Complex::new(particle_refr_index_re, particle_refr_index_im)];
        let orientation: OrientationScheme = OrientationScheme {
            scheme: Scheme::Discrete {
                eulers: vec![Euler::new(0.0, 0.0, 0.0)],
            },
            euler_convention: EulerConvention::ZYZ,
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
        }
    }
}

impl Settings {
    pub fn beam_area_threshold(&self) -> f32 {
        self.wavelength * self.wavelength * self.beam_area_threshold_fac
    }
}

pub fn load_config() -> Result<Settings> {
    let exe_path = env::current_exe().expect("Failed to get current executable path");
    let goad_dir = exe_path
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent())
        .expect("Failed to get target directory.");

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
