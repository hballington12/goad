use anyhow::Result;
use clap::Args;
use clap::Parser;
use config::{Config, Environment, File};
use nalgebra::Complex;
use pyo3::prelude::*;
use serde::Deserialize;
use std::env;
use std::path::PathBuf;

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
/// Minimum Distortion factor for the geometry.
pub const MIN_DISTORTION: f32 = 1e-5;

/// Central configuration for electromagnetic scattering simulations.
/// 
/// **Context**: Electromagnetic scattering simulations require numerous physical
/// and numerical parameters that control beam propagation, material properties,
/// geometric discretization, and result analysis. These settings must be
/// configurable through files, environment variables, and command-line arguments
/// while maintaining physical consistency.
/// 
/// **How it Works**: Combines all simulation parameters into a single structure
/// with hierarchical loading from default configuration files, local overrides,
/// environment variables, and command-line arguments. Includes validation to
/// ensure physical plausibility and numerical stability.
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
    pub distortion: Option<f32>,
    #[serde(default = "default_geom_scale")]
    pub geom_scale: Option<Vec<f32>>,
    #[serde(default = "default_directory")]
    pub directory: PathBuf,
    #[serde(default = "default_fov_factor")]
    pub fov_factor: Option<f32>,
}

fn default_scale_factor() -> f32 {
    1.0
}

fn default_geom_scale() -> Option<Vec<f32>> {
    None
}

fn default_fov_factor() -> Option<f32> {
    None
}

fn default_directory() -> PathBuf {
    // Get current directory or default to a new PathBuf if it fails
    let current_dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::new());

    // Find the next available run number by checking existing directories
    let mut run_number = 1;
    let mut run_dir;

    loop {
        let run_name = format!("run{:05}", run_number);
        run_dir = current_dir.join(&run_name);

        if !run_dir.exists() {
            break;
        }

        run_number += 1;

        // Safety check to prevent infinite loops in extreme cases
        if run_number > 99999 {
            eprintln!("Warning: Exceeded maximum run number. Using timestamp instead.");
            let timestamp = chrono::Local::now().format("%Y%m%d_%H%M%S");
            run_dir = current_dir.join(format!("run_{}", timestamp));
            break;
        }
    }

    run_dir
}

#[pymethods]
impl Settings {
    #[new]
    #[pyo3(signature = (
        wavelength = None,
        beam_power_threshold = None, 
        beam_area_threshold_fac = None, 
        cutoff = None, 
        medium_refr_index_re = None, 
        medium_refr_index_im = None, 
        particle_refr_index_re = None, 
        particle_refr_index_im = None, 
        geom_name = None, 
        max_rec = None, 
        max_tir = None, 
        theta_res = None, 
        phi_res = None, 
        euler = None
    ))]
    fn py_new(
        wavelength: Option<f32>,
        beam_power_threshold: Option<f32>,
        beam_area_threshold_fac: Option<f32>,
        cutoff: Option<f32>,
        medium_refr_index_re: Option<f32>,
        medium_refr_index_im: Option<f32>,
        particle_refr_index_re: Option<f32>,
        particle_refr_index_im: Option<f32>,
        geom_name: Option<String>,
        max_rec: Option<i32>,
        max_tir: Option<i32>,
        theta_res: Option<usize>,
        phi_res: Option<usize>,
        euler: Option<Vec<f32>>,
    ) -> Self {
        // Load default settings from config file
        let mut settings = load_config_with_cli(false).expect("Failed to load config");

        // Override with provided parameters if they exist
        if let Some(w) = wavelength {
            settings.wavelength = w;
        }

        if let Some(bpt) = beam_power_threshold {
            settings.beam_power_threshold = bpt;
        }

        if let Some(batf) = beam_area_threshold_fac {
            settings.beam_area_threshold_fac = batf;
        }

        if let Some(c) = cutoff {
            settings.cutoff = c;
        }

        // Update medium refractive index if either component is provided
        if medium_refr_index_re.is_some() || medium_refr_index_im.is_some() {
            let re = medium_refr_index_re.unwrap_or(settings.medium_refr_index.re);
            let im = medium_refr_index_im.unwrap_or(settings.medium_refr_index.im);
            settings.medium_refr_index = Complex::new(re, im);
        }

        // Update particle refractive index if either component is provided
        if particle_refr_index_re.is_some() || particle_refr_index_im.is_some() {
            let re = particle_refr_index_re
                .unwrap_or(settings.particle_refr_index.first().map_or(1.0, |c| c.re));
            let im = particle_refr_index_im
                .unwrap_or(settings.particle_refr_index.first().map_or(0.0, |c| c.im));
            settings.particle_refr_index = vec![Complex::new(re, im)];
        }

        if let Some(name) = geom_name {
            settings.geom_name = name;
        }

        if let Some(rec) = max_rec {
            settings.max_rec = rec;
        }

        if let Some(tir) = max_tir {
            settings.max_tir = tir;
        }

        // Update binning if either theta_res or phi_res is provided
        if theta_res.is_some() || phi_res.is_some() {
            // Get the current values or extract from the current binning scheme
            let current_theta = match &settings.binning.scheme {
                bins::Scheme::Simple { num_theta, .. } => *num_theta,
                _ => 180, // Default if not a simple scheme
            };

            let current_phi = match &settings.binning.scheme {
                bins::Scheme::Simple { num_phi, .. } => *num_phi,
                _ => 360, // Default if not a simple scheme
            };

            settings.binning = BinningScheme {
                scheme: bins::Scheme::Simple {
                    num_theta: theta_res.unwrap_or(current_theta),
                    num_phi: phi_res.unwrap_or(current_phi),
                },
            };
        }

        // Update orientation if euler is provided
        if let Some(e) = euler {
            if e.len() >= 3 {
                settings.orientation = Orientation {
                    scheme: Scheme::Discrete {
                        eulers: vec![Euler::new(e[0], e[1], e[2])],
                    },
                    euler_convention: settings.orientation.euler_convention,
                };
            }
        }

        settings
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
    /// Computes the absolute beam area threshold for simulation termination.
    /// 
    /// **Context**: Geometric optics becomes invalid for beam cross-sections smaller
    /// than approximately one wavelength squared. This threshold prevents numerical
    /// instabilities and maintains physical validity by terminating sub-wavelength beams.
    /// 
    /// **How it Works**: Scales the wavelength-squared by the user-specified factor
    /// and geometry scaling to get the absolute area threshold in geometry units.
    pub fn beam_area_threshold(&self) -> f32 {
        self.wavelength * self.wavelength * self.beam_area_threshold_fac * self.scale.powi(2)
    }
}

/// Loads the default configuration from the project's default.toml file.
/// 
/// **Context**: Default configurations provide baseline parameter sets that
/// ensure the simulation runs with physically reasonable values. These defaults
/// serve as starting points for user customization.
/// 
/// **How it Works**: Locates the project root directory and loads the default
/// configuration file, applying validation to ensure parameter consistency.
pub fn load_default_config() -> Result<Settings> {
    let goad_dir = retrieve_project_root()?;
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

/// Loads configuration with full hierarchy including CLI argument processing.
/// 
/// **Context**: Production simulations require loading configurations from
/// multiple sources with proper precedence ordering. This provides the complete
/// configuration loading workflow.
/// 
/// **How it Works**: Calls load_config_with_cli with CLI processing enabled.
pub fn load_config() -> Result<Settings> {
    load_config_with_cli(true)
}

/// Loads configuration with optional CLI argument integration.
/// 
/// **Context**: Different use cases require different levels of configuration
/// processing. Library usage may skip CLI processing while command-line tools
/// need full argument parsing and override capability.
/// 
/// **How it Works**: Loads configuration from files and environment variables,
/// optionally applies CLI argument overrides, and validates the final configuration.
pub fn load_config_with_cli(apply_cli_updates: bool) -> Result<Settings> {
    let config_file = get_config_file()?;

    let settings: Config = Config::builder()
        .add_source(File::from(config_file).required(true))
        .add_source(Environment::with_prefix("goad"))
        .build()
        .unwrap_or_else(|err| {
            eprintln!("Error loading configuration: {}", err);
            std::process::exit(1);
        });

    // println!("config: {:#?}", settings);

    let mut config: Settings = settings.try_deserialize().unwrap_or_else(|err| {
        eprintln!("Error deserializing configuration: {}", err);
        std::process::exit(1);
    });

    if apply_cli_updates {
        update_settings_from_cli(&mut config);
    }

    validate_config(&config);

    // println!("{:#?}", config);

    Ok(config)
}

/// Updates configuration settings from parsed command-line arguments.
/// 
/// **Context**: Command-line arguments provide the highest precedence in the
/// configuration hierarchy, allowing users to override file-based settings
/// for individual simulation runs.
/// 
/// **How it Works**: Parses CLI arguments and selectively updates configuration
/// fields, handling complex nested structures like orientation schemes and
/// binning configurations with validation.
fn update_settings_from_cli(config: &mut Settings) {
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

    // Store the Euler convention to use (default or user-specified)
    let euler_convention = args.orientation.euler.unwrap_or(DEFAULT_EULER_ORDER);

    // Handle orientation schemes
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
        // If only the convention is specified but no orientation scheme,
        // just update the convention on the existing scheme
        config.orientation.euler_convention = convention;
    }

    // Handle binning scheme
    if let Some(custom_path) = &args.binning.custom {
        // Custom binning scheme from file takes precedence over other binning options
        config.binning = BinningScheme {
            scheme: bins::Scheme::Custom {
                bins: vec![], // Empty vector, will be filled from file at runtime
                file: Some(custom_path.clone()),
            },
        };
    } else if let Some(simple_bins) = &args.binning.simple {
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

    // Handle output directory if specified
    if let Some(dir) = args.dir {
        config.directory = dir;
    }

    // Distortion
    if let Some(distortion) = args.material.distortion {
        config.distortion = Some(distortion);
    }

    // Field of view factor
    if let Some(fov_factor) = args.fov_factor {
        config.fov_factor = Some(fov_factor);
    }

    if let Some(geom_scale) = args.material.geom_scale {
        if geom_scale.len() != 3 {
            panic!("Geometry scale must have exactly 3 values (x, y, z)");
        } else {
            config.geom_scale = Some(geom_scale);
        }
    }
}

fn get_config_file() -> Result<PathBuf, anyhow::Error> {
    let current_dir_config = std::env::current_dir()
        .map(|dir| dir.join("local.toml"))
        .unwrap();
    let config_file = if current_dir_config.exists() {
        // println!(
        //     "Using current directory configuration: {:?}",
        //     current_dir_config
        // );
        current_dir_config
    } else {
        // then check local config file, then use default
        let goad_dir = retrieve_project_root()?;
        let default_config_file = goad_dir.join("config/default.toml");
        let local_config = goad_dir.join("config/local.toml");
        // println!("current_dir_config: {:?}", current_dir_config);

        if local_config.exists() {
            println!("Using local configuration: {:?}", local_config);
            local_config
        } else {
            println!("Using default configuration: {:?}", default_config_file);
            default_config_file
        }
    };
    Ok(config_file)
}

/// Locates the project root directory for configuration file access.
/// 
/// **Context**: Configuration files are stored relative to the project root,
/// but the executable may be located in different directories depending on
/// the deployment method (cargo run, installed binary, development build).
/// 
/// **How it Works**: Uses multiple fallback strategies to locate the project
/// root: CARGO_MANIFEST_DIR for development, GOAD_ROOT_DIR for deployment,
/// and directory traversal to find the nearest config directory.
fn retrieve_project_root() -> Result<std::path::PathBuf> {
    if let Ok(manifest_dir) = env::var("CARGO_MANIFEST_DIR") {
        // When running through cargo (e.g. cargo run, cargo test)
        Ok(std::path::PathBuf::from(manifest_dir))
    } else if let Ok(path) = env::var("GOAD_ROOT_DIR") {
        // Allow explicit configuration via environment variable
        Ok(std::path::PathBuf::from(path))
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
            Ok(current_dir)
        } else {
            Err(anyhow::anyhow!("Could not find project root directory"))
        }
    }
}

/// Validates configuration parameters for physical plausibility.
/// 
/// **Context**: Invalid configuration parameters can cause simulation failures,
/// numerical instabilities, or unphysical results. Validation catches common
/// parameter errors before expensive computations begin.
/// 
/// **How it Works**: Applies assertion checks for critical parameters like
/// positive wavelengths and reasonable threshold values.
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
    \x1b[36mgoad --discrete=\"-30.0,20.0,1.0 -40.0,13.0,12.1\" --interval \\\x1b[0m
    \x1b[36m     --theta 0 1 10 2 180 --phi 0 2 180\x1b[0m

    \x1b[32m# Run inside a medium other than air\x1b[0m
    \x1b[36mgoad --ri0 1.5+0.0i\x1b[0m

    \x1b[32m# Run with multiple shapes with different refractive indices\x1b[0m
    \x1b[36mgoad --ri 1.31+0.0i 1.5+0.1i --geo geometries.obj\x1b[0m
    
    \x1b[32m# Save output to a specific directory\x1b[0m
    \x1b[36mgoad --dir /path/to/output\x1b[0m
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

    /// Output directory for simulation results.
    /// If not specified, a directory in the format 'run00001' will be created automatically.
    #[arg(long)]
    pub dir: Option<PathBuf>,

    /// Set the field of view truncation factor for diffraction of beams.
    /// Beams outside an angle of lambda/d * fov_factor will be truncated,
    /// where d is the maximum dimension of the aperture
    #[arg(long)]
    pub fov_factor: Option<f32>,
}

/// Command-line arguments controlling beam propagation behavior.
/// 
/// **Context**: Beam propagation parameters directly affect simulation accuracy,
/// computational cost, and convergence properties. These parameters need
/// accessible command-line control for parameter studies and optimization.
/// 
/// **How it Works**: Groups propagation-related CLI arguments with descriptive
/// help text and appropriate validation constraints.
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

/// Command-line arguments for material properties and geometry specification.
/// 
/// **Context**: Material properties and geometric parameters define the physical
/// system being simulated. These fundamental parameters require easy modification
/// for comparing different materials and particle shapes.
/// 
/// **How it Works**: Groups material and geometry CLI arguments with appropriate
/// parsing for complex numbers and file paths.
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

    /// Distortion factor for the geometry.
    /// Applies distortion sampled from a Gaussian distribution.
    /// Default: sigma = 0.0 (no distortion).
    /// Sigma is the standard deviation of the facet theta tilt (in radians).
    #[arg(long)]
    pub distortion: Option<f32>,

    /// Geometry scale factors for each axis (x, y, z).
    /// Format: "x y z" (e.g., "1.0 1.0 1.0").
    /// Default: "1.0 1.0 1.0" (no scaling).
    #[arg(long, value_parser, num_args = 1..=3, value_delimiter = ' ')]
    pub geom_scale: Option<Vec<f32>>,
}

/// Command-line arguments for particle orientation specification.
/// 
/// **Context**: Particle orientation significantly affects scattering patterns.
/// Different applications require fixed orientations, random orientation
/// averaging, or discrete orientation sets for systematic studies.
/// 
/// **How it Works**: Provides mutually exclusive orientation schemes with
/// appropriate parsing for Euler angles and orientation conventions.
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

    /// Specify Euler angle convention for orientation.
    /// Valid values: XYZ, XZY, YXZ, YZX, ZXY, ZYX, etc.
    /// Default: ZYZ
    #[arg(long, value_parser = parse_euler_convention)]
    pub euler: Option<EulerConvention>,
}

/// Command-line arguments for angular binning scheme configuration.
/// 
/// **Context**: Far-field scattering analysis requires discretization of the
/// observation sphere into angular bins. Different applications need different
/// angular resolutions and sampling patterns.
/// 
/// **How it Works**: Provides multiple binning schemes with validation for
/// interval specifications and custom file formats.
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

    /// Path to custom binning scheme file.
    /// Contains a list of (theta, phi) bin pairs in TOML format.
    /// Overrides other binning parameters.
    #[arg(long)]
    pub custom: Option<String>,
}

/// Parses command-line Euler angle specification into structured format.
/// 
/// **Context**: Euler angles from command-line arguments arrive as comma-separated
/// strings but need validation and conversion to structured types for use
/// in rotation calculations.
/// 
/// **How it Works**: Splits the input string on commas, validates that exactly
/// three angles are provided, and parses each component as a floating-point value.
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

/// Parses variable-spacing interval specification for angular binning.
/// 
/// **Context**: Variable-spacing angular grids require specification of both
/// breakpoints and spacing values. The command-line format must be intuitive
/// while supporting arbitrary numbers of intervals.
/// 
/// **How it Works**: Interprets alternating step/position values to build
/// position and spacing vectors, validating monotonicity and positive step sizes.
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

/// Parses Euler angle convention string into enumerated type.
/// 
/// **Context**: Euler angle conventions affect rotation sequence interpretation.
/// Command-line specification needs validation against supported conventions
/// with clear error messages for invalid inputs.
/// 
/// **How it Works**: Case-insensitive string matching against all supported
/// Euler conventions with descriptive error messages listing valid options.
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
        _ => Err(format!("Invalid Euler convention: '{}'. Valid values are: XYZ, XZY, YXZ, YZX, ZXY, ZYX, XYX, XZX, YXY, YZY, ZXZ, ZYZ", s)),
    }
}
