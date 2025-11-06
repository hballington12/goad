pub mod cli;
pub mod constants;
pub mod loading;
pub mod validation;

use nalgebra::Complex;
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use crate::bins;
use crate::diff::Mapping;
use crate::orientation::Euler;
use crate::{bins::BinningScheme, orientation::*};

/// Configuration for output file generation
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
pub struct OutputConfig {
    /// Enable writing of results.dat summary file
    pub results_summary: bool,
    /// Enable writing of settings.json configuration file
    pub settings_json: bool,
    /// Enable writing of powers.json power distribution file
    pub powers_json: bool,
    /// Enable writing of params.json derived parameters file
    pub params_json: bool,
    /// Enable writing of 2D Mueller matrix files
    pub mueller_2d: bool,
    /// Enable writing of 1D integrated Mueller matrix files
    pub mueller_1d: bool,
    /// Enable writing of specific Mueller components
    pub mueller_components: MuellerComponentConfig,
}

/// Configuration for Mueller matrix component outputs
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
pub struct MuellerComponentConfig {
    /// Enable total scattering component output
    pub total: bool,
    /// Enable beam component output
    pub beam: bool,
    /// Enable external diffraction component output
    pub external: bool,
}

// Re-export constants, defaults, and loading functions for backward compatibility
pub use self::constants::*;
pub use self::loading::{load_config, load_config_with_cli, load_default_config};

/// Runtime configuration for the application.
#[pyclass]
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
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
    /// Problem scaling factor - scales the entire problem (geometry, wavelength, and beam area thresholds)
    #[serde(default = "constants::default_scale_factor")]
    pub scale: f32,
    pub distortion: Option<f32>,
    /// Per-axis geometry scaling [x, y, z] - scales only the geometry in each dimension
    #[serde(default = "constants::default_geom_scale")]
    pub geom_scale: Option<Vec<f32>>,
    #[serde(default = "constants::default_directory")]
    pub directory: PathBuf,
    #[serde(default = "constants::default_fov_factor")]
    pub fov_factor: Option<f32>,
    pub mapping: Mapping,
    #[serde(default = "constants::default_output_config")]
    pub output: OutputConfig,
    pub coherence: bool,
    /// Suppress progress bars and status messages
    #[serde(default = "constants::default_quiet")]
    pub quiet: bool,
}

#[pymethods]
impl Settings {
    #[new]
    #[pyo3(signature = (
        geom_path,
        wavelength = DEFAULT_WAVELENGTH,
        particle_refr_index_re = DEFAULT_PARTICLE_REFR_INDEX_RE,
        particle_refr_index_im = DEFAULT_PARTICLE_REFR_INDEX_IM,
        medium_refr_index_re = DEFAULT_MEDIUM_REFR_INDEX_RE,
        medium_refr_index_im = DEFAULT_MEDIUM_REFR_INDEX_IM,
        orientation = None,
        binning = None,
        beam_power_threshold = DEFAULT_BEAM_POWER_THRESHOLD,
        beam_area_threshold_fac = DEFAULT_BEAM_AREA_THRESHOLD_FAC,
        cutoff = DEFAULT_CUTOFF,
        max_rec = DEFAULT_MAX_REC,
        max_tir = DEFAULT_MAX_TIR,
        scale = 1.0,
        directory = "goad_run",
        mapping = DEFAULT_MAPPING,
        coherence = DEFAULT_COHERENCE,
        quiet = DEFAULT_QUIET,
    ))]
    fn py_new(
        geom_path: String,
        wavelength: f32,
        particle_refr_index_re: f32,
        particle_refr_index_im: f32,
        medium_refr_index_re: f32,
        medium_refr_index_im: f32,
        orientation: Option<Orientation>,
        binning: Option<BinningScheme>,
        beam_power_threshold: f32,
        beam_area_threshold_fac: f32,
        cutoff: f32,
        max_rec: i32,
        max_tir: i32,
        scale: f32,
        directory: &str,
        mapping: Mapping,
        coherence: bool,
        quiet: bool,
    ) -> PyResult<Self> {
        // Input validation
        if wavelength <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Wavelength must be positive, got: {}",
                wavelength
            )));
        }

        if !std::path::Path::new(&geom_path).exists() {
            return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
                "Geometry file not found: {}",
                geom_path
            )));
        }

        if cutoff < 0.0 || cutoff > 1.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Cutoff must be between 0 and 1, got: {}",
                cutoff
            )));
        }

        if max_rec < 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "max_rec must be non-negative, got: {}",
                max_rec
            )));
        }

        if max_tir < 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "max_tir must be non-negative, got: {}",
                max_tir
            )));
        }
        // Create default orientation if none provided (single random orientation)
        let orientation = orientation.unwrap_or_else(|| Orientation {
            scheme: Scheme::Uniform { num_orients: 1 },
            euler_convention: DEFAULT_EULER_ORDER,
        });

        // Create default binning if none provided (interval binning with high resolution)
        let binning = binning.unwrap_or_else(|| BinningScheme {
            scheme: bins::Scheme::Interval {
                thetas: vec![0.0, 5.0, 175.0, 179.0, 180.0],
                theta_spacings: vec![0.1, 2.0, 0.5, 0.1],
                phis: vec![0.0, 360.0],
                phi_spacings: vec![7.5],
            },
        });

        let mut settings = Settings {
            wavelength,
            beam_power_threshold,
            beam_area_threshold_fac,
            cutoff,
            medium_refr_index: Complex::new(medium_refr_index_re, medium_refr_index_im),
            particle_refr_index: vec![Complex::new(particle_refr_index_re, particle_refr_index_im)],
            orientation,
            geom_name: geom_path,
            max_rec,
            max_tir,
            binning,
            seed: None,
            scale,
            distortion: None,
            geom_scale: None,
            directory: PathBuf::from(directory),
            fov_factor: None,
            mapping,
            output: constants::default_output_config(),
            coherence,
            quiet,
        };

        validation::validate_config(&mut settings);

        Ok(settings)
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

    /// Set the full orientation object
    #[setter]
    fn set_orientation(&mut self, orientation: Orientation) {
        self.orientation = orientation;
    }

    /// Get the full orientation object
    #[getter]
    fn get_orientation(&self) -> Orientation {
        self.orientation.clone()
    }

    /// Set the geometry file path
    #[setter]
    fn set_geom_path(&mut self, geom_path: String) {
        self.geom_name = geom_path;
    }

    /// Get the geometry file path
    #[getter]
    fn get_geom_path(&self) -> String {
        self.geom_name.clone()
    }

    /// Set the wavelength
    #[setter]
    fn set_wavelength(&mut self, wavelength: f32) {
        self.wavelength = wavelength;
    }

    /// Get the wavelength
    #[getter]
    fn get_wavelength(&self) -> f32 {
        self.wavelength
    }

    /// Set the particle refractive index (real part)
    #[setter]
    fn set_particle_refr_index_re(&mut self, re: f32) {
        if !self.particle_refr_index.is_empty() {
            self.particle_refr_index[0].re = re;
        }
    }

    /// Get the particle refractive index (real part)
    #[getter]
    fn get_particle_refr_index_re(&self) -> f32 {
        if !self.particle_refr_index.is_empty() {
            self.particle_refr_index[0].re
        } else {
            0.0
        }
    }

    /// Set the particle refractive index (imaginary part)
    #[setter]
    fn set_particle_refr_index_im(&mut self, im: f32) {
        if !self.particle_refr_index.is_empty() {
            self.particle_refr_index[0].im = im;
        }
    }

    /// Get the particle refractive index (imaginary part)
    #[getter]
    fn get_particle_refr_index_im(&self) -> f32 {
        if !self.particle_refr_index.is_empty() {
            self.particle_refr_index[0].im
        } else {
            0.0
        }
    }

    /// Set the medium refractive index (real part)
    #[setter]
    fn set_medium_refr_index_re(&mut self, re: f32) {
        self.medium_refr_index.re = re;
    }

    /// Get the medium refractive index (real part)
    #[getter]
    fn get_medium_refr_index_re(&self) -> f32 {
        self.medium_refr_index.re
    }

    /// Set the medium refractive index (imaginary part)
    #[setter]
    fn set_medium_refr_index_im(&mut self, im: f32) {
        self.medium_refr_index.im = im;
    }

    /// Get the medium refractive index (imaginary part)
    #[getter]
    fn get_medium_refr_index_im(&self) -> f32 {
        self.medium_refr_index.im
    }

    /// Set the beam power threshold
    #[setter]
    fn set_beam_power_threshold(&mut self, threshold: f32) {
        self.beam_power_threshold = threshold;
    }

    /// Get the beam power threshold
    #[getter]
    fn get_beam_power_threshold(&self) -> f32 {
        self.beam_power_threshold
    }

    /// Set the cutoff
    #[setter]
    fn set_cutoff(&mut self, cutoff: f32) {
        self.cutoff = cutoff;
    }

    /// Get the cutoff
    #[getter]
    fn get_cutoff(&self) -> f32 {
        self.cutoff
    }

    /// Set the max recursion depth
    #[setter]
    fn set_max_rec(&mut self, max_rec: i32) {
        self.max_rec = max_rec;
    }

    /// Get the max recursion depth
    #[getter]
    fn get_max_rec(&self) -> i32 {
        self.max_rec
    }

    /// Set the max TIR bounces
    #[setter]
    fn set_max_tir(&mut self, max_tir: i32) {
        self.max_tir = max_tir;
    }

    /// Get the max TIR bounces
    #[getter]
    fn get_max_tir(&self) -> i32 {
        self.max_tir
    }

    /// Set the binning scheme
    #[setter]
    fn set_binning(&mut self, binning: BinningScheme) {
        self.binning = binning;
    }

    /// Get the binning scheme
    #[getter]
    fn get_binning(&self) -> BinningScheme {
        self.binning.clone()
    }

    /// Set the per-axis geometry scaling [x, y, z]
    #[setter]
    fn set_geom_scale(&mut self, geom_scale: Option<Vec<f32>>) {
        self.geom_scale = geom_scale;
    }

    /// Get the per-axis geometry scaling [x, y, z]
    #[getter]
    fn get_geom_scale(&self) -> Option<Vec<f32>> {
        self.geom_scale.clone()
    }

    /// Set the seed for random number generation
    #[setter]
    fn set_seed(&mut self, seed: Option<u64>) {
        self.seed = seed;
    }

    /// Get the seed for random number generation
    #[getter]
    fn get_seed(&self) -> Option<u64> {
        self.seed
    }

    /// Set the distortion factor
    #[setter]
    fn set_distortion(&mut self, distortion: Option<f32>) {
        self.distortion = distortion;
    }

    /// Get the distortion factor
    #[getter]
    fn get_distortion(&self) -> Option<f32> {
        self.distortion
    }

    /// Set the field of view factor
    #[setter]
    fn set_fov_factor(&mut self, fov_factor: Option<f32>) {
        self.fov_factor = fov_factor;
    }

    /// Get the field of view factor
    #[getter]
    fn get_fov_factor(&self) -> Option<f32> {
        self.fov_factor
    }

    /// Set quiet mode (suppress progress bars)
    #[setter]
    fn set_quiet(&mut self, quiet: bool) {
        self.quiet = quiet;
    }

    /// Get quiet mode
    #[getter]
    fn get_quiet(&self) -> bool {
        self.quiet
    }
}

impl Settings {
    pub fn beam_area_threshold(&self) -> f32 {
        self.wavelength * self.wavelength * self.beam_area_threshold_fac * self.scale.powi(2)
    }
}
