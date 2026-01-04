use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rand::SeedableRng;

use crate::geom::Geom;
use crate::multiproblem::init_result;
use crate::orientation::{OrientationSampler, Scheme};
use crate::params::Param;
use crate::problem::init_geom;
use crate::result::Results;
use crate::settings::Settings;

use super::{Convergence, ConvergenceTracker, MAX_CONVERGENCE_ORIENTATIONS};

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl Convergence {
    #[new]
    #[pyo3(signature = (settings, geoms = None))]
    fn py_new(settings: Settings, geoms: Option<Vec<Geom>>) -> PyResult<Self> {
        let mut geoms = match geoms {
            Some(g) => g,
            None => Geom::load(&settings.geom_name).map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "Failed to load geometry file '{}': {}\n\
                    Hint: This may be caused by degenerate faces (zero cross product), \
                    faces that are too small, or non-planar geometry. \
                    Please check and fix the geometry file.",
                    settings.geom_name, e
                ))
            })?,
        };

        for geom in geoms.iter_mut() {
            init_geom(&settings, geom);
        }

        let template = init_result(&settings);

        // Create sampler based on scheme setting
        let sampler = match &settings.orientation.scheme {
            Scheme::Uniform { .. } => OrientationSampler::uniform(settings.seed),
            Scheme::Discrete { eulers } => OrientationSampler::discrete(eulers.clone()),
            Scheme::Sobol { .. } => OrientationSampler::sobol(settings.seed),
            Scheme::Halton { .. } => OrientationSampler::halton(),
        };

        let rng = if let Some(seed) = settings.seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };

        Ok(Self {
            geoms,
            settings,
            max_orientations: MAX_CONVERGENCE_ORIENTATIONS,
            targets: Vec::new(),
            tracker: ConvergenceTracker::new(&template),
            sampler,
            rng,
            log_file: None,
        })
    }

    /// Solve the multi-orientation scattering problem using work-stealing.
    /// Periodically checks for Python signals (Ctrl-C) and interrupts if needed.
    #[pyo3(name = "solve")]
    pub fn py_solve(&mut self, py: Python) -> PyResult<()> {
        self.solve_with_interrupt(|| py.check_signals().is_err())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Access the current mean results (live during solve).
    #[getter]
    pub fn get_mean(&self) -> Results {
        self.tracker.mean()
    }

    /// Access the current standard error of the mean (live during solve).
    #[getter]
    pub fn get_sem(&self) -> Results {
        self.tracker.sem()
    }

    /// Get the max orientations (safety cap).
    #[getter]
    pub fn get_count(&self) -> usize {
        self.tracker.count()
    }

    #[getter]
    pub fn get_max_orientations(&self) -> usize {
        self.max_orientations
    }

    /// Set the max orientations (safety cap).
    #[setter]
    pub fn set_max_orientations(&mut self, max_orientations: usize) {
        self.max_orientations = max_orientations;
    }

    /// Add a convergence target for a parameter (Python API).
    /// Solver terminates when ALL targets are satisfied.
    #[pyo3(name = "add_target")]
    pub fn py_add_target(&mut self, param: Param, relative_error: f32) {
        self.add_target(param, relative_error);
    }

    /// Clear all convergence targets (Python API).
    #[pyo3(name = "clear_targets")]
    pub fn py_clear_targets(&mut self) {
        self.clear_targets();
    }

    /// Reset the solver to initial state.
    pub fn py_reset(&mut self) -> PyResult<()> {
        self.reset();
        Ok(())
    }

    /// Reset the orientation sampler (for reproducibility).
    pub fn py_reset_sampler(&mut self) -> PyResult<()> {
        self.reset_sampler();
        Ok(())
    }

    /// Save simulation results to disk.
    ///
    /// Writes Mueller matrices, parameters, and other output files to the
    /// specified directory (or the directory configured in settings).
    ///
    /// Args:
    ///     directory: Optional output directory path. If not provided, uses
    ///                the directory from settings.
    #[pyo3(signature = (directory=None))]
    pub fn save(&self, directory: Option<String>) -> PyResult<()> {
        let mut result = self.mean();
        result.mueller_to_1d();
        let _ = result.compute_params(self.settings.wavelength);

        let settings = if let Some(dir) = directory {
            let mut s = self.settings.clone();
            s.directory = std::path::PathBuf::from(dir);
            s
        } else {
            self.settings.clone()
        };

        let output_manager = crate::output::OutputManager::new(&settings, &result);
        output_manager.write_all().map_err(|e| {
            pyo3::exceptions::PyIOError::new_err(format!("Failed to save results: {}", e))
        })?;
        Ok(())
    }
}
