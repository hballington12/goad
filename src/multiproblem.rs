// use std::time::Instant;

use crate::{
    convergence::Convergeable,
    geom::Geom,
    orientation::{Euler, Orientations},
    output,
    problem::{self, Problem},
    result::Results,
    settings::Settings,
    zones::Zones,
};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use nalgebra::Complex;
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::time::Duration;

// ============================================================================
// Helper functions for solver initialization
// ============================================================================

/// Loads settings from config if not provided.
pub fn load_settings_or_default(settings: Option<Settings>) -> Settings {
    settings.unwrap_or_else(|| crate::settings::load_config().expect("Failed to load config"))
}

/// Loads and initializes geometries. Used by MultiProblem and Convergence.
pub fn load_and_init_geoms(
    geoms: Option<Vec<Geom>>,
    settings: &Settings,
) -> anyhow::Result<Vec<Geom>> {
    let mut geoms = match geoms {
        Some(g) => g,
        None => Geom::load(&settings.geom_name).map_err(|e| {
            anyhow::anyhow!(
                "Failed to load geometry file '{}': {}\n\
                Hint: This may be caused by degenerate faces (zero cross product), \
                faces that are too small, or non-planar geometry. \
                Please check and fix the geometry file.",
                settings.geom_name,
                e
            )
        })?,
    };

    for geom in geoms.iter_mut() {
        problem::init_geom(settings, geom);
    }

    Ok(geoms)
}

/// Initializes zones and creates an empty Results struct.
pub fn init_result(settings: &Settings) -> Results {
    let zones = Zones::from_configs(&settings.zones);
    Results::new_with_zones(zones)
}

// ============================================================================

/// Multi-orientation light scattering simulation for a single geometry.
///
/// Computes orientation-averaged scattering properties by running multiple
/// single-orientation simulations and averaging the results. Supports both
/// random and systematic orientation sampling schemes. Results include
/// Mueller matrices, cross-sections, and derived optical parameters.
///
/// # Examples
/// ```python
/// import goad_py as goad
///
/// # Create orientation scheme and settings
/// orientations = goad.create_uniform_orientation(100)
/// settings = goad.Settings("particle.obj", orientation=orientations)
///
/// # Run multi-orientation simulation
/// mp = goad.MultiProblem(settings)
/// mp.py_solve()
///
/// # Access averaged results
/// results = mp.results
/// print(f"Scattering cross-section: {results.scat_cross}")
/// ```
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
#[derive(Debug)]
pub struct MultiProblem {
    pub geoms: Vec<Geom>,
    pub orientations: Orientations,
    pub settings: Settings, // runtime settings
    pub result: Results,    // averaged result of the problems
}

impl MultiProblem {
    /// Creates a new `MultiProblem` from optional `Geom` and `Settings`.
    /// If settings not provided, loads from config file.
    /// If geoms not provided, load from file
    pub fn new(geoms: Option<Vec<Geom>>, settings: Option<Settings>) -> anyhow::Result<Self> {
        let settings = load_settings_or_default(settings);
        let geoms = load_and_init_geoms(geoms, &settings)?;
        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let result = init_result(&settings);

        Ok(Self {
            geoms,
            orientations,
            settings,
            result,
        })
    }

    /// Regenerates the orientations for the problem.
    /// Useful for rerunning a random orientation problem with no seed set.
    pub fn regenerate_orientations(&mut self) {
        self.orientations =
            Orientations::generate(&self.settings.orientation.scheme, self.settings.seed);
    }

    /// Resets a `MultiOrientProblem` to its initial state.
    pub fn reset(&mut self) {
        self.result = init_result(&self.settings);
        self.regenerate_orientations();
    }

    /// Solves a `MultiOrientProblem` by averaging over the problems.
    pub fn solve(&mut self) {
        // Initialize file-based logging (avoids conflicts with indicatif progress bar)
        if let Err(e) = crate::filelog::init(&self.settings.directory) {
            log::warn!("Could not initialize file logging: {}", e);
        }

        let n = self.orientations.num_orientations;

        // Initialize progress display only if not in quiet mode
        let (status_pb, pb, info_pb) = if !self.settings.quiet {
            let m = MultiProgress::new();

            // Status spinner (top) - shows current phase
            let status_pb = m.add(ProgressBar::new_spinner());
            status_pb
                .set_style(ProgressStyle::with_template("{spinner:.cyan} Status: {msg}").unwrap());
            status_pb.enable_steady_tick(Duration::from_millis(100));

            // Main progress bar for orientations
            let pb = m.add(ProgressBar::new(n as u64));
            pb.set_style(
                ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.green/blue}] {pos:>5}/{len:5} {msg} | ETA: {eta_precise}",
                )
                .unwrap()
                .progress_chars("█▇▆▅▄▃▂▁")
            );
            pb.set_message("Computing orientations");

            // Info display (bottom) - shows additional context
            let info_pb = m.add(ProgressBar::new_spinner());
            info_pb.set_style(ProgressStyle::with_template("{msg}").unwrap());
            info_pb.enable_steady_tick(Duration::from_millis(500));

            // Phase 1: Initialization
            status_pb.set_message("Initializing geometry and solver...");
            info_pb.set_message(format!(
                "Geometry: {} | Orientations: {}",
                self.settings.geom_name, n
            ));

            (status_pb, pb, info_pb)
        } else {
            // In quiet mode, create hidden progress bars
            (
                ProgressBar::hidden(),
                ProgressBar::hidden(),
                ProgressBar::hidden(),
            )
        };

        // init a set of base problems that can be reset
        let problems_base: Vec<Problem> = self
            .geoms
            .iter()
            .map(|geom| {
                Problem::new(Some(geom.clone()), Some(self.settings.clone()))
                    .expect("Failed to create Problem")
            })
            .collect();
        let num_problems = problems_base.iter().len();
        // let problem_base = Problem::new(Some(self.geoms.clone()), Some(self.settings.clone()));

        // Phase 2: Main computation
        status_pb.set_message("Running orientation averaging...");
        info_pb.set_message(format!("Processing {} orientations in parallel", n));

        // Solve for each orientation and reduce results on the fly
        self.result = self
            .orientations
            .eulers
            .par_iter()
            .map(|(a, b, g)| {
                let mut rng = if let Some(seed) = self.settings.seed {
                    rand::rngs::StdRng::seed_from_u64(seed)
                } else {
                    rand::rngs::StdRng::from_rng(&mut rand::rng())
                };
                let problem_idx = rng.random_range(0..num_problems);
                // choose a random problem from the base set
                let mut problem = problems_base[problem_idx].clone();
                let euler = Euler::new(*a, *b, *g);

                if let Err(err) = problem.run(Some(&euler)) {
                    log::error!("Error running problem (will skip this iteration): {}", err);
                }

                pb.inc(1);
                problem.result
            })
            .reduce(
                || self.result.zero_like(),
                |accum, item| self.reduce_results(accum, item),
            );

        // Phase 3: Post-processing
        pb.finish_with_message("Orientations complete");
        status_pb.set_message("Post-processing results...");

        // Normalize results by the number of orientations
        info_pb.set_message("Normalizing by orientation count...");
        self.normalize_results(self.orientations.num_orientations as f32);

        // Compute 1D integration
        info_pb.set_message("Computing 1D integrated Mueller matrices...");
        self.result.mueller_to_1d();

        // Compute derived parameters
        info_pb.set_message("Computing scattering parameters...");
        let _ = self.result.compute_params(self.settings.wavelength);

        // Phase 4: Complete
        status_pb.finish_with_message("✓ Computation complete");
        info_pb.finish_with_message(format!(
            "Power ratio: {:.3} | Results ready for output",
            self.result.powers.output / self.result.powers.input.max(1e-10)
        ));
    }

    /// Combines two Results objects by adding their fields
    fn reduce_results(&self, mut acc: Results, item: Results) -> Results {
        // Add powers
        acc.powers += item.powers;

        // Add zones (uses zone arithmetic)
        for (acc_zone, item_zone) in acc.zones.iter_mut().zip(item.zones.iter()) {
            for (a, i) in acc_zone.field_2d.iter_mut().zip(item_zone.field_2d.iter()) {
                // Amplitude matrices (complex)
                a.ampl_total += i.ampl_total;
                a.ampl_beam += i.ampl_beam;
                a.ampl_ext += i.ampl_ext;
                // Mueller matrices (real)
                a.mueller_total += i.mueller_total;
                a.mueller_beam += i.mueller_beam;
                a.mueller_ext += i.mueller_ext;
            }
        }

        acc
    }

    /// Normalizes the results by dividing by the number of orientations
    fn normalize_results(&mut self, num_orientations: f32) {
        // Powers
        self.result.powers /= num_orientations;

        // Normalize all zones
        let div_c = Complex::from(num_orientations);
        for zone in self.result.zones.iter_mut() {
            for field in zone.field_2d.iter_mut() {
                // Amplitude Matrices - divide by complex representation
                field.ampl_total /= div_c;
                field.ampl_beam /= div_c;
                field.ampl_ext /= div_c;

                // Mueller Matrices - divide by real value
                field.mueller_total /= num_orientations;
                field.mueller_beam /= num_orientations;
                field.mueller_ext /= num_orientations;
            }
        }
    }

    pub fn writeup(&self) {
        // Create progress display for writing phase
        let m = MultiProgress::new();

        let status_pb = m.add(ProgressBar::new_spinner());
        status_pb
            .set_style(ProgressStyle::with_template("{spinner:.cyan} Writing: {msg}").unwrap());
        status_pb.enable_steady_tick(Duration::from_millis(100));
        status_pb.set_message("Preparing output files...");

        // Use the new unified output system
        let output_manager = output::OutputManager::new(&self.settings, &self.result);
        let _ = output_manager.write_all();

        status_pb.finish_with_message(format!(
            "✓ Output written to {}",
            self.settings.directory.display()
        ));
    }
}

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl MultiProblem {
    #[new]
    #[pyo3(signature = (settings, geoms = None))]
    fn py_new(settings: Settings, geoms: Option<Vec<Geom>>) -> PyResult<Self> {
        // Load geometries from file if not provided
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
            problem::init_geom(&settings, geom);
        }
        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let result = init_result(&settings);

        Ok(Self {
            geoms,
            orientations,
            settings,
            result,
        })
    }

    /// Solve the multi-orientation scattering problem.
    ///
    /// Computes scattering properties averaged over all orientations using
    /// parallel processing. The Global Interpreter Lock (GIL) is released
    /// during computation to allow concurrent Python operations.
    ///
    /// # Returns
    /// PyResult<()> - Success or error if computation fails
    #[pyo3(name = "solve")]
    pub fn py_solve(&mut self, py: Python) -> PyResult<()> {
        py.detach(|| {
            self.solve();
        });
        Ok(())
    }

    /// Access the orientation-averaged simulation results.
    ///
    /// Returns the complete Results object containing Mueller matrices,
    /// amplitude matrices, power distributions, and derived parameters
    /// averaged over all orientations.
    ///
    /// # Returns
    /// Results - Complete scattering simulation results
    #[getter]
    pub fn get_results(&self) -> Results {
        self.result.clone()
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
        if let Some(dir) = directory {
            let mut settings = self.settings.clone();
            settings.directory = std::path::PathBuf::from(dir);
            let output_manager = crate::output::OutputManager::new(&settings, &self.result);
            output_manager.write_all().map_err(|e| {
                pyo3::exceptions::PyIOError::new_err(format!("Failed to save results: {}", e))
            })?;
        } else {
            self.writeup();
        }
        Ok(())
    }

    /// Reset the multiproblem to initial state
    pub fn py_reset(&mut self) -> PyResult<()> {
        self.reset();
        Ok(())
    }

    /// Regenerate orientations (useful for random schemes)
    pub fn py_regenerate_orientations(&mut self) -> PyResult<()> {
        self.regenerate_orientations();
        Ok(())
    }

    /// Get the number of orientations
    #[getter]
    pub fn get_num_orientations(&self) -> usize {
        self.orientations.num_orientations
    }
}
