use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

use crossbeam_deque::{Injector, Steal};

use crate::{
    geom::Geom,
    orientation::{Euler, Orientations},
    problem::{self, Problem},
    result::Results,
    settings::Settings,
};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use pyo3::prelude::*;
use rand::{Rng, SeedableRng};
use std::time::Duration;

/// Trait for types that can be tracked for convergence.
/// Provides operations needed for online mean/variance computation.
pub trait Convergeable: Clone + Sized {
    /// Create a zero/empty version with the same structure
    fn zero_like(&self) -> Self;

    /// Weighted addition: combines self (with weight w1) and other (with weight w2).
    /// For simple quantities: (self * w1 + other * w2) / (w1 + w2)
    /// For derived quantities like asymmetry: uses appropriate weighting (e.g., by ScatCross)
    fn weighted_add(&self, other: &Self, self_weight: f32, other_weight: f32) -> Self;

    /// Element-wise multiplication (for computing x²)
    fn mul_elem(&self, other: &Self) -> Self;

    /// Element-wise division
    fn div_elem(&self, other: &Self) -> Self;

    /// Element-wise addition
    fn add_elem(&self, other: &Self) -> Self;

    /// Element-wise subtraction
    fn sub_elem(&self, other: &Self) -> Self;

    /// Scale by a scalar
    fn scale(&self, scalar: f32) -> Self;

    /// Element-wise square root (for SEM computation)
    fn sqrt_elem(&self) -> Self;

    /// Returns a pre-weighted version for convergence tracking.
    /// e.g., asymmetry becomes asymmetry * scat_cross
    fn to_weighted(&self) -> Self;

    /// Returns the weights for each field.
    /// e.g., asymmetry weight is scat_cross, powers weight is 1.0
    fn weights(&self) -> Self;
}

/// Tracks running statistics for convergence using Welford's online algorithm.
/// Matches the Python implementation in goad/convergence/convergable.py exactly.
/// Computes mean and standard error of the mean (SEM) incrementally.
pub struct ConvergenceTracker<T: Convergeable> {
    i: usize, // iteration counter
    m: T,     // running weighted sum (stores value*weight accumulated via Welford)
    s: T,     // sum of squared deltas (for variance)
    w: T,     // running mean weight
}

impl<T: Convergeable> ConvergenceTracker<T> {
    /// Create a new tracker using a template for structure
    pub fn new(template: &T) -> Self {
        Self {
            i: 0,
            m: template.zero_like(),
            s: template.zero_like(),
            w: template.zero_like(),
        }
    }

    /// Update with a new result.
    /// Internally computes weighted value and weight via to_weighted()/weights(),
    /// then applies Welford's algorithm matching Python exactly.
    pub fn update(&mut self, result: &T) {
        self.i += 1;

        // Get pre-weighted value and weights from result
        let value = result.to_weighted();
        let weight = result.weights();

        if self.i == 1 {
            self.m = value;
            self.w = weight;
            // s stays zero
        } else {
            // delta = value - m_old
            let delta = value.sub_elem(&self.m);

            // m = m_old + delta / i
            self.m = self.m.add_elem(&delta.scale(1.0 / self.i as f32));

            // s = s + delta^2 * (i-1)/i
            let delta_sq = delta.mul_elem(&delta);
            let factor = (self.i - 1) as f32 / self.i as f32;
            self.s = self.s.add_elem(&delta_sq.scale(factor));

            // w = w_old + (weight - w_old) / i
            let dw = weight.sub_elem(&self.w);
            self.w = self.w.add_elem(&dw.scale(1.0 / self.i as f32));
        }
    }

    /// Get the current count
    pub fn count(&self) -> usize {
        self.i
    }

    /// Get the running mean: m / w
    pub fn mean(&self) -> T {
        if self.i == 0 {
            return self.m.zero_like();
        }
        self.m.div_elem(&self.w)
    }

    /// Get the standard error of the mean: sqrt(s / (i-1)^2 / w^2)
    pub fn sem(&self) -> T {
        if self.i < 2 {
            return self.m.zero_like();
        }
        // SEM = sqrt(s / (i-1)^2 / w^2)
        let n_minus_1 = (self.i - 1) as f32;
        let w_sq = self.w.mul_elem(&self.w);
        self.s
            .scale(1.0 / (n_minus_1 * n_minus_1))
            .div_elem(&w_sq)
            .sqrt_elem()
    }
}

/// A task representing a single orientation to be computed.
#[derive(Clone)]
struct OrientationTask {
    euler: Euler,
    problem_idx: usize,
}

/// Work-stealing based multi-orientation solver with convergence support.
///
/// Uses crossbeam-deque for work distribution:
/// - Master thread owns the Injector and performs live reduction
/// - Worker threads steal tasks and send results back via channel
/// - Convergence checked after each result (currently: 100 orientations)
#[pyclass]
#[derive(Debug)]
pub struct Convergence {
    pub geoms: Vec<Geom>,
    pub orientations: Orientations,
    pub settings: Settings,
    pub result: Results,
    pub error: Results, // tracks orientation averaging error
    pub convergence_target: usize,
}

impl Convergence {
    /// Creates a new Convergence solver from geometries and settings.
    pub fn new(geoms: Option<Vec<Geom>>, settings: Option<Settings>) -> anyhow::Result<Self> {
        let settings = settings
            .unwrap_or_else(|| crate::settings::load_config().expect("Failed to load config"));

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
            problem::init_geom(&settings, geom);
        }

        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let bins = &settings.binning.scheme.generate();
        let result = Results::new_empty(bins);
        let error = Results::new_empty(bins);

        Ok(Self {
            geoms,
            orientations,
            settings,
            result,
            error,
            convergence_target: 100, // dummy convergence: stop after 100 orientations
        })
    }

    /// Regenerates the orientations for the problem.
    pub fn regenerate_orientations(&mut self) {
        self.orientations =
            Orientations::generate(&self.settings.orientation.scheme, self.settings.seed);
    }

    /// Resets the solver to its initial state.
    pub fn reset(&mut self) {
        let bins: Vec<_> = self.result.field_2d.iter().map(|f| f.bin).collect();
        self.result = Results::new_empty(&bins);
        self.error = Results::new_empty(&bins);
        self.regenerate_orientations();
    }

    /// Solves using work-stealing parallelism.
    ///
    /// Architecture:
    /// - Master thread: owns Injector, receives results, performs reduction
    /// - Worker threads: steal from Injector, compute, send results via channel
    pub fn solve(&mut self) {
        let num_workers = std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(4)
            .saturating_sub(1) // reserve 1 for master
            .max(1);

        let n = self.orientations.num_orientations;
        let target = self.convergence_target.min(n);

        // Progress bars
        let (status_pb, pb, info_pb) = self.setup_progress_bars(target);
        status_pb.set_message("Initializing work-stealing solver...");

        // Prepare base problems (cloned per worker later)
        let problems_base: Vec<Problem> = self
            .geoms
            .iter()
            .map(|geom| Problem::new(Some(geom.clone()), Some(self.settings.clone())))
            .collect();
        let num_problems = problems_base.len();

        // Create the injector (global task queue)
        let injector: Injector<OrientationTask> = Injector::new();

        // Pre-populate tasks from the orientation pool
        let mut rng = if let Some(seed) = self.settings.seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };

        for (a, b, g) in self.orientations.eulers.iter() {
            let task = OrientationTask {
                euler: Euler::new(*a, *b, *g),
                problem_idx: rng.random_range(0..num_problems),
            };
            injector.push(task);
        }

        // Channel for results: workers send, master receives
        let (tx, rx): (Sender<Results>, Receiver<Results>) = mpsc::channel();

        // Spawn worker threads
        status_pb.set_message("Spawning worker threads...");
        info_pb.set_message(format!(
            "Workers: {} | Target: {} orientations",
            num_workers, target
        ));

        let injector_ref = &injector;
        let problems_ref = &problems_base;

        thread::scope(|s| {
            // Spawn workers
            for _ in 0..num_workers {
                let tx = tx.clone();
                s.spawn(move || {
                    Self::worker_loop(injector_ref, problems_ref, tx);
                });
            }

            // Drop the original sender so rx knows when all workers are done
            drop(tx);

            // Master reduction loop with convergence tracking
            status_pb.set_message("Running orientation averaging...");
            let mut tracker = ConvergenceTracker::new(&self.result);

            while tracker.count() < target {
                match rx.recv() {
                    Ok(result) => {
                        tracker.update(&result);
                        pb.inc(1);
                    }
                    Err(_) => {
                        // Channel closed, no more results coming
                        break;
                    }
                }
            }

            // Extract mean and SEM from tracker
            self.result = tracker.mean();
            self.error = tracker.sem();
        });

        // Post-processing
        pb.finish_with_message("Orientations complete");
        status_pb.set_message("Post-processing results...");

        info_pb.set_message("Computing 1D integrated Mueller matrices...");
        self.result.mueller_to_1d(&self.settings.binning.scheme);

        info_pb.set_message("Computing scattering parameters...");
        let _ = self.result.compute_params(self.settings.wavelength);

        status_pb.finish_with_message("✓ Computation complete");
        info_pb.finish_with_message(format!(
            "Power ratio: {:.3} | Results ready for output",
            self.result.powers.output / self.result.powers.input.max(1e-10)
        ));
    }

    /// Worker loop: steal tasks, compute, send results.
    fn worker_loop(
        injector: &Injector<OrientationTask>,
        problems_base: &[Problem],
        tx: Sender<Results>,
    ) {
        loop {
            // Try to steal a task
            match injector.steal() {
                Steal::Success(task) => {
                    // Clone the base problem for this task
                    let mut problem = problems_base[task.problem_idx].clone();

                    // Run the simulation
                    if let Err(err) = problem.run(Some(&task.euler)) {
                        eprintln!("Error running problem (will skip this iteration): {}", err);
                    }

                    // Send result back to master (ignore send errors on shutdown)
                    let _ = tx.send(problem.result);
                }
                Steal::Empty => {
                    // No more work, exit
                    break;
                }
                Steal::Retry => {
                    // Contention, try again
                    std::hint::spin_loop();
                }
            }
        }
    }

    /// Sets up progress bars for the solve.
    fn setup_progress_bars(&self, target: usize) -> (ProgressBar, ProgressBar, ProgressBar) {
        if !self.settings.quiet {
            let m = MultiProgress::new();

            let status_pb = m.add(ProgressBar::new_spinner());
            status_pb
                .set_style(ProgressStyle::with_template("{spinner:.cyan} Status: {msg}").unwrap());
            status_pb.enable_steady_tick(Duration::from_millis(100));

            let pb = m.add(ProgressBar::new(target as u64));
            pb.set_style(
                ProgressStyle::with_template(
                    "{spinner:.green} [{elapsed_precise}] [{bar:40.green/blue}] {pos:>5}/{len:5} {msg} | ETA: {eta_precise}",
                )
                .unwrap()
                .progress_chars("█▇▆▅▄▃▂▁"),
            );
            pb.set_message("Computing orientations");

            let info_pb = m.add(ProgressBar::new_spinner());
            info_pb.set_style(ProgressStyle::with_template("{msg}").unwrap());
            info_pb.enable_steady_tick(Duration::from_millis(500));
            info_pb.set_message(format!(
                "Geometry: {} | Target: {} orientations",
                self.settings.geom_name, target
            ));

            (status_pb, pb, info_pb)
        } else {
            (
                ProgressBar::hidden(),
                ProgressBar::hidden(),
                ProgressBar::hidden(),
            )
        }
    }
}

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
            problem::init_geom(&settings, geom);
        }

        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let bins = &settings.binning.scheme.generate();
        let result = Results::new_empty(bins);
        let error = Results::new_empty(bins);

        Ok(Self {
            geoms,
            orientations,
            settings,
            result,
            error,
            convergence_target: 100,
        })
    }

    /// Solve the multi-orientation scattering problem using work-stealing.
    #[pyo3(name = "solve")]
    pub fn py_solve(&mut self, py: Python) -> PyResult<()> {
        py.detach(|| {
            self.solve();
        });
        Ok(())
    }

    /// Access the orientation-averaged simulation results.
    #[getter]
    pub fn get_results(&self) -> Results {
        self.result.clone()
    }

    /// Access the orientation averaging error.
    #[getter]
    pub fn get_error(&self) -> Results {
        self.error.clone()
    }

    /// Get the number of orientations.
    #[getter]
    pub fn get_num_orientations(&self) -> usize {
        self.orientations.num_orientations
    }

    /// Get the convergence target.
    #[getter]
    pub fn get_convergence_target(&self) -> usize {
        self.convergence_target
    }

    /// Set the convergence target.
    #[setter]
    pub fn set_convergence_target(&mut self, target: usize) {
        self.convergence_target = target;
    }

    /// Reset the solver to initial state.
    pub fn py_reset(&mut self) -> PyResult<()> {
        self.reset();
        Ok(())
    }

    /// Regenerate orientations (useful for random schemes).
    pub fn py_regenerate_orientations(&mut self) -> PyResult<()> {
        self.regenerate_orientations();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_convergence_creation() {
        // Smoke test - full comparison in tests/convergence_tests.rs
    }
}
