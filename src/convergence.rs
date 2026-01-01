mod progress;

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::Arc;
use std::thread;
use std::time::Duration;

use crossbeam_deque::{Injector, Steal};
use rand::rngs::StdRng;

use crate::{
    geom::Geom,
    multiproblem::{init_result, load_and_init_geoms, load_settings_or_default},
    orientation::{Euler, OrientationSampler},
    params::Param,
    problem::{init_geom, Problem},
    result::{GOComponent, Results},
    settings::Settings,
};
use progress::ConvergenceProgress;
use pyo3::prelude::*;
use rand::{Rng, SeedableRng};

const MAX_CONVERGENCE_ORIENTATIONS: usize = 100_000;

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
#[derive(Debug)]
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

/// A convergence target for a specific parameter.
#[derive(Clone, Debug)]
pub struct ParamConvergenceTarget {
    pub param: Param,
    pub relative_error: f32,
}

impl ParamConvergenceTarget {
    pub fn new(param: Param, relative_error: f32) -> Self {
        Self {
            param,
            relative_error,
        }
    }
}

/// A task representing a single orientation to be computed.
#[derive(Clone)]
struct OrientationTask {
    euler: Euler,
    problem_idx: usize,
}

use crate::settings::constants::MIN_ORIENTATIONS;

#[pyclass]
pub struct Convergence {
    pub geoms: Vec<Geom>,
    pub settings: Settings,
    pub max_orientations: usize,
    pub targets: Vec<ParamConvergenceTarget>,
    tracker: ConvergenceTracker<Results>,
    sampler: OrientationSampler,
    rng: StdRng,
}

impl Convergence {
    /// Creates a new Convergence solver from geometries and settings.
    pub fn new(geoms: Option<Vec<Geom>>, settings: Option<Settings>) -> anyhow::Result<Self> {
        let settings = load_settings_or_default(settings);
        let geoms = load_and_init_geoms(geoms, &settings)?;
        let result = init_result(&settings);
        let rng = if let Some(seed) = settings.seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };

        // Convergence always uses uniform random sampling (infinite supply)
        let sampler = OrientationSampler::uniform(settings.seed);

        Ok(Self {
            geoms,
            settings,
            max_orientations: MAX_CONVERGENCE_ORIENTATIONS, // safety cap
            targets: Vec::new(),
            tracker: ConvergenceTracker::new(&result),
            sampler,
            rng,
        })
    }

    /// Add a convergence target for a parameter.
    /// Solver will terminate when ALL targets are satisfied.
    pub fn add_target(&mut self, param: Param, relative_error: f32) {
        self.targets
            .push(ParamConvergenceTarget::new(param, relative_error));
    }

    /// Clear all convergence targets.
    pub fn clear_targets(&mut self) {
        self.targets.clear();
    }

    /// Get the number of orientations computed so far.
    pub fn count(&self) -> usize {
        self.tracker.count()
    }

    /// Get the current mean results (live during solve).
    pub fn mean(&self) -> Results {
        self.tracker.mean()
    }

    /// Get the current standard error of the mean (live during solve).
    pub fn sem(&self) -> Results {
        self.tracker.sem()
    }

    /// Check if all convergence targets are satisfied.
    fn is_converged(&self) -> bool {
        // Need minimum orientations for stable SEM
        if self.tracker.count() < MIN_ORIENTATIONS {
            return false;
        }

        // No targets means use max_orientations only
        if self.targets.is_empty() {
            return false;
        }

        let mean = self.tracker.mean();
        let sem = self.tracker.sem();

        self.targets.iter().all(|t| {
            let mean_val = mean.params.get(&t.param, &GOComponent::Total);
            let sem_val = sem.params.get(&t.param, &GOComponent::Total);

            match (mean_val, sem_val) {
                (Some(m), Some(s)) if m.abs() > 1e-10 => (s / m.abs()) < t.relative_error,
                _ => false, // can't check, not converged
            }
        })
    }

    /// Solves using work-stealing parallelism (non-interruptible version).
    pub fn solve(&mut self) -> anyhow::Result<()> {
        self.solve_with_interrupt(|| false)
    }

    /// Resets the sampler to its initial state.
    pub fn reset_sampler(&mut self) {
        self.sampler.reset();
    }

    /// Resets the solver to its initial state.
    pub fn reset(&mut self) {
        let bins = self.settings.binning.scheme.generate();
        let template = Results::new_empty(&bins);
        self.tracker = ConvergenceTracker::new(&template);
        self.reset_sampler();
    }

    // ========================================================================
    // Helper methods for solve_with_interrupt
    // ========================================================================

    /// Determines the number of worker threads to use.
    /// Reserves 1 thread for the master (reduction), minimum 1 worker.
    fn num_workers() -> usize {
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or_else(|e| {
                eprintln!(
                    "Warning: Could not determine available parallelism ({}), defaulting to 4",
                    e
                );
                4
            })
            .saturating_sub(1)
            .max(1)
    }

    /// Creates base problems from geometries.
    fn create_base_problems(&self) -> Vec<Problem> {
        self.geoms
            .iter()
            .map(|geom| {
                Problem::new(Some(geom.clone()), Some(self.settings.clone()))
                    .expect("Failed to create Problem")
            })
            .collect()
    }

    /// Runs the work-stealing solver with worker threads.
    fn run<F>(&mut self, mut check_interrupt: F)
    where
        F: FnMut() -> bool,
    {
        let num_workers = Self::num_workers();
        let progress = ConvergenceProgress::new(self.targets.len(), self.max_orientations);
        let problems_base = self.create_base_problems();
        let injector: Injector<OrientationTask> = Injector::new();

        // Initial task queue fill
        let buffer_size = (num_workers * 2).min(self.max_orientations);
        self.fill_initial_tasks(&injector, buffer_size);

        // Channel for results: workers send, master receives
        let (tx, rx): (Sender<Results>, Receiver<Results>) = mpsc::channel();

        // Shutdown flag for workers
        let done = Arc::new(AtomicBool::new(false));

        let injector_ref = &injector;
        let problems_ref = &problems_base;

        thread::scope(|s| {
            // Spawn workers
            for _ in 0..num_workers {
                let tx = tx.clone();
                let done = Arc::clone(&done);
                s.spawn(move || {
                    Self::worker_loop(injector_ref, problems_ref, tx, &done);
                });
            }

            // Drop the original sender so rx knows when all workers are done
            drop(tx);

            // Master reduction loop with convergence tracking
            let mut converged = false;
            let mut interrupted = false;

            progress.set_running();

            while self.tracker.count() < self.max_orientations && !converged && !interrupted {
                // Use timeout so we can periodically check for interrupts
                match rx.recv_timeout(Duration::from_millis(100)) {
                    Ok(result) => self.update(&progress, &injector, &mut converged, result),
                    Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                        // Check for interrupt (e.g., Ctrl-C from Python)
                        if check_interrupt() {
                            interrupted = true;
                        }
                    }
                    Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                        // Channel closed, no more results coming
                        break;
                    }
                }
            }

            // Signal workers to exit and set status to FINALISING
            done.store(true, Ordering::Relaxed);
            progress.set_finalising();
        });

        progress.finish();
    }

    fn fill_initial_tasks(&mut self, injector: &Injector<OrientationTask>, buffer_size: usize) {
        for _ in 0..buffer_size {
            self.push_task(injector);
        }
    }

    fn push_task(&mut self, injector: &Injector<OrientationTask>) {
        if let Some(task) = self.sample_next_problem() {
            injector.push(task);
        }
    }

    fn update(
        &mut self,
        progress: &ConvergenceProgress,
        injector: &Injector<OrientationTask>,
        converged: &mut bool,
        result: Results,
    ) {
        self.tracker.update(&result);
        let count = self.tracker.count();
        progress.update_info(count);
        // Update per-target progress bars
        if count >= MIN_ORIENTATIONS {
            for (i, target) in self.targets.iter().enumerate() {
                self.update_target(progress, i, target);
            }
        }
        // Check convergence periodically (every orientation after minimum)
        if self.tracker.count() >= MIN_ORIENTATIONS {
            *converged = self.is_converged();
        }
        // Replenish task queue if not converged and more orientations available
        if !*converged && self.tracker.count() < self.max_orientations {
            self.push_task(injector);
        }
    }

    fn sample_next_problem(&mut self) -> Option<OrientationTask> {
        self.sampler.next().map(|euler| OrientationTask {
            euler,
            problem_idx: self.rng.random_range(0..self.geoms.len()),
        })
    }

    // ========================================================================

    /// Solves using work-stealing parallelism.
    ///
    /// Architecture:
    /// - Master thread: owns Injector, receives results, performs reduction
    /// - Worker threads: steal from Injector, compute, send results via channel
    ///
    /// Termination: stops when all convergence targets are satisfied,
    /// or when max_orientations is reached (whichever comes first).
    ///
    /// The optional `check_interrupt` closure is called periodically to allow
    /// signal handling (e.g., Ctrl-C from Python). Return `true` to interrupt.
    pub fn solve_with_interrupt<F>(&mut self, check_interrupt: F) -> anyhow::Result<()>
    where
        F: FnMut() -> bool,
    {
        // Validation
        if self.targets.is_empty() {
            anyhow::bail!("No convergence targets set. Use add_target() before solving.");
        }

        // Run the worker pool
        self.run(check_interrupt);

        // Print final status
        if self.is_converged() {
            println!("Converged after {} orientations", self.tracker.count());
        } else {
            println!(
                "Completed {} orientations (max reached or interrupted)",
                self.tracker.count()
            );
        }

        Ok(())
    }

    /// Worker loop: steal tasks, compute, send results.
    fn worker_loop(
        injector: &Injector<OrientationTask>,
        problems_base: &[Problem],
        tx: Sender<Results>,
        done: &AtomicBool,
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
                    // Check if we should exit
                    if done.load(Ordering::Relaxed) {
                        break;
                    }
                    // Otherwise wait for more work
                    std::hint::spin_loop();
                }
                Steal::Retry => {
                    // Contention, try again
                    std::hint::spin_loop();
                }
            }
        }
    }

    fn update_target(
        &self,
        progress: &ConvergenceProgress,
        i: usize,
        target: &ParamConvergenceTarget,
    ) {
        let Some(mean_val) = self
            .tracker
            .mean()
            .params
            .get(&target.param, &GOComponent::Total)
        else {
            return;
        };
        let Some(sem_val) = self
            .tracker
            .sem()
            .params
            .get(&target.param, &GOComponent::Total)
        else {
            return;
        };

        progress.update_target(
            i,
            target.param.clone(),
            mean_val,
            sem_val,
            target.relative_error,
        );
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
            init_geom(&settings, geom);
        }

        let bins = &settings.binning.scheme.generate();
        let template = Results::new_empty(bins);

        // Convergence always uses uniform random sampling (infinite supply)
        let sampler = OrientationSampler::uniform(settings.seed);

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
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_convergence_creation() {
        // Smoke test - full comparison in tests/convergence_tests.rs
    }
}
