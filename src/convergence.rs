use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::Arc;
use std::thread;

use chrono::Local;
use crossbeam_deque::{Injector, Steal};

use crate::{
    geom::Geom,
    orientation::{Euler, OrientationSampler, UniformSampler},
    params::Param,
    problem::{self, Problem},
    result::{GOComponent, Results},
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

/// Work-stealing based multi-orientation solver with convergence support.
///
/// Uses crossbeam-deque for work distribution:
/// - Master thread owns the Injector and performs live reduction
/// - Worker threads steal tasks and send results back via channel
/// - Convergence checked after each result (currently: 100 orientations)
use crate::settings::constants::MIN_ORIENTATIONS;

/// Check if all convergence targets are satisfied (standalone version for disjoint borrows).
fn is_converged_check(
    tracker: &ConvergenceTracker<Results>,
    targets: &[ParamConvergenceTarget],
) -> bool {
    if tracker.count() < MIN_ORIENTATIONS {
        return false;
    }
    if targets.is_empty() {
        return false;
    }

    let mean = tracker.mean();
    let sem = tracker.sem();

    targets.iter().all(|t| {
        let mean_val = match t.param {
            Param::Asymmetry => mean.params.asymmetry(&GOComponent::Total),
            Param::Albedo => mean.params.albedo(&GOComponent::Total),
            Param::ScatCross => mean.params.scatt_cross(&GOComponent::Total),
            Param::ExtCross => mean.params.ext_cross(&GOComponent::Total),
        };
        let sem_val = match t.param {
            Param::Asymmetry => sem.params.asymmetry(&GOComponent::Total),
            Param::Albedo => sem.params.albedo(&GOComponent::Total),
            Param::ScatCross => sem.params.scatt_cross(&GOComponent::Total),
            Param::ExtCross => sem.params.ext_cross(&GOComponent::Total),
        };

        match (mean_val, sem_val) {
            (Some(m), Some(s)) if m.abs() > 1e-10 => (s / m.abs()) < t.relative_error,
            _ => false,
        }
    })
}

#[pyclass]
pub struct Convergence {
    pub geoms: Vec<Geom>,
    pub settings: Settings,
    pub max_orientations: usize,
    pub targets: Vec<ParamConvergenceTarget>,
    tracker: ConvergenceTracker<Results>,
    sampler: UniformSampler,
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

        let bins = &settings.binning.scheme.generate();
        let template = Results::new_empty(bins);
        let sampler = UniformSampler::new(settings.seed);

        Ok(Self {
            geoms,
            settings,
            max_orientations: 100_000, // safety cap
            targets: Vec::new(),
            tracker: ConvergenceTracker::new(&template),
            sampler,
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
            let mean_val = match t.param {
                Param::Asymmetry => mean.params.asymmetry(&GOComponent::Total),
                Param::Albedo => mean.params.albedo(&GOComponent::Total),
                Param::ScatCross => mean.params.scatt_cross(&GOComponent::Total),
                Param::ExtCross => mean.params.ext_cross(&GOComponent::Total),
            };
            let sem_val = match t.param {
                Param::Asymmetry => sem.params.asymmetry(&GOComponent::Total),
                Param::Albedo => sem.params.albedo(&GOComponent::Total),
                Param::ScatCross => sem.params.scatt_cross(&GOComponent::Total),
                Param::ExtCross => sem.params.ext_cross(&GOComponent::Total),
            };

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
    pub fn solve_with_interrupt<F>(&mut self, mut check_interrupt: F) -> anyhow::Result<()>
    where
        F: FnMut() -> bool,
    {
        if self.targets.is_empty() {
            anyhow::bail!("No convergence targets set. Use add_target() before solving.");
        }

        let num_workers = std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or_else(|e| {
                eprintln!(
                    "Warning: Could not determine available parallelism ({}), defaulting to 4",
                    e
                );
                4
            })
            .saturating_sub(1) // reserve 1 for master thread doing reduction
            .max(1);

        let max_target = self.max_orientations;

        // Progress display
        let m = MultiProgress::new();
        let title_pb = m.add(ProgressBar::new_spinner());
        title_pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} GOAD: [Convergence]  [Elapsed: {elapsed}]  {msg}  [{prefix}]",
            )
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
        );
        title_pb.set_message("[Status: \x1b[33mINITIALISING\x1b[0m]");
        title_pb.set_prefix(Local::now().format("%Y-%m-%d %H:%M:%S").to_string());
        title_pb.enable_steady_tick(Duration::from_millis(100));

        let info_pb = m.add(ProgressBar::new_spinner());
        info_pb.set_style(ProgressStyle::with_template("  {msg}").unwrap());

        // Create progress bars for each target
        let target_pbs: Vec<ProgressBar> = self
            .targets
            .iter()
            .map(|_| {
                let pb = m.add(ProgressBar::new(100));
                pb.set_style(
                    ProgressStyle::with_template("  {msg} [{bar:20.green/dim}] {pos:>3}%")
                        .unwrap()
                        .progress_chars("█▓░"),
                );
                pb
            })
            .collect();

        let start_time = std::time::Instant::now();

        // Prepare base problems (cloned per worker later)
        let problems_base: Vec<Problem> = self
            .geoms
            .iter()
            .map(|geom| Problem::new(Some(geom.clone()), Some(self.settings.clone())))
            .collect();
        let num_problems = problems_base.len();

        // Create the injector (global task queue)
        let injector: Injector<OrientationTask> = Injector::new();

        // RNG for selecting problem index
        let mut rng = if let Some(seed) = self.settings.seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };

        // Disjoint borrows of self fields to avoid borrow conflicts in thread::scope
        let sampler = &mut self.sampler;
        let tracker = &mut self.tracker;
        let targets = &self.targets;

        let mut tasks_pushed = 0;

        // Initial fill: 2 tasks per worker
        let buffer_size = (num_workers * 2).min(max_target);
        for _ in 0..buffer_size {
            if let Some(euler) = sampler.next() {
                let task = OrientationTask {
                    euler,
                    problem_idx: rng.random_range(0..num_problems),
                };
                injector.push(task);
                tasks_pushed += 1;
            }
        }

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

            // Set status to RUNNING
            title_pb.set_message("[Status: \x1b[32mRUNNING\x1b[0m]");

            while tracker.count() < max_target && !converged && !interrupted {
                // Use timeout so we can periodically check for interrupts
                match rx.recv_timeout(Duration::from_millis(100)) {
                    Ok(result) => {
                        tracker.update(&result);

                        // Update info bar
                        let count = tracker.count();
                        let elapsed = start_time.elapsed().as_secs_f64();
                        let sec_per_orient = if count > 0 {
                            elapsed / count as f64
                        } else {
                            0.0
                        };
                        let min_color = if count >= MIN_ORIENTATIONS {
                            "\x1b[32m" // green
                        } else {
                            "\x1b[31m" // red
                        };
                        info_pb.set_message(format!(
                            "[Orientations: {} ({}{}{}|{})] [{:.3} sec/orientation]",
                            count,
                            min_color,
                            MIN_ORIENTATIONS,
                            "\x1b[0m",
                            max_target,
                            sec_per_orient
                        ));
                        title_pb.set_prefix(Local::now().format("%Y-%m-%d %H:%M:%S").to_string());

                        // Update per-target progress bars
                        if count >= MIN_ORIENTATIONS {
                            let mean = tracker.mean();
                            let sem = tracker.sem();

                            for (i, target) in targets.iter().enumerate() {
                                let mean_val = match target.param {
                                    Param::Asymmetry => mean.params.asymmetry(&GOComponent::Total),
                                    Param::Albedo => mean.params.albedo(&GOComponent::Total),
                                    Param::ScatCross => {
                                        mean.params.scatt_cross(&GOComponent::Total)
                                    }
                                    Param::ExtCross => mean.params.ext_cross(&GOComponent::Total),
                                };
                                let sem_val = match target.param {
                                    Param::Asymmetry => sem.params.asymmetry(&GOComponent::Total),
                                    Param::Albedo => sem.params.albedo(&GOComponent::Total),
                                    Param::ScatCross => sem.params.scatt_cross(&GOComponent::Total),
                                    Param::ExtCross => sem.params.ext_cross(&GOComponent::Total),
                                };

                                if let (Some(m), Some(s)) = (mean_val, sem_val) {
                                    let current_rel_sem =
                                        if m.abs() > 1e-10 { s / m.abs() } else { 0.0 };
                                    let target_rel_sem = target.relative_error;

                                    // Progress with sqrt scaling, capped at 100%
                                    let progress = if current_rel_sem > 1e-10 {
                                        ((target_rel_sem / current_rel_sem).sqrt()).min(1.0)
                                    } else {
                                        1.0
                                    };

                                    let param_name = format!("{:?}", target.param);
                                    target_pbs[i].set_message(format!(
                                        "{:<9} {:>10.4e} ± {:<10.4e} [{:>5.2}% / {:>5.2}%]",
                                        param_name,
                                        m,
                                        s,
                                        current_rel_sem * 100.0,
                                        target_rel_sem * 100.0
                                    ));
                                    target_pbs[i].set_position((progress * 100.0) as u64);
                                }
                            }
                        }

                        // Check convergence periodically (every orientation after minimum)
                        if tracker.count() >= MIN_ORIENTATIONS {
                            converged = is_converged_check(tracker, targets);
                        }

                        // Replenish task queue if not converged and more orientations available
                        if !converged && tasks_pushed < max_target {
                            if let Some(euler) = sampler.next() {
                                let task = OrientationTask {
                                    euler,
                                    problem_idx: rng.random_range(0..num_problems),
                                };
                                injector.push(task);
                                tasks_pushed += 1;
                            }
                        }
                    }
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
            title_pb.set_message("[Status: \x1b[33mFINALISING\x1b[0m]");
        });

        // Finish progress bars
        title_pb.finish_and_clear();
        for pb in &target_pbs {
            pb.finish_and_clear();
        }
        info_pb.finish_and_clear();

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

        let bins = &settings.binning.scheme.generate();
        let template = Results::new_empty(bins);
        let sampler = UniformSampler::new(settings.seed);

        Ok(Self {
            geoms,
            settings,
            max_orientations: 100_000,
            targets: Vec::new(),
            tracker: ConvergenceTracker::new(&template),
            sampler,
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
