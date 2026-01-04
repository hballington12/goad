mod convergeable;
mod progress;
mod python;

pub use convergeable::{Convergeable, ConvergenceTracker};
use log::{error, info, warn};

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

use crossbeam_deque::{Injector, Steal};
use rand::rngs::StdRng;

use crate::{
    geom::Geom,
    multiproblem::{init_result, load_and_init_geoms, load_settings_or_default},
    orientation::{Euler, OrientationSampler, Scheme},
    output,
    params::Param,
    problem::Problem,
    result::{GOComponent, Results},
    settings::Settings,
};
use progress::ConvergenceProgress;
use pyo3::pyclass;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rand::{Rng, SeedableRng};

const MAX_CONVERGENCE_ORIENTATIONS: usize = 100_000;

/// Minimum batch size to avoid excessive overhead.
const MIN_BATCH_SIZE: usize = 2;

/// Maximum batch size to maintain load balancing and convergence responsiveness.
const MAX_BATCH_SIZE: usize = 1024;

/// Number of samples to run during prognosis for timing measurements.
const PROGNOSIS_SAMPLES: usize = 3;

/// Fallback batch size if prognosis fails or produces invalid results.
const FALLBACK_BATCH_SIZE: usize = 10;

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

#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
pub struct Convergence {
    pub geoms: Vec<Geom>,
    pub settings: Settings,
    pub max_orientations: usize,
    pub targets: Vec<ParamConvergenceTarget>,
    tracker: ConvergenceTracker<Results>,
    sampler: OrientationSampler,
    rng: StdRng,
    /// Optional path to log file for mean values during convergence
    pub log_file: Option<PathBuf>,
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

        // Create sampler based on scheme setting
        let sampler = match &settings.orientation.scheme {
            Scheme::Uniform { .. } => OrientationSampler::uniform(settings.seed),
            Scheme::Discrete { eulers } => OrientationSampler::discrete(eulers.clone()),
            Scheme::Sobol { .. } => OrientationSampler::sobol(settings.seed),
            Scheme::Halton { .. } => OrientationSampler::halton(),
        };

        Ok(Self {
            geoms,
            settings,
            max_orientations: MAX_CONVERGENCE_ORIENTATIONS, // safety cap
            targets: Vec::new(),
            tracker: ConvergenceTracker::new(&result),
            sampler,
            rng,
            log_file: None,
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

    /// Set the log file path for logging mean values during convergence.
    /// The file will be created/truncated when solve() is called.
    pub fn set_log_file(&mut self, path: impl Into<PathBuf>) {
        self.log_file = Some(path.into());
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

    /// Write results to output files.
    pub fn writeup(&self) {
        let mut result = self.mean();
        // Recompute params from averaged Mueller matrices
        result.mueller_to_1d();
        let _ = result.compute_params(self.settings.wavelength);
        let output_manager = output::OutputManager::new(&self.settings, &result);
        let _ = output_manager.write_all();
        info!("Output written to {}", self.settings.directory.display());
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
        let template = init_result(&self.settings);
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
                warn!(
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

    /// Runs a prognosis to determine optimal batch size based on actual timing.
    ///
    /// Measures:
    /// - compute_time: average time to run one orientation
    /// - merge_time: average time to merge one tracker into the master
    ///
    /// The optimal batch size ensures the master can keep up with all workers:
    ///   batch_size = ceil(num_workers * merge_time / compute_time)
    ///
    /// Returns the computed batch size, clamped to [MIN_BATCH_SIZE, MAX_BATCH_SIZE].
    fn run_prognosis(&mut self, problems_base: &[Problem], num_workers: usize) -> usize {
        let result_template = init_result(&self.settings);

        // Measure compute time: run PROGNOSIS_SAMPLES orientations
        let mut compute_times = Vec::with_capacity(PROGNOSIS_SAMPLES);

        for _ in 0..PROGNOSIS_SAMPLES {
            let Some(task) = self.sample_next_problem() else {
                break;
            };

            let mut problem = problems_base[task.problem_idx].clone();
            let start = Instant::now();
            if problem.run(Some(&task.euler)).is_ok() {
                compute_times.push(start.elapsed());
            }
        }

        if compute_times.is_empty() {
            warn!("Prognosis: no successful compute samples, using fallback batch size");
            return FALLBACK_BATCH_SIZE;
        }

        let avg_compute_time = compute_times.iter().sum::<Duration>() / compute_times.len() as u32;

        // Measure merge time: create sample trackers and merge them
        let mut merge_times = Vec::with_capacity(PROGNOSIS_SAMPLES);
        let mut dummy_tracker = ConvergenceTracker::new(&result_template);

        for _ in 0..PROGNOSIS_SAMPLES {
            // Create a tracker with one sample (simulates what workers send)
            let mut sample_tracker = ConvergenceTracker::new(&result_template);
            let Some(task) = self.sample_next_problem() else {
                break;
            };

            let mut problem = problems_base[task.problem_idx].clone();
            if problem.run(Some(&task.euler)).is_ok() {
                sample_tracker.update(&problem.result);

                let start = Instant::now();
                dummy_tracker.merge(&sample_tracker);
                merge_times.push(start.elapsed());
            }
        }

        if merge_times.is_empty() {
            warn!("Prognosis: no successful merge samples, using fallback batch size");
            return FALLBACK_BATCH_SIZE;
        }

        let avg_merge_time = merge_times.iter().sum::<Duration>() / merge_times.len() as u32;

        // Calculate optimal batch size:
        // Workers produce at rate: num_workers / compute_time
        // Master can handle at rate: 1 / merge_time
        // To balance: batch_size * (1 / merge_time) >= num_workers / compute_time
        // => batch_size >= num_workers * merge_time / compute_time
        let batch_size = if avg_compute_time.as_nanos() > 0 {
            let ratio = (num_workers as f64 * avg_merge_time.as_nanos() as f64)
                / avg_compute_time.as_nanos() as f64;
            // Add 20% headroom to keep master comfortably ahead
            (ratio * 1.2).ceil() as usize
        } else {
            warn!(
                "Prognosis: compute time was zero (too fast to measure), using fallback batch size"
            );
            return FALLBACK_BATCH_SIZE;
        };

        let clamped = batch_size.clamp(MIN_BATCH_SIZE, MAX_BATCH_SIZE);

        info!(
            "Prognosis: compute={:?}, merge={:?}, workers={}, optimal_batch={}",
            avg_compute_time, avg_merge_time, num_workers, clamped
        );

        clamped
    }

    /// Runs the work-stealing solver with worker threads.
    ///
    /// Workers accumulate results locally into a ConvergenceTracker, then send
    /// the tracker to the master for merging. This reduces master thread overhead
    /// when many workers are active.
    fn run<F>(&mut self, mut check_interrupt: F)
    where
        F: FnMut() -> bool,
    {
        let num_workers = Self::num_workers();
        let progress = ConvergenceProgress::new(self.targets.len(), self.max_orientations);
        let problems_base = self.create_base_problems();

        // Initialize log file if configured
        let mut log_file: Option<File> = self.log_file.as_ref().and_then(|path| {
            match File::create(path) {
                Ok(mut f) => {
                    // Write header
                    let _ = writeln!(f, "count,lidar_ratio,lidar_ratio_sem,relative_sem_pct");
                    Some(f)
                }
                Err(e) => {
                    warn!("Failed to create log file {:?}: {}", path, e);
                    None
                }
            }
        });

        // Run prognosis to determine optimal batch size based on actual timing
        let batch_size = self
            .run_prognosis(&problems_base, num_workers)
            .min(self.max_orientations);
        let injector: Injector<OrientationTask> = Injector::new();

        // Template for workers to create their local trackers
        let result_template = init_result(&self.settings);

        // Initial task queue fill - enough for all workers to have a full batch
        let buffer_size = (num_workers * batch_size * 2).min(self.max_orientations);
        self.fill_initial_tasks(&injector, buffer_size);

        // Channel for batched results: workers send trackers, master merges
        let (tx, rx): (
            Sender<ConvergenceTracker<Results>>,
            Receiver<ConvergenceTracker<Results>>,
        ) = mpsc::channel();

        // Shutdown flag for workers
        let done = Arc::new(AtomicBool::new(false));

        let injector_ref = &injector;
        let problems_ref = &problems_base;
        let template_ref = &result_template;

        thread::scope(|s| {
            // Spawn workers
            for _ in 0..num_workers {
                let tx = tx.clone();
                let done = Arc::clone(&done);
                s.spawn(move || {
                    Self::worker_loop_batched(
                        injector_ref,
                        problems_ref,
                        template_ref,
                        batch_size,
                        tx,
                        &done,
                    );
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
                    Ok(batch_tracker) => self.update_from_batch(
                        &progress,
                        &injector,
                        &mut converged,
                        batch_tracker,
                        &mut log_file,
                    ),
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

    /// Update from a batch of results (merged tracker from a worker).
    fn update_from_batch(
        &mut self,
        progress: &ConvergenceProgress,
        injector: &Injector<OrientationTask>,
        converged: &mut bool,
        batch_tracker: ConvergenceTracker<Results>,
        log_file: &mut Option<File>,
    ) {
        let batch_count = batch_tracker.count();
        self.tracker.merge(&batch_tracker);
        let count = self.tracker.count();
        progress.update_info(count);

        // Update per-target progress bars and log mean values
        if count >= MIN_ORIENTATIONS {
            let mean_results = self.tracker.mean();
            let sem_results = self.tracker.sem();
            for (i, target) in self.targets.iter().enumerate() {
                self.update_target(progress, i, target, &mean_results, &sem_results);
            }

            // Log mean values to file if configured
            if let Some(ref mut file) = log_file {
                let lidar = mean_results
                    .params
                    .lidar_ratio(&GOComponent::Total)
                    .unwrap_or(0.0);
                let lidar_sem = sem_results
                    .params
                    .lidar_ratio(&GOComponent::Total)
                    .unwrap_or(0.0);
                let relative_sem = if lidar.abs() > 1e-10 {
                    (lidar_sem / lidar.abs()) * 100.0
                } else {
                    0.0
                };
                let _ = writeln!(
                    file,
                    "{},{:.6},{:.6},{:.4}",
                    count, lidar, lidar_sem, relative_sem
                );
            }
        }

        // Check convergence periodically (every batch after minimum)
        if self.tracker.count() >= MIN_ORIENTATIONS {
            *converged = self.is_converged();
        }

        // Replenish task queue with as many tasks as we just processed
        if !*converged && self.tracker.count() < self.max_orientations {
            for _ in 0..batch_count {
                self.push_task(injector);
            }
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
        // Initialize file-based logging (avoids conflicts with indicatif progress bar)
        if let Err(e) = crate::filelog::init(&self.settings.directory) {
            warn!("Could not initialize file logging: {}", e);
        }

        // Validation
        if self.targets.is_empty() {
            anyhow::bail!("No convergence targets set. Use add_target() before solving.");
        }

        // Run the worker pool
        self.run(check_interrupt);

        // Print final status
        if self.is_converged() {
            info!("Converged after {} orientations", self.tracker.count());
        } else {
            info!(
                "Did not converge after {} orientations (max reached or interrupted)",
                self.tracker.count()
            );
        }

        Ok(())
    }

    /// Worker loop with batching: steal batch_size tasks, process them, send tracker to master.
    ///
    /// This reduces master thread overhead by having workers perform local
    /// reduction using Welford's algorithm, then sending the accumulated
    /// statistics for merging.
    fn worker_loop_batched(
        injector: &Injector<OrientationTask>,
        problems_base: &[Problem],
        result_template: &Results,
        batch_size: usize,
        tx: Sender<ConvergenceTracker<Results>>,
        done: &AtomicBool,
    ) {
        loop {
            // Try to steal a batch of tasks
            let tasks = Self::steal_batch(injector, batch_size, done);

            if tasks.is_empty() {
                // No tasks and done signal received
                if done.load(Ordering::Relaxed) {
                    break;
                }
                // Otherwise keep waiting
                std::hint::spin_loop();
                continue;
            }

            // Process all stolen tasks, accumulating into local tracker
            let mut local_tracker = ConvergenceTracker::new(result_template);

            for task in tasks {
                let mut problem = problems_base[task.problem_idx].clone();

                if let Err(err) = problem.run(Some(&task.euler)) {
                    error!("Error running problem (will skip this iteration): {}", err);
                    continue;
                }

                local_tracker.update(&problem.result);
            }

            // Send batch to master (if we accumulated any results)
            if local_tracker.count() > 0 {
                if tx.send(local_tracker).is_err() {
                    // Channel closed, master is done
                    break;
                }
            }
        }
    }

    /// Steal up to `batch_size` tasks from the injector.
    /// Returns early if done signal is set and no tasks are available.
    fn steal_batch(
        injector: &Injector<OrientationTask>,
        batch_size: usize,
        done: &AtomicBool,
    ) -> Vec<OrientationTask> {
        let mut tasks = Vec::with_capacity(batch_size);

        while tasks.len() < batch_size {
            match injector.steal() {
                Steal::Success(task) => {
                    tasks.push(task);
                }
                Steal::Empty => {
                    // No more tasks available right now
                    if !tasks.is_empty() {
                        // Return what we have (partial batch)
                        break;
                    }
                    // Nothing stolen yet - check if we should exit
                    if done.load(Ordering::Relaxed) {
                        break;
                    }
                    // Wait a bit for more tasks
                    std::hint::spin_loop();
                }
                Steal::Retry => {
                    // Contention, try again immediately
                    std::hint::spin_loop();
                }
            }
        }

        tasks
    }

    fn update_target(
        &self,
        progress: &ConvergenceProgress,
        i: usize,
        target: &ParamConvergenceTarget,
        mean: &Results,
        sem: &Results,
    ) {
        let Some(mean_val) = mean.params.get(&target.param, &GOComponent::Total) else {
            return;
        };
        let Some(sem_val) = sem.params.get(&target.param, &GOComponent::Total) else {
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
