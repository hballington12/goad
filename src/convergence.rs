mod convergeable;
mod progress;
mod python;

pub use convergeable::{Convergeable, ConvergenceTracker};
use log::{error, info, warn};

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
    output,
    params::Param,
    problem::Problem,
    result::{GOComponent, Results},
    settings::Settings,
};
use progress::ConvergenceProgress;
use pyo3::pyclass;
use pyo3_stub_gen::derive::*;
use rand::{Rng, SeedableRng};

const MAX_CONVERGENCE_ORIENTATIONS: usize = 100_000;

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

#[gen_stub_pyclass]
#[pyclass(module = "goad._goad")]
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
            let mean_results = self.tracker.mean();
            let sem_results = self.tracker.sem();
            for (i, target) in self.targets.iter().enumerate() {
                self.update_target(progress, i, target, &mean_results, &sem_results);
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
            info!("Converged after {} orientations", self.tracker.count());
        } else {
            info!(
                "Did not converge after {} orientations (max reached or interrupted)",
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
                        error!("Error running problem (will skip this iteration): {}", err);
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
