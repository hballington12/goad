//! Wall-time instrumentation for GOAD simulations.
//!
//! Three layers:
//! - [`StageTimings`]: per-orientation stage durations recorded by a `Problem`.
//! - [`TimingAccumulator`]: lock-free aggregation of stage timings across
//!   worker threads (relaxed atomic adds, negligible overhead per orientation).
//! - [`Timings`]: the final report, serialized into `results.json` and
//!   printable as a tree-style breakdown on the console.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

use serde::Serialize;

/// Stage durations for a single orientation solve, recorded by a `Problem`.
#[derive(Debug, Clone, Default)]
pub struct StageTimings {
    pub near_field: Duration,
    pub far_field: Duration,
    pub mueller_1d: Duration,
    pub params: Duration,
}

impl StageTimings {
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}

/// Aggregates per-orientation stage timings across worker threads.
#[derive(Debug, Default)]
pub struct TimingAccumulator {
    near_field_ns: AtomicU64,
    far_field_ns: AtomicU64,
    mueller_1d_ns: AtomicU64,
    params_ns: AtomicU64,
    total_ns: AtomicU64,
    orientations: AtomicU64,
}

impl TimingAccumulator {
    /// Adds one orientation's stage timings plus its total run wall time.
    pub fn add(&self, stages: &StageTimings, run_total: Duration) {
        let add = |counter: &AtomicU64, d: Duration| {
            counter.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
        };
        add(&self.near_field_ns, stages.near_field);
        add(&self.far_field_ns, stages.far_field);
        add(&self.mueller_1d_ns, stages.mueller_1d);
        add(&self.params_ns, stages.params);
        add(&self.total_ns, run_total);
        self.orientations.fetch_add(1, Ordering::Relaxed);
    }

    pub fn snapshot(&self) -> WorkerTimings {
        let secs = |counter: &AtomicU64| counter.load(Ordering::Relaxed) as f64 * 1e-9;
        WorkerTimings {
            near_field: secs(&self.near_field_ns),
            far_field: secs(&self.far_field_ns),
            mueller_1d: secs(&self.mueller_1d_ns),
            params: secs(&self.params_ns),
            total: secs(&self.total_ns),
            orientations: self.orientations.load(Ordering::Relaxed),
        }
    }
}

/// Stage times summed over all orientations across all worker threads,
/// in seconds. These are accumulated thread times, not wall times.
#[derive(Debug, Clone, Default, Serialize)]
pub struct WorkerTimings {
    pub near_field: f64,
    pub far_field: f64,
    pub mueller_1d: f64,
    pub params: f64,
    /// Full per-orientation run time (includes init/orientation/illumination
    /// overhead not covered by the named stages).
    pub total: f64,
    pub orientations: u64,
}

impl WorkerTimings {
    fn other(&self) -> f64 {
        (self.total - self.near_field - self.far_field - self.mueller_1d - self.params).max(0.0)
    }
}

/// Wall-clock times for the outer pipeline stages, in seconds.
#[derive(Debug, Clone, Default, Serialize)]
pub struct WallTimings {
    /// Geometry and solver setup.
    pub init: f64,
    /// Convergence-mode batch-size prognosis (absent for fixed orientations).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub prognosis: Option<f64>,
    /// The parallel orientation loop.
    pub solve: f64,
    /// Normalization, 1D integration, and parameter computation.
    pub post_process: f64,
    /// Output file writing (excluding results.json itself).
    pub file_io: f64,
}

/// Complete timing report for a simulation run.
#[derive(Debug, Clone, Default, Serialize)]
pub struct Timings {
    /// Worker threads used for the orientation loop.
    pub threads: usize,
    pub wall: WallTimings,
    pub worker: WorkerTimings,
}

impl Timings {
    /// Total wall time across all measured pipeline stages, in seconds.
    pub fn total_wall(&self) -> f64 {
        self.wall.init
            + self.wall.prognosis.unwrap_or(0.0)
            + self.wall.solve
            + self.wall.post_process
            + self.wall.file_io
    }

    /// Fraction of available thread time spent doing orientation work
    /// during the solve loop. 1.0 means perfect scaling. Can slightly
    /// exceed 1.0 because the far-field solve is itself parallel: stage
    /// clocks of concurrent orientations overlap when rayon steals
    /// far-field chunks across orientations.
    pub fn parallel_efficiency(&self) -> Option<f64> {
        let capacity = self.wall.solve * self.threads as f64;
        (capacity > 0.0).then(|| self.worker.total / capacity)
    }

    /// Renders the tree-style console breakdown.
    pub fn report(&self) -> String {
        let mut out = String::new();
        let line = |out: &mut String, indent: usize, label: &str, secs: f64| {
            out.push_str(&format!(
                "{:indent$}{:<width$}{:>12.4} s\n",
                "",
                label,
                secs,
                indent = indent,
                width = 28 - indent,
            ));
        };

        out.push_str("Timing Breakdown\n");
        out.push_str("----------------------------------------\n");
        line(&mut out, 0, "Total wall time:", self.total_wall());
        line(&mut out, 2, "Initialization:", self.wall.init);
        if let Some(prognosis) = self.wall.prognosis {
            line(&mut out, 2, "Prognosis:", prognosis);
        }
        line(&mut out, 2, "Solve:", self.wall.solve);
        out.push_str(&format!(
            "    ({} orientations on {} threads — worker times below are summed)\n",
            self.worker.orientations, self.threads
        ));
        line(&mut out, 4, "Near-field:", self.worker.near_field);
        line(&mut out, 4, "Far-field:", self.worker.far_field);
        line(&mut out, 4, "1D integration:", self.worker.mueller_1d);
        line(&mut out, 4, "Parameters:", self.worker.params);
        line(&mut out, 4, "Other:", self.worker.other());
        if self.worker.orientations > 0 {
            out.push_str(&format!(
                "    Mean per orientation: {:.4} ms\n",
                self.worker.total * 1000.0 / self.worker.orientations as f64
            ));
        }
        if let Some(eff) = self.parallel_efficiency() {
            out.push_str(&format!("    Parallel efficiency:  {:.1}%\n", eff * 100.0));
        }
        line(&mut out, 2, "Post-processing:", self.wall.post_process);
        line(&mut out, 2, "File I/O:", self.wall.file_io);
        out
    }

    pub fn print(&self) {
        print!("{}", self.report());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accumulator_sums_across_adds() {
        let acc = TimingAccumulator::default();
        let stages = StageTimings {
            near_field: Duration::from_millis(10),
            far_field: Duration::from_millis(20),
            mueller_1d: Duration::from_millis(1),
            params: Duration::from_millis(2),
        };
        acc.add(&stages, Duration::from_millis(40));
        acc.add(&stages, Duration::from_millis(40));

        let worker = acc.snapshot();
        assert_eq!(worker.orientations, 2);
        assert!((worker.near_field - 0.020).abs() < 1e-9);
        assert!((worker.far_field - 0.040).abs() < 1e-9);
        assert!((worker.total - 0.080).abs() < 1e-9);
        assert!((worker.other() - 0.014).abs() < 1e-9);
    }

    #[test]
    fn efficiency_and_total() {
        let timings = Timings {
            threads: 4,
            wall: WallTimings {
                init: 1.0,
                prognosis: None,
                solve: 10.0,
                post_process: 0.5,
                file_io: 0.5,
            },
            worker: WorkerTimings {
                total: 30.0,
                orientations: 100,
                ..Default::default()
            },
        };
        assert!((timings.total_wall() - 12.0).abs() < 1e-9);
        assert!((timings.parallel_efficiency().unwrap() - 0.75).abs() < 1e-9);
    }
}
