use chrono::Local;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::time::Duration;

use crate::params::Param;
use crate::settings::constants::MIN_ORIENTATIONS;

pub(crate) struct ConvergenceProgress {
    title_pb: ProgressBar,
    info_pb: ProgressBar,
    target_pbs: Vec<ProgressBar>,
    start_time: std::time::Instant,
    max_target: usize,
}

impl ConvergenceProgress {
    pub fn new(num_targets: usize, max_target: usize) -> Self {
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

        let target_pbs: Vec<ProgressBar> = (0..num_targets)
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

        Self {
            title_pb,
            info_pb,
            target_pbs,
            start_time: std::time::Instant::now(),
            max_target,
        }
    }

    pub fn set_running(&self) {
        self.title_pb
            .set_message("[Status: \x1b[32mRUNNING\x1b[0m]");
    }

    pub fn set_finalising(&self) {
        self.title_pb
            .set_message("[Status: \x1b[33mFINALISING\x1b[0m]");
    }

    pub fn update_info(&self, count: usize) {
        let elapsed = self.start_time.elapsed().as_secs_f64();
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
        self.info_pb.set_message(format!(
            "[Orientations: {} ({}{}{}|{})] [{:.3} sec/orientation]",
            count, min_color, MIN_ORIENTATIONS, "\x1b[0m", self.max_target, sec_per_orient
        ));
        self.title_pb
            .set_prefix(Local::now().format("%Y-%m-%d %H:%M:%S").to_string());
    }

    pub fn update_target(
        &self,
        index: usize,
        param: Param,
        mean_val: f32,
        sem_val: f32,
        target_rel_sem: f32,
    ) {
        let current_rel_sem = if mean_val.abs() > 1e-10 {
            sem_val / mean_val.abs()
        } else {
            0.0
        };

        // Progress with sqrt scaling, capped at 100%
        let progress = if current_rel_sem > 1e-10 {
            ((target_rel_sem / current_rel_sem).sqrt()).min(1.0)
        } else {
            1.0
        };

        let param_name = format!("{:?}", param);
        self.target_pbs[index].set_message(format!(
            "{:<22} {:>10.4e} ± {:<10.4e} [{:>5.2}% / {:>5.2}%]",
            param_name,
            mean_val,
            sem_val,
            current_rel_sem * 100.0,
            target_rel_sem * 100.0
        ));
        self.target_pbs[index].set_position((progress * 100.0) as u64);
    }

    pub fn finish(&self) {
        self.title_pb.finish_and_clear();
        for pb in &self.target_pbs {
            pb.finish_and_clear();
        }
        self.info_pb.finish_and_clear();
    }
}
