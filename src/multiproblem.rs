use std::time::Instant;

use crate::{
    bins::{generate_bins, Scheme},
    geom::Geom,
    orientation::{Euler, Orientations},
    output,
    problem::{self, Problem},
    result::{self, Results},
    settings::Settings,
};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use macroquad::prelude::*;
use nalgebra::Complex;
use rayon::prelude::*;

/// A problem for a single geometry with multiple orientations.
#[derive(Debug)] // Added Default derive
pub struct MultiProblem {
    pub geom: Geom,
    pub orientations: Orientations,
    pub settings: Settings, // runtime settings
    pub result: Results,    // averaged result of the problems
}

impl MultiProblem {
    /// Creates a new `MultiOrientProblem` from a `settings: Settings` configuration.
    pub fn new(settings: Settings) -> Self {
        let mut geom = Geom::from_file(&settings.geom_name).unwrap();

        problem::init_geom(&settings, &mut geom);

        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let bins = generate_bins(&settings.binning.scheme);
        let result = Results::new_empty(&bins);

        Self {
            geom,
            orientations,
            settings,
            result,
        }
    }

    /// Resets a `MultiOrientProblem` to its initial state.
    pub fn reset(&mut self) {
        self.result = Results::new_empty(&self.result.bins);
    }

    /// Solves a `MultiOrientProblem` by averaging over the problems.
    pub fn solve(&mut self) {
        let start = Instant::now();
        println!("Solving problem...");

        // init a base problem that can be reset
        let problem_base = Problem::new(self.geom.clone(), Some(self.settings.clone()));

        let m = MultiProgress::new();
        let n = self.orientations.num_orientations;
        let pb = m.add(ProgressBar::new(n as u64));
        pb.set_style(
            ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] {bar:40.green/blue} {pos:>5}/{len:5} {msg} ETA: {eta_precise}",
            )
            .unwrap()
            .progress_chars("█▇▆▅▄▃▂▁")
        );
        pb.set_message("orientation".to_string());

        // Solve for each orientation and reduce results on the fly
        self.result = self
            .orientations
            .eulers
            .par_iter()
            .map(|(a, b, g)| {
                let mut problem = problem_base.clone();
                let euler = Euler::new(*a, *b, *g);

                problem.run(Some(&euler)); // run the problem with an euler rotation

                pb.inc(1);
                problem.result
            })
            .reduce(
                || Results::new_empty(&self.result.bins),
                |accum, item| self.reduce_results(accum, item),
            );

        // Normalize results by the number of orientations
        self.normalize_results(self.orientations.num_orientations as f32);

        let end = Instant::now();
        let duration = end.duration_since(start);
        let time_per_orientation = duration / self.orientations.num_orientations as u32;

        println!(
            "Time taken: {:.2?}, Time per orientation: {:.2?}",
            duration, time_per_orientation
        );

        // try compute 1d mueller
        match self.settings.binning.scheme {
            Scheme::Custom { .. } => {} // 1d mueller not supported for custom bins
            _ => {
                match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller) {
                    Ok((theta, mueller_1d)) => {
                        self.result.bins_1d = Some(theta);
                        self.result.mueller_1d = Some(mueller_1d);

                        // compute params
                        let _ = self.result.compute_params(self.settings.wavelength);
                    }
                    Err(..) => {}
                };
                match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller_beam) {
                    Ok((theta, mueller_1d_beam)) => {
                        self.result.bins_1d = Some(theta);
                        self.result.mueller_1d_beam = Some(mueller_1d_beam);
                    }
                    Err(e) => {
                        println!("Failed to compute 1d mueller (beam): {}", e);
                    }
                };
                match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller_ext) {
                    Ok((theta, mueller_1d_ext)) => {
                        self.result.bins_1d = Some(theta);
                        self.result.mueller_1d_ext = Some(mueller_1d_ext);
                    }
                    Err(e) => {
                        println!("Failed to compute 1d mueller (ext): {}", e);
                    }
                };
            }
        }

        println!("Results:");
        self.result.print();
    }

    /// Combines two Results objects by adding their fields
    fn reduce_results(&self, mut acc: Results, item: Results) -> Results {
        // Add Mueller matrix elements
        for (a, i) in acc.mueller.iter_mut().zip(item.mueller.iter()) {
            *a += i;
        }

        // Add powers
        acc.powers += item.powers;

        // Add amplitude matrices if they exist
        for (a, i) in acc.ampl.iter_mut().zip(item.ampl.iter()) {
            *a += i;
        }

        for (a, i) in acc.ampl_beam.iter_mut().zip(item.ampl_beam.iter()) {
            *a += i;
        }

        for (a, i) in acc.ampl_ext.iter_mut().zip(item.ampl_ext.iter()) {
            *a += i;
        }

        for (a, i) in acc.mueller_beam.iter_mut().zip(item.mueller_beam.iter()) {
            *a += i;
        }

        for (a, i) in acc.mueller_ext.iter_mut().zip(item.mueller_ext.iter()) {
            *a += i;
        }

        acc
    }

    /// Normalizes the results by dividing by the number of orientations
    fn normalize_results(&mut self, num_orientations: f32) {
        // Normalize powers
        self.result.powers /= num_orientations;

        for ampl in self.result.ampl.iter_mut() {
            *ampl /= Complex::new(num_orientations, 0.0);
        }

        for ampl in self.result.ampl_beam.iter_mut() {
            *ampl /= Complex::new(num_orientations, 0.0);
        }

        for ampl in self.result.ampl_ext.iter_mut() {
            *ampl /= Complex::new(num_orientations, 0.0);
        }

        for mut row in self.result.mueller.outer_iter_mut() {
            for val in row.iter_mut() {
                *val /= num_orientations;
            }
        }

        for mut row in self.result.mueller_beam.outer_iter_mut() {
            for val in row.iter_mut() {
                *val /= num_orientations;
            }
        }

        for mut row in self.result.mueller_ext.outer_iter_mut() {
            for val in row.iter_mut() {
                *val /= num_orientations;
            }
        }
    }

    pub fn writeup(&self) {
        // Write 2D mueller matrices
        let _ = output::write_mueller(
            &self.result.bins,
            &self.result.mueller,
            "",
            &self.settings.directory,
        );
        let _ = output::write_mueller(
            &self.result.bins,
            &self.result.mueller_beam,
            "_beam",
            &self.settings.directory,
        );
        let _ = output::write_mueller(
            &self.result.bins,
            &self.result.mueller_ext,
            "_ext",
            &self.settings.directory,
        );

        // Write generic results
        let _ = output::write_result(&self.result, &self.settings.directory);

        // (Try to) write 1D mueller matrices
        match self.result.mueller_1d {
            Some(ref mueller_1d) => {
                let _ = output::write_mueller_1d(
                    &self.result.bins_1d.as_ref().unwrap(),
                    mueller_1d,
                    "",
                    &self.settings.directory,
                );
            }
            None => {
                println!("Failed to write 1D mueller matrix");
            }
        }
        match self.result.mueller_1d_beam {
            Some(ref mueller_1d) => {
                let _ = output::write_mueller_1d(
                    &self.result.bins_1d.as_ref().unwrap(),
                    mueller_1d,
                    "_beam",
                    &self.settings.directory,
                );
            }
            None => {
                println!("Failed to write 1D mueller matrix (beam)");
            }
        }
        match self.result.mueller_1d_ext {
            Some(ref mueller_1d) => {
                let _ = output::write_mueller_1d(
                    &self.result.bins_1d.as_ref().unwrap(),
                    mueller_1d,
                    "_ext",
                    &self.settings.directory,
                );
            }
            None => {
                println!("Failed to write 1D mueller matrix (ext)");
            }
        }
    }
}
