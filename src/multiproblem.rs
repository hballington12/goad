//! Multi-orientation simulation orchestration and ensemble averaging.
//!
//! This module manages electromagnetic scattering simulations across multiple
//! particle orientations to compute ensemble-averaged scattering properties.
//! It provides parallel computation, progress tracking, and proper statistical
//! averaging for both discrete orientation sets and random orientation sampling.
//!
//! The multi-orientation system provides:
//! - Parallel orientation processing with rayon
//! - Progress tracking for long-running calculations
//! - Statistical averaging with proper normalization
//! - Result aggregation across orientations
//! - Comprehensive output for ensemble analysis
//! - Memory-efficient streaming computation
//!
//! # Key Features
//!
//! - [`MultiProblem`]: Main orchestrator for multi-orientation simulations
//! - Parallel execution with progress bars
//! - On-the-fly result reduction for memory efficiency
//! - Automatic 1D integration for symmetric cases
//! - Complete output file generation
//! - Performance timing and analysis

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
use pyo3::prelude::*;
use rayon::prelude::*;

/// Multi-orientation scattering simulation with orientation averaging.
/// 
/// **Context**: Many scattering applications require orientation-averaged results
/// to model ensembles of randomly oriented particles or systematic orientation
/// studies. This requires running multiple single-orientation simulations and
/// properly averaging the results with parallel computation for efficiency.
/// 
/// **How it Works**: Manages a set of orientations, runs parallel simulations
/// for each orientation using the same base geometry, and aggregates results
/// through proper averaging. Handles progress tracking, result normalization,
/// and output file generation for the ensemble-averaged scattering patterns.
#[pyclass]
#[derive(Debug)] // Added Default derive
pub struct MultiProblem {
    pub geom: Geom,
    pub orientations: Orientations,
    pub settings: Settings, // runtime settings
    pub result: Results,    // averaged result of the problems
}

impl MultiProblem {
    /// Creates a multi-orientation simulation from configuration settings.
    /// 
    /// **Context**: Multi-orientation simulations require initialization of
    /// the base geometry, orientation sets, and result storage before
    /// simulation execution.
    /// 
    /// **How it Works**: Loads geometry from file, applies initial transformations,
    /// generates the orientation set from the specified scheme, and allocates
    /// result storage for the expected angular bins.
    pub fn new(geom: Option<Geom>, settings: Option<Settings>) -> Self {
        let settings = settings.unwrap_or_else(|| crate::settings::load_config().expect("Failed to load config"));
        let mut geom = geom.unwrap_or_else(|| Geom::from_file(&settings.geom_name).expect("Failed to load geometry"));

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

    /// Regenerates orientation set for repeated random sampling.
    /// 
    /// **Context**: Statistical convergence studies may require multiple
    /// independent realizations of random orientation sets to assess
    /// convergence and estimate uncertainties.
    /// 
    /// **How it Works**: Generates a new orientation set using the same
    /// scheme parameters but different random values if no seed is set.
    pub fn regenerate_orientations(&mut self) {
        self.orientations =
            Orientations::generate(&self.settings.orientation.scheme, self.settings.seed);
    }

    /// Resets simulation to initial state for rerunning.
    /// 
    /// **Context**: Parameter studies or convergence analysis may require
    /// resetting the simulation state while preserving the configuration.
    /// 
    /// **How it Works**: Clears accumulated results and optionally regenerates
    /// orientations to prepare for a fresh simulation run.
    pub fn reset(&mut self) {
        self.result = Results::new_empty(&self.result.bins);
        self.regenerate_orientations();
    }

    /// Executes parallel multi-orientation simulation with progress tracking.
    /// 
    /// **Context**: Multi-orientation simulations involve substantial computational
    /// work that benefits from parallel execution. Progress tracking provides
    /// user feedback for long-running calculations.
    /// 
    /// **How it Works**: Creates a base problem template, runs parallel simulations
    /// for each orientation using rayon, aggregates results on-the-fly using
    /// reduction, normalizes by orientation count, and post-processes to compute
    /// 1D results and integral parameters.
    pub fn solve(&mut self) {
        let start = Instant::now();
        println!("Solving problem...");

        // init a base problem that can be reset
        let problem_base = Problem::new(Some(self.geom.clone()), Some(self.settings.clone()));

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

    /// Aggregates simulation results from individual orientations.
    /// 
    /// **Context**: Parallel reduction requires a function to combine results
    /// from different orientation simulations. This must handle all result
    /// components including amplitude matrices, Mueller matrices, and power data.
    /// 
    /// **How it Works**: Element-wise addition of all result arrays and
    /// structures, accumulating contributions from each orientation simulation
    /// for later normalization by orientation count.
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

    /// Normalizes accumulated results by the number of orientations.
    /// 
    /// **Context**: Orientation averaging requires dividing accumulated
    /// results by the number of orientations to obtain proper ensemble
    /// averages. This applies to all result components.
    /// 
    /// **How it Works**: Divides all accumulated amplitude matrices, Mueller
    /// matrices, and power values by the orientation count to get ensemble
    /// averages representing the orientation-averaged scattering properties.
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

    /// Writes all simulation results to output files.
    /// 
    /// **Context**: Multi-orientation results include multiple data types
    /// (2D and 1D Mueller matrices, amplitude data, summary statistics)
    /// that need to be written in appropriate formats for analysis.
    /// 
    /// **How it Works**: Calls output functions for each result type,
    /// handling both successful writes and cases where certain results
    /// (like 1D data) may not be available.
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

#[pymethods]
impl MultiProblem {
    #[new]
    #[pyo3(signature = (settings = None, geom = None))]
    fn py_new(settings: Option<Settings>, geom: Option<Geom>) -> Self {
        MultiProblem::new(geom, settings)
    }

    /// Python wrapper for solve method
    pub fn py_solve(&mut self) -> PyResult<()> {
        self.solve();
        Ok(())
    }

    /// Get the results object (same pattern as Problem)
    #[getter]
    pub fn get_results(&self) -> Results {
        self.result.clone()
    }

    /// Python wrapper for writeup method
    pub fn py_writeup(&self) -> PyResult<()> {
        self.writeup();
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
