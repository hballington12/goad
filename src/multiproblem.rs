// use std::time::Instant;

use std::io::Write;
use std::{fs::File, io::BufWriter};

use crate::result::{Mueller, MuellerMatrix};
use crate::{
    bins::generate_bins,
    geom::Geom,
    orientation::{Euler, Orientations},
    output,
    problem::{self, Problem},
    result::{self, Results},
    settings::Settings,
};
use anyhow::Result;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
#[cfg(feature = "macroquad")]
use macroquad::prelude::*;
use nalgebra::Complex;
use pyo3::prelude::*;
use rayon::prelude::*;

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
#[pyclass]
#[derive(Debug)] // Added Default derive
pub struct MultiProblem {
    pub geom: Geom,
    pub orientations: Orientations,
    pub settings: Settings, // runtime settings
    pub result: Results,    // averaged result of the problems
}

impl MultiProblem {
    /// Creates a new `MultiProblem` from optional `Geom` and `Settings`.
    /// If settings not provided, loads from config file.
    /// If geom not provided, loads from file using settings.geom_name.
    pub fn new(geom: Option<Geom>, settings: Option<Settings>) -> Self {
        let settings = settings
            .unwrap_or_else(|| crate::settings::load_config().expect("Failed to load config"));
        let mut geom = geom.unwrap_or_else(|| {
            Geom::from_file(&settings.geom_name).expect("Failed to load geometry")
        });

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

    /// Regenerates the orientations for the problem.
    /// Useful for rerunning a random orientation problem with no seed set.
    pub fn regenerate_orientations(&mut self) {
        self.orientations =
            Orientations::generate(&self.settings.orientation.scheme, self.settings.seed);
    }

    /// Resets a `MultiOrientProblem` to its initial state.
    pub fn reset(&mut self) {
        self.result = Results::new_empty(
            &self
                .result
                .field_2d
                .iter()
                .map(|f| f.bin)
                .collect::<Vec<_>>(),
        );
        self.regenerate_orientations();
    }

    /// Solves a `MultiOrientProblem` by averaging over the problems.
    pub fn solve(&mut self) {
        // let start = Instant::now();
        // println!("Solving problem...");

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
                || {
                    let bins = &self.result.bins();
                    Results::new_empty(bins)
                },
                |accum, item| self.reduce_results(accum, item),
            );

        // Normalize results by the number of orientations
        self.normalize_results(self.orientations.num_orientations as f32);
        self.result.try_mueller_to_1d(&self.settings.binning.scheme);
        let _ = self.result.compute_params(self.settings.wavelength);

        // let end = Instant::now();
        // let duration = end.duration_since(start);
        // let time_per_orientation = duration / self.orientations.num_orientations as u32;

        // println!(
        //     "Time taken: {:.2?}, Time per orientation: {:.2?}",
        //     duration, time_per_orientation
        // );

        // // try compute 1d mueller
        // match self.settings.binning.scheme {
        //     Scheme::Custom { .. } => {} // 1d mueller not supported for custom bins
        //     _ => {
        //         match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller) {
        //             Ok((theta, mueller_1d)) => {
        //                 self.result.bins_1d = Some(theta);
        //                 self.result.mueller_1d = Some(mueller_1d);

        //                 // compute params
        //                 let _ = self.result.compute_params(self.settings.wavelength);
        //             }
        //             Err(..) => {}
        //         };
        //         match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller_beam) {
        //             Ok((theta, mueller_1d_beam)) => {
        //                 self.result.bins_1d = Some(theta);
        //                 self.result.mueller_1d_beam = Some(mueller_1d_beam);
        //             }
        //             Err(e) => {
        //                 println!("Failed to compute 1d mueller (beam): {}", e);
        //             }
        //         };
        //         match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller_ext) {
        //             Ok((theta, mueller_1d_ext)) => {
        //                 self.result.bins_1d = Some(theta);
        //                 self.result.mueller_1d_ext = Some(mueller_1d_ext);
        //             }
        //             Err(e) => {
        //                 println!("Failed to compute 1d mueller (ext): {}", e);
        //             }
        //         };
        //     }
        // }

        println!("Results:");
        self.result.print();
    }

    /// Combines two Results objects by adding their fields
    fn reduce_results(&self, mut acc: Results, item: Results) -> Results {
        // Add powers
        acc.powers += item.powers;

        // Add Mueller matrix elements
        for (a, i) in acc.field_2d.iter_mut().zip(item.field_2d.into_iter()) {
            // Handle Mueller matrices
            a.mueller_total = match (a.mueller_total, i.mueller_total) {
                (Some(a_val), Some(i_val)) => Some(a_val + i_val),
                (Some(a_val), None) => Some(a_val),
                (None, Some(i_val)) => Some(i_val),
                (None, None) => None,
            };

            a.mueller_beam = match (a.mueller_beam, i.mueller_beam) {
                (Some(a_val), Some(i_val)) => Some(a_val + i_val),
                (Some(a_val), None) => Some(a_val),
                (None, Some(i_val)) => Some(i_val),
                (None, None) => None,
            };

            a.mueller_ext = match (a.mueller_ext, i.mueller_ext) {
                (Some(a_val), Some(i_val)) => Some(a_val + i_val),
                (Some(a_val), None) => Some(a_val),
                (None, Some(i_val)) => Some(i_val),
                (None, None) => None,
            };
        }

        acc
    }

    /// Normalizes the results by dividing by the number of orientations
    fn normalize_results(&mut self, num_orientations: f32) {
        // Powers
        self.result.powers /= num_orientations;
        println!("powers: {:?}", self.result.powers);

        for field in self.result.field_2d.iter_mut() {
            // Amplitude Matrices - divide by complex representation
            let div_c = Complex::from(num_orientations);
            if let Some(ampl) = field.ampl_total.as_mut() {
                *ampl /= div_c;
            }
            if let Some(ampl) = field.ampl_beam.as_mut() {
                *ampl /= div_c;
            }
            if let Some(ampl) = field.ampl_ext.as_mut() {
                *ampl /= div_c;
            }

            // Mueller Matrices - divide by real value
            if let Some(mueller) = field.mueller_total.as_mut() {
                *mueller /= num_orientations;
            }
            if let Some(mueller) = field.mueller_beam.as_mut() {
                *mueller /= num_orientations;
            }
            if let Some(mueller) = field.mueller_ext.as_mut() {
                *mueller /= num_orientations;
            }
        }
    }

    pub fn writeup(&self) {
        // Helper closure to write Mueller matrices to files
        let write_mueller_grid =
            |file_suffix: &str,
             mueller_getter: &dyn Fn(&result::ScattResult2D) -> Option<Mueller>|
             -> Result<()> {
                let file_name = format!("mueller_scatgrid_{}", file_suffix);
                let path = output::output_path(Some(&self.settings.directory), &file_name)?;
                let file = File::create(&path)?;
                let mut writer = BufWriter::new(file);

                for result in &self.result.field_2d {
                    let bin = result.bin;
                    write!(writer, "{} {} ", bin.theta_bin.center, bin.phi_bin.center)?;

                    if let Some(mueller) = mueller_getter(result) {
                        for element in mueller.to_vec() {
                            write!(writer, "{} ", element)?;
                        }
                    }
                    writeln!(writer)?;
                }
                Ok(())
            };

        // Write all three Mueller matrix types
        let _ = write_mueller_grid("", &|r| r.mueller_total);
        let _ = write_mueller_grid("beam", &|r| r.mueller_beam);
        let _ = write_mueller_grid("ext", &|r| r.mueller_ext);

        // Write generic results
        let _ = output::write_result(&self.result, &self.settings.directory);

        // Write 1d results (todo)
        // todo!("Implement 1d results writing")
    }
}

#[pymethods]
impl MultiProblem {
    #[new]
    #[pyo3(signature = (settings, geom = None))]
    fn py_new(settings: Settings, geom: Option<Geom>) -> PyResult<Self> {
        // Load geometry from file if not provided
        let mut geom = match geom {
            Some(g) => g,
            None => Geom::from_file(&settings.geom_name).map_err(|e| {
                pyo3::exceptions::PyFileNotFoundError::new_err(format!(
                    "Failed to load geometry file '{}': {}",
                    settings.geom_name, e
                ))
            })?,
        };

        problem::init_geom(&settings, &mut geom);

        let orientations = Orientations::generate(&settings.orientation.scheme, settings.seed);
        let bins = generate_bins(&settings.binning.scheme);
        let result = Results::new_empty(&bins);

        Ok(Self {
            geom,
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
    pub fn py_solve(&mut self, py: Python) -> PyResult<()> {
        py.allow_threads(|| {
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

    /// Python wrapper for writeup method
    pub fn py_writeup(&self) -> PyResult<()> {
        let _ = self.writeup();
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
