//! # convergence.rs
//!
//! This module solves a multi-orientation problem by averaging over each
//! orientation until the desired results have converged within some tolerance.
//! In theory, one can choose any number of outputs to use for awaiting
//! convergence, but in practice, the most useful are:
//! - asymmetry parameter
//! - scattering cross section
//! - single scattering albedo
//! - parts of the $S_{11}$ 1D Mueller matrix

use crate::{
    bins::generate_bins, geom::Geom, multiproblem::MultiProblem, orientation::Orientation,
    orientation::Scheme, problem::init_geom, result::Results, settings::Settings,
};
use ndarray::Array2; // Added for Array2 type in trait

#[derive(Debug)]
pub enum ConvergenceType {
    Asymmetry,
    ScatteringCrossSection,
    SingleScatteringAlbedo,
    Mueller1D,
}

/// Trait for structs that can provide data for convergence checking.
pub trait ConvergenceDataSource {
    fn get_asymmetry_parameter(&self) -> Option<f32>;
    fn get_scattering_cross_section(&self) -> Option<f32>;
    fn get_single_scattering_albedo(&self) -> Option<f32>;
    fn get_mueller_1d(&self) -> Option<&Array2<f32>>;
    // Potentially, more specific methods for parts of Mueller matrix if needed
}

/// Struct to hold data from a problem solved by convergence.
#[derive(Debug)]
pub struct Convergence {
    /// Runtime settings
    pub settings: Settings,
    /// Geometry of the problem
    pub geom: Geom,
    /// Batch of results to check
    // pub batch: T,
    /// Size of the batch
    pub batch_size: usize,
    /// Types of convergence to check
    pub convergence_types: Vec<ConvergenceType>,
    /// Maximum number of iterations
    pub max_iterations: usize,
    /// The averaged result of the problems
    pub result: Results,
}

// Example of how you might implement a method in Convergence
// (This is a conceptual example, actual implementation details might vary)
impl Convergence {
    pub fn new(
        settings: Settings,
        // batch: T,
        batch_size: usize,
        convergence_types: Vec<ConvergenceType>,
        max_iterations: usize,
    ) -> Self {
        let mut geom = Geom::from_file(&settings.geom_name).unwrap();
        init_geom(&settings, &mut geom);
        let bins = generate_bins(&settings.binning.scheme);
        let result = Results::new_empty(&bins);
        Self {
            settings,
            geom,
            // batch,
            batch_size,
            convergence_types,
            max_iterations,
            result,
        }
    }

    // pub fn check_convergence(&self) -> bool {
    //     // Example: Iterate through convergence_types and check values from self.result
    //     for conv_type in &self.convergence_types {
    //         match conv_type {
    //             ConvergenceType::Asymmetry => {
    //                 if let Some(asym) = self.batch.get_asymmetry_parameter() {
    //                     // Perform convergence check for asymmetry
    //                     println!("Checking asymmetry: {}", asym);
    //                     // ... logic to compare with previous values or a threshold ...
    //                 } else {
    //                     // Handle case where asymmetry parameter is not available
    //                     return false; // Or some other error handling
    //                 }
    //             }
    //             ConvergenceType::ScatteringCrossSection => {
    //                 if let Some(scat) = self.batch.get_scattering_cross_section() {
    //                     // Perform convergence check for scattering cross section
    //                     println!("Checking scattering cross section: {}", scat);
    //                     // ... logic ...
    //                 } else {
    //                     return false;
    //                 }
    //             }
    //             // ... handle other ConvergenceType variants ...
    //             _ => {}
    //         }
    //     }
    //     // Return true if all specified types have converged
    //     true // Placeholder
    // }

    pub fn run(&mut self) {
        // set the settings for batches
        let mut batch_settings = self.settings.clone();
        // set the orientation scheme to be random uniform with the batch size
        batch_settings.orientation = Orientation {
            scheme: Scheme::Uniform {
                num_orients: self.batch_size,
            },
            euler_convention: self.settings.orientation.euler_convention,
        };
        batch_settings.seed = None; // Seed must be None to allow randomness across batches

        // Create a new multiproblem with settings according to batch size
        let mut multiproblem = MultiProblem::new(batch_settings);

        for i in 0..self.max_iterations {
            // Run the simulation or computation
            // Update self.batch with new results
            println!("Running simulation... Iteration: {}", i);

            // Ensure the multiproblem is reset for each iteration
            multiproblem.reset();

            multiproblem.solve();

            // // Check convergence
            // if self.check_convergence() {
            //     println!("Convergence achieved!");
            //     break;
            // } else {
            //     println!("Not converged yet, continuing...");
            // }
        }
    }
}
