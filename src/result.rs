//! Result types for GOAD scattering simulations.
//!
//! This module contains the complete output from a light scattering simulation,
//! including Mueller matrices, amplitude matrices, and derived parameters.
//!
//! # Module Structure
//!
//! - [`component`]: GOComponent enum for distinguishing scattering contributions
//! - [`mueller`]: Mueller matrix types and traits
//! - [`scatt_result`]: Scattering result types for 1D and 2D distributions
//! - [`results`]: Main Results struct
//! - [`integrate`]: Integration helpers
//! - [`python`]: Python bindings

mod component;
mod integrate;
mod mueller;
mod python;
mod results;
mod scatt_result;

// Re-export main types at module root
pub use component::GOComponent;
pub use integrate::integrate_theta_weighted_component;
pub use mueller::{Ampl, ApproxEq, Mueller, MuellerMatrix};
pub use python::collect_mueller;
pub use results::Results;
pub use scatt_result::{ScattResult, ScattResult1D, ScattResult2D, ScatteringBin};
