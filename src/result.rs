//! Result types and Mueller matrix handling for GOAD light scattering simulations.
//!
//! This module contains:
//! - `GOComponent` - Enum for total/beam/ext-diff components
//! - `Mueller`, `Ampl` - Matrix types for Mueller and amplitude matrices
//! - `ScattResult` - Generic scattering result type (1D or 2D)
//! - `Results` - Complete simulation results with zones integration
//! - Helper functions for integration

mod component;
mod integrate;
mod mueller;
mod python;
mod results;
mod scatt_result;

// Re-export public types
pub use component::GOComponent;
pub use integrate::integrate_theta_weighted_component;
pub use mueller::{Ampl, ApproxEq, Mueller, MuellerMatrix};
pub use results::Results;
pub use scatt_result::{ScattResult, ScattResult1D, ScattResult2D, ScatteringBin};
