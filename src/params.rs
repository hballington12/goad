//! Integral scattering parameters for comparative analysis.
//!
//! This module defines the storage and computation of key integral parameters
//! that characterize electromagnetic scattering behavior. These parameters
//! provide compact representations of complex angular scattering patterns
//! for comparison with analytical theories, experimental data, and radiative
//! transfer models.
//!
//! The parameter system provides:
//! - Scattering and extinction cross sections
//! - Single-scattering albedo
//! - Asymmetry parameter for angular distribution characterization
//! - Optional parameter storage for various analysis scenarios
//! - Integration with result processing pipelines
//!
//! # Key Parameters
//!
//! - **Scattering cross section**: Total scattered power normalized by incident intensity
//! - **Extinction cross section**: Total removed power (scattered + absorbed)
//! - **Albedo**: Ratio of scattering to extinction (probability of scattering)
//! - **Asymmetry parameter**: Average cosine of scattering angle

/// Integral scattering parameters derived from angular distributions.
/// 
/// **Context**: While full angular scattering patterns provide detailed information,
/// many applications require integrated parameters that characterize the overall
/// scattering behavior. These parameters enable comparison with analytical theories,
/// radiative transfer models, and experimental measurements.
/// 
/// **How it Works**: Stores key integral parameters as optional values since they
/// may not be computable for all simulation configurations (e.g., 2D patterns
/// without azimuthal symmetry). Parameters include the asymmetry parameter for
/// angular distribution characterization, cross sections for scattering strength,
/// and albedo for the scattering-to-extinction ratio.
#[derive(Debug, PartialEq, Clone)]
pub struct Params {
    pub asymettry: Option<f32>,
    pub scat_cross: Option<f32>,
    pub ext_cross: Option<f32>,
    pub albedo: Option<f32>,
}

impl Params {
    /// Creates empty parameter storage.
    /// 
    /// **Context**: Parameters are computed after simulation completion from
    /// the angular scattering results, starting with no values.
    /// 
    /// **How it Works**: Initializes all parameters as None, to be populated
    /// by post-processing calculations.
    pub fn new() -> Self {
        Self {
            asymettry: None,
            scat_cross: None,
            ext_cross: None,
            albedo: None,
        }
    }
}
