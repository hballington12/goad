//! Scattering result types for 1D and 2D angular distributions.

use std::fmt::Debug;

use crate::bins::{AngleBin, SolidAngleBin};

use super::mueller::{Ampl, Mueller};

/// Trait for different types of scattering bins (1D or 2D).
pub trait ScatteringBin: Clone + Debug {
    /// Get the theta center value.
    fn theta_center(&self) -> f32;

    /// Get the theta bin.
    fn theta_bin(&self) -> &AngleBin;

    /// Check if this bin has the same theta as another.
    fn same_theta(&self, other: &Self) -> bool {
        self.theta_bin() == other.theta_bin()
    }
}

impl ScatteringBin for SolidAngleBin {
    fn theta_center(&self) -> f32 {
        self.theta_bin.center
    }

    fn theta_bin(&self) -> &AngleBin {
        &self.theta_bin
    }
}

impl ScatteringBin for AngleBin {
    fn theta_center(&self) -> f32 {
        self.center
    }

    fn theta_bin(&self) -> &AngleBin {
        self
    }
}

/// A generic far-field scattering result that can be 1D or 2D.
#[derive(Debug, Clone)]
pub struct ScattResult<B: ScatteringBin> {
    pub bin: B,
    pub ampl_total: Ampl,
    pub ampl_beam: Ampl,
    pub ampl_ext: Ampl,
    pub mueller_total: Mueller,
    pub mueller_beam: Mueller,
    pub mueller_ext: Mueller,
}

impl<B: ScatteringBin> ScattResult<B> {
    /// Creates a new empty ScattResult.
    pub fn new(bin: B) -> Self {
        Self {
            bin,
            ampl_total: Ampl::zeros(),
            ampl_beam: Ampl::zeros(),
            ampl_ext: Ampl::zeros(),
            mueller_total: Mueller::zeros(),
            mueller_beam: Mueller::zeros(),
            mueller_ext: Mueller::zeros(),
        }
    }
}

/// Type alias for 2D scattering results (full solid angle).
pub type ScattResult2D = ScattResult<SolidAngleBin>;

/// Type alias for 1D scattering results (theta only).
pub type ScattResult1D = ScattResult<AngleBin>;
