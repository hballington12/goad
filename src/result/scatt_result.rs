use std::fmt::Debug;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use crate::bins::AngleBin;
use crate::bins::SolidAngleBin;
use crate::convergence::Convergeable;
use nalgebra::Complex;
use rand_distr::num_traits::Pow;

use super::mueller::{Ampl, Mueller};

/// Trait for different types of scattering bins (1D or 2D)
pub trait ScatteringBin: Clone + Debug {
    /// Get the theta center value
    fn theta_center(&self) -> f32;

    /// Get the theta bin
    fn theta_bin(&self) -> &AngleBin;

    /// Check if this bin has the same theta as another
    fn same_theta(&self, other: &Self) -> bool {
        self.theta_bin() == other.theta_bin()
    }
}

impl ScatteringBin for SolidAngleBin {
    fn theta_center(&self) -> f32 {
        self.theta.center
    }

    fn theta_bin(&self) -> &AngleBin {
        &self.theta
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

impl<B: ScatteringBin> Pow<f32> for ScattResult<B> {
    type Output = Self;

    fn pow(self, rhs: f32) -> Self {
        Self {
            bin: self.bin,
            ampl_total: self.ampl_total.map(|c| c.powf(rhs)),
            ampl_beam: self.ampl_beam.map(|c| c.powf(rhs)),
            ampl_ext: self.ampl_ext.map(|c| c.powf(rhs)),
            mueller_total: self.mueller_total.map(|m| m.powf(rhs)),
            mueller_beam: self.mueller_beam.map(|m| m.powf(rhs)),
            mueller_ext: self.mueller_ext.map(|m| m.powf(rhs)),
        }
    }
}

impl<B: ScatteringBin> Mul<f32> for ScattResult<B> {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        Self {
            bin: self.bin,
            ampl_total: self.ampl_total * Complex::from(rhs),
            ampl_beam: self.ampl_beam * Complex::from(rhs),
            ampl_ext: self.ampl_ext * Complex::from(rhs),
            mueller_total: self.mueller_total * rhs,
            mueller_beam: self.mueller_beam * rhs,
            mueller_ext: self.mueller_ext * rhs,
        }
    }
}

impl<B: ScatteringBin> Mul for ScattResult<B> {
    type Output = ScattResult<B>;

    fn mul(self, other: ScattResult<B>) -> Self::Output {
        ScattResult {
            bin: self.bin,
            ampl_total: self.ampl_total * other.ampl_total,
            ampl_beam: self.ampl_beam * other.ampl_beam,
            ampl_ext: self.ampl_ext * other.ampl_ext,
            mueller_total: self.mueller_total * other.mueller_total,
            mueller_beam: self.mueller_beam * other.mueller_beam,
            mueller_ext: self.mueller_ext * other.mueller_ext,
        }
    }
}

impl<B: ScatteringBin> Add for ScattResult<B> {
    type Output = ScattResult<B>;

    fn add(self, other: ScattResult<B>) -> Self::Output {
        ScattResult {
            bin: self.bin,
            ampl_total: self.ampl_total + other.ampl_total,
            ampl_beam: self.ampl_beam + other.ampl_beam,
            ampl_ext: self.ampl_ext + other.ampl_ext,
            mueller_total: self.mueller_total + other.mueller_total,
            mueller_beam: self.mueller_beam + other.mueller_beam,
            mueller_ext: self.mueller_ext + other.mueller_ext,
        }
    }
}

impl<B: ScatteringBin> Sub for ScattResult<B> {
    type Output = ScattResult<B>;

    fn sub(self, other: ScattResult<B>) -> Self::Output {
        ScattResult {
            bin: self.bin,
            ampl_total: self.ampl_total - other.ampl_total,
            ampl_beam: self.ampl_beam - other.ampl_beam,
            ampl_ext: self.ampl_ext - other.ampl_ext,
            mueller_total: self.mueller_total - other.mueller_total,
            mueller_beam: self.mueller_beam - other.mueller_beam,
            mueller_ext: self.mueller_ext - other.mueller_ext,
        }
    }
}

impl<B: ScatteringBin> Div<f32> for ScattResult<B> {
    type Output = Self;

    fn div(self, other: f32) -> Self {
        Self {
            bin: self.bin,
            ampl_total: self.ampl_total / Complex::from(other),
            ampl_beam: self.ampl_beam / Complex::from(other),
            ampl_ext: self.ampl_ext / Complex::from(other),
            mueller_total: self.mueller_total / other,
            mueller_beam: self.mueller_beam / other,
            mueller_ext: self.mueller_ext / other,
        }
    }
}

impl<B: ScatteringBin> Div for ScattResult<B> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Self {
            bin: self.bin,
            ampl_total: self.ampl_total.component_div(&other.ampl_total),
            ampl_beam: self.ampl_beam.component_div(&other.ampl_beam),
            ampl_ext: self.ampl_ext.component_div(&other.ampl_ext),
            mueller_total: self.mueller_total.component_div(&other.mueller_total),
            mueller_beam: self.mueller_beam.component_div(&other.mueller_beam),
            mueller_ext: self.mueller_ext.component_div(&other.mueller_ext),
        }
    }
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

    /// Returns a ScattResult with all values set to 1.0 (for weights)
    pub fn ones_like(&self) -> Self {
        Self {
            bin: self.bin.clone(),
            ampl_total: Ampl::from_element(Complex::new(1.0, 0.0)),
            ampl_beam: Ampl::from_element(Complex::new(1.0, 0.0)),
            ampl_ext: Ampl::from_element(Complex::new(1.0, 0.0)),
            mueller_total: Mueller::from_element(1.0),
            mueller_beam: Mueller::from_element(1.0),
            mueller_ext: Mueller::from_element(1.0),
        }
    }
}

impl<B: ScatteringBin> Convergeable for ScattResult<B> {
    fn zero_like(&self) -> Self {
        ScattResult::new(self.bin.clone())
    }

    fn weighted_add(&self, other: &Self, w1: f32, w2: f32) -> Self {
        // Simple weighted average by count
        let total = w1 + w2;
        (self.clone() * w1 + other.clone() * w2) / total
    }

    fn mul_elem(&self, other: &Self) -> Self {
        self.clone() * other.clone()
    }

    fn div_elem(&self, other: &Self) -> Self {
        self.clone() / other.clone()
    }

    fn add_elem(&self, other: &Self) -> Self {
        self.clone() + other.clone()
    }

    fn sub_elem(&self, other: &Self) -> Self {
        self.clone() - other.clone()
    }

    fn scale(&self, scalar: f32) -> Self {
        self.clone() * scalar
    }

    fn sqrt_elem(&self) -> Self {
        self.clone().pow(0.5)
    }

    fn to_weighted(&self) -> Self {
        // ScattResult doesn't need special weighting
        self.clone()
    }

    fn weights(&self) -> Self {
        // All weights are 1.0
        self.ones_like()
    }
}

/// Type alias for 2D scattering results (full solid angle)
pub type ScattResult2D = ScattResult<SolidAngleBin>;

/// Type alias for 1D scattering results (theta only)
pub type ScattResult1D = ScattResult<AngleBin>;
