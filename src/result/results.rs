use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use crate::bins::Scheme;
use crate::bins::SolidAngleBin;
use crate::convergence::Convergeable;
use crate::params::Params;
use crate::powers::Powers;
use crate::zones::{Zone, ZoneType, Zones};
use itertools::Itertools;
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rand_distr::num_traits::Pow;

use super::component::GOComponent;
use super::scatt_result::{ScattResult1D, ScattResult2D};

/// Complete results from a GOAD light scattering simulation.
///
/// Contains all computed scattering data including Mueller matrices,
/// amplitude matrices, power distributions, and derived parameters.
/// Supports both 2D angular distributions and 1D integrated results.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
#[derive(Debug, Clone)]
pub struct Results {
    pub zones: Zones,
    pub powers: Powers,
    pub params: Params,
}

impl AddAssign for Results {
    fn add_assign(&mut self, other: Self) {
        *self = self.clone() + other;
    }
}

impl Pow<f32> for Results {
    type Output = Self;

    fn pow(self, rhs: f32) -> Self {
        Self {
            zones: self.zones.pow(rhs),
            powers: self.powers.pow(rhs),
            params: self.params.pow(rhs),
        }
    }
}

impl Mul<f32> for Results {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        Self {
            zones: self.zones * rhs,
            powers: self.powers * rhs,
            params: self.params * rhs,
        }
    }
}

impl Mul for Results {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            zones: self.zones * other.zones,
            powers: self.powers * other.powers,
            params: self.params * other.params,
        }
    }
}

impl Add for Results {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            zones: self.zones + other.zones,
            powers: self.powers + other.powers,
            params: self.params + other.params,
        }
    }
}

impl Sub for Results {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            zones: self.zones - other.zones,
            powers: self.powers - other.powers,
            params: self.params - other.params,
        }
    }
}

impl Div<f32> for Results {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        Self {
            zones: self.zones / rhs,
            powers: self.powers / rhs,
            params: self.params / rhs,
        }
    }
}

impl Div for Results {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Self {
            zones: self.zones / other.zones,
            powers: self.powers.div_elem(&other.powers),
            params: self.params.div_elem(&other.params),
        }
    }
}

impl Convergeable for Results {
    fn zero_like(&self) -> Self {
        Results::new_with_zones(self.zones.zero_like())
    }

    fn weighted_add(&self, other: &Self, w1: f32, w2: f32) -> Self {
        Self {
            zones: self.zones.weighted_add(&other.zones, w1, w2),
            powers: self.powers.weighted_add(&other.powers, w1, w2),
            params: self.params.weighted_add(&other.params, w1, w2),
        }
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
        Results::to_weighted(self)
    }

    fn weights(&self) -> Self {
        Results::weights(self)
    }
}

impl Results {
    /// Returns an owned vector of solid angle bins
    pub fn bins(&self) -> Vec<SolidAngleBin> {
        self.zones
            .full_zone()
            .map(|z| z.bins.clone())
            .unwrap_or_default()
    }

    /// Creates a new `Result` with empty mueller and amplitude matrix
    pub fn new_empty(bins: &[SolidAngleBin]) -> Self {
        let field_2d: Vec<ScattResult2D> =
            bins.iter().map(|&bin| ScattResult2D::new(bin)).collect();
        let full_zone = Zone::new(
            ZoneType::Full,
            bins.to_vec(),
            field_2d,
            None, // field_1d populated later by mueller_to_1d
        );
        Self {
            zones: Zones::new(vec![full_zone]),
            powers: Powers::new(),
            params: Params::new(),
        }
    }

    /// Create a new empty Results with initialized zones.
    pub fn new_with_zones(zones: Zones) -> Self {
        Self {
            zones,
            powers: Powers::new(),
            params: Params::new(),
        }
    }

    pub fn mueller_to_1d(&mut self) {
        for zone in self.zones.iter_mut() {
            // Skip Custom schemes (no regular grid for integration)
            match &zone.scheme {
                Scheme::Custom { .. } => continue,
                Scheme::Simple { .. } | Scheme::Interval { .. } => {}
            }

            // Group by theta using chunk_by (leveraging sorted property)
            let theta_groups: Vec<Vec<&ScattResult2D>> = zone
                .field_2d
                .iter()
                .chunk_by(|result| result.bin.theta)
                .into_iter()
                .map(|(_, group)| group.collect())
                .collect();

            // Rectangular integration over phi for each theta
            let field_1d: Vec<ScattResult1D> = theta_groups
                .into_iter()
                .map(|group| Self::integrate_over_phi(group))
                .collect();

            zone.field_1d = Some(field_1d);
        }
    }

    /// Integrates Mueller matrices over phi using rectangular rule
    /// Weighted by phi bin width in radians
    fn integrate_over_phi(phi_group: Vec<&ScattResult2D>) -> ScattResult1D {
        // All results in group have same theta bin
        let theta_bin = phi_group[0].bin.theta;
        let mut result = ScattResult1D::new(theta_bin);

        for phi_result in phi_group {
            // Convert phi width to radians to match theta integration units
            let phi_width_rad = phi_result.bin.phi.width().to_radians();

            // Integrate Mueller (weighted by phi bin width in radians)
            result.mueller_total += phi_result.mueller_total * phi_width_rad;
            result.mueller_beam += phi_result.mueller_beam * phi_width_rad;
            result.mueller_ext += phi_result.mueller_ext * phi_width_rad;
        }

        // Return integrated values without normalization to preserve 2π factor
        result
    }

    /// Computes and sets the parameters for each zone, then aggregates to global params.
    pub fn compute_params(&mut self, wavelength: f32) {
        let absorbed = self.powers.absorbed;

        // Compute params for each zone
        // Order matters: Forward must come before Backward for lidar ratio
        for zone in self.zones.iter_mut() {
            zone.compute_params(wavelength, absorbed, &self.params);

            // Aggregate zone params to global params for backwards compatibility
            self.params.merge(&zone.params);
        }
    }

    /// Returns a weighted version of Results for convergence tracking.
    /// - asymmetry becomes asymmetry * scat_cross
    /// - albedo becomes albedo * ext_cross
    /// - powers and other params stay the same (weight = 1)
    pub fn to_weighted(&self) -> Self {
        let mut result = self.clone();
        result.params = self.params.to_weighted();
        result
    }

    /// Returns a Results struct containing the weights for each field.
    /// - asymmetry slot contains scat_cross
    /// - albedo slot contains ext_cross
    /// - powers fields are all 1.0
    /// - scat_cross and ext_cross slots contain 1.0
    pub fn weights(&self) -> Self {
        Self {
            zones: self.zones.ones_like(),
            powers: Powers::ones(),
            params: self.params.weights(),
        }
    }

    pub fn print(&self) {
        println!("Powers: {:?}", self.powers);

        // Print parameters for each component
        for component in [GOComponent::Total, GOComponent::Beam, GOComponent::ExtDiff] {
            let comp_str = match component {
                GOComponent::Total => "Total",
                GOComponent::Beam => "Beam",
                GOComponent::ExtDiff => "ExtDiff",
            };

            if let Some(val) = self.params.asymmetry(&component) {
                println!("{} Asymmetry: {}", comp_str, val);
            }
            if let Some(val) = self.params.scatt_cross(&component) {
                println!("{} Scat Cross: {}", comp_str, val);
            }
            if let Some(val) = self.params.ext_cross(&component) {
                println!("{} Ext Cross: {}", comp_str, val);
            }
            if let Some(val) = self.params.albedo(&component) {
                println!("{} Albedo: {}", comp_str, val);
            }
        }
        if let Some(val) = self.params.scatt_cross(&GOComponent::Beam) {
            println!(
                "Beam Scat Cross / Output power: {}",
                val / self.powers.output
            );
        }
    }
}
