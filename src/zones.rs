//! Zone-based binning system for scattering calculations.
//!
//! Zones replace the single global binning scheme with a flexible list of
//! angular regions, each with its own binning configuration, results, and
//! computed parameters.

use log::info;
use pyo3::prelude::*;
use rand_distr::num_traits::Pow;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

use crate::bins::{Scheme, SolidAngleBin};
use crate::convergence::Convergeable;
use crate::params::Params;
use crate::result::{ScattResult1D, ScattResult2D};

/// The type of zone, which determines what parameters can be computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[pyclass(eq)]
pub enum ZoneType {
    /// Full 0-180 degree theta coverage. Computes: asymmetry, scattering cross-section.
    Full,
    /// Forward scattering zone. Computes: extinction cross-section (optical theorem).
    Forward,
    /// Backscatter zone. Computes: lidar ratio, backscatter cross-section.
    Backward,
    /// Custom angular range. Parameters depend on coverage.
    Custom,
}

impl ZoneType {
    /// Infer zone type from theta range.
    /// If theta spans 0-180 (within tolerance), it's Full; otherwise Custom.
    pub fn infer_from_scheme(scheme: &Scheme) -> Self {
        let (theta_min, theta_max) = scheme.theta_range();
        const TOL: f32 = 0.01;

        if (theta_min.abs() < TOL) && ((theta_max - 180.0).abs() < TOL) {
            ZoneType::Full
        } else {
            ZoneType::Custom
        }
    }
}

/// Configuration for a zone, as specified in TOML or via CLI.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ZoneConfig {
    /// Optional user-provided label for the zone.
    #[serde(default)]
    pub label: Option<String>,
    /// The binning scheme for this zone.
    pub scheme: Scheme,
}

impl ZoneConfig {
    /// Create a new zone config with the given scheme.
    pub fn new(scheme: Scheme) -> Self {
        Self {
            label: None,
            scheme,
        }
    }

    /// Create a new zone config with a label.
    pub fn with_label(label: impl Into<String>, scheme: Scheme) -> Self {
        Self {
            label: Some(label.into()),
            scheme,
        }
    }
}

/// A zone represents a region of the scattering sphere with its own binning,
/// results, and computed parameters.
#[derive(Debug, Clone)]
pub struct Zone {
    /// Optional user-provided label.
    pub label: Option<String>,
    /// The type of zone (Full, Forward, Backward, Custom).
    pub zone_type: ZoneType,
    /// The binning scheme for this zone.
    pub scheme: Scheme,
    /// Generated bins for this zone.
    pub bins: Vec<SolidAngleBin>,
    /// 2D scattering results (one per bin).
    pub field_2d: Vec<ScattResult2D>,
    /// 1D scattering results (integrated over phi), if applicable.
    pub field_1d: Option<Vec<ScattResult1D>>,
    /// Zone-specific computed parameters.
    pub params: Params,
}

impl Zone {
    /// Create a new zone with explicit fields.
    pub fn new(
        zone_type: ZoneType,
        bins: Vec<SolidAngleBin>,
        field_2d: Vec<ScattResult2D>,
        field_1d: Option<Vec<ScattResult1D>>,
    ) -> Self {
        Self {
            label: None,
            zone_type,
            scheme: Scheme::Custom {
                bins: vec![],
                file: None,
            },
            bins,
            field_2d,
            field_1d,
            params: Params::new(),
        }
    }

    /// Create a new zone from a config, generating bins and initializing empty results.
    pub fn from_config(config: &ZoneConfig) -> Self {
        let zone_type = ZoneType::infer_from_scheme(&config.scheme);
        let bins = config.scheme.generate();

        let label_str = config.label.as_deref().unwrap_or("<unnamed>");
        info!(
            "Processing zone '{}': {:?} ({} bins)",
            label_str,
            zone_type,
            bins.len()
        );

        let field_2d = bins.iter().map(|&bin| ScattResult2D::new(bin)).collect();

        Self {
            label: config.label.clone(),
            zone_type,
            scheme: config.scheme.clone(),
            bins,
            field_2d,
            field_1d: None,
            params: Params::new(),
        }
    }

    /// Create a forward scattering zone (single bin at theta≈0).
    /// Uses theta=0.01 to match legacy behavior and avoid singularity at exact zero.
    pub fn forward() -> Self {
        let scheme = Scheme::Custom {
            bins: vec![[[0.01, 0.01], [0.0, 0.0]]],
            file: None,
        };
        let bins = scheme.generate();
        let field_2d = bins.iter().map(|&bin| ScattResult2D::new(bin)).collect();

        info!("Processing zone 'forward': Forward (1 bin)");

        Self {
            label: Some("forward".to_string()),
            zone_type: ZoneType::Forward,
            scheme,
            bins,
            field_2d,
            field_1d: None,
            params: Params::new(),
        }
    }

    /// Create a backscatter zone (single bin at theta=180).
    pub fn backward() -> Self {
        let scheme = Scheme::Custom {
            bins: vec![[[180.0, 180.0], [0.0, 0.0]]],
            file: None,
        };
        let bins = scheme.generate();
        let field_2d = bins.iter().map(|&bin| ScattResult2D::new(bin)).collect();

        info!("Processing zone 'backward': Backward (1 bin)");

        Self {
            label: Some("backward".to_string()),
            zone_type: ZoneType::Backward,
            scheme,
            bins,
            field_2d,
            field_1d: None,
            params: Params::new(),
        }
    }

    /// Get a display name for this zone.
    pub fn display_name(&self) -> String {
        self.label
            .clone()
            .unwrap_or_else(|| format!("{:?}", self.zone_type).to_lowercase())
    }

    /// Reset the zone's results to empty state.
    pub fn reset(&mut self) {
        self.field_2d = self
            .bins
            .iter()
            .map(|&bin| ScattResult2D::new(bin))
            .collect();
        self.field_1d = None;
        self.params = Params::new();
    }

    /// Returns a Zone with all values set to 1.0 (for weights).
    pub fn ones_like(&self) -> Self {
        Self {
            label: self.label.clone(),
            zone_type: self.zone_type,
            scheme: self.scheme.clone(),
            bins: self.bins.clone(),
            field_2d: self.field_2d.iter().map(|f| f.ones_like()).collect(),
            field_1d: self
                .field_1d
                .as_ref()
                .map(|f| f.iter().map(|x| x.ones_like()).collect()),
            params: self.params.weights(),
        }
    }
}

/// A collection of zones for a simulation.
#[derive(Debug, Clone)]
pub struct Zones {
    zones: Vec<Zone>,
}

impl Zones {
    /// Create a new Zones collection from zone configs.
    /// Automatically adds Forward and Backward zones.
    pub fn from_configs(configs: &[ZoneConfig]) -> Self {
        let mut zones: Vec<Zone> = configs.iter().map(Zone::from_config).collect();

        // Add forward and backward zones
        zones.push(Zone::forward());
        zones.push(Zone::backward());

        Self { zones }
    }

    /// Create a Zones collection from a vector of zones.
    pub fn new(zones: Vec<Zone>) -> Self {
        Self { zones }
    }

    /// Create an empty Zones collection.
    pub fn empty() -> Self {
        Self { zones: Vec::new() }
    }

    /// Get all zones.
    pub fn all(&self) -> &[Zone] {
        &self.zones
    }

    /// Get all zones mutably.
    pub fn all_mut(&mut self) -> &mut [Zone] {
        &mut self.zones
    }

    /// Get a zone by label.
    pub fn get(&self, label: &str) -> Option<&Zone> {
        self.zones
            .iter()
            .find(|z| z.label.as_deref() == Some(label))
    }

    /// Get a zone by label mutably.
    pub fn get_mut(&mut self, label: &str) -> Option<&mut Zone> {
        self.zones
            .iter_mut()
            .find(|z| z.label.as_deref() == Some(label))
    }

    /// Get all zones of a specific type.
    pub fn by_type(&self, zone_type: ZoneType) -> Vec<&Zone> {
        self.zones
            .iter()
            .filter(|z| z.zone_type == zone_type)
            .collect()
    }

    /// Get the first Full zone, if any.
    pub fn full_zone(&self) -> Option<&Zone> {
        self.zones.iter().find(|z| z.zone_type == ZoneType::Full)
    }

    /// Get the first Full zone mutably, if any.
    pub fn full_zone_mut(&mut self) -> Option<&mut Zone> {
        self.zones
            .iter_mut()
            .find(|z| z.zone_type == ZoneType::Full)
    }

    /// Get the forward zone.
    pub fn forward_zone(&self) -> Option<&Zone> {
        self.zones.iter().find(|z| z.zone_type == ZoneType::Forward)
    }

    /// Get the backward zone.
    pub fn backward_zone(&self) -> Option<&Zone> {
        self.zones
            .iter()
            .find(|z| z.zone_type == ZoneType::Backward)
    }

    /// Get the number of zones.
    pub fn len(&self) -> usize {
        self.zones.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.zones.is_empty()
    }

    /// Iterate over zones.
    pub fn iter(&self) -> impl Iterator<Item = &Zone> {
        self.zones.iter()
    }

    /// Iterate over zones mutably.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Zone> {
        self.zones.iter_mut()
    }

    /// Reset all zones to empty state.
    pub fn reset(&mut self) {
        for zone in &mut self.zones {
            zone.reset();
        }
    }

    /// Returns a Zones collection with all values set to 1.0 (for weights).
    pub fn ones_like(&self) -> Self {
        Self {
            zones: self.zones.iter().map(|z| z.ones_like()).collect(),
        }
    }
}

impl IntoIterator for Zones {
    type Item = Zone;
    type IntoIter = std::vec::IntoIter<Zone>;

    fn into_iter(self) -> Self::IntoIter {
        self.zones.into_iter()
    }
}

impl<'a> IntoIterator for &'a Zones {
    type Item = &'a Zone;
    type IntoIter = std::slice::Iter<'a, Zone>;

    fn into_iter(self) -> Self::IntoIter {
        self.zones.iter()
    }
}

impl<'a> IntoIterator for &'a mut Zones {
    type Item = &'a mut Zone;
    type IntoIter = std::slice::IterMut<'a, Zone>;

    fn into_iter(self) -> Self::IntoIter {
        self.zones.iter_mut()
    }
}

// ============================================================================
// Arithmetic operations for Zone
// ============================================================================

impl Add for Zone {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a + b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(f1), Some(f2)) => Some(f1.into_iter().zip(f2).map(|(a, b)| a + b).collect()),
            (Some(f1), None) => Some(f1),
            (None, Some(f2)) => Some(f2),
            (None, None) => None,
        };
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params + other.params,
        }
    }
}

impl Sub for Zone {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a - b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(f1), Some(f2)) => Some(f1.into_iter().zip(f2).map(|(a, b)| a - b).collect()),
            (Some(f1), None) => Some(f1),
            (None, Some(f2)) => Some(f2),
            (None, None) => None,
        };
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params - other.params,
        }
    }
}

impl Mul for Zone {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a * b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(f1), Some(f2)) => Some(f1.into_iter().zip(f2).map(|(a, b)| a * b).collect()),
            (Some(f1), None) => Some(f1),
            (None, Some(f2)) => Some(f2),
            (None, None) => None,
        };
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params * other.params,
        }
    }
}

impl Mul<f32> for Zone {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        let field_2d = self.field_2d.into_iter().map(|f| f * rhs).collect();
        let field_1d = self
            .field_1d
            .map(|f| f.into_iter().map(|x| x * rhs).collect());
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params * rhs,
        }
    }
}

impl Div for Zone {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a / b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(f1), Some(f2)) => Some(f1.into_iter().zip(f2).map(|(a, b)| a / b).collect()),
            (Some(f1), None) => Some(f1),
            (None, Some(_)) => None,
            (None, None) => None,
        };
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params.div_elem(&other.params),
        }
    }
}

impl Div<f32> for Zone {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        let field_2d = self.field_2d.into_iter().map(|f| f / rhs).collect();
        let field_1d = self
            .field_1d
            .map(|f| f.into_iter().map(|x| x / rhs).collect());
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params / rhs,
        }
    }
}

impl Pow<f32> for Zone {
    type Output = Self;

    fn pow(self, rhs: f32) -> Self {
        let field_2d = self.field_2d.into_iter().map(|f| f.pow(rhs)).collect();
        let field_1d = self
            .field_1d
            .map(|f| f.into_iter().map(|x| x.pow(rhs)).collect());
        Self {
            label: self.label,
            zone_type: self.zone_type,
            scheme: self.scheme,
            bins: self.bins,
            field_2d,
            field_1d,
            params: self.params.pow(rhs),
        }
    }
}

// ============================================================================
// Arithmetic operations for Zones
// ============================================================================

impl Add for Zones {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let zones = self
            .zones
            .into_iter()
            .zip(other.zones)
            .map(|(a, b)| a + b)
            .collect();
        Self { zones }
    }
}

impl Sub for Zones {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let zones = self
            .zones
            .into_iter()
            .zip(other.zones)
            .map(|(a, b)| a - b)
            .collect();
        Self { zones }
    }
}

impl Mul for Zones {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let zones = self
            .zones
            .into_iter()
            .zip(other.zones)
            .map(|(a, b)| a * b)
            .collect();
        Self { zones }
    }
}

impl Mul<f32> for Zones {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        let zones = self.zones.into_iter().map(|z| z * rhs).collect();
        Self { zones }
    }
}

impl Div for Zones {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let zones = self
            .zones
            .into_iter()
            .zip(other.zones)
            .map(|(a, b)| a / b)
            .collect();
        Self { zones }
    }
}

impl Div<f32> for Zones {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        let zones = self.zones.into_iter().map(|z| z / rhs).collect();
        Self { zones }
    }
}

impl Pow<f32> for Zones {
    type Output = Self;

    fn pow(self, rhs: f32) -> Self {
        let zones = self.zones.into_iter().map(|z| z.pow(rhs)).collect();
        Self { zones }
    }
}

// ============================================================================
// Convergeable implementations for Zone and Zones
// ============================================================================

impl Convergeable for Zone {
    fn zero_like(&self) -> Self {
        Self {
            label: self.label.clone(),
            zone_type: self.zone_type,
            scheme: self.scheme.clone(),
            bins: self.bins.clone(),
            field_2d: self
                .bins
                .iter()
                .map(|&bin| ScattResult2D::new(bin))
                .collect(),
            field_1d: None,
            params: self.params.zero_like(),
        }
    }

    fn weighted_add(&self, other: &Self, w1: f32, w2: f32) -> Self {
        let field_2d = self
            .field_2d
            .iter()
            .zip(other.field_2d.iter())
            .map(|(a, b)| a.weighted_add(b, w1, w2))
            .collect();
        let field_1d = match (&self.field_1d, &other.field_1d) {
            (Some(f1), Some(f2)) => Some(
                f1.iter()
                    .zip(f2.iter())
                    .map(|(a, b)| a.weighted_add(b, w1, w2))
                    .collect(),
            ),
            (Some(f1), None) => Some(f1.clone()),
            (None, Some(f2)) => Some(f2.clone()),
            (None, None) => None,
        };
        Self {
            label: self.label.clone(),
            zone_type: self.zone_type,
            scheme: self.scheme.clone(),
            bins: self.bins.clone(),
            field_2d,
            field_1d,
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
        Pow::pow(self.clone(), 0.5)
    }

    fn to_weighted(&self) -> Self {
        // For zones, each field's to_weighted is delegated
        Self {
            label: self.label.clone(),
            zone_type: self.zone_type,
            scheme: self.scheme.clone(),
            bins: self.bins.clone(),
            field_2d: self.field_2d.iter().map(|f| f.to_weighted()).collect(),
            field_1d: self
                .field_1d
                .as_ref()
                .map(|fs| fs.iter().map(|f| f.to_weighted()).collect()),
            params: self.params.to_weighted(),
        }
    }

    fn weights(&self) -> Self {
        // For zones, each field's weights is delegated
        Self {
            label: self.label.clone(),
            zone_type: self.zone_type,
            scheme: self.scheme.clone(),
            bins: self.bins.clone(),
            field_2d: self.field_2d.iter().map(|f| f.weights()).collect(),
            field_1d: self
                .field_1d
                .as_ref()
                .map(|fs| fs.iter().map(|f| f.weights()).collect()),
            params: self.params.weights(),
        }
    }
}

impl Convergeable for Zones {
    fn zero_like(&self) -> Self {
        Self {
            zones: self.zones.iter().map(|z| z.zero_like()).collect(),
        }
    }

    fn weighted_add(&self, other: &Self, w1: f32, w2: f32) -> Self {
        let zones = self
            .zones
            .iter()
            .zip(other.zones.iter())
            .map(|(a, b)| a.weighted_add(b, w1, w2))
            .collect();
        Self { zones }
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
        Pow::pow(self.clone(), 0.5)
    }

    fn to_weighted(&self) -> Self {
        Self {
            zones: self.zones.iter().map(|z| z.to_weighted()).collect(),
        }
    }

    fn weights(&self) -> Self {
        Self {
            zones: self.zones.iter().map(|z| z.weights()).collect(),
        }
    }
}
