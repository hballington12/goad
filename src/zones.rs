//! Zone-based binning system for scattering calculations.
//!
//! Zones replace the single global binning scheme with a flexible list of
//! angular regions, each with its own binning configuration, results, and
//! computed parameters.

use log::{info, warn};
use numpy::IntoPyArray;
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rand_distr::num_traits::Pow;

use crate::result::MuellerMatrix;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;
use std::ops::{Add, Div, Mul, Sub};

use crate::bins::{BinningScheme, Scheme, SolidAngleBin};
use crate::convergence::Convergeable;
use crate::params::{Param, Params};
use crate::result::{
    integrate_theta_weighted_component, GOComponent, ScattResult1D, ScattResult2D,
};
use crate::settings::constants::ZONE_THETA_OFFSET;

/// The type of zone, which determines what parameters can be computed.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass_enum)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[pyclass(module = "goad._goad", eq)]
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
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
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

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl ZoneConfig {
    /// Create a new zone configuration.
    ///
    /// Args:
    ///     binning: The binning scheme for this zone
    ///     label: Optional label for the zone
    #[new]
    #[pyo3(signature = (binning, label=None))]
    fn py_new(binning: BinningScheme, label: Option<String>) -> Self {
        Self {
            label,
            scheme: binning.scheme,
        }
    }

    /// Get the zone label
    #[getter]
    fn get_label(&self) -> Option<String> {
        self.label.clone()
    }

    /// Set the zone label
    #[setter]
    fn set_label(&mut self, label: Option<String>) {
        self.label = label;
    }

    /// Get the binning scheme
    #[getter]
    fn get_binning(&self) -> BinningScheme {
        BinningScheme {
            scheme: self.scheme.clone(),
        }
    }

    /// Set the binning scheme
    #[setter]
    fn set_binning(&mut self, binning: BinningScheme) {
        self.scheme = binning.scheme;
    }

    fn __repr__(&self) -> String {
        let label_str = self
            .label
            .as_ref()
            .map(|l| format!("'{}'", l))
            .unwrap_or_else(|| "None".to_string());
        format!("ZoneConfig(label={}, binning=...)", label_str)
    }
}

/// A zone represents a region of the scattering sphere with its own binning,
/// results, and computed parameters.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[derive(Debug, Clone)]
#[pyclass(module = "goad._goad")]
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
    /// Uses theta=ZONE_THETA_OFFSET to avoid singularity at exact zero.
    pub fn forward() -> Self {
        let scheme = Scheme::Custom {
            bins: vec![[[ZONE_THETA_OFFSET, ZONE_THETA_OFFSET], [0.0, 0.0]]],
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

    /// Create a backscatter zone (single bin at theta≈180).
    /// Uses theta=180-ZONE_THETA_OFFSET to avoid singularity at exact 180.
    pub fn backward() -> Self {
        let theta = 180.0 - ZONE_THETA_OFFSET;
        let scheme = Scheme::Custom {
            bins: vec![[[theta, theta], [0.0, 0.0]]],
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

    /// Compute zone-specific parameters based on zone type.
    ///
    /// - Full: scatt_cross, ext_cross, asymmetry, albedo
    /// - Forward: ext_cross_optical_theorem
    /// - Backward: backscatter_cross, depolarization_ratio, lidar_ratio
    /// - Custom: nothing (for now)
    pub fn compute_params(&mut self, wavelength: f32, absorbed: f32, global_params: &Params) {
        let k = 2.0 * PI / wavelength;

        match self.zone_type {
            ZoneType::Full => {
                self.compute_full_params(k, absorbed);
            }
            ZoneType::Forward => {
                self.compute_forward_params(k);
            }
            ZoneType::Backward => {
                self.compute_backward_params(k, global_params);
            }
            ZoneType::Custom => {
                // No params computed for custom zones (yet)
            }
        }
    }

    fn compute_full_params(&mut self, k: f32, absorbed: f32) {
        let Some(field_1d) = &self.field_1d else {
            return;
        };

        for component in [GOComponent::Total, GOComponent::Beam, GOComponent::ExtDiff] {
            let scatt = integrate_theta_weighted_component(field_1d, component, |theta, s11| {
                theta.sin() * s11 / k.powi(2)
            });
            let asymmetry_scatt =
                integrate_theta_weighted_component(field_1d, component, |theta, s11| {
                    theta.sin() * theta.cos() * s11 / k.powi(2)
                });

            self.params.set_param(Param::ScatCross, component, scatt);
            self.params
                .set_param(Param::Asymmetry, component, asymmetry_scatt / scatt);

            // ext_cross and albedo only for Total component
            if component == GOComponent::Total {
                let ext = scatt + absorbed;
                self.params.set_param(Param::ExtCross, component, ext);
                self.params.set_param(Param::Albedo, component, scatt / ext);
            }
        }
    }

    fn compute_forward_params(&mut self, k: f32) {
        let Some(field_fs) = self.field_2d.first() else {
            return;
        };

        let s2 = field_fs.ampl_total[(0, 0)];
        let ext_cross = s2.im * 4.0 * PI / k.powi(2);
        self.params
            .set_param(Param::ExtCrossOpticalTheorem, GOComponent::Total, ext_cross);
    }

    fn compute_backward_params(&mut self, k: f32, global_params: &Params) {
        let Some(field_bs) = self.field_2d.first() else {
            return;
        };

        for (component, mueller) in [
            (GOComponent::Total, field_bs.mueller_total),
            (GOComponent::Beam, field_bs.mueller_beam),
            (GOComponent::ExtDiff, field_bs.mueller_ext),
        ] {
            let s11 = mueller[(0, 0)];
            let s22 = mueller[(1, 1)];
            let bs_cross = s11 * 4.0 * PI / k.powi(2);

            self.params
                .set_param(Param::BackscatterCross, component, bs_cross);

            // Lidar ratio = ext_cross_optical_theorem / backscatter_cross
            if let Some(ext_cross_ot) = global_params.ext_cross_optical_theorem(&component) {
                if bs_cross > 1e-10 {
                    self.params
                        .set_param(Param::LidarRatio, component, ext_cross_ot / bs_cross);
                }
            } else if component == GOComponent::Total {
                warn!(
                    "Cannot compute lidar ratio: ext_cross_optical_theorem not available. \
                     Ensure Forward zone is processed before Backward zone."
                );
            }

            // Depolarization ratio = (S11 - S22) / (S11 + S22)
            let s11_plus_s22 = s11 + s22;
            if s11_plus_s22.abs() > 1e-10 {
                let depol = (s11 - s22) / s11_plus_s22;
                self.params
                    .set_param(Param::DepolarizationRatio, component, depol);
                self.params
                    .set_param(Param::BackscatterS11S22, component, s11_plus_s22);
            }
        }
    }
}

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl Zone {
    /// Get the zone label
    #[getter]
    pub fn get_label(&self) -> Option<String> {
        self.label.clone()
    }

    /// Get the zone type
    #[getter]
    pub fn get_zone_type(&self) -> ZoneType {
        self.zone_type
    }

    /// Get the display name for this zone
    #[getter]
    pub fn get_name(&self) -> String {
        self.display_name()
    }

    /// Get the number of bins in this zone
    #[getter]
    pub fn get_num_bins(&self) -> usize {
        self.bins.len()
    }

    /// Get the bins as a numpy array of shape (n_bins, 2) with columns [theta, phi]
    #[getter]
    pub fn get_bins<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        let bins: Vec<f32> = self
            .bins
            .iter()
            .flat_map(|bin| vec![bin.theta.center, bin.phi.center])
            .collect();
        ndarray::Array2::from_shape_vec((bins.len() / 2, 2), bins)
            .unwrap()
            .into_pyarray(py)
    }

    /// Get the Mueller matrix as a numpy array of shape (n_bins, 16)
    #[getter]
    pub fn get_mueller<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        let muellers: Vec<f32> = self
            .field_2d
            .iter()
            .flat_map(|r| r.mueller_total.to_vec())
            .collect();
        ndarray::Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
            .unwrap()
            .into_pyarray(py)
    }

    /// Get the 1D Mueller matrix as a numpy array (if available)
    #[getter]
    pub fn get_mueller_1d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, numpy::PyArray2<f32>>> {
        self.field_1d.as_ref().map(|field_1d| {
            let muellers: Vec<f32> = field_1d
                .iter()
                .flat_map(|r| r.mueller_total.to_vec())
                .collect();
            ndarray::Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                .unwrap()
                .into_pyarray(py)
        })
    }

    /// Get the 1D theta bins as a numpy array (if available)
    #[getter]
    pub fn get_bins_1d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, numpy::PyArray1<f32>>> {
        self.field_1d.as_ref().map(|field_1d| {
            let bins: Vec<f32> = field_1d.iter().map(|r| r.bin.center).collect();
            ndarray::Array1::from_vec(bins).into_pyarray(py)
        })
    }

    /// Get zone-specific parameters as a dict
    #[getter]
    pub fn get_params(&self) -> PyResult<Py<PyAny>> {
        use pyo3::types::PyDict;
        Python::attach(|py| {
            let dict = PyDict::new(py);
            // Add all available params
            if let Some(v) = self.params.asymmetry(&crate::result::GOComponent::Total) {
                dict.set_item("asymmetry", v)?;
            }
            if let Some(v) = self.params.scatt_cross(&crate::result::GOComponent::Total) {
                dict.set_item("scatt_cross", v)?;
            }
            if let Some(v) = self.params.ext_cross(&crate::result::GOComponent::Total) {
                dict.set_item("ext_cross", v)?;
            }
            if let Some(v) = self.params.albedo(&crate::result::GOComponent::Total) {
                dict.set_item("albedo", v)?;
            }
            if let Some(v) = self
                .params
                .ext_cross_optical_theorem(&crate::result::GOComponent::Total)
            {
                dict.set_item("ext_cross_optical_theorem", v)?;
            }
            if let Some(v) = self
                .params
                .backscatter_cross(&crate::result::GOComponent::Total)
            {
                dict.set_item("backscatter_cross", v)?;
            }
            if let Some(v) = self.params.lidar_ratio(&crate::result::GOComponent::Total) {
                dict.set_item("lidar_ratio", v)?;
            }
            if let Some(v) = self
                .params
                .depolarization_ratio(&crate::result::GOComponent::Total)
            {
                dict.set_item("depolarization_ratio", v)?;
            }
            Ok(dict.into())
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "Zone(name='{}', type={:?}, bins={})",
            self.display_name(),
            self.zone_type,
            self.bins.len()
        )
    }
}

/// A collection of zones for a simulation.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[derive(Debug, Clone)]
#[pyclass(module = "goad._goad")]
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

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl Zones {
    /// Get the number of zones
    fn __len__(&self) -> usize {
        self.zones.len()
    }

    /// Get a zone by index
    fn __getitem__(&self, index: isize) -> PyResult<Zone> {
        let len = self.zones.len() as isize;
        let idx = if index < 0 { len + index } else { index };
        if idx < 0 || idx >= len {
            return Err(pyo3::exceptions::PyIndexError::new_err(
                "zone index out of range",
            ));
        }
        Ok(self.zones[idx as usize].clone())
    }

    /// Get a zone by label
    #[pyo3(name = "get")]
    pub fn py_get(&self, label: &str) -> Option<Zone> {
        self.get(label).cloned()
    }

    /// Get a zone by type (returns first matching zone)
    pub fn get_by_type(&self, zone_type: ZoneType) -> Option<Zone> {
        self.zones
            .iter()
            .find(|z| z.zone_type == zone_type)
            .cloned()
    }

    /// Get the full zone (convenience method)
    #[getter]
    pub fn full(&self) -> Option<Zone> {
        self.full_zone().cloned()
    }

    /// Get the forward zone (convenience method)
    #[getter]
    pub fn forward(&self) -> Option<Zone> {
        self.zones
            .iter()
            .find(|z| z.zone_type == ZoneType::Forward)
            .cloned()
    }

    /// Get the backward zone (convenience method)
    #[getter]
    pub fn backward(&self) -> Option<Zone> {
        self.backward_zone().cloned()
    }

    /// Get all zones as a list
    #[getter]
    pub fn all_zones(&self) -> Vec<Zone> {
        self.zones.clone()
    }

    fn __repr__(&self) -> String {
        let zone_names: Vec<String> = self.zones.iter().map(|z| z.display_name()).collect();
        format!("Zones([{}])", zone_names.join(", "))
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<ZonesIterator>> {
        let iter = ZonesIterator {
            zones: slf.zones.clone(),
            index: 0,
        };
        Py::new(slf.py(), iter)
    }
}

/// Iterator for Zones in Python
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
pub struct ZonesIterator {
    zones: Vec<Zone>,
    index: usize,
}

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl ZonesIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Zone> {
        if slf.index < slf.zones.len() {
            let zone = slf.zones[slf.index].clone();
            slf.index += 1;
            Some(zone)
        } else {
            None
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
