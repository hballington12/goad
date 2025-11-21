use std::f32::consts::PI;
use std::fmt::Debug;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use crate::bins::AngleBin;
use crate::bins::Scheme;
use crate::bins::SolidAngleBin;
use crate::params::Param;
use crate::params::Params;
use crate::powers::Powers;
use itertools::Itertools;
use nalgebra::Matrix4;
use nalgebra::{Complex, Matrix2};
use ndarray::Array2;
use numpy::PyArrayMethods;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rand_distr::num_traits::Pow;
use serde::Serialize;

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

#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Serialize)]
pub enum GOComponent {
    Total,
    Beam,
    ExtDiff,
}

impl GOComponent {
    /// Returns the file extension for the given GOComponent.
    pub fn file_extension(&self) -> &'static str {
        match self {
            GOComponent::Total => "",
            GOComponent::Beam => "beam",
            GOComponent::ExtDiff => "ext",
        }
    }
}

type Ampl = Matrix2<Complex<f32>>;
pub type Mueller = Matrix4<f32>;

/// Trait for approximate equality comparison with tolerance
pub trait ApproxEq {
    /// Check if two values are approximately equal within the given tolerance
    fn approx_eq(&self, other: &Self, tolerance: f32) -> bool;
}

impl ApproxEq for Ampl {
    /// Check if two amplitude matrices are approximately equal within tolerance
    ///
    /// Compares both real and imaginary parts of each complex element.
    /// Returns true if all corresponding elements differ by less than the tolerance.
    fn approx_eq(&self, other: &Self, tolerance: f32) -> bool {
        for i in 0..2 {
            for j in 0..2 {
                let a = self[(i, j)];
                let b = other[(i, j)];

                // Check both real and imaginary parts
                if (a.re - b.re).abs() > tolerance || (a.im - b.im).abs() > tolerance {
                    return false;
                }
            }
        }
        true
    }
}

pub trait MuellerMatrix {
    fn s11(&self) -> f32;
    fn s12(&self) -> f32;
    fn s13(&self) -> f32;
    fn s14(&self) -> f32;
    fn s21(&self) -> f32;
    fn s22(&self) -> f32;
    fn s23(&self) -> f32;
    fn s24(&self) -> f32;
    fn s31(&self) -> f32;
    fn s32(&self) -> f32;
    fn s33(&self) -> f32;
    fn s34(&self) -> f32;
    fn s41(&self) -> f32;
    fn s42(&self) -> f32;
    fn s43(&self) -> f32;
    fn s44(&self) -> f32;
    fn to_vec(&self) -> Vec<f32>;
}

impl MuellerMatrix for Mueller {
    fn s11(&self) -> f32 {
        self[(0, 0)]
    }
    fn s12(&self) -> f32 {
        self[(0, 1)]
    }
    fn s13(&self) -> f32 {
        self[(0, 2)]
    }
    fn s14(&self) -> f32 {
        self[(0, 3)]
    }
    fn s21(&self) -> f32 {
        self[(1, 0)]
    }
    fn s22(&self) -> f32 {
        self[(1, 1)]
    }
    fn s23(&self) -> f32 {
        self[(1, 2)]
    }
    fn s24(&self) -> f32 {
        self[(1, 3)]
    }
    fn s31(&self) -> f32 {
        self[(2, 0)]
    }
    fn s32(&self) -> f32 {
        self[(2, 1)]
    }
    fn s33(&self) -> f32 {
        self[(2, 2)]
    }
    fn s34(&self) -> f32 {
        self[(2, 3)]
    }
    fn s41(&self) -> f32 {
        self[(3, 0)]
    }
    fn s42(&self) -> f32 {
        self[(3, 1)]
    }
    fn s43(&self) -> f32 {
        self[(3, 2)]
    }
    fn s44(&self) -> f32 {
        self[(3, 3)]
    }
    /// Returns the Mueller matrix as a vector of its elements.
    fn to_vec(&self) -> Vec<f32> {
        vec![
            self.s11(),
            self.s12(),
            self.s13(),
            self.s14(),
            self.s21(),
            self.s22(),
            self.s23(),
            self.s24(),
            self.s31(),
            self.s32(),
            self.s33(),
            self.s34(),
            self.s41(),
            self.s42(),
            self.s43(),
            self.s44(),
        ]
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

/// Type alias for 2D scattering results (full solid angle)
pub type ScattResult2D = ScattResult<SolidAngleBin>;

/// Type alias for 1D scattering results (theta only)
pub type ScattResult1D = ScattResult<AngleBin>;
/// Complete results from a GOAD light scattering simulation.
///
/// Contains all computed scattering data including Mueller matrices,
/// amplitude matrices, power distributions, and derived parameters.
/// Supports both 2D angular distributions and 1D integrated results.
#[pyclass]
#[derive(Debug, Clone)]
pub struct Results {
    pub field_2d: Vec<ScattResult2D>,
    pub field_1d: Option<Vec<ScattResult1D>>,
    pub powers: Powers,
    pub params: Params,
}

impl AddAssign for Results {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            field_2d: self
                .field_2d
                .clone()
                .into_iter()
                .zip(other.field_2d)
                .map(|(a, b)| a + b)
                .collect(),
            field_1d: match (self.field_1d.clone(), other.field_1d) {
                (Some(field_1d), Some(other_field_1d)) => Some(
                    field_1d
                        .into_iter()
                        .zip(other_field_1d)
                        .map(|(a, b)| a + b)
                        .collect(),
                ),
                (Some(field_1d), None) => Some(field_1d),
                (None, Some(other_field_1d)) => Some(other_field_1d),
                (None, None) => None,
            },
            powers: self.powers.clone() + other.powers,
            params: self.params.clone() + other.params,
        };
    }
}

impl Pow<f32> for Results {
    type Output = Self;

    fn pow(self, rhs: f32) -> Self {
        let field_1d = match self.field_1d {
            Some(field_1d) => Some(field_1d.into_iter().map(|a| a.pow(rhs)).collect()),
            None => None,
        };
        Self {
            field_2d: self.field_2d.into_iter().map(|a| a.pow(rhs)).collect(),
            field_1d,
            powers: self.powers.pow(rhs),
            params: self.params.pow(rhs),
        }
    }
}

impl Mul<f32> for Results {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        let field_1d = match self.field_1d {
            Some(field_1d) => Some(field_1d.into_iter().map(|a| a * rhs).collect()),
            None => None,
        };
        Self {
            field_2d: self.field_2d.into_iter().map(|a| a * rhs).collect(),
            field_1d,
            powers: self.powers * rhs,
            params: self.params * rhs,
        }
    }
}

impl Mul for Results {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a * b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(field_1d), Some(other_field_1d)) => Some(
                field_1d
                    .into_iter()
                    .zip(other_field_1d)
                    .map(|(a, b)| a * b)
                    .collect(),
            ),
            (Some(field_1d), None) => Some(field_1d),
            (None, Some(other_field_1d)) => Some(other_field_1d),
            (None, None) => None,
        };
        let powers = self.powers * other.powers;
        let params = self.params * other.params;
        Self {
            field_2d,
            field_1d,
            powers,
            params,
        }
    }
}

impl Add for Results {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a + b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(field_1d), Some(other_field_1d)) => Some(
                field_1d
                    .into_iter()
                    .zip(other_field_1d)
                    .map(|(a, b)| a + b)
                    .collect(),
            ),
            (Some(field_1d), None) => Some(field_1d),
            (None, Some(other_field_1d)) => Some(other_field_1d),
            (None, None) => None,
        };
        let powers = self.powers + other.powers;
        let params = self.params + other.params;
        Self {
            field_2d,
            field_1d,
            powers,
            params,
        }
    }
}

impl Sub for Results {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let field_2d = self
            .field_2d
            .into_iter()
            .zip(other.field_2d)
            .map(|(a, b)| a - b)
            .collect();
        let field_1d = match (self.field_1d, other.field_1d) {
            (Some(field_1d), Some(other_field_1d)) => Some(
                field_1d
                    .into_iter()
                    .zip(other_field_1d)
                    .map(|(a, b)| a - b)
                    .collect(),
            ),
            (Some(field_1d), None) => Some(field_1d),
            (None, Some(other_field_1d)) => Some(other_field_1d),
            (None, None) => None,
        };
        let powers = self.powers - other.powers;
        let params = self.params - other.params;
        Self {
            field_2d,
            field_1d,
            powers,
            params,
        }
    }
}

impl Div<f32> for Results {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        let field_1d = match self.field_1d {
            Some(field_1d) => Some(field_1d.into_iter().map(|a| a / rhs).collect()),
            None => None,
        };
        Self {
            field_2d: self.field_2d.into_iter().map(|a| a / rhs).collect(),
            field_1d,
            powers: self.powers / rhs,
            params: self.params / rhs,
        }
    }
}

impl Results {
    /// Returns an owned vector of solid angle bins
    pub fn bins(&self) -> Vec<SolidAngleBin> {
        self.field_2d.iter().map(|a| a.bin.clone()).collect()
    }

    /// Writes some stuff to a file

    /// Creates a new `Result` with empty mueller and amplitude matrix
    pub fn new_empty(bins: &[SolidAngleBin]) -> Self {
        let field = bins.iter().map(|&bin| ScattResult2D::new(bin)).collect();
        Self {
            field_2d: field,
            powers: Powers::new(),
            field_1d: None,
            params: Params::new(),
        }
    }

    pub fn mueller_to_1d(&mut self, binning_scheme: &crate::bins::Scheme) {
        // Step 1: Check scheme compatibility
        match binning_scheme {
            Scheme::Custom { .. } => {
                return;
            }
            Scheme::Simple { .. } | Scheme::Interval { .. } => {}
        }

        // Step 2: Group by theta using chunk_by (leveraging sorted property)
        let theta_groups: Vec<Vec<&ScattResult2D>> = self
            .field_2d
            .iter()
            .chunk_by(|result| result.bin.theta_bin)
            .into_iter()
            .map(|(_, group)| group.collect())
            .collect();

        // Step 3: Rectangular integration over phi for each theta
        let field_1d: Vec<ScattResult1D> = theta_groups
            .into_iter()
            .map(|group| Self::integrate_over_phi(group))
            .collect();

        // Step 4: Update struct
        self.field_1d = Some(field_1d);
    }

    /// Integrates Mueller matrices over phi using rectangular rule
    /// Weighted by phi bin width in radians
    fn integrate_over_phi(phi_group: Vec<&ScattResult2D>) -> ScattResult1D {
        // All results in group have same theta bin
        let theta_bin = phi_group[0].bin.theta_bin;
        let mut result = ScattResult1D::new(theta_bin);

        for phi_result in phi_group {
            // Convert phi width to radians to match theta integration units
            let phi_width_rad = phi_result.bin.phi_bin.width().to_radians();

            // Integrate Mueller (weighted by phi bin width in radians)
            result.mueller_total += phi_result.mueller_total * phi_width_rad;
            result.mueller_beam += phi_result.mueller_beam * phi_width_rad;
            result.mueller_ext += phi_result.mueller_ext * phi_width_rad;
        }

        // Return integrated values without normalization to preserve 2π factor
        result
    }

    /// Computes and sets the parameters of the result
    pub fn compute_params(&mut self, wavelength: f32) {
        // Total field
        if let Some(field_1d) = &self.field_1d {
            let k = 2.0 * PI / wavelength;
            let scatt_total =
                integrate_theta_weighted_component(field_1d, GOComponent::Total, |theta, s11| {
                    theta.sin() * s11 / k.powi(2)
                });
            let asymmetry_scatt_total =
                integrate_theta_weighted_component(field_1d, GOComponent::Total, |theta, s11| {
                    theta.sin() * theta.cos() * s11 / k.powi(2)
                });

            self.params
                .set_param(Param::ScatCross, GOComponent::Total, scatt_total);
            self.params.set_param(
                Param::AsymmetryScatt,
                GOComponent::Total,
                asymmetry_scatt_total,
            );
            let ext_total = scatt_total + self.powers.absorbed;
            self.params
                .set_param(Param::ExtCross, GOComponent::Total, ext_total);

            // Beam field
            let scatt_beam =
                integrate_theta_weighted_component(field_1d, GOComponent::Beam, |theta, s11| {
                    theta.sin() * s11 / k.powi(2)
                });
            let asymmetry_scatt_beam =
                integrate_theta_weighted_component(field_1d, GOComponent::Beam, |theta, s11| {
                    theta.sin() * theta.cos() * s11 / k.powi(2)
                });
            self.params
                .set_param(Param::ScatCross, GOComponent::Beam, scatt_beam);
            self.params.set_param(
                Param::AsymmetryScatt,
                GOComponent::Beam,
                asymmetry_scatt_beam,
            );

            // Ext field
            let scatt_ext =
                integrate_theta_weighted_component(field_1d, GOComponent::ExtDiff, |theta, s11| {
                    theta.sin() * s11 / k.powi(2)
                });
            let asymmetry_scatt_ext =
                integrate_theta_weighted_component(field_1d, GOComponent::Beam, |theta, s11| {
                    theta.sin() * theta.cos() * s11 / k.powi(2)
                });
            self.params
                .set_param(Param::ScatCross, GOComponent::ExtDiff, scatt_ext);
            self.params.set_param(
                Param::AsymmetryScatt,
                GOComponent::ExtDiff,
                asymmetry_scatt_ext,
            );
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

#[pymethods]
impl Results {
    fn __add__(&self, other: &Results) -> Results {
        self.clone() + other.clone()
    }

    fn __sub__(&self, other: &Results) -> Results {
        self.clone() - other.clone()
    }

    fn __pow__(&self, exponent: f32, _modulo: Option<u32>) -> Results {
        self.clone().pow(exponent)
    }

    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Results> {
        // Try to extract as Results first
        if let Ok(other_results) = other.extract::<Results>() {
            return Ok(self.clone() * other_results);
        }

        // Try to extract as f32
        if let Ok(scalar) = other.extract::<f32>() {
            return Ok(self.clone() * scalar);
        }

        // If neither works, return an error
        Err(pyo3::exceptions::PyTypeError::new_err(
            "unsupported operand type(s) for *: 'Results' and the provided type",
        ))
    }

    fn __truediv__(&self, rhs: f32) -> Results {
        self.clone() / rhs
    }
    /// Get the bins as a numpy array of shape (n_bins, 2) with columns [theta, phi]
    #[getter]
    pub fn get_bins<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let bins: Vec<f32> = self
            .bins()
            .iter()
            .flat_map(|bin| vec![bin.theta_bin.center, bin.phi_bin.center])
            .collect();

        Array2::from_shape_vec((bins.len() / 2, 2), bins)
            .unwrap()
            .into_pyarray(py)
    }

    /// Get the 1D bins (theta values) as a numpy array
    #[getter]
    pub fn get_bins_1d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        self.field_1d.as_ref().map(|field_1d| {
            let bins: Vec<f32> = field_1d.iter().map(|result| result.bin.center).collect();
            Array2::from_shape_vec((bins.len(), 1), bins)
                .unwrap()
                .into_pyarray(py)
        })
    }

    /// Get the Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let muellers: Vec<f32> = self
            .field_2d
            .iter()
            .flat_map(|r| r.mueller_total.to_vec())
            .collect();

        Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
            .unwrap()
            .into_pyarray(py)
    }

    /// Set the Mueller matrix from a numpy array
    #[setter]
    pub fn set_mueller(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        // unsafe view in numpy array memory without bounds checking
        let array_view = unsafe { array.as_array() };

        for (i, field) in self.field_2d.iter_mut().enumerate() {
            let row = array_view.row(i);
            let slice = row.as_slice().unwrap();
            field.mueller_total = Mueller::from_row_slice(slice);
        }
    }

    /// Get the beam Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_beam<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let muellers: Vec<f32> = self
            .field_2d
            .iter()
            .flat_map(|r| r.mueller_beam.to_vec())
            .collect();
        Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
            .unwrap()
            .into_pyarray(py)
    }

    /// Set the Mueller beam matrix from a numpy array
    #[setter]
    pub fn set_mueller_beam(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        // unsafe view in numpy array memory without bounds checking
        let array_view = unsafe { array.as_array() };

        for (i, field) in self.field_2d.iter_mut().enumerate() {
            let row = array_view.row(i);
            let slice = row.as_slice().unwrap();
            field.mueller_beam = Mueller::from_row_slice(slice);
        }
    }

    /// Get the external diffraction Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_ext<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let muellers: Vec<f32> = self
            .field_2d
            .iter()
            .flat_map(|r| r.mueller_ext.to_vec())
            .collect();
        Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
            .unwrap()
            .into_pyarray(py)
    }

    /// Set the Mueller ext matrix from a numpy array
    #[setter]
    pub fn set_mueller_ext(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        // unsafe view in numpy array memory without bounds checking
        let array_view = unsafe { array.as_array() };

        for (i, field) in self.field_2d.iter_mut().enumerate() {
            let row = array_view.row(i);
            let slice = row.as_slice().unwrap();
            field.mueller_ext = Mueller::from_row_slice(slice);
        }
    }

    /// Get the 1D Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_1d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<f32> = field_1d
                .iter()
                .flat_map(|r| r.mueller_total.to_vec())
                .collect();
            Some(
                Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                    .unwrap()
                    .into_pyarray(py),
            )
        } else {
            None
        }
    }

    /// Set the 1D Mueller matrix from a numpy array
    #[setter]
    pub fn set_mueller_1d(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        let array_view = unsafe { array.as_array() };

        if let Some(ref mut field_1d) = self.field_1d {
            for (i, field) in field_1d.iter_mut().enumerate() {
                let row = array_view.row(i);
                let slice = row.as_slice().unwrap();
                field.mueller_total = Mueller::from_row_slice(slice);
            }
        }
    }

    /// Get the 1D beam Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_1d_beam<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<f32> = field_1d
                .iter()
                .flat_map(|r| r.mueller_beam.to_vec())
                .collect();
            Some(
                Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                    .unwrap()
                    .into_pyarray(py),
            )
        } else {
            None
        }
    }

    /// Set the 1D Mueller beam matrix from a numpy array
    #[setter]
    pub fn set_mueller_1d_beam(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        let array_view = unsafe { array.as_array() };

        if let Some(ref mut field_1d) = self.field_1d {
            for (i, field) in field_1d.iter_mut().enumerate() {
                let row = array_view.row(i);
                let slice = row.as_slice().unwrap();
                field.mueller_beam = Mueller::from_row_slice(slice);
            }
        }
    }

    /// Get the 1D external diffraction Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_1d_ext<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<f32> = field_1d
                .iter()
                .flat_map(|r| r.mueller_ext.to_vec())
                .collect();
            Some(
                Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                    .unwrap()
                    .into_pyarray(py),
            )
        } else {
            None
        }
    }

    /// Set the 1D Mueller beam matrix from a numpy array
    #[setter]
    pub fn set_mueller_1d_ext(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        let array_view = unsafe { array.as_array() };

        if let Some(ref mut field_1d) = self.field_1d {
            for (i, field) in field_1d.iter_mut().enumerate() {
                let row = array_view.row(i);
                let slice = row.as_slice().unwrap();
                field.mueller_ext = Mueller::from_row_slice(slice);
            }
        }
    }

    /// Get the asymmetry parameter
    #[getter]
    pub fn get_asymmetry(&self) -> Option<f32> {
        self.params.asymmetry(&GOComponent::Total)
    }

    /// Get the scattering cross section
    #[getter]
    pub fn get_scat_cross(&self) -> Option<f32> {
        self.params.scatt_cross(&GOComponent::Total)
    }

    /// Get the extinction cross section
    #[getter]
    pub fn get_ext_cross(&self) -> Option<f32> {
        self.params.ext_cross(&GOComponent::Total)
    }

    /// Get the albedo
    #[getter]
    pub fn get_albedo(&self) -> Option<f32> {
        self.params.albedo(&GOComponent::Total)
    }

    /// Get the powers as a dictionary
    #[getter]
    pub fn get_powers(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("input", self.powers.input)?;
            dict.set_item("output", self.powers.output)?;
            dict.set_item("absorbed", self.powers.absorbed)?;
            dict.set_item("trnc_ref", self.powers.trnc_ref)?;
            dict.set_item("trnc_rec", self.powers.trnc_rec)?;
            dict.set_item("trnc_clip", self.powers.trnc_clip)?;
            dict.set_item("trnc_energy", self.powers.trnc_energy)?;
            dict.set_item("clip_err", self.powers.clip_err)?;
            dict.set_item("trnc_area", self.powers.trnc_area)?;
            dict.set_item("trnc_cop", self.powers.trnc_cop)?;
            dict.set_item("ext_diff", self.powers.ext_diff)?;
            dict.set_item("missing", self.powers.missing())?;
            Ok(dict.into())
        })
    }

    /// Set the powers from a dictionary
    #[setter]
    pub fn set_powers(&mut self, dict: &Bound<'_, pyo3::types::PyDict>) -> PyResult<()> {
        if let Some(val) = dict.get_item("input")? {
            self.powers.input = val.extract()?;
        }
        if let Some(val) = dict.get_item("output")? {
            self.powers.output = val.extract()?;
        }
        if let Some(val) = dict.get_item("absorbed")? {
            self.powers.absorbed = val.extract()?;
        }
        if let Some(val) = dict.get_item("trnc_ref")? {
            self.powers.trnc_ref = val.extract()?;
        }
        if let Some(val) = dict.get_item("trnc_rec")? {
            self.powers.trnc_rec = val.extract()?;
        }
        if let Some(val) = dict.get_item("trnc_clip")? {
            self.powers.trnc_clip = val.extract()?;
        }
        if let Some(val) = dict.get_item("trnc_energy")? {
            self.powers.trnc_energy = val.extract()?;
        }
        if let Some(val) = dict.get_item("clip_err")? {
            self.powers.clip_err = val.extract()?;
        }
        if let Some(val) = dict.get_item("trnc_area")? {
            self.powers.trnc_area = val.extract()?;
        }
        if let Some(val) = dict.get_item("trnc_cop")? {
            self.powers.trnc_cop = val.extract()?;
        }
        if let Some(val) = dict.get_item("ext_diff")? {
            self.powers.ext_diff = val.extract()?;
        }
        // Note: "missing" is computed, not stored, so we skip it
        Ok(())
    }
}

/// Helper function to integrate over theta for a specific component with custom weighting
fn integrate_theta_weighted_component<F>(
    field_1d: &[ScattResult1D],
    component: GOComponent,
    weight_fn: F,
) -> f32
where
    F: Fn(f32, f32) -> f32, // (theta_radians, s11_value) -> weighted_value
{
    let sum: f32 = field_1d
        .iter()
        .map(|result| {
            let mueller = match component {
                GOComponent::Total => result.mueller_total,
                GOComponent::Beam => result.mueller_beam,
                GOComponent::ExtDiff => result.mueller_ext,
            };
            let s11 = mueller.s11();
            let theta_rad = result.bin.center.to_radians();
            let bin_width = result.bin.width().to_radians();
            weight_fn(theta_rad, s11) * bin_width
        })
        .sum();

    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Complex;

    #[test]
    fn test_ampl_approx_eq() {
        // Create two similar amplitude matrices
        let ampl1 = Ampl::new(
            Complex::new(1.0, 2.0),
            Complex::new(3.0, 4.0),
            Complex::new(5.0, 6.0),
            Complex::new(7.0, 8.0),
        );

        let ampl2 = Ampl::new(
            Complex::new(1.001, 2.001),
            Complex::new(3.001, 4.001),
            Complex::new(5.001, 6.001),
            Complex::new(7.001, 8.001),
        );

        // Should be equal within tolerance
        assert!(ampl1.approx_eq(&ampl2, 0.01));

        // Should not be equal with stricter tolerance
        assert!(!ampl1.approx_eq(&ampl2, 0.0001));

        // Test with exact equality
        assert!(ampl1.approx_eq(&ampl1, 0.0));
    }
}
