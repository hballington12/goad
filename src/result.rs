use std::f32::consts::PI;
use std::fmt::Debug;

use crate::bins::AngleBin;
use crate::bins::Scheme;
use crate::bins::SolidAngleBin;
use crate::params::Params;
use crate::powers::Powers;
use itertools::Itertools;
use nalgebra::Matrix4;
use nalgebra::{Complex, Matrix2};
use pyo3::prelude::*;

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

#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, serde::Serialize)]
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

    /// Computes the parameters of the result
    pub fn compute_params(&mut self, wavelength: f32) {
        // Compute all 4 parameters for Total
        self.compute_scat_cross(wavelength, GOComponent::Total);
        self.compute_asymmetry(wavelength, GOComponent::Total);
        self.compute_ext_cross(GOComponent::Total);
        self.compute_albedo(GOComponent::Total);

        // Compute scat_cross and asymmetry for Beam
        self.compute_scat_cross(wavelength, GOComponent::Beam);
        self.compute_asymmetry(wavelength, GOComponent::Beam);

        // Compute scat_cross and asymmetry for ExtDiff
        self.compute_scat_cross(wavelength, GOComponent::ExtDiff);
        self.compute_asymmetry(wavelength, GOComponent::ExtDiff);
    }

    pub fn compute_asymmetry(&mut self, wavelength: f32, component: GOComponent) {
        if let Some(field_1d) = &self.field_1d {
            if let Some(scatt) = self.params.scat_cross.get(&component) {
                let k = 2.0 * PI / wavelength;
                let asymmetry =
                    integrate_theta_weighted_component(field_1d, component, |theta, s11| {
                        theta.sin() * theta.cos() * s11 / (scatt * k.powi(2))
                    });

                self.params.asymmetry.insert(component, asymmetry);
            }
        }
    }

    /// Computes the scattering cross section from the 1D Mueller matrix
    pub fn compute_scat_cross(&mut self, wavelength: f32, component: GOComponent) {
        if let Some(field_1d) = &self.field_1d {
            let k = 2.0 * PI / wavelength;
            let scat_cross =
                integrate_theta_weighted_component(field_1d, component, |theta, s11| {
                    theta.sin() * s11 / k.powi(2)
                });

            self.params.scat_cross.insert(component, scat_cross);
        }
    }

    /// Computes the extinction cross section from the scattering cross section and absorbed power
    pub fn compute_ext_cross(&mut self, component: GOComponent) {
        if let Some(scat) = self.params.scat_cross.get(&component) {
            // For Total component, add absorbed power; for others, ext = scat (no absorption in partial components)
            let ext = match component {
                GOComponent::Total => scat + self.powers.absorbed,
                GOComponent::Beam | GOComponent::ExtDiff => *scat,
            };
            self.params.ext_cross.insert(component, ext);
        }
    }

    /// Computes the albedo from the scattering and extinction cross sections
    pub fn compute_albedo(&mut self, component: GOComponent) {
        if let (Some(scat), Some(ext)) = (
            self.params.scat_cross.get(&component),
            self.params.ext_cross.get(&component),
        ) {
            if *ext > 0.0 {
                self.params.albedo.insert(component, scat / ext);
            }
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

            if let Some(val) = self.params.asymmetry.get(&component) {
                println!("{} Asymmetry: {}", comp_str, val);
            }
            if let Some(val) = self.params.scat_cross.get(&component) {
                println!("{} Scat Cross: {}", comp_str, val);
            }
            if let Some(val) = self.params.ext_cross.get(&component) {
                println!("{} Ext Cross: {}", comp_str, val);
            }
            if let Some(val) = self.params.albedo.get(&component) {
                println!("{} Albedo: {}", comp_str, val);
            }
        }
        if let Some(val) = self.params.scat_cross.get(&GOComponent::Beam) {
            println!(
                "Beam Scat Cross / Output power: {}",
                val / self.powers.output
            );
        }
    }
}

#[pymethods]
impl Results {
    /// Get the bins as a list of tuples (returns bin centers for backwards compatibility)
    #[getter]
    pub fn get_bins(&self) -> Vec<(f32, f32)> {
        self.bins()
            .iter()
            .map(|bin| (bin.theta_bin.center, bin.phi_bin.center))
            .collect()
    }

    /// Get the 1D bins (theta values)
    #[getter]
    pub fn get_bins_1d(&self) -> Option<Vec<f32>> {
        self.field_1d
            .as_ref()
            .map(|field_1d| field_1d.iter().map(|result| result.bin.center).collect())
    }

    /// Get the Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().map(|r| r.mueller_total).collect();
        crate::problem::collect_mueller(&muellers)
    }

    /// Get the beam Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_beam(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().map(|r| r.mueller_beam).collect();
        crate::problem::collect_mueller(&muellers)
    }

    /// Get the external diffraction Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_ext(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().map(|r| r.mueller_ext).collect();
        crate::problem::collect_mueller(&muellers)
    }

    /// Get the 1D Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().map(|r| r.mueller_total).collect();
            crate::problem::collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D beam Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d_beam(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().map(|r| r.mueller_beam).collect();
            crate::problem::collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D external diffraction Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d_ext(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().map(|r| r.mueller_ext).collect();
            crate::problem::collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the asymmetry parameter
    #[getter]
    pub fn get_asymmetry(&self) -> Option<f32> {
        self.params.asymmetry()
    }

    /// Get the scattering cross section
    #[getter]
    pub fn get_scat_cross(&self) -> Option<f32> {
        self.params.scat_cross()
    }

    /// Get the extinction cross section
    #[getter]
    pub fn get_ext_cross(&self) -> Option<f32> {
        self.params.ext_cross()
    }

    /// Get the albedo
    #[getter]
    pub fn get_albedo(&self) -> Option<f32> {
        self.params.albedo()
    }

    /// Get all parameters as a dictionary
    #[getter]
    pub fn get_params(&self) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            let dict = pyo3::types::PyDict::new(py);

            // Add backwards-compatible top-level keys for Total component
            dict.set_item("asymmetry", self.params.asymmetry())?;
            dict.set_item("scat_cross", self.params.scat_cross())?;
            dict.set_item("ext_cross", self.params.ext_cross())?;
            dict.set_item("albedo", self.params.albedo())?;

            // Add component-specific dictionaries
            for (comp_name, component) in [
                ("total", GOComponent::Total),
                ("beam", GOComponent::Beam),
                ("ext_diff", GOComponent::ExtDiff),
            ] {
                let comp_dict = pyo3::types::PyDict::new(py);

                if let Some(val) = self.params.asymmetry.get(&component) {
                    comp_dict.set_item("asymmetry", val)?;
                }
                if let Some(val) = self.params.scat_cross.get(&component) {
                    comp_dict.set_item("scat_cross", val)?;
                }
                if let Some(val) = self.params.ext_cross.get(&component) {
                    comp_dict.set_item("ext_cross", val)?;
                }
                if let Some(val) = self.params.albedo.get(&component) {
                    comp_dict.set_item("albedo", val)?;
                }

                dict.set_item(comp_name, comp_dict)?;
            }

            Ok(dict.into())
        })
    }

    /// Get the powers as a dictionary
    #[getter]
    pub fn get_powers(&self) -> PyResult<PyObject> {
        Python::with_gil(|py| {
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
