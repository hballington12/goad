//! Results struct containing complete simulation output.

use std::f32::consts::PI;

use itertools::Itertools;
use pyo3::prelude::*;

use crate::bins::{Scheme, SolidAngleBin};
use crate::params::Params;
use crate::powers::Powers;

use super::component::GOComponent;
use super::integrate::integrate_theta_weighted_component;
use super::scatt_result::{ScattResult1D, ScattResult2D};

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
    /// Returns an owned vector of solid angle bins.
    pub fn bins(&self) -> Vec<SolidAngleBin> {
        self.field_2d.iter().map(|a| a.bin.clone()).collect()
    }

    /// Creates a new `Result` with empty mueller and amplitude matrix.
    pub fn new_empty(bins: &[SolidAngleBin]) -> Self {
        let field = bins.iter().map(|&bin| ScattResult2D::new(bin)).collect();
        Self {
            field_2d: field,
            powers: Powers::new(),
            field_1d: None,
            params: Params::new(),
        }
    }

    /// Converts 2D Mueller matrices to 1D by integrating over phi.
    pub fn mueller_to_1d(&mut self, binning_scheme: &Scheme) {
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

    /// Integrates Mueller matrices over phi using rectangular rule.
    /// Weighted by phi bin width in radians.
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

    /// Computes the parameters of the result.
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

    /// Computes the asymmetry parameter for a given component.
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

    /// Computes the scattering cross section from the 1D Mueller matrix.
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

    /// Computes the extinction cross section from the scattering cross section and absorbed power.
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

    /// Computes the albedo from the scattering and extinction cross sections.
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

    /// Prints a summary of the results.
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
