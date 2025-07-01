//! Simulation result aggregation and analysis.
//!
//! This module handles the collection, processing, and analysis of electromagnetic
//! scattering simulation results. It manages multiple data types including amplitude
//! matrices, Mueller matrices, power conservation data, and derived integral
//! parameters for comprehensive scattering characterization.
//!
//! The result system provides:
//! - Multi-component result storage (amplitude, Mueller, power)
//! - 2D to 1D integration for azimuthally symmetric cases
//! - Integral parameter computation (cross sections, asymmetry, albedo)
//! - Power conservation analysis and validation
//! - Result combination for multi-orientation averaging
//! - Comprehensive output formatting
//!
//! # Key Components
//!
//! - [`Results`]: Main result container with all simulation outputs
//! - Integration utilities for parameter computation
//! - Mueller matrix manipulation and analysis
//! - Cross section and albedo calculations
//! - Statistical analysis for ensemble results

use std::f32::consts::PI;

use crate::output;
use crate::params::Params;
use crate::powers::Powers;
use anyhow::anyhow;
use anyhow::Result;
use itertools::Itertools;
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2};
use ndarray::{s, Array1, Array2, Axis};
use pyo3::prelude::*;

/// Aggregated results from electromagnetic scattering simulations.
/// 
/// **Context**: Scattering simulations produce multiple types of output data including
/// power tracking, Mueller matrices at different angular bins, amplitude matrices
/// for polarization analysis, and derived parameters like cross sections. This
/// structure consolidates all results for analysis and [`crate::output`].
/// 
/// **How it Works**: Stores raw amplitude matrices and derived Mueller matrices for
/// total scattering, beam-only contributions, and external diffraction. Optionally
/// includes phi-integrated 1D results for azimuthally symmetric cases. The [`crate::params::Params`]
/// field contains derived integral parameters computed from the angular distributions.
#[pyclass]
#[derive(Debug, Clone)]
pub struct Results {
    pub powers: Powers,           // [`crate::powers::Powers`] conservation data
    pub bins: Vec<(f32, f32)>,    // Angular bins (theta, phi) for 2D results
    pub mueller: Array2<f32>,
    pub mueller_beam: Array2<f32>,
    pub mueller_ext: Array2<f32>,
    pub ampl: Vec<Matrix2<Complex<f32>>>,
    pub ampl_beam: Vec<Matrix2<Complex<f32>>>,
    pub ampl_ext: Vec<Matrix2<Complex<f32>>>,
    pub bins_1d: Option<Vec<f32>>,
    pub mueller_1d: Option<Array2<f32>>,
    pub mueller_1d_beam: Option<Array2<f32>>,
    pub mueller_1d_ext: Option<Array2<f32>>,
    pub params: Params,
}

impl Results {
    /// Creates empty result storage for the specified angular bins.
    /// 
    /// **Context**: Results structures must be pre-allocated with the correct
    /// dimensions based on the angular binning scheme before simulation begins.
    /// 
    /// **How it Works**: Allocates zero-filled arrays for Mueller matrices
    /// (16 elements per bin) and amplitude matrices (2x2 complex per bin).
    pub fn new_empty(bins: &[(f32, f32)]) -> Self {
        let mueller = Array2::<f32>::zeros((bins.len(), 16));
        let mueller_beam = mueller.clone();
        let mueller_ext = mueller.clone();
        let ampl = vec![Matrix2::<Complex<f32>>::zeros(); bins.len()];
        let ampl_beam = ampl.clone();
        let ampl_ext = ampl.clone();
        Self {
            powers: Powers::new(),
            bins: bins.to_vec(),
            mueller,
            mueller_beam,
            mueller_ext,
            ampl,
            ampl_beam,
            ampl_ext,
            bins_1d: None,
            mueller_1d: None,
            mueller_1d_beam: None,
            mueller_1d_ext: None,
            params: Params::new(),
        }
    }

    /// Attempts phi integration to produce 1D theta-dependent results.
    /// 
    /// **Context**: Many scattering problems exhibit azimuthal symmetry, allowing
    /// reduction from 2D (theta, phi) to 1D (theta) results through phi integration.
    /// This simplifies analysis and visualization.
    /// 
    /// **How it Works**: Calls the integration function and stores the 1D results
    /// if successful, enabling computation of integral parameters.
    pub fn try_mueller_to_1d(&mut self) -> std::result::Result<(), anyhow::Error> {
        match try_mueller_to_1d(&self.bins, &self.mueller) {
            Ok((theta, mueller_1d)) => {
                self.bins_1d = Some(theta);
                self.mueller_1d = Some(mueller_1d);

                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    /// Computes integral scattering parameters from angular distributions.
    /// 
    /// **Context**: Key scattering parameters like cross sections, asymmetry parameter,
    /// and albedo are computed by integrating the Mueller matrix over solid angle.
    /// These parameters enable comparison with analytical theories and experiments.
    /// 
    /// **How it Works**: Sequentially computes scattering cross section, asymmetry
    /// parameter, extinction cross section, and single-scattering albedo from the
    /// 1D Mueller matrix and power absorption data.
    pub fn compute_params(&mut self, wavelength: f32) -> std::result::Result<(), anyhow::Error> {
        self.compute_scat_cross(wavelength);
        self.compute_asymmetry(wavelength);
        self.compute_ext_cross();
        self.compute_albedo();

        Ok(())
    }

    /// Computes the asymmetry parameter from the phase function.
    /// 
    /// **Context**: The asymmetry parameter quantifies the average cosine of the
    /// scattering angle, indicating forward (positive) or backward (negative)
    /// scattering preference. This parameter is essential for radiative transfer.
    /// 
    /// **How it Works**: Integrates P11 * cos(theta) * sin(theta) over theta,
    /// normalized by the scattering cross section.
    pub fn compute_asymmetry(&mut self, wavelength: f32) {
        if let (Some(theta), Some(mueller_1d), Some(scatt)) =
            (&self.bins_1d, &self.mueller_1d, self.params.scat_cross)
        {
            self.params.asymettry = Some(compute_asymmetry(
                theta,
                mueller_1d,
                2.0 * PI / wavelength,
                scatt,
            ));
        }
    }

    /// Computes the scattering cross section from angular distribution.
    /// 
    /// **Context**: The scattering cross section quantifies the total scattered
    /// power relative to incident intensity. This fundamental parameter enables
    /// comparison with Mie theory and experimental measurements.
    /// 
    /// **How it Works**: Integrates the P11 Mueller element over solid angle
    /// using the trapezoidal rule, accounting for the sin(theta) Jacobian.
    pub fn compute_scat_cross(&mut self, wavelength: f32) {
        if let (Some(theta), Some(mueller_1d)) = (&self.bins_1d, &self.mueller_1d) {
            self.params.scat_cross =
                Some(compute_scat_cross(theta, mueller_1d, 2.0 * PI / wavelength));
        }
    }

    /// Computes the extinction cross section from optical theorem.
    /// 
    /// **Context**: The extinction cross section represents total power removed
    /// from the incident beam through both scattering and absorption. The optical
    /// theorem relates this to the forward scattering amplitude.
    /// 
    /// **How it Works**: Sums the scattering cross section and absorbed power
    /// to get total extinction following energy conservation.
    pub fn compute_ext_cross(&mut self) {
        match self.params.scat_cross {
            Some(scat) => {
                self.params.ext_cross = Some(scat + self.powers.absorbed);
            }
            None => {
                self.params.ext_cross = None;
            }
        }
    }

    /// Computes the single-scattering albedo.
    /// 
    /// **Context**: The albedo quantifies the fraction of extinction due to
    /// scattering versus absorption. This parameter is crucial for radiative
    /// transfer in atmospheric and oceanic applications.
    /// 
    /// **How it Works**: Divides scattering cross section by extinction cross
    /// section to get the scattering probability.
    pub fn compute_albedo(&mut self) {
        if let (Some(scat), Some(ext)) = (self.params.scat_cross, self.params.ext_cross) {
            self.params.albedo = Some(scat / ext);
        }
    }

    pub fn print(&self) {
        println!("Powers: {:?}", self.powers);
        println!("Asymmetry: {:?}", self.params.asymettry);
        println!("Scat Cross: {:?}", self.params.scat_cross);
        println!("Ext Cross: {:?}", self.params.ext_cross);
        println!("Albedo: {:?}", self.params.albedo);
    }
}

#[pymethods]
impl Results {
    /// Get the bins as a list of tuples
    #[getter]
    pub fn get_bins(&self) -> Vec<(f32, f32)> {
        self.bins.clone()
    }

    /// Get the 1D bins (theta values)
    #[getter]
    pub fn get_bins_1d(&self) -> Option<Vec<f32>> {
        self.bins_1d.clone()
    }

    /// Get the Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller(&self) -> Vec<Vec<f32>> {
        crate::problem::collect_mueller(&self.mueller)
    }

    /// Get the beam Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_beam(&self) -> Vec<Vec<f32>> {
        crate::problem::collect_mueller(&self.mueller_beam)
    }

    /// Get the external diffraction Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_ext(&self) -> Vec<Vec<f32>> {
        crate::problem::collect_mueller(&self.mueller_ext)
    }

    /// Get the 1D Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d(&self) -> Vec<Vec<f32>> {
        if let Some(ref mueller_1d) = self.mueller_1d {
            crate::problem::collect_mueller(mueller_1d)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D beam Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d_beam(&self) -> Vec<Vec<f32>> {
        if let Some(ref mueller_1d_beam) = self.mueller_1d_beam {
            crate::problem::collect_mueller(mueller_1d_beam)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D external diffraction Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d_ext(&self) -> Vec<Vec<f32>> {
        if let Some(ref mueller_1d_ext) = self.mueller_1d_ext {
            crate::problem::collect_mueller(mueller_1d_ext)
        } else {
            Vec::new()
        }
    }

    /// Get the asymmetry parameter
    #[getter]
    pub fn get_asymmetry(&self) -> Option<f32> {
        self.params.asymettry
    }

    /// Get the scattering cross section
    #[getter]
    pub fn get_scat_cross(&self) -> Option<f32> {
        self.params.scat_cross
    }

    /// Get the extinction cross section
    #[getter]
    pub fn get_ext_cross(&self) -> Option<f32> {
        self.params.ext_cross
    }

    /// Get the albedo
    #[getter]
    pub fn get_albedo(&self) -> Option<f32> {
        self.params.albedo
    }

    /// Get all parameters as a dictionary
    #[getter]
    pub fn get_params(&self) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("asymmetry", self.params.asymettry)?;
            dict.set_item("scat_cross", self.params.scat_cross)?;
            dict.set_item("ext_cross", self.params.ext_cross)?;
            dict.set_item("albedo", self.params.albedo)?;
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

/// Integrates 2D Mueller matrix over azimuthal angle to produce 1D results.
/// 
/// **Context**: Azimuthally symmetric scattering patterns allow reduction from
/// 2D (theta, phi) to 1D (theta) representation through phi integration. This
/// simplifies analysis and enables direct comparison with 1D scattering theories.
/// 
/// **How it Works**: Groups Mueller matrix elements by theta value, integrates
/// each group over phi using the trapezoidal rule, and returns theta values
/// with corresponding integrated Mueller matrix elements.
pub fn try_mueller_to_1d(
    bins: &[(f32, f32)],
    mueller: &Array2<f32>,
) -> Result<(Vec<f32>, Array2<f32>)> {
    // Check that the mueller matrix and bins are the same length
    if mueller.len_of(Axis(0)) != bins.len() {
        return Err(anyhow!(
            "Mueller matrix and bins must have the same length. Got {} and {}",
            mueller.len_of(Axis(0)),
            bins.len()
        ));
    }

    // Create indices and sort them by corresponding theta values
    let mueller = mueller.to_owned();
    let mut bins = bins.to_owned();
    let mut indices: Vec<usize> = (0..bins.len()).collect();
    indices.sort_by(|&i, &j| bins[i].0.partial_cmp(&bins[j].0).unwrap());

    // Sort bins according to the sorted indices
    bins = indices.iter().map(|&i| bins[i]).collect();

    // Create a new sorted mueller matrix using the same indices
    let mut sorted_mueller = Array2::<f32>::zeros(mueller.dim());
    for (new_idx, &old_idx) in indices.iter().enumerate() {
        sorted_mueller
            .slice_mut(s![new_idx, ..])
            .assign(&mueller.slice(s![old_idx, ..]));
    }

    // zip the bins and mueller matrix
    let combined: Vec<_> = bins
        .iter()
        .zip(mueller.outer_iter())
        .map(|(bin, row)| (*bin, row.to_owned()))
        .collect();

    // group the combined Vec by theta
    let grouped: Vec<Vec<_>> = combined
        .into_iter()
        .chunk_by(|((key, _), _)| *key)
        .into_iter()
        .map(|(_, group)| group.map(|x| x).collect())
        .collect();

    let mut thetas = Vec::new();
    let mut mueller_1d = Array2::<f32>::zeros((grouped.len(), 16));

    // loop over vectors at each theta
    for (i, muellers) in grouped.iter().enumerate() {
        // Unzip the theta, phi, and mueller values
        let thetas_phi: Vec<_> = muellers
            .iter()
            .map(|((theta, phi), _)| (*theta, *phi))
            .collect();
        let mueller_phi: Vec<_> = muellers.iter().map(|(_, mueller)| mueller).collect();
        let mut mueller_1d_row = Array1::<f32>::zeros(16);

        // loop over the mueller values at each phi
        for j in 0..16 {
            // Create 1D arrays for x and y, where x is phi and y is 1 of the 16 mueller values
            let y = Array1::from(mueller_phi.iter().map(|row| row[j]).collect::<Vec<_>>());
            let phi_values: Vec<_> = thetas_phi.iter().map(|(_, phi)| phi.to_radians()).collect();

            // Check if phi values are sorted in ascending order (probably dont need this)
            for i in 1..phi_values.len() {
                if phi_values[i] < phi_values[i - 1] {
                    return Err(anyhow!(
                        "Phi values must be sorted in ascending order for integration"
                    ));
                }
            }

            let x = Array1::from(phi_values);
            mueller_1d_row[j] = output::integrate_trapezoidal(&x, &y, |_, y| y);
            // integrate over phi
        }

        // Assign the theta and mueller values to the final arrays
        thetas.push(thetas_phi[0].0);
        mueller_1d
            .slice_mut(s![i, ..])
            .assign(&mueller_1d_row.slice(s![..]));
    }

    Ok((thetas, mueller_1d))
}

/// Computes asymmetry parameter from 1D Mueller matrix.
/// 
/// **Context**: The asymmetry parameter (g) characterizes the angular distribution
/// of scattered light, with g=0 for isotropic scattering, g>0 for forward
/// scattering, and g<0 for backscattering. Essential for radiative transfer models.
/// 
/// **How it Works**: Integrates P11(theta) * sin(theta) * cos(theta) over theta,
/// normalized by the scattering cross section and wavenumber squared.
pub fn compute_asymmetry(theta: &[f32], mueller_1d: &Array2<f32>, waveno: f32, scatt: f32) -> f32 {
    // get first column of mueller matrix
    let y = mueller_1d.slice(s![.., 0]).to_owned();

    let x = Array1::from(theta.iter().map(|&t| t.to_radians()).collect::<Vec<f32>>());

    // integrate p11 sin(theta) cos(theta) / scatt cross / k^2
    let asymmetry = output::integrate_trapezoidal(&x, &y, |x, y| {
        x.sin() * x.cos() * y / scatt / waveno.powi(2)
    });

    asymmetry
}

/// Computes scattering cross section from 1D Mueller matrix.
/// 
/// **Context**: The scattering cross section quantifies the effective area for
/// electromagnetic scattering, enabling comparison between particles of different
/// sizes and shapes. Normalized by wavenumber squared for dimensionless form.
/// 
/// **How it Works**: Integrates P11(theta) * sin(theta) over theta using
/// trapezoidal rule, accounting for solid angle element and wavenumber normalization.
pub fn compute_scat_cross(theta: &[f32], mueller_1d: &Array2<f32>, waveno: f32) -> f32 {
    // get first column of mueller matrix
    let y = mueller_1d.slice(s![.., 0]).to_owned();

    // convert theta values from degrees to radians for integration
    let x = Array1::from(theta.iter().map(|&t| t.to_radians()).collect::<Vec<f32>>());

    // integrate p11 sin(theta) / k^2
    let scat_cross = output::integrate_trapezoidal(&x, &y, |x, y| x.sin() * y / waveno.powi(2));

    scat_cross
}
