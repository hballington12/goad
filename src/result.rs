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

/// Aggregated results from electromagnetic scattering simulations.
/// 
/// **Context**: Scattering simulations produce multiple types of output data including
/// power tracking, Mueller matrices at different angular bins, amplitude matrices
/// for polarization analysis, and derived parameters like cross sections. This
/// structure consolidates all results for analysis and output.
/// 
/// **How it Works**: Stores raw amplitude matrices and derived Mueller matrices for
/// total scattering, beam-only contributions, and external diffraction. Optionally
/// includes phi-integrated 1D results for azimuthally symmetric cases. The params
/// field contains derived integral parameters computed from the angular distributions.
#[derive(Debug, Clone)]
pub struct Results {
    pub powers: Powers,
    pub bins: Vec<(f32, f32)>,
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
