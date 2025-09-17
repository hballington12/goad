use std::f32::consts::PI;
use std::fmt::Debug;

use crate::bins::AngleBin;
use crate::bins::SolidAngleBin;
use crate::params::Params;
use crate::powers::Powers;
#[cfg(feature = "macroquad")]
use macroquad::prelude::*;
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

// type Mueller = Matrix4<f32>;
// type Ampl = Matrix2<Complex<f32>>;

#[derive(Debug, Clone, PartialEq, Hash, Eq)]
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

// /// A Mueller matrix
// #[derive(Debug, Clone)]
// pub struct Mueller {
//     pub meta: ScattResultMeta,
//     pub matrix: Matrix4<f32>,
// }

// /// An amplitude matrix
// #[derive(Debug, Clone)]
// pub struct Ampl {
//     pub meta: ScattResultMeta,
//     pub matrix: Matrix2<Complex<f32>>,
// }

pub type Ampl = Matrix2<Complex<f32>>;
pub type Mueller = Matrix4<f32>;

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

pub trait AmplMatrix {
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
    fn to_mueller(&self) -> Mueller;
}

impl AmplMatrix for Ampl {
    fn s11(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj()
            + self[(0, 1)] * self[(0, 1)].conj()
            + self[(1, 0)] * self[(1, 0)].conj()
            + self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s12(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj() - self[(0, 1)] * self[(0, 1)].conj()
            + self[(1, 0)] * self[(1, 0)].conj()
            - self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s13(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() + self[(1, 1)] * self[(1, 0)].conj()).re
    }
    fn s14(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() - self[(1, 1)] * self[(1, 0)].conj()).im
    }
    fn s21(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj() + self[(0, 1)] * self[(0, 1)].conj()
            - self[(1, 0)] * self[(1, 0)].conj()
            - self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s22(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj()
            - self[(0, 1)] * self[(0, 1)].conj()
            - self[(1, 0)] * self[(1, 0)].conj()
            + self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s23(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() - self[(1, 1)] * self[(1, 0)].conj()).re
    }
    fn s24(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() + self[(1, 1)] * self[(1, 0)].conj()).im
    }
    fn s31(&self) -> f32 {
        (self[(0, 0)] * self[(1, 0)].conj() + self[(1, 1)] * self[(0, 1)].conj()).re
    }
    fn s32(&self) -> f32 {
        (self[(0, 0)] * self[(1, 0)].conj() - self[(1, 1)] * self[(0, 1)].conj()).re
    }
    fn s33(&self) -> f32 {
        (self[(0, 0)] * self[(1, 1)].conj() + self[(0, 1)] * self[(1, 0)].conj()).re
    }
    fn s34(&self) -> f32 {
        (self[(0, 0)] * self[(1, 1)].conj() + self[(0, 1)] * self[(1, 0)].conj()).im
    }
    fn s41(&self) -> f32 {
        (self[(1, 0)] * self[(0, 0)].conj() + self[(1, 1)] * self[(0, 1)].conj()).im
    }
    fn s42(&self) -> f32 {
        (self[(1, 0)] * self[(0, 0)].conj() - self[(1, 1)] * self[(0, 1)].conj()).im
    }
    fn s43(&self) -> f32 {
        (self[(1, 1)] * self[(0, 0)].conj() - self[(0, 1)] * self[(1, 0)].conj()).im
    }
    fn s44(&self) -> f32 {
        (self[(1, 1)] * self[(0, 0)].conj() - self[(0, 1)] * self[(1, 0)].conj()).re
    }
    fn to_mueller(&self) -> Mueller {
        Mueller::new(
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
        )
    }
}

/// A generic far-field scattering result that can be 1D or 2D.
#[derive(Debug, Clone)]
pub struct ScattResult<B: ScatteringBin> {
    pub bin: B,
    pub ampl_total: Option<Ampl>,
    pub ampl_beam: Option<Ampl>,
    pub ampl_ext: Option<Ampl>,
    pub mueller_total: Option<Mueller>,
    pub mueller_beam: Option<Mueller>,
    pub mueller_ext: Option<Mueller>,
}

impl<B: ScatteringBin> ScattResult<B> {
    /// Creates a new empty ScattResult.
    pub fn new_empty(bin: B) -> Self {
        Self {
            bin,
            ampl_total: None,
            ampl_beam: None,
            ampl_ext: None,
            mueller_total: None,
            mueller_beam: None,
            mueller_ext: None,
        }
    }
}

/// Type alias for 2D scattering results (full solid angle)
pub type ScattResult2D = ScattResult<SolidAngleBin>;

/// Type alias for 1D scattering results (theta only)
pub type ScattResult1D = ScattResult<AngleBin>;

// /// A mueller matrix
// #[derive(Debug, Clone, Default)]
// pub struct Mueller {
//     pub s11: f32,
//     pub s12: f32,
//     pub s13: f32,
//     pub s14: f32,
//     pub s21: f32,
//     pub s22: f32,
//     pub s23: f32,
//     pub s24: f32,
//     pub s31: f32,
//     pub s32: f32,
//     pub s33: f32,
//     pub s34: f32,
//     pub s41: f32,
//     pub s42: f32,
//     pub s43: f32,
//     pub s44: f32,
// }

// impl Mueller {
//     // /// Creates a blank Mueller matrix.
//     // pub fn new_empty() -> Self {
//     //     Self {
//     //         s11: 0.0,
//     //         s12: 0.0,
//     //         s13: 0.0,
//     //         s14: 0.0,
//     //         s21: 0.0,
//     //         s22: 0.0,
//     //         s23: 0.0,
//     //         s24: 0.0,
//     //         s31: 0.0,
//     //         s32: 0.0,
//     //         s33: 0.0,
//     //         s34: 0.0,
//     //         s41: 0.0,
//     //         s42: 0.0,
//     //         s43: 0.0,
//     //         s44: 0.0,
//     //     }
//     // }

//     pub fn s11(&self) -> f32 {
//         self.matrix[(0, 0)]
//     }

//     pub fn s12(&self) -> f32 {
//         self.matrix[(0, 1)]
//     }

//     pub fn s13(&self) -> f32 {
//         self.matrix[(0, 2)]
//     }

//     pub fn s14(&self) -> f32 {
//         self.matrix[(0, 3)]
//     }

//     pub fn s21(&self) -> f32 {
//         self.matrix[(1, 0)]
//     }

//     pub fn s22(&self) -> f32 {
//         self.matrix[(1, 1)]
//     }

//     pub fn s23(&self) -> f32 {
//         self.matrix[(1, 2)]
//     }

//     pub fn s24(&self) -> f32 {
//         self.matrix[(1, 3)]
//     }

//     pub fn s31(&self) -> f32 {
//         self.matrix[(2, 0)]
//     }

//     pub fn s32(&self) -> f32 {
//         self.matrix[(2, 1)]
//     }

//     pub fn s33(&self) -> f32 {
//         self.matrix[(2, 2)]
//     }

//     pub fn s34(&self) -> f32 {
//         self.matrix[(2, 3)]
//     }

//     pub fn s41(&self) -> f32 {
//         self.matrix[(3, 0)]
//     }

//     pub fn s42(&self) -> f32 {
//         self.matrix[(3, 1)]
//     }

//     pub fn s43(&self) -> f32 {
//         self.matrix[(3, 2)]
//     }

//     pub fn s44(&self) -> f32 {
//         self.matrix[(3, 3)]
//     }

//     /// Returns the Mueller matrix as a vector of its elements.
//     pub fn to_vec(&self) -> Vec<f32> {
//         vec![
//             self.s11(),
//             self.s12(),
//             self.s13(),
//             self.s14(),
//             self.s21(),
//             self.s22(),
//             self.s23(),
//             self.s24(),
//             self.s31(),
//             self.s32(),
//             self.s33(),
//             self.s34(),
//             self.s41(),
//             self.s42(),
//             self.s43(),
//             self.s44(),
//         ]
//     }
// }

// impl Add for Mueller {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         Self {
//             s11: self.s11 + other.s11,
//             s12: self.s12 + other.s12,
//             s13: self.s13 + other.s13,
//             s14: self.s14 + other.s14,
//             s21: self.s21 + other.s21,
//             s22: self.s22 + other.s22,
//             s23: self.s23 + other.s23,
//             s24: self.s24 + other.s24,
//             s31: self.s31 + other.s31,
//             s32: self.s32 + other.s32,
//             s33: self.s33 + other.s33,
//             s34: self.s34 + other.s34,
//             s41: self.s41 + other.s41,
//             s42: self.s42 + other.s42,
//             s43: self.s43 + other.s43,
//             s44: self.s44 + other.s44,
//         }
//     }
// }

// impl AddAssign for Mueller {
//     fn add_assign(&mut self, other: Self) {
//         *self = Self {
//             s11: self.s11 + other.s11,
//             s12: self.s12 + other.s12,
//             s13: self.s13 + other.s13,
//             s14: self.s14 + other.s14,
//             s21: self.s21 + other.s21,
//             s22: self.s22 + other.s22,
//             s23: self.s23 + other.s23,
//             s24: self.s24 + other.s24,
//             s31: self.s31 + other.s31,
//             s32: self.s32 + other.s32,
//             s33: self.s33 + other.s33,
//             s34: self.s34 + other.s34,
//             s41: self.s41 + other.s41,
//             s42: self.s42 + other.s42,
//             s43: self.s43 + other.s43,
//             s44: self.s44 + other.s44,
//         };
//     }
// }

// impl DivAssign<f32> for Mueller {
//     fn div_assign(&mut self, fac: f32) {
//         self.s11 /= fac;
//         self.s12 /= fac;
//         self.s13 /= fac;
//         self.s14 /= fac;
//         self.s21 /= fac;
//         self.s22 /= fac;
//         self.s23 /= fac;
//         self.s24 /= fac;
//         self.s31 /= fac;
//         self.s32 /= fac;
//         self.s33 /= fac;
//         self.s34 /= fac;
//         self.s41 /= fac;
//         self.s42 /= fac;
//         self.s43 /= fac;
//         self.s44 /= fac;
//     }
// }

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
    // pub bins: Vec<(AngleBin, AngleBin)>,
    // pub mueller: Vec<Mueller>,
    // pub mueller_beam: Vec<Mueller>,
    // pub mueller_ext: Vec<Mueller>,
    // pub ampl: Vec<Matrix2<Complex<f32>>>,
    // pub ampl_beam: Vec<Matrix2<Complex<f32>>>,
    // pub ampl_ext: Vec<Matrix2<Complex<f32>>>,
    // pub bins_1d: Option<Vec<f32>>,
    // pub mueller_1d: Option<Vec<Mueller>>,
    // pub mueller_1d_beam: Option<Vec<Mueller>>,
    // pub mueller_1d_ext: Option<Vec<Mueller>>,
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
        let field = bins
            .iter()
            .map(|&bin| ScattResult2D::new_empty(bin))
            .collect();
        Self {
            field_2d: field,
            powers: Powers::new(),
            field_1d: None,
            params: Params::new(),
        }
    }

    pub fn try_mueller_to_1d(&mut self, binning_scheme: &crate::bins::Scheme) {
        use crate::bins::Scheme;
        use itertools::Itertools;

        // Step 1: Check scheme compatibility
        match binning_scheme {
            Scheme::Custom { .. } => return, // Skip custom binning
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
    /// Weighted by phi bin width
    fn integrate_over_phi(phi_group: Vec<&ScattResult2D>) -> ScattResult1D {
        // All results in group have same theta bin
        let theta_bin = phi_group[0].bin.theta_bin;

        // Calculate total phi width and weighted sums
        let mut total_phi_width = 0.0;
        let mut mueller_total_sum = Mueller::zeros();
        let mut mueller_beam_sum = Mueller::zeros();
        let mut mueller_ext_sum = Mueller::zeros();
        let mut has_total = false;
        let mut has_beam = false;
        let mut has_ext = false;

        for result in phi_group {
            let phi_width = result.bin.phi_bin.width();
            total_phi_width += phi_width;

            // Integrate Mueller matrices (weighted by phi bin width)
            if let Some(mueller) = result.mueller_total {
                mueller_total_sum += mueller * phi_width;
                has_total = true;
            }
            if let Some(mueller) = result.mueller_beam {
                mueller_beam_sum += mueller * phi_width;
                has_beam = true;
            }
            if let Some(mueller) = result.mueller_ext {
                mueller_ext_sum += mueller * phi_width;
                has_ext = true;
            }
        }

        // Normalize by total phi width (complete the integration)
        ScattResult1D {
            bin: theta_bin,
            ampl_total: None, // Amplitudes not integrated
            ampl_beam: None,
            ampl_ext: None,
            mueller_total: if has_total {
                Some(mueller_total_sum / total_phi_width)
            } else {
                None
            },
            mueller_beam: if has_beam {
                Some(mueller_beam_sum / total_phi_width)
            } else {
                None
            },
            mueller_ext: if has_ext {
                Some(mueller_ext_sum / total_phi_width)
            } else {
                None
            },
        }
    }

    /// Computes the parameters of the result
    pub fn compute_params(&mut self, wavelength: f32) -> std::result::Result<(), anyhow::Error> {
        self.compute_scat_cross(wavelength);
        self.compute_asymmetry(wavelength);
        self.compute_ext_cross();
        self.compute_albedo();

        Ok(())
    }

    pub fn compute_asymmetry(&mut self, wavelength: f32) {
        if let (Some(field_1d), Some(scatt)) = (&self.field_1d, self.params.scat_cross) {
            let k = 2.0 * PI / wavelength;
            self.params.asymettry = Some(integrate_theta_weighted(field_1d, |theta, s11| {
                theta.sin() * theta.cos() * s11 / (scatt * k.powi(2))
            }));
        }
    }

    /// Computes the scattering cross section from the 1D Mueller matrix
    pub fn compute_scat_cross(&mut self, wavelength: f32) {
        if let Some(field_1d) = &self.field_1d {
            let k = 2.0 * PI / wavelength;
            self.params.scat_cross = Some(integrate_theta_weighted(field_1d, |theta, s11| {
                theta.sin() * s11 / k.powi(2)
            }));
        }
    }

    /// Computes the extinction cross section from the scattering cross section and absorbed power
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

    /// Computes the albedo from the scattering and extinction cross sections
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
        let muellers: Vec<Mueller> = self
            .field_2d
            .iter()
            .filter_map(|r| r.mueller_total)
            .collect();
        crate::problem::collect_mueller(&muellers)
    }

    /// Get the beam Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_beam(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self
            .field_2d
            .iter()
            .filter_map(|r| r.mueller_beam)
            .collect();
        crate::problem::collect_mueller(&muellers)
    }

    /// Get the external diffraction Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_ext(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().filter_map(|r| r.mueller_ext).collect();
        crate::problem::collect_mueller(&muellers)
    }

    /// Get the 1D Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().filter_map(|r| r.mueller_total).collect();
            crate::problem::collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D beam Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d_beam(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().filter_map(|r| r.mueller_beam).collect();
            crate::problem::collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D external diffraction Mueller matrix as a list of lists
    #[getter]
    pub fn get_mueller_1d_ext(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().filter_map(|r| r.mueller_ext).collect();
            crate::problem::collect_mueller(&muellers)
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

/// Helper function to integrate over theta with custom weighting
/// Takes field_1d results and applies a weight function to (theta_radians, s11_value)
fn integrate_theta_weighted<F>(field_1d: &[ScattResult1D], weight_fn: F) -> f32
where
    F: Fn(f32, f32) -> f32, // (theta_radians, s11_value) -> weighted_value
{
    field_1d
        .iter()
        .filter_map(|result| {
            result.mueller_total.map(|mueller| {
                let theta_rad = result.bin.center.to_radians();
                let s11 = mueller.s11();
                let bin_width = result.bin.width().to_radians();
                weight_fn(theta_rad, s11) * bin_width
            })
        })
        .sum()
}
