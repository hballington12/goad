use std::f32::consts::PI;
use std::fmt::Debug;

use crate::bins::AngleBin;
use crate::bins::SolidAngleBin;
use crate::params::Params;
use crate::powers::Powers;
use anyhow::anyhow;
use anyhow::Result;
#[cfg(feature = "macroquad")]
use macroquad::prelude::*;
use nalgebra::Matrix4;
use nalgebra::{Complex, Matrix2};
use ndarray::Array1;
use pyo3::prelude::*;

/// Trait for different types of scattering bins (1D or 2D)
pub trait ScatteringBin: Clone + Debug {
    /// Get the theta center value
    fn theta_center(&self) -> f32;

    /// Get the theta bin
    fn theta_bin(&self) -> &AngleBin;
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

    /// Computes the parameters of the result
    pub fn compute_params(&mut self, wavelength: f32) -> std::result::Result<(), anyhow::Error> {
        self.compute_scat_cross(wavelength);
        self.compute_asymmetry(wavelength);
        self.compute_ext_cross();
        self.compute_albedo();

        Ok(())
    }

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

    /// Computes the scattering cross section from the 1D Mueller matrix
    pub fn compute_scat_cross(&mut self, wavelength: f32) {
        if let (Some(theta), Some(mueller_1d)) = (&self.bins_1d, &self.mueller_1d) {
            self.params.scat_cross =
                Some(compute_scat_cross(theta, mueller_1d, 2.0 * PI / wavelength));
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

/// Integrate over phi (second bin of the tuple) to get the 1D Mueller matrix
/// Uses the trapezoidal rule
/// Returns a tuple of the theta bins and the 1D Mueller matrix
/// NOTE: Assumes phi is ordered
pub fn try_mueller_to_1d(
    bins: &[(AngleBin, AngleBin)],
    mueller: &[Mueller],
) -> Result<(Vec<f32>, Vec<Mueller>)> {
    // Check that the mueller matrix and bins are the same length
    if mueller.len() != bins.len() {
        return Err(anyhow!(
            "Mueller matrix and bins must have the same length. Got {} and {}",
            mueller.len(),
            bins.len()
        ));
    }

    // // Create indices and sort them by corresponding theta values
    // let mueller = mueller.to_owned();
    // let mut bins = bins.to_owned();
    // let mut indices: Vec<usize> = (0..bins.len()).collect();
    // indices.sort_by(|&i, &j| bins[i].0.center.partial_cmp(&bins[j].0.center).unwrap());

    // // Sort bins according to the sorted indices
    // bins = indices.iter().map(|&i| bins[i]).collect();

    // // Create a new sorted mueller matrix using the same indices
    // let mut sorted_mueller = Array2::<f32>::zeros(mueller.dim());
    // for (new_idx, &old_idx) in indices.iter().enumerate() {
    //     sorted_mueller
    //         .slice_mut(s![new_idx, ..])
    //         .assign(&mueller.slice(s![old_idx, ..]));
    // }

    // // zip the bins and mueller matrix
    // let combined: Vec<_> = bins
    //     .iter()
    //     .zip(mueller.outer_iter())
    //     .map(|(bin, row)| (*bin, row.to_owned()))
    //     .collect();

    // // group the combined Vec by theta center
    // let grouped: Vec<Vec<_>> = combined
    //     .into_iter()
    //     .chunk_by(|((theta_bin, _), _)| theta_bin.center)
    //     .into_iter()
    //     .map(|(_, group)| group.map(|x| x).collect())
    //     .collect();

    // let mut thetas = Vec::new();
    // let mut mueller_1d = Array2::<f32>::zeros((grouped.len(), 16));

    // // loop over vectors at each theta
    // for (i, muellers) in grouped.iter().enumerate() {
    //     // Unzip the theta, phi, and mueller values
    //     let thetas_phi: Vec<_> = muellers
    //         .iter()
    //         .map(|((theta_bin, phi_bin), _)| (theta_bin.center, phi_bin.center))
    //         .collect();
    //     let mueller_phi: Vec<_> = muellers.iter().map(|(_, mueller)| mueller).collect();
    //     let mut mueller_1d_row = Array1::<f32>::zeros(16);

    //     // loop over the mueller values at each phi
    //     for j in 0..16 {
    //         // Create 1D arrays for x and y, where x is phi and y is 1 of the 16 mueller values
    //         let y = Array1::from(mueller_phi.iter().map(|row| row[j]).collect::<Vec<_>>());
    //         let phi_values: Vec<_> = thetas_phi.iter().map(|(_, phi)| phi.to_radians()).collect();

    //         // Check if phi values are sorted in ascending order (probably dont need this)
    //         for i in 1..phi_values.len() {
    //             if phi_values[i] < phi_values[i - 1] {
    //                 return Err(anyhow!(
    //                     "Phi values must be sorted in ascending order for integration"
    //                 ));
    //             }
    //         }

    //         let x = Array1::from(phi_values);
    //         mueller_1d_row[j] = output::integrate_trapezoidal(&x, &y, |_, y| y);
    //         // integrate over phi
    //     }

    //     // Assign the theta and mueller values to the final arrays
    //     thetas.push(thetas_phi[0].0);
    //     mueller_1d
    //         .slice_mut(s![i, ..])
    //         .assign(&mueller_1d_row.slice(s![..]));
    // }
    !todo!();

    // Ok((thetas, mueller_1d))
}

/// Integrate the first mueller element over theta to get the asymmetry parameter
pub fn compute_asymmetry(theta: &[f32], mueller_1d: &[Mueller], waveno: f32, scatt: f32) -> f32 {
    // get first column of mueller matrix
    let y = mueller_1d.into_iter().map(|m| m.s11).collect::<Vec<f32>>();

    let x = Array1::from(theta.iter().map(|&t| t.to_radians()).collect::<Vec<f32>>());

    // integrate p11 sin(theta) cos(theta) / scatt cross / k^2
    // let asymmetry = output::integrate_trapezoidal(&x, &y, |x, y| {
    //     x.sin() * x.cos() * y / scatt / waveno.powi(2)
    // });
    !todo!();

    // asymmetry
}

/// Integrate the first mueller element over theta to get the scattering cross section
pub fn compute_scat_cross(theta: &[f32], mueller_1d: &[Mueller], waveno: f32) -> f32 {
    // get first column of mueller matrix
    let y = mueller_1d.into_iter().map(|m| m.s11).collect::<Vec<f32>>();

    // convert theta values from degrees to radians for integration
    let x = Array1::from(theta.iter().map(|&t| t.to_radians()).collect::<Vec<f32>>());

    // integrate p11 sin(theta) / k^2
    // let scat_cross = output::integrate_trapezoidal(&x, &y, |x, y| x.sin() * y / waveno.powi(2));
    !todo!();

    // scat_cross
}
