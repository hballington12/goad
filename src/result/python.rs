//! Python bindings for Results.

use pyo3::prelude::*;

use super::component::GOComponent;
use super::mueller::Mueller;
use super::results::Results;

#[pymethods]
impl Results {
    /// Get the bins as a list of tuples (returns bin centers for backwards compatibility).
    #[getter]
    pub fn get_bins(&self) -> Vec<(f32, f32)> {
        self.bins()
            .iter()
            .map(|bin| (bin.theta_bin.center, bin.phi_bin.center))
            .collect()
    }

    /// Get the 1D bins (theta values).
    #[getter]
    pub fn get_bins_1d(&self) -> Option<Vec<f32>> {
        self.field_1d
            .as_ref()
            .map(|field_1d| field_1d.iter().map(|result| result.bin.center).collect())
    }

    /// Get the Mueller matrix as a list of lists.
    #[getter]
    pub fn get_mueller(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().map(|r| r.mueller_total).collect();
        collect_mueller(&muellers)
    }

    /// Get the beam Mueller matrix as a list of lists.
    #[getter]
    pub fn get_mueller_beam(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().map(|r| r.mueller_beam).collect();
        collect_mueller(&muellers)
    }

    /// Get the external diffraction Mueller matrix as a list of lists.
    #[getter]
    pub fn get_mueller_ext(&self) -> Vec<Vec<f32>> {
        let muellers: Vec<Mueller> = self.field_2d.iter().map(|r| r.mueller_ext).collect();
        collect_mueller(&muellers)
    }

    /// Get the 1D Mueller matrix as a list of lists.
    #[getter]
    pub fn get_mueller_1d(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().map(|r| r.mueller_total).collect();
            collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D beam Mueller matrix as a list of lists.
    #[getter]
    pub fn get_mueller_1d_beam(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().map(|r| r.mueller_beam).collect();
            collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the 1D external diffraction Mueller matrix as a list of lists.
    #[getter]
    pub fn get_mueller_1d_ext(&self) -> Vec<Vec<f32>> {
        if let Some(ref field_1d) = self.field_1d {
            let muellers: Vec<Mueller> = field_1d.iter().map(|r| r.mueller_ext).collect();
            collect_mueller(&muellers)
        } else {
            Vec::new()
        }
    }

    /// Get the asymmetry parameter.
    #[getter]
    pub fn get_asymmetry(&self) -> Option<f32> {
        self.params.asymmetry()
    }

    /// Get the scattering cross section.
    #[getter]
    pub fn get_scat_cross(&self) -> Option<f32> {
        self.params.scat_cross()
    }

    /// Get the extinction cross section.
    #[getter]
    pub fn get_ext_cross(&self) -> Option<f32> {
        self.params.ext_cross()
    }

    /// Get the albedo.
    #[getter]
    pub fn get_albedo(&self) -> Option<f32> {
        self.params.albedo()
    }

    /// Get all parameters as a dictionary.
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

    /// Get the powers as a dictionary.
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

/// Helper function to collect Mueller matrices into a list of lists.
pub fn collect_mueller(muellers: &[Mueller]) -> Vec<Vec<f32>> {
    use super::mueller::MuellerMatrix;
    muellers.iter().map(|m| m.to_vec()).collect()
}
