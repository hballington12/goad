use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

use crate::zones::{Zone, ZoneType, Zones};

use super::mueller::{Mueller, MuellerMatrix};
use super::results::Results;

// Note: gen_stub_pymethods is not used here because the pymethods are in a separate
// file from the pyclass definition. The stub generation will capture the pyclass
// attributes from results.rs.
#[pymethods]
impl Results {
    fn __add__(&self, other: &Results) -> Results {
        self.clone() + other.clone()
    }

    fn __sub__(&self, other: &Results) -> Results {
        self.clone() - other.clone()
    }

    fn __pow__(&self, exponent: f32, _modulo: Option<u32>) -> Results {
        use rand_distr::num_traits::Pow;
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

    /// Get the zones collection
    #[getter]
    pub fn get_zones(&self) -> Zones {
        self.zones.clone()
    }

    /// Get a specific zone by label
    pub fn get_zone(&self, label: &str) -> Option<Zone> {
        self.zones.get(label).cloned()
    }

    /// Get a zone by type (returns first matching)
    pub fn get_zone_by_type(&self, zone_type: ZoneType) -> Option<Zone> {
        self.zones
            .iter()
            .find(|z| z.zone_type == zone_type)
            .cloned()
    }

    /// Get the full zone (convenience method)
    #[getter]
    pub fn get_full_zone(&self) -> Option<Zone> {
        self.zones.full_zone().cloned()
    }

    /// Get the forward zone (convenience method)
    #[getter]
    pub fn get_forward_zone(&self) -> Option<Zone> {
        self.zones
            .iter()
            .find(|z| z.zone_type == ZoneType::Forward)
            .cloned()
    }

    /// Get the backward zone (convenience method)
    #[getter]
    pub fn get_backward_zone(&self) -> Option<Zone> {
        self.zones.backward_zone().cloned()
    }

    /// Get the bins as a numpy array of shape (n_bins, 2) with columns [theta, phi]
    #[getter]
    pub fn get_bins<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let bins: Vec<f32> = self
            .bins()
            .iter()
            .flat_map(|bin| vec![bin.theta.center, bin.phi.center])
            .collect();

        Array2::from_shape_vec((bins.len() / 2, 2), bins)
            .unwrap()
            .into_pyarray(py)
    }

    /// Get the 1D bins (theta values) as a numpy array
    #[getter]
    pub fn get_bins_1d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        self.zones.full_zone().and_then(|zone| {
            zone.field_1d.as_ref().map(|field_1d| {
                let bins: Vec<f32> = field_1d.iter().map(|result| result.bin.center).collect();
                Array2::from_shape_vec((bins.len(), 1), bins)
                    .unwrap()
                    .into_pyarray(py)
            })
        })
    }

    /// Get the Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let field_2d = self
            .zones
            .full_zone()
            .map(|z| &z.field_2d[..])
            .unwrap_or(&[]);
        let muellers: Vec<f32> = field_2d
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
        let array_view = unsafe { array.as_array() };
        if let Some(zone) = self.zones.full_zone_mut() {
            for (i, field) in zone.field_2d.iter_mut().enumerate() {
                let row = array_view.row(i);
                let slice = row.as_slice().unwrap();
                field.mueller_total = Mueller::from_row_slice(slice);
            }
        }
    }

    /// Get the beam Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_beam<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let field_2d = self
            .zones
            .full_zone()
            .map(|z| &z.field_2d[..])
            .unwrap_or(&[]);
        let muellers: Vec<f32> = field_2d
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
        let array_view = unsafe { array.as_array() };
        if let Some(zone) = self.zones.full_zone_mut() {
            for (i, field) in zone.field_2d.iter_mut().enumerate() {
                let row = array_view.row(i);
                let slice = row.as_slice().unwrap();
                field.mueller_beam = Mueller::from_row_slice(slice);
            }
        }
    }

    /// Get the external diffraction Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_ext<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let field_2d = self
            .zones
            .full_zone()
            .map(|z| &z.field_2d[..])
            .unwrap_or(&[]);
        let muellers: Vec<f32> = field_2d
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
        let array_view = unsafe { array.as_array() };
        if let Some(zone) = self.zones.full_zone_mut() {
            for (i, field) in zone.field_2d.iter_mut().enumerate() {
                let row = array_view.row(i);
                let slice = row.as_slice().unwrap();
                field.mueller_ext = Mueller::from_row_slice(slice);
            }
        }
    }

    /// Get the 1D Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_1d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        self.zones.full_zone().and_then(|zone| {
            zone.field_1d.as_ref().map(|field_1d| {
                let muellers: Vec<f32> = field_1d
                    .iter()
                    .flat_map(|r| r.mueller_total.to_vec())
                    .collect();
                Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                    .unwrap()
                    .into_pyarray(py)
            })
        })
    }

    /// Set the 1D Mueller matrix from a numpy array
    #[setter]
    pub fn set_mueller_1d(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        let array_view = unsafe { array.as_array() };
        if let Some(zone) = self.zones.full_zone_mut() {
            if let Some(ref mut field_1d) = zone.field_1d {
                for (i, field) in field_1d.iter_mut().enumerate() {
                    let row = array_view.row(i);
                    let slice = row.as_slice().unwrap();
                    field.mueller_total = Mueller::from_row_slice(slice);
                }
            }
        }
    }

    /// Get the 1D beam Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_1d_beam<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        self.zones.full_zone().and_then(|zone| {
            zone.field_1d.as_ref().map(|field_1d| {
                let muellers: Vec<f32> = field_1d
                    .iter()
                    .flat_map(|r| r.mueller_beam.to_vec())
                    .collect();
                Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                    .unwrap()
                    .into_pyarray(py)
            })
        })
    }

    /// Set the 1D Mueller beam matrix from a numpy array
    #[setter]
    pub fn set_mueller_1d_beam(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        let array_view = unsafe { array.as_array() };
        if let Some(zone) = self.zones.full_zone_mut() {
            if let Some(ref mut field_1d) = zone.field_1d {
                for (i, field) in field_1d.iter_mut().enumerate() {
                    let row = array_view.row(i);
                    let slice = row.as_slice().unwrap();
                    field.mueller_beam = Mueller::from_row_slice(slice);
                }
            }
        }
    }

    /// Get the 1D external diffraction Mueller matrix as a numpy array
    #[getter]
    pub fn get_mueller_1d_ext<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f32>>> {
        self.zones.full_zone().and_then(|zone| {
            zone.field_1d.as_ref().map(|field_1d| {
                let muellers: Vec<f32> = field_1d
                    .iter()
                    .flat_map(|r| r.mueller_ext.to_vec())
                    .collect();
                Array2::from_shape_vec((muellers.len() / 16, 16), muellers)
                    .unwrap()
                    .into_pyarray(py)
            })
        })
    }

    /// Set the 1D Mueller ext matrix from a numpy array
    #[setter]
    pub fn set_mueller_1d_ext(&mut self, array: &Bound<'_, PyArray2<f32>>) {
        let array_view = unsafe { array.as_array() };
        if let Some(zone) = self.zones.full_zone_mut() {
            if let Some(ref mut field_1d) = zone.field_1d {
                for (i, field) in field_1d.iter_mut().enumerate() {
                    let row = array_view.row(i);
                    let slice = row.as_slice().unwrap();
                    field.mueller_ext = Mueller::from_row_slice(slice);
                }
            }
        }
    }

    /// Get the asymmetry parameter
    #[getter]
    pub fn get_asymmetry(&self) -> Option<f32> {
        use super::component::GOComponent;
        self.params.asymmetry(&GOComponent::Total)
    }

    /// Get the scattering cross section
    #[getter]
    pub fn get_scat_cross(&self) -> Option<f32> {
        use super::component::GOComponent;
        self.params.scatt_cross(&GOComponent::Total)
    }

    /// Get the extinction cross section
    #[getter]
    pub fn get_ext_cross(&self) -> Option<f32> {
        use super::component::GOComponent;
        self.params.ext_cross(&GOComponent::Total)
    }

    /// Get the albedo
    #[getter]
    pub fn get_albedo(&self) -> Option<f32> {
        use super::component::GOComponent;
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
