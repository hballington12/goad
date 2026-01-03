use crate::{
    bins::BinningScheme,
    convergence::Convergence,
    diff::Mapping,
    geom::Geom,
    geom::Shape,
    multiproblem::MultiProblem,
    orientation::{Euler, EulerConvention, Orientation, Scheme},
    problem::Problem,
    result::Results,
    settings::Settings,
    zones::{Zone, ZoneConfig, ZoneType, Zones, ZonesIterator},
};
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::{define_stub_info_gatherer, derive::*};

/// Formats the sum of two numbers as string.
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// Create a uniform orientation scheme with specified number of orientations
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
fn uniform_orientation(num_orients: usize) -> PyResult<Scheme> {
    Ok(Scheme::Uniform { num_orients })
}

/// Create a discrete orientation scheme from a list of Euler angles
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
fn discrete_orientation(eulers: Vec<Euler>) -> PyResult<Scheme> {
    Ok(Scheme::Discrete { eulers })
}

/// Create a Sobol quasi-random orientation scheme (faster convergence)
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
fn sobol_orientation(num_orients: usize) -> PyResult<Scheme> {
    Ok(Scheme::Sobol { num_orients })
}

/// Create a Halton quasi-random orientation scheme (faster convergence)
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
fn halton_orientation(num_orients: usize) -> PyResult<Scheme> {
    Ok(Scheme::Halton { num_orients })
}

/// Create an Orientation with uniform scheme and default convention
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
#[pyo3(signature = (num_orients, euler_convention = None))]
fn create_uniform_orientation(
    num_orients: usize,
    euler_convention: Option<EulerConvention>,
) -> PyResult<Orientation> {
    Ok(Orientation {
        scheme: Scheme::Uniform { num_orients },
        euler_convention: euler_convention.unwrap_or(EulerConvention::ZYZ),
    })
}

/// Create an Orientation with discrete scheme and default convention
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
#[pyo3(signature = (eulers, euler_convention = None))]
fn create_discrete_orientation(
    eulers: Vec<Euler>,
    euler_convention: Option<EulerConvention>,
) -> PyResult<Orientation> {
    Ok(Orientation {
        scheme: Scheme::Discrete { eulers },
        euler_convention: euler_convention.unwrap_or(EulerConvention::ZYZ),
    })
}

/// Create an Orientation with Sobol quasi-random scheme (faster convergence)
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
#[pyo3(signature = (num_orients, euler_convention = None))]
fn create_sobol_orientation(
    num_orients: usize,
    euler_convention: Option<EulerConvention>,
) -> PyResult<Orientation> {
    Ok(Orientation {
        scheme: Scheme::Sobol { num_orients },
        euler_convention: euler_convention.unwrap_or(EulerConvention::ZYZ),
    })
}

/// Create an Orientation with Halton quasi-random scheme (faster convergence)
#[cfg_attr(feature = "stub-gen", gen_stub_pyfunction(module = "goad._goad"))]
#[pyfunction]
#[pyo3(signature = (num_orients, euler_convention = None))]
fn create_halton_orientation(
    num_orients: usize,
    euler_convention: Option<EulerConvention>,
) -> PyResult<Orientation> {
    Ok(Orientation {
        scheme: Scheme::Halton { num_orients },
        euler_convention: euler_convention.unwrap_or(EulerConvention::ZYZ),
    })
}

// Gather stub info from all annotated items
#[cfg(feature = "stub-gen")]
define_stub_info_gatherer!(stub_info);

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "_goad")]
fn _goad_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;

    // Core classes
    m.add_class::<Shape>()?;
    m.add_class::<Geom>()?;
    m.add_class::<Settings>()?;
    m.add_class::<Problem>()?;
    m.add_class::<MultiProblem>()?;
    m.add_class::<Results>()?;
    m.add_class::<BinningScheme>()?;

    // Orientation classes
    m.add_class::<Euler>()?;
    m.add_class::<EulerConvention>()?;
    m.add_class::<Orientation>()?;
    m.add_class::<Scheme>()?;

    // Mapping enum
    m.add_class::<Mapping>()?;

    // Param enum (for convergence targets)
    m.add_class::<crate::params::Param>()?;

    // Convergence solver
    m.add_class::<Convergence>()?;

    // Zone classes
    m.add_class::<ZoneType>()?;
    m.add_class::<ZoneConfig>()?;
    m.add_class::<Zone>()?;
    m.add_class::<Zones>()?;
    m.add_class::<ZonesIterator>()?;

    // Helper functions for orientations
    m.add_function(wrap_pyfunction!(uniform_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(discrete_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(sobol_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(halton_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(create_uniform_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(create_discrete_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(create_sobol_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(create_halton_orientation, m)?)?;

    Ok(())
}
