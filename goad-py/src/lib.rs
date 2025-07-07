use goad::{
    self,
    bins::BinningScheme,
    geom::Geom,
    geom::Shape,
    multiproblem::MultiProblem,
    orientation::{Euler, EulerConvention, Orientation, Scheme},
    problem::Problem,
    result::Results,
    settings::Settings,
};
use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// Create a uniform orientation scheme with specified number of orientations
#[pyfunction]
fn uniform_orientation(num_orients: usize) -> PyResult<Scheme> {
    Ok(Scheme::Uniform { num_orients })
}

/// Create a discrete orientation scheme from a list of Euler angles
#[pyfunction]
fn discrete_orientation(eulers: Vec<Euler>) -> PyResult<Scheme> {
    Ok(Scheme::Discrete { eulers })
}

/// Create an Orientation with uniform scheme and default convention
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

/// A Python module implemented in Rust.
#[pymodule]
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

    // Helper functions for orientations
    m.add_function(wrap_pyfunction!(uniform_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(discrete_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(create_uniform_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(create_discrete_orientation, m)?)?;

    Ok(())
}
