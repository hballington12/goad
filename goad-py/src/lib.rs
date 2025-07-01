use goad::{self, geom::Geom, geom::Shape, problem::Problem, result::Results, settings::Settings};
use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
#[pyfunction]
fn goad_py_add() -> PyResult<()> {
    let input = 1;
    let output = goad::helpers::add_one(input);
    println!("{} + 1 = {}", input, output);
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn goad_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(goad_py_add, m)?)?;
    m.add_class::<Shape>()?;
    m.add_class::<Geom>()?;
    m.add_class::<Settings>()?;
    m.add_class::<Problem>()?;
    m.add_class::<Results>()?;
    Ok(())
}
