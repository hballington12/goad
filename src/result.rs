use crate::{output, powers::Powers};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2};
use ndarray::Array2;

#[derive(Debug, PartialEq, Clone)]
pub struct Result {
    pub powers: Powers,
    pub bins: Vec<(f32, f32)>,
    pub mueller: Array2<f32>,
    pub ampl: Vec<Matrix2<Complex<f32>>>,
    pub bins_1d: Option<Vec<f32>>,
    pub mueller_1d: Option<Array2<f32>>,
}

impl Result {
    /// Creates a new `Result` with empty mueller and amplitude matrix
    pub fn new_empty(bins: Vec<(f32, f32)>) -> Self {
        let mueller = Array2::<f32>::zeros((bins.len(), 16));
        let ampl = vec![Matrix2::<Complex<f32>>::zeros(); bins.len()];
        Self {
            powers: Powers::new(),
            bins,
            mueller,
            ampl,
            bins_1d: None,
            mueller_1d: None,
        }
    }

    pub fn try_mueller_to_1d(&mut self) -> std::result::Result<(), anyhow::Error> {
        match output::try_mueller_to_1d(&self.bins, &self.mueller) {
            Ok((theta, mueller_1d)) => {
                self.bins_1d = Some(theta);
                self.mueller_1d = Some(mueller_1d);

                Ok(())
            }
            Err(e) => Err(e),
        }
    }
}
