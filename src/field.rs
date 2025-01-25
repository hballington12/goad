use crate::config;
use anyhow::Result;
use core::fmt;
use std::fmt::{Debug, Error};

use nalgebra::{Complex, ComplexField, Matrix2, Vector3};

#[cfg(test)]
mod tests {

    use nalgebra::Complex;

    use super::*;

    #[test]
    fn identity_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let field = Field::new_identity(e_perp, prop).unwrap();
        assert!((field.intensity() - 1.0).abs() < config::COLINEAR_THRESHOLD);
    }

    #[test]
    fn null_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let mut field = Field::new_identity(e_perp, prop).unwrap();
        field.ampl *= Complex::ZERO;
        assert!((field.intensity()).abs() < config::COLINEAR_THRESHOLD);
    }
}

#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct Field {
    pub ampl: Matrix2<Complex<f32>>,
    pub e_perp: Vector3<f32>,
    pub e_par: Vector3<f32>,
}

impl Field {
    /// Creates a new unit electric field with the given input perpendicular
    /// and propagation vectors.
    pub fn new_identity(e_perp: Vector3<f32>, prop: Vector3<f32>) -> Result<Self> {
        #[cfg(debug_assertions)]
        let norm_e_perp_diff = e_perp.norm() - 1.0;
        if norm_e_perp_diff.abs() >= config::COLINEAR_THRESHOLD {
            return Err(anyhow::anyhow!("e-perp is not normalised: {:?}", e_perp));
        }

        let norm_prop_diff = prop.norm() - 1.0;
        if norm_prop_diff.abs() >= config::COLINEAR_THRESHOLD {
            return Err(anyhow::anyhow!(
                "propagation vector is not normalised: {:?}",
                prop
            ));
        }

        let dot_product = e_perp.dot(&prop);
        if dot_product.abs() >= config::COLINEAR_THRESHOLD {
            return Err(anyhow::anyhow!(
                "e-perp and propagation vector are not perpendicular, e_perp is: {:?}, prop is: {:?}, dot product is: {:?}",
                e_perp,
                prop,
                dot_product
            ));
        }

        let field = Self {
            ampl: Matrix2::identity(),
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        };

        Ok(field)
    }

    /// Creates an electric field with the given input perpendicular field
    /// vector, propagation vector, and amplitude matrix.
    pub fn new(
        e_perp: Vector3<f32>,
        prop: Vector3<f32>,
        ampl: Matrix2<Complex<f32>>,
    ) -> Result<Self> {
        #[cfg(debug_assertions)]
        let norm_e_perp_diff = e_perp.norm() - 1.0;
        if norm_e_perp_diff.abs() >= config::COLINEAR_THRESHOLD {
            return Err(anyhow::anyhow!("e-perp is not normalised: {:?}", e_perp));
        }

        let norm_prop_diff = prop.norm() - 1.0;
        if norm_prop_diff.abs() >= config::COLINEAR_THRESHOLD {
            return Err(anyhow::anyhow!(
                "propagation vector is not normalised: {:?}",
                prop
            ));
        }

        let dot_product = e_perp.dot(&prop);
        if dot_product.abs() >= config::COLINEAR_THRESHOLD {
            return Err(anyhow::anyhow!(
                "e-perp and propagation vector are not perpendicular, e_perp is: {:?}, prop is: {:?}, dot product is: {:?}",
                e_perp,
                prop,
                dot_product
            ));
        }

        let field = Self {
            ampl,
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        };

        Ok(field)
    }

    /// Returns the 2x2 rotation matrix for rotating an amplitude matrix
    /// about the propagation vector `prop` from
    /// the plane perpendicular to `e_perp_in` to the plane perpendicular to
    /// `e_perp_out`.
    pub fn rotation_matrix(
        e_perp_in: Vector3<f32>,
        e_perp_out: Vector3<f32>,
        prop: Vector3<f32>,
    ) -> Matrix2<f32> {
        let dot1 = e_perp_out.dot(&e_perp_in);
        let evo2 = prop.cross(&e_perp_in).normalize();
        let dot2 = e_perp_out.dot(&evo2);

        println!("dot1: {}", dot1);
        println!("evo2: {}", evo2);

        let result = Matrix2::new(dot1, -dot2, dot2, dot1);
        let det = result.determinant();

        result / det.abs().sqrt()
    }

    /// Returns the field intensity.
    pub fn intensity(&self) -> f32 {
        Self::ampl_intensity(&self.ampl)
    }

    pub fn ampl_intensity(ampl: &Matrix2<Complex<f32>>) -> f32 {
        0.5 * ampl.norm_squared()
    }
}
