#[cfg(debug_assertions)]
use crate::settings;

use anyhow::Result;
use std::fmt::Debug;

use nalgebra::{Complex, Matrix2, RealField, Vector3};

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn identity_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let field = Field::new_identity(e_perp, prop).unwrap();
        assert!((field.intensity() - 1.0).abs() < settings::COLINEAR_THRESHOLD);
    }

    #[test]
    fn null_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let mut field = Field::new_identity(e_perp, prop).unwrap();
        field.ampl *= Complex::ZERO;
        assert!((field.intensity()).abs() < settings::COLINEAR_THRESHOLD);
    }
}

/// Electromagnetic field representation for polarized wave propagation.
/// 
/// **Context**: Electromagnetic waves in scattering simulations require tracking
/// both field amplitude and polarization state. The field must be decomposed
/// into perpendicular and parallel components relative to the scattering plane
/// for proper application of Fresnel equations and Mueller matrix calculations.
/// 
/// **How it Works**: Uses a 2x2 complex amplitude matrix to represent the field
/// in a coordinate system defined by perpendicular (e_perp) and parallel (e_par)
/// polarization vectors. The amplitude matrix encodes both field strength and
/// phase relationships between polarization components. The coordinate system
/// maintains orthogonality with the propagation direction.
#[derive(Debug, Clone, PartialEq)] // Added Default derive
pub struct Field {
    pub ampl: Matrix2<Complex<f32>>,
    pub e_perp: Vector3<f32>,
    pub e_par: Vector3<f32>,
}

impl Field {
    /// Creates a unit electromagnetic field with identity amplitude matrix.
    /// 
    /// **Context**: Initial incident beams in scattering simulations typically
    /// start with unit field strength and identity amplitude matrices before
    /// surface interactions modify the field through reflection and refraction.
    /// 
    /// **How it Works**: Establishes orthogonal polarization basis vectors from
    /// the perpendicular field vector and propagation direction. Validates vector
    /// normalization and orthogonality in debug mode, then creates field with
    /// identity amplitude matrix representing unit field strength.
    /// 
    /// # Example
    /// ```rust
    /// let field = Field::new_identity(e_perp, prop)?;
    /// ```
    pub fn new_identity(e_perp: Vector3<f32>, prop: Vector3<f32>) -> Result<Self> {
        #[cfg(debug_assertions)]
        {
            let norm_e_perp_diff = e_perp.norm() - 1.0;
            if norm_e_perp_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!("e-perp is not normalised: {:?}", e_perp));
            }

            let norm_prop_diff = prop.norm() - 1.0;
            if norm_prop_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                    "propagation vector is not normalised: {:?}",
                    prop
                ));
            }

            let dot_product = e_perp.dot(&prop);
            if dot_product.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                "e-perp and propagation vector are not perpendicular, e_perp is: {:?}, prop is: {:?}, dot product is: {:?}",
                e_perp,
                prop,
                dot_product
            ));
            }
        }

        let field = Self {
            ampl: Matrix2::identity(),
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        };

        Ok(field)
    }

    /// Creates an electromagnetic field with specified amplitude matrix.
    /// 
    /// **Context**: Surface interactions and beam propagation modify field
    /// amplitudes through Fresnel coefficients, phase accumulation, and
    /// coordinate transformations. This constructor allows creation of
    /// fields with arbitrary amplitude matrices.
    /// 
    /// **How it Works**: Similar to new_identity but accepts a pre-computed
    /// amplitude matrix. Validates vector properties and establishes the
    /// orthogonal polarization coordinate system.
    /// 
    /// # Example
    /// ```rust
    /// Field::new(e_perp, prop, refl_ampl)?
    /// ```
    pub fn new(
        e_perp: Vector3<f32>,
        prop: Vector3<f32>,
        ampl: Matrix2<Complex<f32>>,
    ) -> Result<Self> {
        #[cfg(debug_assertions)]
        {
            let norm_e_perp_diff = e_perp.norm() - 1.0;
            if norm_e_perp_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!("e-perp is not normalised: {:?}", e_perp));
            }

            let norm_prop_diff = prop.norm() - 1.0;
            if norm_prop_diff.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                    "propagation vector is not normalised: {:?}",
                    prop
                ));
            }

            let dot_product = e_perp.dot(&prop);
            if dot_product.abs() >= settings::COLINEAR_THRESHOLD {
                return Err(anyhow::anyhow!(
                "e-perp and propagation vector are not perpendicular, e_perp is: {:?}, prop is: {:?}, dot product is: {:?}",
                e_perp,
                prop,
                dot_product
            ));
            }
        }

        let field = Self {
            ampl,
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        };

        Ok(field)
    }

    /// Computes rotation matrix for polarization coordinate transformations.
    /// 
    /// **Context**: When electromagnetic fields interact with surfaces at
    /// different orientations, the polarization coordinate system changes.
    /// Field amplitudes must be rotated to maintain physical consistency
    /// between incident and scattered coordinate systems.
    /// 
    /// **How it Works**: Constructs rotation matrix from dot products between
    /// input and output perpendicular vectors and their cross products with
    /// the propagation direction. Normalizes by determinant to ensure proper
    /// orthogonal transformation.
    /// 
    /// # Example
    /// ```rust
    /// let rot = Field::rotation_matrix(beam.field.e_perp, e_perp, beam.prop)
    ///     .map(|x| nalgebra::Complex::new(x, 0.0));
    /// ```
    pub fn rotation_matrix<T: RealField + std::marker::Copy>(
        e_perp_in: Vector3<T>,
        e_perp_out: Vector3<T>,
        prop: Vector3<T>,
    ) -> Matrix2<T> {
        let dot1 = e_perp_out.dot(&e_perp_in);
        let evo2 = prop.cross(&e_perp_in).normalize();
        let dot2 = e_perp_out.dot(&evo2);

        let result = Matrix2::new(dot1, -dot2, dot2.clone(), dot1.clone());
        let det = result.determinant();

        result / det.abs().sqrt()
    }

    /// Returns the electromagnetic field intensity.
    /// 
    /// **Context**: Field intensity determines beam power when multiplied by
    /// cross-sectional area and refractive index. Intensity calculations are
    /// essential for power conservation verification and threshold-based
    /// beam termination.
    /// 
    /// **How it Works**: Calls the static amplitude intensity method.
    pub fn intensity(&self) -> f32 {
        Self::ampl_intensity(&self.ampl)
    }

    /// Computes intensity from a complex amplitude matrix.
    /// 
    /// **Context**: Field intensity represents the time-averaged electromagnetic
    /// energy density. For complex amplitude matrices, this requires computing
    /// the squared magnitude of all matrix elements.
    /// 
    /// **How it Works**: Uses the Frobenius norm squared of the amplitude matrix
    /// scaled by 0.5 to get the proper intensity in electromagnetic units.
    pub fn ampl_intensity(ampl: &Matrix2<Complex<f32>>) -> f32 {
        0.5 * ampl.norm_squared()
    }
}
