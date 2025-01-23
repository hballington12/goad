use nalgebra::{Complex, ComplexField, Matrix2, Vector3};

#[cfg(test)]
mod tests {

    use nalgebra::Complex;

    use super::*;

    #[test]
    fn identity_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let field = Field::new_identity(e_perp, prop);
        assert!((field.intensity() - 1.0).abs() < 0.01);
    }

    #[test]
    fn null_ampl() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();
        let mut field = Field::new_identity(e_perp, prop);
        field.ampl *= Complex::ZERO;
        assert!((field.intensity()).abs() < 0.01);
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
    pub fn new_identity(e_perp: Vector3<f32>, prop: Vector3<f32>) -> Self {
        #[cfg(debug_assertions)]
        debug_assert!(
            (e_perp.norm() - 1.0).abs() < 0.01,
            "e-perp is not normalised: {}",
            e_perp
        );
        debug_assert!(
            (prop.norm() - 1.0).abs() < 0.01,
            "propagation vector is not normalised: {}",
            prop
        );
        debug_assert!(
            e_perp.dot(&prop) < 0.01,
            "e-perp and propagation vector are not perpendicular, dot product is: {}",
            e_perp.dot(&prop)
        );

        let field = Self {
            ampl: Matrix2::identity(),
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        };

        field
    }

    /// Creates an electric field with the given input perpendicular field
    /// vector, propagation vector, and amplitude matrix.
    pub fn new(e_perp: Vector3<f32>, prop: Vector3<f32>, ampl: Matrix2<Complex<f32>>) -> Self {
        #[cfg(debug_assertions)]
        debug_assert!(
            (e_perp.norm() - 1.0).abs() < 0.01,
            "e-perp is not normalised: {}",
            e_perp
        );
        debug_assert!(
            (prop.norm() - 1.0).abs() < 0.01,
            "propagation vector is not normalised: {}",
            prop
        );
        debug_assert!(
            e_perp.dot(&prop) < 0.01,
            "e-perp and propagation vector are not perpendicular, dot product is: {}",
            e_perp.dot(&prop)
        );

        let field = Self {
            ampl,
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        };

        field
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
