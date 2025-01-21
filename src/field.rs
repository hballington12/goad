use nalgebra::{Complex, Matrix2, Vector3};

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
            "e-perp is not normalised"
        );
        debug_assert!(
            (prop.norm() - 1.0).abs() < 0.01,
            "propagation vector is not normalised"
        );
        debug_assert!(
            e_perp.dot(&prop) < 0.01,
            "e-perp and propagation vector are not perpendicular"
        );

        Self {
            ampl: Matrix2::identity(),
            e_perp,
            e_par: e_perp.cross(&prop).normalize(),
        }
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

        Matrix2::new(dot1, -dot2, dot2, dot1).normalize()
    }
}
