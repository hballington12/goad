use crate::result::{ApproxEq, Mueller};
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
    fn field_partial_eq_with_tolerance() {
        let e_perp = Vector3::x();
        let prop = Vector3::z();

        // Create two nearly identical fields
        let field1 = Field::new_identity(e_perp, prop).unwrap();

        // Create a slightly different field
        let e_perp2 = Vector3::new(1.0 + 1e-6, 0.0, 0.0).normalize();
        let mut field2 = Field::new_identity(e_perp2, prop).unwrap();
        field2.phase = 1e-6;

        // Should be equal within tolerance (1e-5)
        assert_eq!(field1, field2);

        // Create a field with larger differences
        let e_perp3 = Vector3::new(1.0 + 1e-4, 0.0, 0.0).normalize();
        let mut field3 = Field::new_identity(e_perp3, prop).unwrap();
        field3.phase = 1e-4;

        // Should not be equal (difference exceeds tolerance)
        assert_ne!(field1, field3);
    }
}

pub type Ampl = Matrix2<Complex<f32>>;

pub trait AmplMatrix {
    fn s11(&self) -> f32;
    fn s12(&self) -> f32;
    fn s13(&self) -> f32;
    fn s14(&self) -> f32;
    fn s21(&self) -> f32;
    fn s22(&self) -> f32;
    fn s23(&self) -> f32;
    fn s24(&self) -> f32;
    fn s31(&self) -> f32;
    fn s32(&self) -> f32;
    fn s33(&self) -> f32;
    fn s34(&self) -> f32;
    fn s41(&self) -> f32;
    fn s42(&self) -> f32;
    fn s43(&self) -> f32;
    fn s44(&self) -> f32;
    fn to_mueller(&self) -> Mueller;
    fn zeros() -> Self;
    fn identity() -> Self;
    fn is_valid(&self) -> bool;
}

impl AmplMatrix for Ampl {
    fn s11(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj()
            + self[(0, 1)] * self[(0, 1)].conj()
            + self[(1, 0)] * self[(1, 0)].conj()
            + self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s12(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj() - self[(0, 1)] * self[(0, 1)].conj()
            + self[(1, 0)] * self[(1, 0)].conj()
            - self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s13(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() + self[(1, 1)] * self[(1, 0)].conj()).re
    }
    fn s14(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() - self[(1, 1)] * self[(1, 0)].conj()).im
    }
    fn s21(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj() + self[(0, 1)] * self[(0, 1)].conj()
            - self[(1, 0)] * self[(1, 0)].conj()
            - self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s22(&self) -> f32 {
        0.5 * (self[(0, 0)] * self[(0, 0)].conj()
            - self[(0, 1)] * self[(0, 1)].conj()
            - self[(1, 0)] * self[(1, 0)].conj()
            + self[(1, 1)] * self[(1, 1)].conj())
        .re
    }
    fn s23(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() - self[(1, 1)] * self[(1, 0)].conj()).re
    }
    fn s24(&self) -> f32 {
        (self[(0, 0)] * self[(0, 1)].conj() + self[(1, 1)] * self[(1, 0)].conj()).im
    }
    fn s31(&self) -> f32 {
        (self[(0, 0)] * self[(1, 0)].conj() + self[(1, 1)] * self[(0, 1)].conj()).re
    }
    fn s32(&self) -> f32 {
        (self[(0, 0)] * self[(1, 0)].conj() - self[(1, 1)] * self[(0, 1)].conj()).re
    }
    fn s33(&self) -> f32 {
        (self[(0, 0)] * self[(1, 1)].conj() + self[(0, 1)] * self[(1, 0)].conj()).re
    }
    fn s34(&self) -> f32 {
        (self[(0, 0)] * self[(1, 1)].conj() + self[(0, 1)] * self[(1, 0)].conj()).im
    }
    fn s41(&self) -> f32 {
        (self[(1, 0)] * self[(0, 0)].conj() + self[(1, 1)] * self[(0, 1)].conj()).im
    }
    fn s42(&self) -> f32 {
        (self[(1, 0)] * self[(0, 0)].conj() - self[(1, 1)] * self[(0, 1)].conj()).im
    }
    fn s43(&self) -> f32 {
        (self[(1, 1)] * self[(0, 0)].conj() - self[(0, 1)] * self[(1, 0)].conj()).im
    }
    fn s44(&self) -> f32 {
        (self[(1, 1)] * self[(0, 0)].conj() - self[(0, 1)] * self[(1, 0)].conj()).re
    }
    fn to_mueller(&self) -> Mueller {
        Mueller::new(
            self.s11(),
            self.s12(),
            self.s13(),
            self.s14(),
            self.s21(),
            self.s22(),
            self.s23(),
            self.s24(),
            self.s31(),
            self.s32(),
            self.s33(),
            self.s34(),
            self.s41(),
            self.s42(),
            self.s43(),
            self.s44(),
        )
    }
    fn zeros() -> Self {
        Self::zeros()
    }
    fn identity() -> Self {
        Self::identity()
    }
    fn is_valid(&self) -> bool {
        self.iter().all(|c| c.re.is_finite() && c.im.is_finite())
    }
}

/// Essentially represents a plane wave field with phase, amplitude matrix, and corresponding electric field vectors.
#[derive(Debug, Clone)]
pub struct Field {
    ampl: Ampl,
    e_perp: Vector3<f32>,
    phase: f32,
    prop: Vector3<f32>,
}

impl PartialEq for Field {
    /// Compare two Field structs with a tolerance of 1e-5
    fn eq(&self, other: &Self) -> bool {
        const TOLERANCE: f32 = 1e-5;

        // Check amplitude matrices using ApproxEq trait
        self.ampl.approx_eq(&other.ampl, TOLERANCE)
            // Check vectors component-wise with tolerance
            && (self.e_perp.x - other.e_perp.x).abs() < TOLERANCE
            && (self.e_perp.y - other.e_perp.y).abs() < TOLERANCE
            && (self.e_perp.z - other.e_perp.z).abs() < TOLERANCE
            && (self.prop.x - other.prop.x).abs() < TOLERANCE
            && (self.prop.y - other.prop.y).abs() < TOLERANCE
            && (self.prop.z - other.prop.z).abs() < TOLERANCE
            // Check phase with tolerance
            && (self.phase - other.phase).abs() < TOLERANCE
    }
}

impl Field {
    /// Multiply the field amplitude by a complex-valued rotation matrix
    pub fn matmul(&mut self, rot: &Matrix2<Complex<f32>>) {
        self.ampl = rot * self.ampl;
    }

    /// Multiply the field amplitude by a real-valued factor
    pub fn mul(&mut self, fac: f32) {
        self.ampl *= Complex::new(fac, 0.0);
    }

    /// Increment the phase by a real value
    pub fn wind(&mut self, arg: f32) {
        self.phase += arg;
    }

    // Getters

    /// Returns the propagation vector of the field
    pub fn prop(&self) -> Vector3<f32> {
        self.prop
    }

    /// Returns the amplitude of the field
    pub fn ampl(&self) -> Ampl {
        let phase = self.phase;
        self.ampl_wo_phase() * Complex::new(phase.cos(), phase.sin())
    }

    /// Returns the amplitude of the field without phase due to path length
    pub fn ampl_wo_phase(&self) -> Ampl {
        self.ampl
    }

    /// Returns the phase of the field
    pub fn phase(&self) -> f32 {
        self.phase
    }

    /// Returns the perpendicular component of the electric field
    pub fn e_perp(&self) -> Vector3<f32> {
        self.e_perp
    }

    /// Returns the parallel component of the electric field
    pub fn e_par(&self) -> Vector3<f32> {
        self.e_perp.cross(&self.prop).normalize()
    }

    // Setters

    pub fn set_ampl(&mut self, ampl0: Ampl) {
        self.ampl = ampl0;
    }

    pub fn set_phase(&mut self, phase: f32) {
        self.phase = phase;
    }

    pub fn set_e_perp(&mut self, e_perp: &Vector3<f32>) {
        self.e_perp = *e_perp;
    }

    pub fn set_prop(&mut self, prop: Vector3<f32>) {
        self.prop = prop;
    }

    /// Returns a rotation matrix for rotating from the current plane perpendicular to the plane perpendicular to `e_perp`.
    pub fn get_rotation_matrix(&self, e_perp: &Vector3<f32>) -> Matrix2<Complex<f32>> {
        Field::rotation_matrix(self.e_perp, *e_perp, self.prop)
            .map(|x| nalgebra::Complex::new(x, 0.0))
    }

    pub fn new_from_e_perp(&self, e_perp: &Vector3<f32>) -> Self {
        let rot = self.get_rotation_matrix(e_perp);
        let mut field = self.clone();
        field.set_e_perp(e_perp);
        field.matmul(&rot);
        field
    }

    /// Creates a new unit electric field with the given input perpendicular
    /// and propagation vectors.
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
            phase: 0.0,
            prop,
        };

        Ok(field)
    }

    /// Creates an electric field with the given input perpendicular field
    /// vector, propagation vector, and amplitude matrix.
    pub fn new(e_perp: Vector3<f32>, prop: Vector3<f32>, ampl0: Ampl, phase: f32) -> Result<Self> {
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
            ampl: ampl0,
            e_perp,
            phase,
            prop,
        };

        Ok(field)
    }

    /// Returns the 2x2 rotation matrix for rotating an amplitude matrix
    /// about the propagation vector `prop` from
    /// the plane perpendicular to `e_perp_in` to the plane perpendicular to
    /// `e_perp_out`.
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

    /// Returns the field intensity.
    pub fn intensity(&self) -> f32 {
        0.5 * self.ampl.norm_squared()
    }

    pub fn ampl_intensity(ampl: &Matrix2<Complex<f32>>) -> f32 {
        0.5 * ampl.norm_squared()
    }
}
