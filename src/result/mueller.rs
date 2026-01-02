use nalgebra::Matrix4;
use nalgebra::{Complex, Matrix2};

pub type Ampl = Matrix2<Complex<f32>>;
pub type Mueller = Matrix4<f32>;

/// Trait for approximate equality comparison with tolerance
pub trait ApproxEq {
    /// Check if two values are approximately equal within the given tolerance
    fn approx_eq(&self, other: &Self, tolerance: f32) -> bool;
}

impl ApproxEq for Ampl {
    /// Check if two amplitude matrices are approximately equal within tolerance
    ///
    /// Compares both real and imaginary parts of each complex element.
    /// Returns true if all corresponding elements differ by less than the tolerance.
    fn approx_eq(&self, other: &Self, tolerance: f32) -> bool {
        for i in 0..2 {
            for j in 0..2 {
                let a = self[(i, j)];
                let b = other[(i, j)];

                // Check both real and imaginary parts
                if (a.re - b.re).abs() > tolerance || (a.im - b.im).abs() > tolerance {
                    return false;
                }
            }
        }
        true
    }
}

pub trait MuellerMatrix {
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
    fn to_vec(&self) -> Vec<f32>;
}

impl MuellerMatrix for Mueller {
    fn s11(&self) -> f32 {
        self[(0, 0)]
    }
    fn s12(&self) -> f32 {
        self[(0, 1)]
    }
    fn s13(&self) -> f32 {
        self[(0, 2)]
    }
    fn s14(&self) -> f32 {
        self[(0, 3)]
    }
    fn s21(&self) -> f32 {
        self[(1, 0)]
    }
    fn s22(&self) -> f32 {
        self[(1, 1)]
    }
    fn s23(&self) -> f32 {
        self[(1, 2)]
    }
    fn s24(&self) -> f32 {
        self[(1, 3)]
    }
    fn s31(&self) -> f32 {
        self[(2, 0)]
    }
    fn s32(&self) -> f32 {
        self[(2, 1)]
    }
    fn s33(&self) -> f32 {
        self[(2, 2)]
    }
    fn s34(&self) -> f32 {
        self[(2, 3)]
    }
    fn s41(&self) -> f32 {
        self[(3, 0)]
    }
    fn s42(&self) -> f32 {
        self[(3, 1)]
    }
    fn s43(&self) -> f32 {
        self[(3, 2)]
    }
    fn s44(&self) -> f32 {
        self[(3, 3)]
    }
    /// Returns the Mueller matrix as a vector of its elements.
    fn to_vec(&self) -> Vec<f32> {
        vec![
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
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Complex;

    #[test]
    fn test_ampl_approx_eq() {
        // Create two similar amplitude matrices
        let ampl1 = Ampl::new(
            Complex::new(1.0, 2.0),
            Complex::new(3.0, 4.0),
            Complex::new(5.0, 6.0),
            Complex::new(7.0, 8.0),
        );

        let ampl2 = Ampl::new(
            Complex::new(1.001, 2.001),
            Complex::new(3.001, 4.001),
            Complex::new(5.001, 6.001),
            Complex::new(7.001, 8.001),
        );

        // Should be equal within tolerance
        assert!(ampl1.approx_eq(&ampl2, 0.01));

        // Should not be equal with stricter tolerance
        assert!(!ampl1.approx_eq(&ampl2, 0.0001));

        // Test with exact equality
        assert!(ampl1.approx_eq(&ampl1, 0.0));
    }
}
