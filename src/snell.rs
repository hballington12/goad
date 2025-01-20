use crate::geom::RefrIndex;
use std::f32::consts::PI;

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn normal_incidence_same_media() {
        let theta_i = 0.0;
        let m1 = RefrIndex {
            real: 1.0,
            imag: 0.0,
        };
        let m2 = m1;
        let theta_t = get_theta_t(theta_i, m1, m2);
        assert!(theta_i - theta_t < 0.01)
    }

    #[test]
    fn normal_incidence() {
        let theta_i = 0.0;
        let m1 = RefrIndex {
            real: 1.0,
            imag: 0.0,
        };
        let m2 = RefrIndex {
            real: 1.31,
            imag: 0.0,
        };
        let theta_t = get_theta_t(theta_i, m1, m2);
        let abs_difference = (theta_i - theta_t).abs();
        assert!(abs_difference < f32::EPSILON)
    }

    #[test]
    fn angle30_incidence() {
        let theta_i = 30.0 * PI / 180.0;
        let m1 = RefrIndex {
            real: 1.0,
            imag: 0.0,
        };
        let m2 = RefrIndex {
            real: 1.31,
            imag: 0.0,
        };
        let theta_t = get_theta_t(theta_i, m1, m2);
        let abs_difference = (theta_t - 0.3916126).abs();
        assert!(abs_difference < 0.001)
    }

    #[test]
    fn absorbing_test() {
        let theta_i = 1.17773;
        let m1 = RefrIndex {
            real: 1.0,
            imag: 0.0,
        };
        let m2 = RefrIndex {
            real: 1.5,
            imag: 0.1,
        };
        let theta_t = get_theta_t(theta_i, m1, m2);
        let abs_difference = (theta_t - 0.662387).abs();
        assert!(abs_difference < 0.001)
    }
}

/// Returns the transmitted angle according to Snell's Law.
/// Port from Fortran code rt_c.f90, Macke 1996.
/// All angles are in radians.
fn get_theta_t(theta_i: f32, m1: RefrIndex, m2: RefrIndex) -> f32 {
    let k1 = m1.imag / m1.real; // imag(inc) / real(inc)
    let k2 = m2.imag / m2.real; // imag(trans) / real(trans)
    let krel = (k2 - k1) / (1.0 + k1 * k2);
    let nrel = m2.real / m1.real * (1.0 + k1 * k2) / (1.0 + k1 * k1);

    let ref1 = nrel * nrel;
    let ref2 = krel * krel;
    let ref3 = (1.0 + ref2) * (1.0 + ref2);
    let ref6 = ref1 * ref3 / ((1.0 + krel * k2) * (1.0 + krel * k2));

    let sintiq = (theta_i).sin().powi(2); // sin²(theta_i)
    let ref4 = 1.0 - (1.0 - ref2) / ref1 / ref3 * sintiq;
    let ref5 = 2.0 * krel / ref1 / ref3 * sintiq;

    let q4 = ref4 * ref4 + ref5 * ref5;
    let q2 = q4.sqrt();

    let test1 = (ref4 / q2).acos() / 2.0;

    let g = test1; // Picking g from test1 (same as Fortran code logic)

    let ref7 = (g.cos() - k2 * g.sin()).powi(2);
    let rnstar = (sintiq + ref6 * q2 * ref7).sqrt();

    (theta_i.sin() / rnstar).asin()
}
