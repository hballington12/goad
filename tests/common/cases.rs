//! Mueller-matrix regression case definitions.
//!
//! This file is intentionally just a list of `mueller_case! { ... }` invocations
//! so it can be re-included verbatim by both:
//!   - the `mueller_regression_tests` integration test target (which produces a
//!     `#[test]` per case), and
//!   - the `regen_mueller_refs` example binary (which iterates the same registry
//!     to rewrite reference files when the physics intentionally changes).
//!
//! Adding a new case = one new `mueller_case! { ... }` block.

use goad::{
    bins,
    orientation::{Euler, EulerConvention, Orientation, Scheme},
    zones::ZoneConfig,
};
use num_complex::Complex32;

mueller_case! {
    name: fixed_hex_30_30_30_mueller_scatgrid,
    settings: |s| {
        s.zones = vec![ZoneConfig::new(bins::Scheme::new_simple(19, 19))];
        s.orientation = Orientation {
            scheme: Scheme::Discrete { eulers: vec![Euler::new(30.0, 30.0, 30.0)] },
            euler_convention: EulerConvention::ZYZ,
        };
    }
}

mueller_case! {
    name: fixed_hex_30_20_20_mueller_scatgrid,
    settings: |s| {
        s.zones = vec![ZoneConfig::new(bins::Scheme::new_simple(19, 19))];
        s.orientation = Orientation {
            scheme: Scheme::Discrete { eulers: vec![Euler::new(30.0, 20.0, 20.0)] },
            euler_convention: EulerConvention::ZYZ,
        };
    },
    refr_index: vec![Complex32::new(1.3117, 0.1)],
}
