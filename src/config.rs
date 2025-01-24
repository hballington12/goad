use std::f32::consts::PI;

use nalgebra::Complex;

pub const BEAM_AREA_THRESHOLD: f32 = 0.01; // minimum area for new beam to propagate
pub const BEAM_POWER_THRESHOLD: f32 = 0.0001; // minimum power for new beam to propagate
pub const CLIP_TOLERANCE: f32 = 1e6; // Named constant for tolerance
pub const COLINEAR_THRESHOLD: f32 = 0.001; // Named constant for avoiding taking cross products from colinear vectors
pub const MEDIUM_REFR_INDEX: Complex<f32> = Complex {
    // outer medium refractive index
    re: 1.0,
    im: 0.0,
};
pub const MAX_REC: i32 = 10; // maximum number of beam recursions before truncation
pub const MAX_TIR: i32 = 10; // maximum number of total internal reflections
pub const WAVENO: f32 = 2.0 * PI / 0.532; // wavenumber
