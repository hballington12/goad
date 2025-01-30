use nalgebra::Complex;
use std::f32::consts::PI;

/// Wavelength of the electric field in the same units as the geometry coordinates.
pub const WAVELENGTH: f32 = 0.532;
/// Wavenumber, equal to 2π divided by `WAVELENGTH`.
pub const WAVENUMBER: f32 = 2.0 * PI / WAVELENGTH;
/// Minimum distance for vertices to be considered the same.
pub const VERTEX_MERGE_DISTANCE: f32 = 0.1;
/// Minimum area for new beam to propagate, set by default to `WAVELENGTH`^2 / 4.
pub const BEAM_AREA_THRESHOLD: f32 = WAVELENGTH * WAVELENGTH / 4.0;
/// Minimum power for new beam to propagate.
pub const BEAM_POWER_THRESHOLD: f32 = 0.005;
/// Scaling factor for integer coordinates during clipping.
pub const CLIP_TOLERANCE: f32 = 1e8;
/// Minimum absolute value of the dot product of two vectors to be considered colinear.
pub const COLINEAR_THRESHOLD: f32 = 0.001;
/// Minimum vector length (in geometry units) to be considered non-degenerate.
pub const VEC_LENGTH_THRESHOLD: f32 = 0.01;
/// Minimum distance traversed by ray to intersection. Intersections closer than this are ignored.
pub const RAYCAST_MINIMUM_DISTANCE: f32 = 0.01;
/// Acceptable energy threshold beyond which beam propagation terminates.
pub const TOTAL_POWER_CUTOFF: f32 = 0.99;
/// Surrounding medium refractive index.
pub const MEDIUM_REFR_INDEX: Complex<f32> = Complex { re: 1.0, im: 0.0 };
/// Maximum number of beam recursions before truncation.
pub const MAX_REC: i32 = 5;
/// Maximum number of total internal reflections.
pub const MAX_TIR: i32 = 8;
/// Distance to far-field. 1e3 - 1e5 is a good range for single precision arithmetic.
pub const RADIUS: f32 = 1e4;
/// Tolerance for diffraction computations, used to avoid divide by zero errors.
pub const DIFF_EPSILON: f32 = 1e-3;
