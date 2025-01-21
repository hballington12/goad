use crate::geom::RefrIndex;
pub const BEAM_AREA_THRESHOLD: f32 = 0.01; // minimum area for new beam to propagate
pub const CLIP_TOLERANCE: f32 = 1e6; // Named constant for tolerance
pub const MEDIUM_REFR_INDEX: RefrIndex = RefrIndex {
    // outer medium refractive index
    real: 1.0,
    imag: 0.0,
};
pub const MAX_REC: i32 = 5; // maximum number of beam recursions before truncation
pub const MAX_TIR: i32 = 5; // maximum number of total internal reflections
