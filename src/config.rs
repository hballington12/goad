use crate::geom::RefrIndex;
pub const BEAM_AREA_THRESHOLD: f32 = 0.1; // minimum area for new beam to propagate
pub const CLIP_TOLERANCE: f32 = 1e6; // Named constant for tolerance
pub const MEDIUM_REFR_INDEX: RefrIndex = RefrIndex {
    // outer medium refractive index
    real: 1.0,
    imag: 0.0,
};
