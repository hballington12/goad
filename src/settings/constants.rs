use super::{MuellerComponentConfig, OutputConfig};
use crate::{diff::Mapping, orientation::EulerConvention};
use std::path::PathBuf;

// ================================
// Computational/Physical Constants
// ================================

/// Minimum distance for vertices to be considered the same.
pub const VERTEX_MERGE_DISTANCE: f32 = 0.001;
/// Scaling factor for integer coordinates during clipping.
pub const CLIP_TOLERANCE: f32 = 1e16;
/// Minimum absolute value of the dot product of two vectors to be considered colinear.
pub const COLINEAR_THRESHOLD: f32 = 0.001;
/// Minimum vector length (in geometry units) to be considered non-degenerate.
pub const VEC_LENGTH_THRESHOLD: f32 = 0.001;
/// Minimum distance traversed by ray to intersection. Intersections closer than this are ignored.
pub const RAYCAST_MINIMUM_DISTANCE: f32 = 0.01;
/// Tolerance for diffraction computations, used to avoid divide by zero errors.
pub const DIFF_EPSILON: f32 = 1e-2;
/// Minimum dx or dy in diffraction computation.
pub const DIFF_DMIN: f32 = 1e-5;
/// Tolerance for kxx or kyy in diffraction computation.
pub const KXY_EPSILON: f32 = 1e-3;
/// Small perturbation for propagation distance to reduce errors in diffraction
pub const PROP_PERTURBATION: f32 = 1e-5;
/// Minimum Distortion factor for the geometry.
pub const MIN_DISTORTION: f32 = 1e-5;
/// Threshold for classification of direct forwards or backwards rays
pub const DIRECT_THRESHOLD: f32 = 1e-4;
/// Tolerance for value matching in interval binning
pub const INTERVAL_IGNORE_TOLERANCE: f32 = 0.0001;

// =============================
// Default Values for Python API
// =============================

/// Default wavelength in geometry units (532nm green laser)
pub const DEFAULT_WAVELENGTH: f32 = 0.532;
/// Default beam power threshold for ray termination
pub const DEFAULT_BEAM_POWER_THRESHOLD: f32 = 0.005;
/// Default beam area threshold factor
pub const DEFAULT_BEAM_AREA_THRESHOLD_FAC: f32 = 0.1;
/// Default power cutoff fraction (0-1)
pub const DEFAULT_CUTOFF: f32 = 0.99;
/// Default medium refractive index (vacuum/air)
pub const DEFAULT_MEDIUM_REFR_INDEX_RE: f32 = 1.0;
pub const DEFAULT_MEDIUM_REFR_INDEX_IM: f32 = 0.0;
/// Default particle refractive index (typical glass)
pub const DEFAULT_PARTICLE_REFR_INDEX_RE: f32 = 1.31;
pub const DEFAULT_PARTICLE_REFR_INDEX_IM: f32 = 0.0;
/// Default maximum recursion depth
pub const DEFAULT_MAX_REC: i32 = 10;
/// Default maximum total internal reflections
pub const DEFAULT_MAX_TIR: i32 = 10;
/// Default number of theta bins
pub const DEFAULT_THETA_BINS: usize = 181;
/// Default number of phi bins
pub const DEFAULT_PHI_BINS: usize = 181;
/// Default Euler angle order for the discrete orientation scheme.
pub const DEFAULT_EULER_ORDER: EulerConvention = EulerConvention::ZYZ;
/// Default mapping from near to far-field
pub const DEFAULT_MAPPING: Mapping = Mapping::ApertureDiffraction;
/// Default coherence settings
pub const DEFAULT_COHERENCE: bool = true;
/// Default quiet mode (false = show progress bars)
pub const DEFAULT_QUIET: bool = false;

// =================
// Default Functions
// =================

pub fn default_scale_factor() -> f32 {
    1.0
}

pub fn default_geom_scale() -> Option<Vec<f32>> {
    None
}

pub fn default_fov_factor() -> Option<f32> {
    None
}

pub fn default_directory() -> PathBuf {
    // Get current directory or default to a new PathBuf if it fails
    let current_dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::new());

    // Find the next available run number by checking existing directories
    let mut run_number = 1;
    let mut run_dir;

    loop {
        let run_name = format!("run{:05}", run_number);
        run_dir = current_dir.join(&run_name);

        if !run_dir.exists() {
            break;
        }

        run_number += 1;

        // Safety check to prevent infinite loops in extreme cases
        if run_number > 99999 {
            eprintln!("Warning: Exceeded maximum run number. Using timestamp instead.");
            let timestamp = chrono::Local::now().format("%Y%m%d_%H%M%S");
            run_dir = current_dir.join(format!("run_{}", timestamp));
            break;
        }
    }

    run_dir
}

pub fn default_output_config() -> OutputConfig {
    OutputConfig {
        results_summary: true,
        settings_json: false,
        powers_json: false,
        params_json: false,
        mueller_2d: true,
        mueller_1d: true,
        mueller_components: MuellerComponentConfig {
            total: true,
            beam: true,
            external: true,
        },
    }
}

pub fn default_quiet() -> bool {
    DEFAULT_QUIET
}
