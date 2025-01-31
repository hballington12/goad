use clap::Parser;
use config::{Config, File};
use nalgebra::Complex;
use serde::Deserialize;
use std::f32::consts::PI;
use std::fmt;
use std::sync::Mutex;

// /// Wavelength of the electric field in the same units as the geometry coordinates.
// pub const WAVELENGTH: f32 = 0.532 * 0.1;
// /// Wavenumber, equal to 2π divided by `WAVELENGTH`.
// pub const WAVENUMBER: f32 = 2.0 * PI / WAVELENGTH;
// /// Minimum distance for vertices to be considered the same.
pub const VERTEX_MERGE_DISTANCE: f32 = 0.01;
// /// Minimum area for new beam to propagate, set by default to `WAVELENGTH`^2 / 4.
// pub const BEAM_AREA_THRESHOLD: f32 = WAVELENGTH * WAVELENGTH / 4.0;
// /// Minimum power for new beam to propagate.
// pub const BEAM_POWER_THRESHOLD: f32 = 0.005;
// /// Scaling factor for integer coordinates during clipping.
pub const CLIP_TOLERANCE: f32 = 1e10;
// /// Minimum absolute value of the dot product of two vectors to be considered colinear.
pub const COLINEAR_THRESHOLD: f32 = 0.001;
// /// Minimum vector length (in geometry units) to be considered non-degenerate.
pub const VEC_LENGTH_THRESHOLD: f32 = 0.01;
// /// Minimum distance traversed by ray to intersection. Intersections closer than this are ignored.
pub const RAYCAST_MINIMUM_DISTANCE: f32 = 0.01;
// /// Acceptable energy threshold beyond which beam propagation terminates.
// pub const TOTAL_POWER_CUTOFF: f32 = 0.99;
// /// Surrounding medium refractive index.
// pub const MEDIUM_REFR_INDEX: Complex<f32> = Complex { re: 1.0, im: 0.0 };
// /// Maximum number of beam recursions before truncation.
// pub const MAX_REC: i32 = 10;
// /// Maximum number of total internal reflections.
// pub const MAX_TIR: i32 = 10;
// /// Distance to far-field. 1e3 - 1e5 is a good range for single precision arithmetic.
pub const RADIUS: f32 = 1e4;
// /// Tolerance for diffraction computations, used to avoid divide by zero errors.
pub const DIFF_EPSILON: f32 = 1e-4;
// /// Base far-field binning reolution.
// pub const FAR_FIELD_RESOLUTION: usize = 360;

/// Runtime configuration for the application.
#[derive(Debug, Deserialize, PartialEq)]
pub struct Settings {
    pub wavelength: f32,
    pub vertex_merge_distance: f32,
    pub beam_power_threshold: f32,
    pub beam_area_threshold_fac: f32,
    pub clip_tolerance: f32,
    pub colinear_threshold: f32,
    pub vec_length_threshold: f32,
    pub raycast_minimum_distance: f32,
    pub total_power_cutoff: f32,
    pub medium_refr_index: Complex<f32>,
    pub max_rec: i32,
    pub max_tir: i32,
    pub radius: f32,
    pub diff_epsilon: f32,
    pub far_field_resolution: usize,
}

impl Settings {
    pub fn beam_area_threshold(&self) -> f32 {
        self.wavelength * self.wavelength * self.beam_area_threshold_fac
    }
}

pub fn load_config() -> Settings {
    let mut settings = Config::builder()
        .add_source(File::with_name("config/default"))
        .build()
        .unwrap();

    // println!("settings: {:#?}", settings);

    settings.try_deserialize().unwrap()
}

// impl Default for Settings {
//     fn default() -> Self {
//         Self {
//             wavelength: WAVELENGTH,
//             vertex_merge_distance: VERTEX_MERGE_DISTANCE,
//             beam_power_threshold: BEAM_POWER_THRESHOLD,
//             clip_tolerance: CLIP_TOLERANCE,
//             colinear_threshold: COLINEAR_THRESHOLD,
//             vec_length_threshold: VEC_LENGTH_THRESHOLD,
//             raycast_minimum_distance: RAYCAST_MINIMUM_DISTANCE,
//             total_power_cutoff: TOTAL_POWER_CUTOFF,
//             medium_refr_index: MEDIUM_REFR_INDEX,
//             max_rec: MAX_REC,
//             max_tir: MAX_TIR,
//             radius: RADIUS,
//             diff_epsilon: DIFF_EPSILON,
//             far_field_resolution: FAR_FIELD_RESOLUTION,
//         }
//     }
// }

// #[derive(Parser, Debug)]
// #[command(version, about = "Light Scattering Simulation")]
// pub struct CliArgs {
//     /// Override wavelength
//     #[arg(long)]
//     wavelength: Option<f32>,

//     /// Override vertex merge distance
//     #[arg(long)]
//     vertex_merge_distance: Option<f32>,
// }

// // Global settings instance, wrapped in a Mutex for thread safety
// pub static GLOBAL_SETTINGS: Lazy<Mutex<Settings>> = Lazy::new(|| {
//     // Load settings from a file or elsewhere
//     Mutex::new(load_config())
// });

impl fmt::Display for Settings {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Settings:
  - Wavelength: {:.6}
  - Vertex Merge Distance: {:.6}
  - Beam Power Threshold: {:.6}
  - Clip Tolerance: {:.6}
  - Colinear Threshold: {:.6}
  - Vec Length Threshold: {:.6}
  - Raycast Minimum Distance: {:.6}
  - Total Power Cutoff: {:.6}
  - Medium Refractive Index: {:.6} + {:.6}i
  - Max Rec: {}
  - Max TIR: {}
  - Radius: {:.6}
  - Diff Epsilon: {:.6}
  - Far Field Resolution: {}",
            self.wavelength,
            self.vertex_merge_distance,
            self.beam_power_threshold,
            self.clip_tolerance,
            self.colinear_threshold,
            self.vec_length_threshold,
            self.raycast_minimum_distance,
            self.total_power_cutoff,
            self.medium_refr_index.re,
            self.medium_refr_index.im,
            self.max_rec,
            self.max_tir,
            self.radius,
            self.diff_epsilon,
            self.far_field_resolution
        )
    }
}
