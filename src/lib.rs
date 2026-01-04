//! > **Geometric Optics with Aperture Diffraction**
//!
//! # What is GOAD?
//! - GOAD is a rust crate for simulating light propagation through the use of
//! geometric optics combined with diffraction theory of a plane wave at an
//! aperture.
//! - Most users will likely be interested in running the `goad` binary, which
//! provides a command line interface for running a general problem. To get
//! started, have a look at the [quick start guide][_quickstart].

// Use jemalloc as the global allocator on macOS only
// This avoids memory corruption issues with the macOS system allocator
// On Linux, jemalloc can cause TLS allocation errors when loaded as a Python extension
#[cfg(all(target_os = "macos", not(target_env = "msvc")))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

pub mod _quickstart;
pub mod beam;
pub mod bins;
pub mod clip;
pub mod containment;
pub mod convergence;
pub mod diff;
pub mod diff2;
pub mod distortion;
pub mod field;
pub mod filelog;
pub mod fresnel;
pub mod geom;
pub mod multiproblem;
pub mod orientation;
pub mod output;
pub mod params;
pub mod powers;
pub mod problem;
pub mod python;
#[cfg(feature = "stub-gen")]
pub use python::stub_info;
pub mod result;
pub mod settings;
pub mod snell;
pub mod zones;
