//! > **Geometric Optics with Aperture Diffraction**
//!
//! # What is GOAD?
//! - GOAD is a rust crate for simulating light propagation through the use of
//! geometric optics combined with diffraction theory of a plane wave at an
//! aperture.
//! - Most users will likely be interested in running the `goad` binary, which
//! provides a command line interface for running a general problem. To get
//! started, have a look at the [quick start guide][_quickstart].

pub mod _quickstart;
pub mod beam;
pub mod bins;
pub mod clip;
pub mod containment;
pub mod diff;
pub mod distortion;
pub mod field;
pub mod fresnel;
pub mod geom;
pub mod helpers;
pub mod multiproblem;
pub mod orientation;
pub mod output;
pub mod params;
pub mod powers;
pub mod problem;
pub mod result;
pub mod settings;
pub mod snell;
