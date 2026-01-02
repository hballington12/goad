//! > **Geometric Optics with Aperture Diffraction**
//!

use goad::{
    multiproblem::MultiProblem,
    settings::{self},
};

fn main() {
    env_logger::init();
    let settings = settings::load_config().unwrap();
    let mut multiproblem =
        MultiProblem::new(None, Some(settings)).expect("Failed to create MultiProblem");

    multiproblem.solve();
    let _ = multiproblem.writeup();
}
