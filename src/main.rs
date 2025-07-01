//! > **Geometric Optics with Aperture Diffraction**
//!

use goad::multiproblem::MultiProblem;
use goad::settings::{self};

fn main() {
    let settings = settings::load_config().unwrap();
    let mut multiproblem = MultiProblem::new(None, Some(settings));

    multiproblem.solve();
    multiproblem.writeup();
}
