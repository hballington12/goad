//! > **Geometric Optics with Aperture Diffraction**
//!

use goad::{
    filelog,
    multiproblem::MultiProblem,
    settings::{self},
};

fn main() {
    let settings = settings::load_config().unwrap();

    // Initialize file-based logging to output directory
    if let Err(e) = filelog::init(&settings.directory) {
        eprintln!("Warning: Could not initialize log file: {}", e);
    }

    let mut multiproblem =
        MultiProblem::new(None, Some(settings)).expect("Failed to create MultiProblem");

    multiproblem.solve();
    let _ = multiproblem.writeup();
}
