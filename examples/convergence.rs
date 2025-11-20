//! > **Geometric Optics with Aperture Diffraction**
//!

use goad::{
    convergence::Convergence,
    settings::{self},
};

fn main() {
    let mut settings = settings::load_config().unwrap();
    settings.seed = Some(6);
    let mut convergence =
        Convergence::new(None, Some(settings), 10).expect("Failed to create MultiProblem");

    convergence.run();
}
