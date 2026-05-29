//! > **Geometric Optics with Aperture Diffraction**
//!

use goad::{
    multiproblem::MultiProblem,
    settings::{self, cli},
};

fn main() {
    let settings = settings::load_config().unwrap();
    let geoms = cli::load_geoms().expect("Failed to load geometry");
    let mut multiproblem =
        MultiProblem::new(geoms, Some(settings)).expect("Failed to create MultiProblem");

    multiproblem.solve();
    let _ = multiproblem.writeup();
}
