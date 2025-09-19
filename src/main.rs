//! > **Geometric Optics with Aperture Diffraction**
//!

use env_logger::Builder;
use std::fs::OpenOptions;

use goad::multiproblem::MultiProblem;
use goad::settings::{self};

fn main() {
    let log_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("goad.log")
        .expect("Failed to open log file");

    Builder::from_default_env()
        .target(env_logger::Target::Pipe(Box::new(log_file)))
        .init();

    let settings = settings::load_config().unwrap();
    let mut multiproblem = MultiProblem::new(None, Some(settings));

    multiproblem.solve();
    let _ = multiproblem.writeup();
}
