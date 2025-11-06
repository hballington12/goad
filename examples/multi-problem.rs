use goad::multiproblem::MultiProblem;
use goad::settings::{self};

fn main() {
    let settings = settings::load_config();
    let mut multiproblem =
        MultiProblem::new(None, settings.ok()).expect("Failed to create MultiProblem");

    multiproblem.solve();
    let _ = multiproblem.writeup();
}
