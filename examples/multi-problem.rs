use goad::multiproblem::MultiProblem;
use goad::settings::{self};

fn main() {
    let settings = settings::load_config();
    let mut multiproblem = MultiProblem::new(None, settings.ok());

    multiproblem.solve();
    multiproblem.writeup();
}
