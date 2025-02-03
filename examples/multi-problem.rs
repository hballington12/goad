use goad::problem::MultiProblem;
use goad::settings::{self};

fn main() {
    let settings = settings::load_config();
    let mut multiproblem = MultiProblem::new(settings);

    multiproblem.solve();
    multiproblem.writeup();
}
