use goad::problem::MultiProblem;
use goad::settings::{self};

fn main() {
    let settings = settings::load_config();
    let mut multiproblem = MultiProblem::new(settings.unwrap());

    multiproblem.solve();
    multiproblem.writeup();
}
