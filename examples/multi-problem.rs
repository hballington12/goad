use pbt::problem::MultiProblem;
use pbt::settings::{self};

fn main() {
    let settings = settings::load_config();
    let mut multiproblem = MultiProblem::new(settings);

    multiproblem.solve();
    multiproblem.writeup();
}
