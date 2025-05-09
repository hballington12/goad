use goad::convergence::Convergence;
use goad::settings::{self};

fn main() {
    let settings = settings::load_config();

    let mut convergence = Convergence::new(settings.unwrap(), 50, vec![], 10);

    convergence.run();

    // multiproblem.solve();
    // multiproblem.writeup();
}
