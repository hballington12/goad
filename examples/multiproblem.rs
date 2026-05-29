// --8<-- [start:multiproblem]
fn main() {
    use goad::multiproblem::MultiProblem;
    use goad::settings::{self, cli};

    // Setup and run a multi-orientation problem with default settings
    let settings = settings::load_config().ok();
    let geoms = cli::load_geoms().expect("load geoms");
    let mut multiproblem = MultiProblem::new(geoms, settings).unwrap();
    multiproblem.solve();
}
// --8<-- [end:multiproblem]
