// --8<-- [start:multiproblem]
fn main() {
    use goad::multiproblem::MultiProblem;
    use goad::settings;

    // Setup and run a multi-orientation problem with default settings
    let mut multiproblem = MultiProblem::new(None, settings::load_config().ok()).unwrap();
    multiproblem.solve();
}
// --8<-- [end:multiproblem]
