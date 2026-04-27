use goad::problem::Problem;
use nalgebra::base;

// --8<-- [start:multiproblem]
fn main() {
    use goad::settings;

    // Setup and run a multi-orientation problem with default settings
    let mut base_settings = settings::load_default_config().unwrap();
    base_settings.quiet = true;
    // base_settings.geom_scale = Some(vec![2.622, 2.622, 3.745]);
    base_settings.geom_name = "examples/data/hex.obj".to_string();

    let distortions = vec![0.1];

    for distortion in distortions {
        let mut settings = base_settings.clone();
        settings.distortion = Some(distortion);
        let mut problem = Problem::new(None, Some(settings)).unwrap();
        let _ = problem.run(None);
        let _ = problem.geom.write_obj(format!("file_{}.obj", distortion));
        let _ = problem.writeup();
    }
}
// --8<-- [end:multiproblem]
