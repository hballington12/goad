use goad::cancel::CancelToken;
use goad::geom::Geom;
use goad::problem::Problem;
use goad::settings::DEFAULT_PARTICLE_REFR_INDEX;

// --8<-- [start:multiproblem]
fn main() {
    use goad::settings;

    // Setup and run a multi-orientation problem with default settings
    let mut base_settings = settings::load_default_config().unwrap();
    base_settings.quiet = true;
    // base_settings.geom_scale = Some(vec![2.622, 2.622, 3.745]);

    let distortions = vec![0.1];

    for distortion in distortions {
        let mut settings = base_settings.clone();
        settings.distortion = Some(distortion);
        let geom = Geom::load("examples/data/hex.obj", vec![DEFAULT_PARTICLE_REFR_INDEX])
            .unwrap()
            .remove(0);
        let mut problem = Problem::new(geom, Some(settings)).unwrap();
        let _ = problem.run(None, &CancelToken::noop());
        let _ = problem.geom.write_obj(format!("file_{}.obj", distortion));
        let _ = problem.writeup();
    }
}
// --8<-- [end:multiproblem]
