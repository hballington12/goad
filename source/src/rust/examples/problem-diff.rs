use goad::cancel::CancelToken;
use goad::geom::{self};
use goad::problem::Problem;
use goad::settings::DEFAULT_PARTICLE_REFR_INDEX;

fn main() {
    let geoms =
        geom::Geom::load("./examples/data/concave2.obj", vec![DEFAULT_PARTICLE_REFR_INDEX])
            .unwrap();
    let mut geom = geoms[0].clone();

    geom.shapes[0].refr_index.re = 1.31;
    geom.shapes[0].refr_index.im = 0.001;

    let mut problem = Problem::new(geom, None).unwrap();
    problem.geom.write_obj("hex_hollow_rot.obj").unwrap();

    println!("Resetting problem...");
    problem.reset();
    println!("Initializing problem...");
    problem.init();
    println!("Illuminating problem...");
    problem.illuminate().unwrap();
    println!("Solving problem...");
    problem.solve(&CancelToken::noop()).unwrap();
    println!("Writing up results...");
    problem.writeup();
}
