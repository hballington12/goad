use goad::geom::{self};
use goad::problem::Problem;

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/concave2.obj").unwrap();

    geom.shapes[0].refr_index.re = 1.31;
    geom.shapes[0].refr_index.im = 0.001;

    let mut problem = Problem::new(Some(geom), None);
    problem.geom.write_obj("hex_hollow_rot.obj").unwrap();

    println!("Resetting problem...");
    problem.reset();
    println!("Initializing problem...");
    problem.init();
    println!("Illuminating problem...");
    problem.illuminate();
    println!("Solving problem...");
    problem.solve();
    println!("Writing up results...");
    problem.writeup();
}
