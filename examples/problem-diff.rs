use goad::geom::{self};
use goad::problem::Problem;

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    geom.shapes[0].refr_index.re = 1.31;
    geom.shapes[0].refr_index.im = 0.001;

    let mut problem = Problem::new(geom, None);
    problem.init();
    problem.solve();
    problem.writeup();
}
