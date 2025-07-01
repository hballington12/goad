use goad::geom::{self};
use goad::problem::Problem;

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    geom.shapes[0].refr_index.re = 1.5;
    geom.shapes[0].refr_index.im = 0.0001;

    let mut problem = Problem::new(Some(geom), None);

    problem.solve_near();
}
