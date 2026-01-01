use goad::geom::{self};
use goad::problem::Problem;

fn main() {
    let geoms = geom::Geom::load("./examples/data/hex2.obj").unwrap();
    let mut geom = geoms[0].clone();

    geom.shapes[0].refr_index.re = 1.5;
    geom.shapes[0].refr_index.im = 0.0001;

    let mut problem = Problem::new(Some(geom), None).unwrap();

    problem.solve_near();
}
