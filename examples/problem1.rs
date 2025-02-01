use pbt::geom::{self};
use pbt::problem::Problem;

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    geom.euler_rotate(30.0, 30.0, 0.0);

    geom.shapes[0].refr_index.re = 1.5;
    geom.shapes[0].refr_index.im = 0.0001;

    let mut problem = Problem::new(geom, None);

    problem.solve_near();
}
