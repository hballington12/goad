use pbt::geom::{self};
use pbt::problem::Problem;

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/8col.obj").unwrap();

    // geom.shapes[0].rescale(2.0);

    geom.shapes[0].refr_index.re = 1.31;
    geom.shapes[0].refr_index.im = 0.0000;

    let mut problem = Problem::new(geom);

    problem.solve_near();

    // problem.solve_far_ext_diff();
    // problem.solve_far_outbeams();
    problem.solve_far();
}
