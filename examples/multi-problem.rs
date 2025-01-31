use geo::Orientation;
use pbt::geom::{self};
use pbt::orientation::Orientations;
use pbt::problem::{MultiProblem, Problem};

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    geom.shapes[0].refr_index.re = 1.31;
    geom.shapes[0].refr_index.im = 0.001;

    let alphas = vec![0.0, 0.0];
    let betas = vec![0.0, 0.0];
    let gammas = vec![0.0, 0.0];

    let orientations = Orientations::new_discrete(alphas, betas, gammas).unwrap();
    let mut multiproblem = MultiProblem::new(geom, orientations);

    multiproblem.solve();
    multiproblem.writeup();
}
