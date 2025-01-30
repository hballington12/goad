use nalgebra::{Complex, Point3, Vector3};
use pbt::problem::Problem;
use pbt::{
    beam::Beam,
    geom::{self, Face},
};

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    geom.shapes[0].refr_index.re = 1.5;
    geom.shapes[0].refr_index.im = 0.00001;

    let mut problem = Problem::new(geom);

    problem.solve_near();
}
