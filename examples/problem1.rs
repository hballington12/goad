use nalgebra::{Complex, Point3, Vector3};
use pbt::problem::Problem;
use pbt::{
    beam::Beam,
    geom::{self, Face},
};

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    let projection = Vector3::new(0.0, 0.0, -1.0).normalize();
    let e_perp = Vector3::x(); // choose e_perp along z-axis for now
    let lower_left = vec![-10.0, -10.0];
    let upper_right = vec![10.0, 10.0];
    let clip_vertices = vec![
        Point3::new(upper_right[0], upper_right[1], 10.0),
        Point3::new(upper_right[0], lower_left[1], 10.0),
        Point3::new(lower_left[0], lower_left[1], 10.0),
        Point3::new(lower_left[0], upper_right[1], 10.0),
    ];

    let mut clip = Face::new_simple(clip_vertices, None).unwrap();
    clip.data_mut().area =
        Some((upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]));
    geom.shapes[0].refr_index.re = 1.5;
    geom.shapes[0].refr_index.im = 0.00001;

    let mut problem = Problem::new(
        geom,
        Beam::new_initial(clip, projection, Complex::new(1.00, 0.0), e_perp).unwrap(),
    );

    problem.solve_near();
}
