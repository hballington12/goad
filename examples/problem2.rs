use macroquad::prelude::*;
use nalgebra::{Complex, Point3, Vector3};
use pbt::beam::Beam;
use pbt::clip::Clipping;
use pbt::geom::{self, Face};
use pbt::helpers::draw_face;
use pbt::problem::Problem;

// #[macroquad::main("Testing...")]
// async fn main() {
fn main() {
    let geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    let projection = Vector3::new(0.0, 0.0, -1.0).normalize();
    let e_perp = Vector3::x(); // choose e_perp along x-axis for now

    let lower_left = vec![-10.0, -3.0];
    let upper_right = vec![10.0, 3.0];
    let clip_vertices = vec![
        Point3::new(lower_left[0], upper_right[1], 10.0),
        Point3::new(lower_left[0], lower_left[1], 10.0),
        Point3::new(upper_right[0], lower_left[1], 10.0),
        Point3::new(upper_right[0], upper_right[1], 10.0),
    ];
    let clip = Face::new_simple(clip_vertices, None);

    let mut problem = Problem::new(
        geom,
        Beam::new_initial(clip, projection, Complex::new(1.31, 0.1), e_perp),
    );

    problem.propagate_next();

    // loop {
    //     clear_background(BLACK);

    //     // draw the original
    //     for shape in &problem.geom.shapes {
    //         for face in &shape.faces {
    //             draw_face(face, GREEN, 4.0);
    //         }
    //     }
    //     // // draw the original clip
    //     // draw_face(&clip, RED, 10.0);
    //     // // draw the remapped intersections
    //     // for face in &intersections {
    //     //     draw_face(face, YELLOW, 2.0);
    //     // }
    //     // // draw the remainders
    //     // for face in &remaining {
    //     //     draw_face(face, BLUE, 2.0);
    //     // }

    //     next_frame().await
    // }
}
