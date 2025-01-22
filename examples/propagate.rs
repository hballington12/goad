use macroquad::prelude::*;
use nalgebra::{Complex, Point3, Vector, Vector3};
use pbt::{
    beam::Beam,
    geom::{self, Face},
};
use pbt::{helpers::draw_face, problem::Problem};

#[macroquad::main("Testing...")]
async fn main() {
    let geom = geom::Geom::from_file("./examples/data/cube.obj").unwrap();

    let projection = Vector3::new(0.0, -1.0, 0.0).normalize();
    let e_perp = Vector3::z(); // choose e_perp along z-axis for now

    let lower_left = vec![-10.0, -2.0];
    let upper_right = vec![10.0, 2.0];
    let clip_vertices = vec![
        Point3::new(lower_left[0], 10.0, upper_right[1]),
        Point3::new(lower_left[0], 10.0, lower_left[1]),
        Point3::new(upper_right[0], 10.0, lower_left[1]),
        Point3::new(upper_right[0], 10.0, upper_right[1]),
    ];
    let mut clip = Face::new_simple(clip_vertices, None);
    clip.data_mut().area =
        Some((upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]));
    let initial = clip.clone();

    let mut problem = Problem::new(
        geom,
        Beam::new_initial(clip, projection, Complex::new(1.00, 0.0), e_perp),
    );

    println!(
        "number of beams in beam queue: {:?}",
        problem.beam_queue.len()
    );

    let mut propagation: Option<pbt::beam::BeamPropagation> = None;

    loop {
        clear_background(BLACK);

        // Draw the current propagation if it exists
        if let Some(ref propagation) = propagation {
            problem.draw_propagation(propagation);
        }

        // Check if "Enter" is pressed
        if is_key_pressed(KeyCode::Enter) {
            if let Some(next_propagation) = problem.propagate_next() {
                propagation = Some(next_propagation);
                println!(
                    "number of beams in beam queue: {:?}",
                    problem.beam_queue.len()
                );
            } else {
                // println!("No more beams to propagate.");
            }
        }

        next_frame().await;
    }
}
