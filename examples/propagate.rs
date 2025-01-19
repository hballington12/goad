use macroquad::prelude::*;
use nalgebra::{Point3, Vector3};
use pbt::beam::Beam;
use pbt::clip::Clipping;
use pbt::geom::{self, Face, RefrIndex};
use pbt::helpers::draw_face;
use pbt::problem::Problem;

#[macroquad::main("Testing...")]
async fn main() {
    let geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    let projection = Vector3::new(0.0, 0.0, -1.0);

    let lower_left = vec![-10.0, -3.0];
    let upper_right = vec![10.0, 3.0];
    let clip_vertices = vec![
        Point3::new(lower_left[0], upper_right[1], 10.0),
        Point3::new(lower_left[0], lower_left[1], 10.0),
        Point3::new(upper_right[0], lower_left[1], 10.0),
        Point3::new(upper_right[0], upper_right[1], 10.0),
    ];
    let clip = Face::new_simple(clip_vertices, None);
    let initial = clip.clone();

    let mut problem = Problem::new(
        geom,
        Beam::new_initial(clip, projection, RefrIndex::new(1.31, 0.1)),
    );

    println!(
        "number of beams in beam queue: {:?}",
        problem.beam_queue.len()
    );

    let mut propagation = None;

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
                println!("No more beams to propagate.");
            }
        }

        next_frame().await;
    }
}
