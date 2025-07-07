use goad::problem::Problem;
use goad::{
    beam::Beam,
    geom::{self, Face},
};
use macroquad::prelude::*;
use nalgebra::{Complex, Point3, Vector3};
use std::io::{self, Write};

#[macroquad::main("Testing...")]
async fn main() {
    // let mut geom = geom::Geom::from_file("./examples/data/hex_20_30_30_face.obj").unwrap();
    let mut geom = geom::Geom::from_file("./examples/data/hex.obj").unwrap();
    let mut geom = geom::Geom::from_file("./examples/data/hex_hollow.obj").unwrap();

    // let euler = Euler::new(0.0, 30.0, 0.0);
    // let _ = geom.euler_rotate(euler, EulerConvention::ZYZ);

    let projection = Vector3::new(0.0, -1.0, 0.0).normalize();
    let e_perp = Vector3::x(); // choose e_perp along z-axis for now

    let lower_left = vec![-10.0, -10.0];
    let upper_right = vec![10.0, 10.0];
    // let clip_vertices = vec![
    //     Point3::new(upper_right[0], upper_right[1], 10.0),
    //     Point3::new(upper_right[0], lower_left[1], 10.0),
    //     Point3::new(lower_left[0], lower_left[1], 10.0),
    //     Point3::new(lower_left[0], upper_right[1], 10.0),
    // ];

    let clip_vertices = vec![
        Point3::new(lower_left[0], 10.0, upper_right[1]),
        Point3::new(lower_left[0], 10.0, lower_left[1]),
        Point3::new(upper_right[0], 10.0, lower_left[1]),
        Point3::new(upper_right[0], 10.0, upper_right[1]),
    ];

    let mut clip = Face::new_simple(clip_vertices, None, None).unwrap();
    clip.data_mut().area =
        Some((upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]));
    geom.shapes[0].refr_index.re = 1.5;
    // geom.shapes[0].refr_index.im = 1e-3;
    // geom.shapes[1].refr_index.re = 2.0;
    // geom.shapes[2].refr_index.re = 2.5;
    let wavelength = 0.532;

    let mut problem = Problem::new_with_field(
        geom,
        Beam::new_initial(
            clip,
            projection,
            Complex::new(1.00, 0.0),
            e_perp,
            wavelength,
        )
        .unwrap(),
    );

    println!(
        "initial number of beams in beam queue: {:?}",
        problem.beam_queue.len()
    );

    let mut propagations: Vec<Option<goad::beam::BeamPropagation>> = vec![None; 0];

    println!("running...");

    loop {
        clear_background(BLACK);

        // Draw the current propagation if it exists
        for prop in propagations.iter() {
            if let Some(ref propagation) = prop {
                problem.draw_propagation(propagation);
            }
        }

        // Check if "Enter" is pressed
        if is_key_pressed(KeyCode::Enter) {
            let next_propagation = problem.propagate_next();
            match next_propagation {
                Some(next_propagation) => {
                    propagations.push(Some(next_propagation));
                }
                _ => {
                    println!("No more beams to propagate.");
                }
            }

            // Print the number of beams in the queue without a new line and flush the output
            println!(
                "\rnumber of beams in beam queue: {:?}    ", // Add spaces to overwrite previous text
                problem.beam_queue.len()
            );
            println!(
                "number of beams in outbeams: {:?}    ",
                problem.out_beam_queue.len()
            );
            io::stdout().flush().unwrap();
        } else {
            println!("Press Enter to propagate another beam.");
            // break;
        }

        next_frame().await;
    }
}
