use goad::clip::Clipping;
use goad::geom::{self, Face};
use goad::helpers::draw_face;
use macroquad::prelude::*;
use nalgebra::{Point3, Vector3};

#[macroquad::main("Testing...")]
async fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/multiple.obj").unwrap();

    let projection = Vector3::new(0.0, 0.0, -1.0);

    let lower_left = vec![-19.0, -3.0];
    let upper_right = vec![10.0, 3.0];
    let mut clip_vertices = vec![
        Point3::new(lower_left[0], upper_right[1], 10.0),
        Point3::new(lower_left[0], lower_left[1], 10.0),
        Point3::new(upper_right[0], lower_left[1], 10.0),
        Point3::new(upper_right[0], upper_right[1], 10.0),
    ];
    clip_vertices.reverse();
    let mut clip = Face::new_simple(clip_vertices, None, None).unwrap();

    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    println!("{}", clipping.stats.unwrap());

    let intersections = clipping.intersections;
    let remaining = clipping.remaining;
    let clip = clipping.clip.clone();

    loop {
        clear_background(BLACK);

        // draw the original
        for shape in &geom.shapes {
            for face in &shape.faces {
                draw_face(face, GREEN, 4.0);
            }
        }
        // draw the original clip
        draw_face(&clip, RED, 10.0);
        // draw the remapped intersections
        for intsn in &intersections {
            draw_face(&intsn, YELLOW, 2.0);
        }
        // draw the remainders
        for face in &remaining {
            draw_face(face, BLUE, 2.0);
        }

        next_frame().await
    }
}
