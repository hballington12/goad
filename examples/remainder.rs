use macroquad::prelude::*;
use nalgebra::{Point3, Vector3};
use pbt::clip::{self, Clipping};
use pbt::geom::{self, Face, FaceData};
use pbt::helpers::draw_face;

#[macroquad::main("Testing...")]
async fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/concave1.obj").unwrap();

    let projection = Vector3::new(0.0, 0.0, -1.0);

    let mut clip_vertices = vec![
        Point3::new(-7.0, 7.0, 10.0),
        Point3::new(-7.0, -7.0, 10.0),
        Point3::new(7.0, -7.0, 10.0),
        Point3::new(7.0, 7.0, 10.0),
    ];
    clip_vertices.reverse();
    let mut clip = Face::new_simple(clip_vertices);

    // start function `do_clip` here:
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    clipping.clip();

    println!("number of remaining: {:?}", clipping.remaining.len());

    loop {
        clear_background(BLACK);

        // draw the original
        for face in &clipping.geom.shapes[0].faces {
            draw_face(face, GREEN, 4.0);
        }
        // draw the remapped intersections
        for face in &clipping.intersections {
            draw_face(face, YELLOW, 2.0);
        }
        // draw the original clip
        draw_face(&clipping.clip, RED, 10.0);
        // draw the remainders
        for face in &clipping.remaining {
            draw_face(face, BLUE, 2.0);
        }

        next_frame().await
    }
}
