use macroquad::prelude::*;
use nalgebra::Vector3;
use pbt::clip::Clipping;
use pbt::geom;
use pbt::helpers::draw_face;

#[macroquad::main("Testing...")]
async fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/concave1.obj").unwrap();

    let clip_index = 4; // the index of the face to be used as the clip
    let projection = Vector3::new(-0.3, 0.0, -1.0);
    let mut clip = geom.shapes[0].faces.remove(clip_index); // choose a face be the clip

    // start function `do_clip` here:
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    clipping.clip();

    loop {
        clear_background(BLACK);

        // draw the original
        for face in &clipping.geom.shapes[0].faces {
            draw_face(face, GREEN, 2.0);
        }
        // draw the remapped intersections
        for face in &clipping.intersections {
            draw_face(face, YELLOW, 2.0);
        }
        // draw the original clip
        draw_face(&clipping.clip, RED, 2.0);

        next_frame().await
    }
}
