use goad::clip::Clipping;
use goad::geom;
use goad::helpers::draw_face;
use macroquad::prelude::*;
use nalgebra::Vector3;

#[macroquad::main("Testing...")]
async fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/cube_inside_ico.obj").unwrap();

    let clip_index = 17; // the index of the face to be used as the clip
    let projection = Vector3::new(-0.2, 0.0, -1.0);
    let mut clip = geom.shapes[0].faces.remove(clip_index); // choose a face be the clip

    // start function `do_clip` here:
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    println!("{}", clipping.stats.unwrap());

    loop {
        clear_background(BLACK);

        // draw the original
        // ico
        for face in &clipping.geom.shapes[0].faces {
            draw_face(face, GREEN, 2.0);
        }
        // // cube
        // for face in &clipping.geom.shapes[1].faces {
        //     draw_face(face, MAGENTA, 2.0);
        // }
        // draw the remapped intersections
        for intsn in &clipping.intersections {
            draw_face(&intsn, YELLOW, 2.0);
        }
        // draw the original clip
        draw_face(&clipping.clip, RED, 2.0);

        next_frame().await
    }
}
