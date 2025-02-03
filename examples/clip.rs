use geo_clipper::Clipper;
use goad::geom;
use goad::helpers;
use macroquad::prelude::*;

#[macroquad::main("Testing...")]
async fn main() {
    let shape = &geom::Geom::from_file("./examples/data/concave1.obj")
        .unwrap()
        .shapes[0];

    let face1 = &shape.faces[4];
    let face2 = &shape.faces[7];

    let subject = face1.to_polygon();
    let clip = face2.to_polygon();

    let result = subject.intersection(&clip, 100000.0);

    loop {
        clear_background(BLACK);

        // Draw the clip and the subject
        helpers::draw_multipolygon(&subject, BLUE);
        helpers::draw_multipolygon(&clip, GREEN);

        // Draw the intersection in red
        for polygon in &result {
            helpers::draw_multipolygon(&polygon, RED);
        }

        next_frame().await
    }
}
