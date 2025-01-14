use geo_clipper::Clipper;
use geo_types::{Coord, LineString, Polygon};
use macroquad::prelude::*;
use pbt::geom;
use pbt::helpers;

#[macroquad::main("Testing...")]
async fn main() {
    let shape = &geom::Geom::from_file("./concave1.obj").shapes[0];

    let face1 = &shape.faces[4];
    let face2 = &shape.faces[7];

    let subject = face1.polygon();
    let clip = face2.polygon();

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
