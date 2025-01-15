use geo_clipper::Clipper;
use geo_types::{Coord, LineString, Polygon};
use macroquad::prelude::*;
use pbt::geom::{self, Face};
use pbt::helpers::{draw_face, draw_multipolygon};

#[macroquad::main("Testing...")]
async fn main() {
    let shape = &geom::Geom::from_file("./concave1.obj").unwrap().shapes[0];

    let face2 = &shape.faces[0]; // to be the clip

    let mut exterior = Vec::new();
    for vertex in &face2.vertices {
        exterior.push(Coord {
            x: vertex.x,
            y: vertex.y,
        });
    }
    exterior.reverse();
    let clip = Polygon::new(LineString(exterior), vec![]);

    // let refs = shape.faces.iter().collect();
    // let output = Face::clip(face2, refs).unwrap();
    // let intersections = output.0;
    // let remaining = output.1;
    // println!("intersection polygons: {:?}", intersections.0.len());
    // println!("remaining polygons: {:?}", remaining.0.len());

    loop {
        clear_background(BLACK);

        // // Draw the intersection in red
        // for polygon in &intersections {
        //     // draw_face(polygon, RED);
        //     draw_multipolygon(&polygon, RED);
        // }

        // for polygon in &remaining {
        //     // draw_face(polygon, RED);
        //     draw_multipolygon(&polygon, GREEN);
        // }

        // Draw the clip and the subject
        // draw_multipolygon(&subject, BLUE);
        // draw_multipolygon(&clip, GREEN);
        // draw_multipolygon(&clip_proj, BLUE);

        next_frame().await
    }
}
