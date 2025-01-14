use macroquad::prelude::*;
use nalgebra::{self as na, Isometry3, Point3, Vector3};
use pbt::clip::clip_faces;
use pbt::geom;
use pbt::helpers::draw_face;

#[macroquad::main("Testing...")]
async fn main() {
    let mut shape = geom::Shape::from_file("./octo2.obj");
    println!("{:?}", shape);
    let clip_index = 4; // the index of the face to be used as the clip
    let projection = Vector3::new(-0.3, 0.0, -1.0);

    let model = Isometry3::new(Vector3::zeros(), na::zero()); // do some sort of projection - set to nothing
    let origin = Point3::origin(); // camera location
    let target = Point3::new(projection.x, projection.y, projection.z); // projection direction, defines negative z-axis in new coords

    let up: Vector3<f32> = if projection.cross(&Vector3::y()).norm() < 0.01 {
        Vector3::x()
    } else {
        Vector3::y()
    };

    let view = Isometry3::look_at_rh(&origin, &target, &up);

    let transform = (view * model).to_homogeneous(); // transform to clipping system
    let itransform = transform.try_inverse().unwrap(); // inverse transform

    shape.transform(&transform); // transform to clipping coordinate system

    let mut clip = shape.faces.remove(clip_index); // choose a face be the clip

    // compute remapped intersections
    let mut remapped = clip_faces(&clip, &shape.faces);

    // transform back to original coordinate system
    shape.transform(&itransform);
    remapped
        .iter_mut()
        .for_each(|face| face.project(&itransform));
    clip.project(&itransform);

    loop {
        clear_background(BLACK);

        // draw the original
        for face in &shape.faces {
            draw_face(face, GREEN);
        }
        // draw the remapped intersections
        for face in &remapped {
            draw_face(face, YELLOW);
        }
        // draw the original clip
        draw_face(&clip, RED);

        next_frame().await
    }
}
