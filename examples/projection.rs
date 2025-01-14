use std::time::Duration;

use geo::{HasDimensions, MultiPolygon};
use geo_clipper::Clipper;
use geo_types::{Coord, LineString, Polygon};
use macroquad::prelude::*;
use nalgebra::Perspective3;
use nalgebra::{self as na, Matrix4, Vector4};
use nalgebra::{Isometry3, Point3, Vector3};
use pbt::geom::{self, Face, Plane};
use pbt::helpers::draw_face;
use pbt::helpers::{self, draw_multipolygon};

trait PolygonExtensions {
    fn project(&self, plane: &Plane) -> Face;
}

impl PolygonExtensions for Polygon<f32> {
    /// Projects the xy coordinates of a polygon onto a plane in 3D
    /// Ignores the last vertex, which is a duplicate of the first
    fn project(&self, plane: &Plane) -> Face {
        let mut vertices = Vec::new();

        for coord in self.exterior().0.iter().take(self.exterior().0.len() - 1) {
            // compute z intersection via z = -(ax + by + d)/c
            let z = -(plane.normal.x * coord.x + plane.normal.y * coord.y + plane.offset)
                / plane.normal.z;

            vertices.push(Point3::new(coord.x, coord.y, z));
        }

        Face::new(&vertices)
    }
}

fn clip_polygons(clip_in: &Face, subjects_in: &Vec<Face>) -> Vec<Face> {
    if subjects_in.is_empty() {
        return Vec::new();
    }

    let clip = clip_in.polygon();

    let mut intersections: Vec<Face> = Vec::new(); // the final remapped
    let mut remaining_clips = vec![clip];

    // for subject in subjects {
    for subject_in in subjects_in {
        let subject = subject_in.polygon();
        let mut next_clips = Vec::new();

        for clip in remaining_clips {
            let mut remapped: Vec<Face> = Vec::new();
            let intersection = subject.intersection(&clip, 100000.0);
            let difference = clip.difference(&subject, 100000.0);

            for intsn in intersection {
                let face = intsn.project(&subject_in.plane());
                remapped.push(face);
            }

            intersections.extend(remapped);
            next_clips.extend(difference);
        }
        remaining_clips = next_clips; // Update clips for the next subject
        if remaining_clips.is_empty() {
            break;
        }
    }

    intersections
}

fn project_face(face: &Face, model_view: &Matrix4<f32>) -> Face {
    let projected_vertices: Vec<Point3<f32>> = face
        .vertices
        .iter()
        .map(|point| {
            let vertex4 = Vector4::new(point.x, point.y, point.z, 1.0);
            let projected_vertex = model_view * vertex4;
            Point3::new(projected_vertex.x, projected_vertex.y, projected_vertex.z)
        })
        .collect();

    Face::new(&projected_vertices) // note coordinate system
}

#[macroquad::main("Testing...")]
async fn main() {
    let shape = &geom::Geom::from_file("./hex3.obj").shapes[0];

    // do some sort of projection

    // let model = Isometry3::new(Vector3::y(), na::zero());
    let model = Isometry3::new(Vector3::zeros(), na::zero());

    let origin = Point3::origin();
    let target = Point3::new(0.5, 0.0, -1.0);
    let view = Isometry3::look_at_rh(&origin, &target, &Vector3::y());

    let transform = (view * model).to_homogeneous();

    let itransform = transform.try_inverse().unwrap();

    let mut proj_faces: Vec<Face> = shape
        .faces
        .iter()
        .map(|face| project_face(face, &transform))
        .collect();

    // choose a face be the clip
    let clip = proj_faces.remove(4);

    // compute remapped intersections
    let remapped = clip_polygons(&clip, &proj_faces);

    // use inverse transform to get back to original coordinates
    let reremapped: Vec<Face> = remapped
        .iter()
        .map(|face| project_face(face, &itransform))
        .collect();
    let reclip = project_face(&clip, &itransform);

    loop {
        clear_background(BLACK);

        // draw the original
        for face in &shape.faces {
            draw_face(face, GREEN);
        }

        // draw the remapped intersections
        for face in &reremapped {
            draw_face(face, YELLOW);
        }

        // draw the original clip
        draw_face(&reclip, RED);

        next_frame().await
    }
}
