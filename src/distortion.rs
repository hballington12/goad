use std::collections::HashMap;

use nalgebra::{Matrix3, Vector3};
use rand::Rng;

use crate::geom::{Face, Geom};

impl Geom {
    pub fn distort(&mut self, sigma: f32) {
        // For each shape in geometry:
        for shape in self.shapes.iter_mut() {
            // Prescan to hold a list of which vertices are in which faces
            let vertex_to_faces = build_vertex_to_face_map(shape);

            // Check if the shape can be distorted
            if shape_can_be_distorted(&vertex_to_faces) {
                return;
            }

            // Perturb the normals of the faces
            let perturbed_normals = perturb_normals(sigma, shape);

            // Solve the linear system to get the new vertex positions
            solve_vertices(shape, vertex_to_faces, perturbed_normals);

            // Update the vertex positions in the faces
            update_face_vertices(shape);
        }
    }
}

fn update_face_vertices(shape: &mut crate::geom::Shape) {
    // Update vertex positions in faces
    for face in shape.faces.iter_mut() {
        match face {
            Face::Simple(data) => {
                if let Some(indices) = &data.exterior_indices {
                    println!(
                        "Updating vertex positions in face with {} vertices",
                        indices.len()
                    );
                    for (pos, &index) in indices.iter().enumerate() {
                        if index < shape.vertices.len() {
                            // println!("old position is {:?}", data.exterior[pos]);
                            let vertex = shape.vertices[index];
                            // println!("new position is {:?}", vertex);
                            data.exterior[pos] = vertex;
                        }
                    }
                }
            }
            Face::Complex { .. } => {
                panic!("Complex faces not supported for distortion");
            }
        }
    }
}

fn solve_vertices(
    shape: &mut crate::geom::Shape,
    vertex_to_faces: HashMap<usize, Vec<usize>>,
    perturbed_normals: Vec<nalgebra::Vector3<f32>>,
) {
    // For each vertex in the shape
    // Get the perturbed normals of the faces it belongs to
    // (use the mapping to get the faces it belongs to)
    // Solve the linear system to get the new vertex position
    for (vertex_index, faces) in vertex_to_faces.iter() {
        let mut norms = Vec::new(); // these are the unperturbed normals
        let mut pnorms = Vec::new(); // these are the perturbed normals
        let mut mids = Vec::new(); // these are the midpoints of the faces

        // loop over faces that this vertex belongs to
        for face_index in faces {
            let face = &shape.faces[*face_index];
            norms.push(face.data().normal);
            pnorms.push(perturbed_normals[*face_index]);
            mids.push(face.data().midpoint);
        }

        println!("normals are {:?}", norms);
        println!("perturbed normals are {:?}", pnorms);
        println!("index is {:?}", vertex_index);

        // Solve the linear system to get the new vertex position
        // the solution is the intersection of the planes defined by the normals
        // This is a simple linear system of equations
        // Ax = b, where A is the matrix of normals, x is the new vertex position,
        // and b is the vector of the original vertex position
        let mut a = Matrix3::zeros();
        let mut b = Vector3::zeros();
        for (i, normal) in norms.iter().enumerate() {
            // a[(i, 0)] = normal.x;
            // a[(i, 1)] = normal.y;
            // a[(i, 2)] = normal.z;
            b[i] = mids[i].coords.dot(normal); // tilt is about the centroid
                                               // b[i] = shape.vertices[*vertex_index].coords.dot(normal);
        }
        for (i, normal) in pnorms.iter().enumerate() {
            a[(i, 0)] = normal.x;
            a[(i, 1)] = normal.y;
            a[(i, 2)] = normal.z;
        }
        // Solve the linear system
        let new_vertex = a.try_inverse().expect("could not invert matrix") * b;
        shape.vertices[*vertex_index].coords = new_vertex;
    }
}

fn perturb_normals(
    sigma: f32,
    shape: &mut crate::geom::Shape,
) -> Vec<
    nalgebra::Matrix<
        f32,
        nalgebra::Const<3>,
        nalgebra::Const<1>,
        nalgebra::ArrayStorage<f32, 3, 1>,
    >,
> {
    // Perturb the normal of each face in the shape
    let mut perturbed_normals = Vec::new();
    for face in shape.faces.iter_mut() {
        let mut rng = rand::rng();
        let perturbation = Vector3::new(
            rng.random_range(-sigma..sigma),
            rng.random_range(-sigma..sigma),
            rng.random_range(-sigma..sigma),
        );
        println!("old normal is {:?}", face.data().normal);
        let mut perturbed = face.data().normal + perturbation;
        perturbed.normalize_mut();
        println!("new normal is {:?}", perturbed);
        perturbed_normals.push(perturbed);
    }
    perturbed_normals
}

fn shape_can_be_distorted(vertex_to_faces: &HashMap<usize, Vec<usize>>) -> bool {
    if vertex_to_faces.values().any(|faces| faces.len() != 3) {
        true
    } else {
        println!("Shape has vertices that do not belong to exactly 3 faces. Skipping distortion.");
        false
    }
}

fn build_vertex_to_face_map(shape: &mut crate::geom::Shape) -> HashMap<usize, Vec<usize>> {
    let mut vertex_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (face_index, face) in shape.faces.iter().enumerate() {
        for vertex in &face.data().exterior {
            let vertex_index = shape.vertices.iter().position(|v| v == vertex).unwrap();
            vertex_to_faces
                .entry(vertex_index)
                .or_insert_with(Vec::new)
                .push(face_index);
        }
    }
    vertex_to_faces
}
