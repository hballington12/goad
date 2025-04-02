use std::collections::HashMap;

use nalgebra::{Matrix3, Vector3};
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};

use crate::{
    geom::{Face, Geom},
    settings::MIN_DISTORTION,
};

impl Geom {
    pub fn distort(&mut self, sigma: f32, seed: Option<u64>) {
        if sigma <= MIN_DISTORTION {
            return;
        }

        // For each shape in geometry:
        for shape in self.shapes.iter_mut() {
            // Use the shape aabb to get the bounding box
            if shape.aabb.is_none() {
                shape.set_aabb();
            }
            let aabb = shape.aabb.clone().unwrap();
            let max_dim = (aabb.max - aabb.min).norm();

            // Prescan to hold a list of which vertices are in which faces
            let vertex_to_faces = build_vertex_to_face_map(shape);

            // Check if the shape can be distorted
            if !shape_can_be_distorted(&vertex_to_faces) {
                return;
            }

            // Keep original vertices in case we need to restore them
            let original_vertices = shape.vertices.clone();

            // Try distortion up to 25 times
            let max_attempts = 25;
            let mut attempt = 0;

            loop {
                attempt += 1;

                // Perturb the normals of the faces
                let perturbed_normals = perturb_normals(sigma, shape, seed);

                // Solve the linear system to get the new vertex positions
                solve_vertices(shape, &vertex_to_faces, &perturbed_normals);

                // Update the vertex positions in the faces
                update_face_vertices(shape);

                // Check for self-intersections
                let mut distort_failed = false;
                for face in &shape.faces {
                    if face.data().self_intersects() {
                        distort_failed = true;
                        break;
                    } else if !face.data().is_convex() {
                        distort_failed = true;
                        break;
                    }
                }

                // If no self-intersections or max attempts reached, break the loop
                if !distort_failed || attempt >= max_attempts {
                    if distort_failed && attempt >= max_attempts {
                        println!(
                            "Maximum retries ({}) reached. Reverting to original shape.",
                            max_attempts
                        );
                        // Restore original vertices
                        shape.vertices = original_vertices;
                        update_face_vertices(shape);
                    }
                    break;
                }

                // Reset vertices for next attempt
                shape.vertices = original_vertices.clone();
            }

            // Get new AABB after distortion
            shape.set_aabb();

            // Get new max dimension after distortion
            let new_aabb = shape.aabb.clone().unwrap();
            let new_max_dim = (new_aabb.max - new_aabb.min).norm();

            let rescale_fac = max_dim / new_max_dim;
            shape.rescale(rescale_fac);
        }
    }
}

fn update_face_vertices(shape: &mut crate::geom::Shape) {
    // Update vertex positions in faces
    for face in shape.faces.iter_mut() {
        match face {
            Face::Simple(data) => {
                if let Some(indices) = &data.exterior_indices {
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
    vertex_to_faces: &HashMap<usize, Vec<usize>>,
    perturbed_normals: &Vec<nalgebra::Vector3<f32>>,
) {
    // For each vertex in the shape
    // Get the perturbed normals of the faces it belongs to
    // (use the mapping to get the faces it belongs to)
    // Solve the linear system to get the new vertex position
    for (vertex_index, faces) in vertex_to_faces.iter() {
        let (norms, pnorms, mids) = fetch_face_data(shape, perturbed_normals, faces);

        let new_vertex = solve_linear_system(norms, pnorms, mids);

        shape.vertices[*vertex_index].coords = new_vertex;
    }
}

fn solve_linear_system(
    norms: Vec<Vector3<f32>>,
    pnorms: Vec<Vector3<f32>>,
    mids: Vec<nalgebra::Point3<f32>>,
) -> Vector3<f32> {
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
    new_vertex
}

fn fetch_face_data(
    shape: &mut crate::geom::Shape,
    perturbed_normals: &Vec<Vector3<f32>>,
    faces: &Vec<usize>,
) -> (
    Vec<Vector3<f32>>,
    Vec<Vector3<f32>>,
    Vec<nalgebra::Point3<f32>>,
) {
    let mut norms = Vec::new();
    // these are the unperturbed normals
    let mut pnorms = Vec::new();
    // these are the perturbed normals
    let mut mids = Vec::new();
    // these are the midpoints of the faces

    // loop over faces that this vertex belongs to
    for face_index in faces {
        let face = &shape.faces[*face_index];
        norms.push(face.data().normal);
        pnorms.push(perturbed_normals[*face_index]);
        mids.push(face.data().midpoint);
    }
    (norms, pnorms, mids)
}

fn perturb_normals(
    sigma: f32,
    shape: &mut crate::geom::Shape,
    seed: Option<u64>,
) -> Vec<Vector3<f32>> {
    // this function is being updated with more mathematical rigour
    // we will first compute the theta and phi angle of the normal vector
    // then we will sample a theta and phi angle from a distribution
    // and then we will compute the new normal vector

    // Perturb the normal of each face in the shape
    let mut perturbed_normals = Vec::new();
    for face in shape.faces.iter_mut() {
        // first, we will compute the theta and phi angle of the normal vector
        let normal = face.data().normal;
        // theta is the angle between the normal and the z-axis
        // let theta = normal.z.acos();
        // phi is the angle between the normal and the x-axis
        // let phi = normal.y.atan2(normal.x);
        // now we will sample a theta and phi angle from a distribution
        // we will use a normal distribution with mean 0 and standard deviation sigma for theta
        // and a uniform distribution for phi from 0 to 2pi

        // get the theta and phi distortion angles
        let mut rng = if let Some(seed) = seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };
        let norm_dist = Normal::new(0.0, sigma).unwrap();
        let dtheta = norm_dist.sample(&mut rng);
        let dphi = rng.random_range(0.0..std::f32::consts::PI * 2.0);

        // compute theta and phi angles for original normal
        let theta = normal.z.acos();

        // get the rotation matrix to rotate z axis to the normal vector
        let fac = 1.0 - normal.z; // 1 - cos(theta)
        let x = normal.x;
        let y = normal.y;
        let z = normal.z;

        // rotation matrix using Rodrigues' rotation formula
        let rotation_matrix = Matrix3::new(
            1.0 + fac * (-y * y),
            fac * (-y * x),
            x * theta.sin(),
            fac * (-y * x),
            1.0 + fac * (-x * x),
            y * theta.sin(),
            -x * theta.sin(),
            -y * theta.sin(),
            z, // z == cos(theta)
        );

        // // assert that rotating the z axis to the normal vector gives the normal vector
        // let rotated_z = rotation_matrix * Vector3::z();
        // assert!(
        //     (rotated_z - normal).norm() < 1e-6,
        //     "Rotation matrix is incorrect"
        // );

        // rotate the perturbation vector to the normal vector
        let perturbation = Vector3::new(
            dtheta.sin() * dphi.cos(),
            dtheta.sin() * dphi.sin(),
            dtheta.cos(),
        );
        let new_normal = rotation_matrix * perturbation;

        perturbed_normals.push(new_normal);
    }
    perturbed_normals
}

fn shape_can_be_distorted(vertex_to_faces: &HashMap<usize, Vec<usize>>) -> bool {
    if vertex_to_faces.values().any(|faces| faces.len() != 3) {
        println!("Shape has vertices that do not belong to exactly 3 faces. Skipping distortion.");
        false
    } else {
        true
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
