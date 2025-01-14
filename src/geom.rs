use geo::Point;
use nalgebra as na;
use nalgebra::Point3;
use nalgebra::Vector3;
use std::env;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Sub;

use geo::Area;
use geo_clipper::Clipper;
use geo_clipper::ToOwnedPolygon;
use geo_types::Coord;
use geo_types::LineString;
use geo_types::MultiPolygon;
use geo_types::Polygon;
use macroquad::camera::Projection;

#[cfg(test)]
mod tests {

    use super::*;
    use geo_clipper::Clipper;
    use geo_types::{Coord, LineString, Polygon};

    #[test]
    fn load_hex_shape() {
        let shape = Shape::from_file("hex.obj");
        assert_eq!(shape.num_faces, 8);
        assert_eq!(shape.num_vertices, 12);
        assert_eq!(shape.faces[0].vertices[0].x, 5.0);
        assert_eq!(shape.faces[4].vertices[4].z, 5.0);
        assert_eq!(shape.faces[4].num_vertices, 6);

        let geom = Geom::from_file("hex.obj");
        assert_eq!(geom.num_shapes, 1);
        assert_eq!(geom.shapes[0].num_faces, 8);
        assert_eq!(geom.shapes[0].num_vertices, 12);
        assert_eq!(geom.shapes[0].faces[0].vertices[0].x, 5.0);
        assert_eq!(geom.shapes[0].faces[4].vertices[4].z, 5.0);
        assert_eq!(geom.shapes[0].faces[4].num_vertices, 6);
    }

    #[test]
    #[should_panic]
    fn load_multiple_shape() {
        let _ = Shape::from_file("multiple.obj");
    }

    #[test]
    fn load_multiple_geom() {
        let geom = Geom::from_file("multiple.obj");
        assert_eq!(geom.num_shapes, 2);
        assert_eq!(geom.shapes[0].num_faces, 8);
        assert_eq!(geom.shapes[0].num_vertices, 12);
        assert_eq!(geom.shapes[0].faces[0].vertices[0].x, 5.0);
        assert_eq!(geom.shapes[0].faces[4].vertices[4].z, 5.0);
        assert_eq!(geom.shapes[0].faces[4].num_vertices, 6);

        assert_eq!(geom.shapes[1].num_faces, 8);
        assert_eq!(geom.shapes[1].num_vertices, 12);
        assert_eq!(geom.shapes[1].faces[0].vertices[0].x, -7.479733);
        assert_eq!(geom.shapes[1].faces[4].vertices[4].z, 3.757555);
        assert_eq!(geom.shapes[1].faces[4].num_vertices, 6);
    }

    #[test]
    fn polygon_clip() {
        let shape = &Geom::from_file("hex2.obj").shapes[0];

        let face1 = &shape.faces[4];
        let face2 = &shape.faces[7];

        let mut exterior = Vec::new();
        for vertex in &face1.vertices {
            exterior.push(Coord {
                x: vertex.x,
                y: vertex.y,
            });
        }
        exterior.reverse();
        let subject = Polygon::new(LineString(exterior), vec![]);

        let mut exterior = Vec::new();
        for vertex in &face2.vertices {
            exterior.push(Coord {
                x: vertex.x,
                y: vertex.y,
            });
        }

        let clip = Polygon::new(LineString(exterior), vec![]);

        let result = subject.intersection(&clip, 100000.0);

        assert!(!result.0.is_empty());
    }
}

/// Represents a plane, defined by a normal and an offset value.
/// Each component of the normal corresponds to a, b, c, respectively.
/// The offset value corresponds to d.
/// The plane is then defined by `ax + by + cz + d = 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct Plane {
    pub normal: Vector3<f32>,
    pub offset: f32,
}

/// Represents a facet in a 3D surface mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Face {
    pub vertices: Vec<Point3<f32>>, // List of vertices
    pub normal: Vector3<f32>,       // Normal vector of the facet
    pub midpoint: Point3<f32>,      // Midpoint
    pub num_vertices: usize,        // Number of vertices
}

impl Face {
    pub fn new(vertices: &Vec<Point3<f32>>) -> Self {
        let vertices = vertices.clone();
        let normal = Self::compute_normal(&vertices);
        let midpoint = Self::compute_midpoint(&vertices);
        let num_vertices = vertices.len();

        Self {
            vertices,
            normal,
            midpoint,
            num_vertices,
        }
    }

    /// Compute the normal vector for the face.
    fn compute_normal(vertices: &[Point3<f32>]) -> Vector3<f32> {
        // Using the first three vertices to compute the normal.
        let v1 = &vertices[0];
        let v2 = &vertices[1];
        let v3 = &vertices[2];

        // Compute edge vectors
        let u = v2 - v1;
        let v = v3 - v1;

        // Compute cross product u × v
        let mut normal = u.cross(&v);

        normal.normalize_mut(); // normalise to unit vector

        #[cfg(debug_assertions)]
        {
            assert!(u.dot(&normal) < 0.01, "value: {}", u.dot(&normal));
            assert!(v.dot(&normal) < 0.01, "value: {}", v.dot(&normal));
        }

        normal
    }

    /// Compute the midpoint of the facet.
    fn compute_midpoint(vertices: &[Point3<f32>]) -> Point3<f32> {
        let len = vertices.len() as f32;
        // let mut mid = vertices.iter().copied();
        let mut sum: Point3<f32> = vertices
            .iter()
            .fold(Point3::origin(), |acc, point| acc + point.coords);

        sum /= len;

        sum
    }

    /// Returns the Polygon of a face, which it's 2D projection in the xy plane.
    pub fn polygon(&self) -> Polygon<f32> {
        let mut exterior = Vec::new();
        for vertex in &self.vertices {
            exterior.push(Coord {
                x: vertex.x,
                y: vertex.y,
            });
        }
        // exterior.reverse();
        Polygon::new(LineString(exterior), vec![])
    }

    /// Computes the plane containing the face.
    /// The components of the normal are a, b, and c, and the offset is d,
    /// such that ax + by + cz + d = 0
    pub fn plane(&self) -> Plane {
        Plane {
            normal: self.normal,
            offset: -self.normal.dot(&self.vertices[0].coords),
        }
    }

    /// Computes the z-distance from one facet to another.
    /// This is defined as the dot product of the position vector between
    ///     their centroids and a given projection vector.
    fn z_distance(&self, other: &Face, proj: &Vector3<f32>) -> f32 {
        let vec = &other.midpoint - &self.midpoint;
        vec.dot(&proj)
    }

    /// Computes the maximum z-distance to the vertices of another.
    /// This is defined as the lowest vertex in the subject to the highest
    /// vertex in the other.
    /// This is used to determine if any part of the other is visible along
    /// the projection direction, in which case the result is positive
    fn z_max(&self, other: &Face, proj: &Vector3<f32>) -> f32 {
        let lowest = self
            .vertices
            .iter()
            .map(|v| v.coords.dot(&proj))
            .collect::<Vec<f32>>()
            .into_iter()
            .reduce(f32::min)
            .unwrap();

        let highest = other
            .vertices
            .iter()
            .map(|v| v.coords.dot(&proj))
            .collect::<Vec<f32>>()
            .into_iter()
            .reduce(f32::max)
            .unwrap();

        highest - lowest
    }
}

/// A refractive index
#[derive(Debug, Clone, PartialEq)]
pub struct RefrIndex {
    pub real: f32,
    pub imag: f32,
}

/// Represents a 3D surface mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Shape {
    pub vertices: Vec<Point3<f32>>, // List of all vertices in the mesh
    pub num_vertices: usize,        // Number of vertices in the mesh
    pub faces: Vec<Face>,           // List of all facets in the mesh
    pub num_faces: usize,           // Number of facets in the mesh
    pub refr_index: RefrIndex,      // Refractive index of this shape
}

impl Shape {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            num_vertices: 0,
            faces: Vec::new(),
            num_faces: 0,
            refr_index: RefrIndex {
                real: 1.0,
                imag: 0.0,
            },
        }
    }

    pub fn from_file(filename: &str) -> Shape {
        let (models, _) = tobj::load_obj(filename, &tobj::LoadOptions::default())
            .expect("Failed to OBJ load file");

        let mut geom = Shape::new();

        for (i, m) in models.iter().enumerate() {
            assert!(i == 0, "found more than 1 mesh in OBJ file."); // only accept 1 mesh

            let mesh = &m.mesh;
            for vtx in 0..mesh.positions.len() / 3 {
                geom.add_vertex(Point3::new(
                    mesh.positions[3 * vtx],
                    mesh.positions[3 * vtx + 1],
                    mesh.positions[3 * vtx + 2],
                ));
            }

            let mut next_face = 0;
            for face in 0..mesh.face_arities.len() {
                let end = next_face + mesh.face_arities[face] as usize;
                let face_indices = &mesh.indices[next_face..end];
                let mut vertices = Vec::new();

                for &i in face_indices {
                    let v = geom.vertices[i as usize].xxx();
                    vertices.push(v);
                }

                let vertices = face_indices
                    .iter()
                    .map(|&i| geom.vertices[i as usize])
                    .collect();

                geom.add_face(Face::new(&vertices));

                next_face = end;
            }
        }

        geom
    }

    /// Adds a vertex to the mesh.
    pub fn add_vertex(&mut self, vertex: Point3<f32>) {
        self.vertices.push(vertex);
        self.num_vertices += 1;
    }

    /// Adds a facet to the mesh from a set of vertex indices.
    pub fn add_face(&mut self, face: Face) {
        self.faces.push(face);
        self.num_faces += 1;
    }
}

pub struct Geom {
    pub shapes: Vec<Shape>,
    pub num_shapes: usize,
}

impl Geom {
    pub fn from_file(filename: &str) -> Geom {
        match env::current_dir() {
            Ok(path) => println!("Current directory: {}", path.display()),
            Err(e) => eprintln!("Error getting current directory: {}", e),
        }
        let (models, _) = tobj::load_obj(filename, &tobj::LoadOptions::default())
            .expect("Failed to OBJ load file");

        let mut shapes = Vec::new();
        let mut num_shapes = 0;

        for (_, m) in models.iter().enumerate() {
            let mut shape = Shape::new();

            let mesh = &m.mesh;
            for vtx in 0..mesh.positions.len() / 3 {
                shape.add_vertex(Point3::new(
                    mesh.positions[3 * vtx],
                    mesh.positions[3 * vtx + 1],
                    mesh.positions[3 * vtx + 2],
                ));
            }

            let mut next_face = 0;
            for face in 0..mesh.face_arities.len() {
                let end = next_face + mesh.face_arities[face] as usize;
                let face_indices = &mesh.indices[next_face..end];
                let mut vertices = Vec::new();

                for &i in face_indices {
                    let v = shape.vertices[i as usize].clone();
                    vertices.push(v);
                }

                shape.add_face(Face::new(&vertices));

                next_face = end;
            }
            shapes.push(shape);
            num_shapes += 1;
        }

        Geom { shapes, num_shapes }
    }
}
