use geo_types::Coord;
use geo_types::LineString;
use geo_types::Polygon;
use nalgebra::Matrix4;
use nalgebra::Point3;
use nalgebra::Vector3;
use nalgebra::Vector4;
use std::path::Path;
use tobj;
use tobj::Model;

#[cfg(test)]
mod tests {

    use super::*;
    use geo_clipper::Clipper;
    use geo_types::{Coord, LineString, Polygon};

    #[test]
    fn load_hex_shape() {
        let shape = &Geom::from_file("./examples/data/hex.obj").unwrap().shapes[0];
        assert_eq!(shape.num_faces, 8);
        assert_eq!(shape.num_vertices, 12);
        assert_eq!(shape.faces[0].vertices[0].x, 5.0);
        assert_eq!(shape.faces[4].vertices[4].z, 5.0);
        assert_eq!(shape.faces[4].num_vertices, 6);

        let geom = Geom::from_file("./examples/data/hex.obj").unwrap();
        assert_eq!(geom.num_shapes, 1);
        assert_eq!(geom.shapes[0].num_faces, 8);
        assert_eq!(geom.shapes[0].num_vertices, 12);
        assert_eq!(geom.shapes[0].faces[0].vertices[0].x, 5.0);
        assert_eq!(geom.shapes[0].faces[4].vertices[4].z, 5.0);
        assert_eq!(geom.shapes[0].faces[4].num_vertices, 6);
    }

    #[test]
    fn load_multiple_geom() {
        let geom = Geom::from_file("./examples/data/multiple.obj").unwrap();
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
        let shape = &Geom::from_file("./examples/data/hex2.obj").unwrap().shapes[0];

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
    pub fn new(vertices: Vec<Point3<f32>>) -> Self {
        let vertices = vertices.clone();
        let num_vertices = vertices.len();

        let mut face = Self {
            vertices,
            num_vertices,
            normal: Vector3::zeros(),
            midpoint: Point3::origin(),
        };

        face.set_normal();
        face.set_midpoint();

        face
    }

    /// Compute the normal vector for the face.
    fn set_normal(&mut self) {
        let vertices = &self.vertices;

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

        self.normal = normal;
    }

    /// Compute the midpoint of the facet.
    fn set_midpoint(&mut self) {
        let vertices = &self.vertices;
        let len = vertices.len() as f32;
        // let mut mid = vertices.iter().copied();
        let mut sum: Point3<f32> = vertices
            .iter()
            .fold(Point3::origin(), |acc, point| acc + point.coords);

        sum /= len;

        self.midpoint = sum;
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

    pub fn vert_min(&self, dim: usize) -> Result<f32, &'static str> {
        if dim > 2 {
            return Err("Dimension must be 0, 1, or 2");
        }

        let min = self
            .vertices
            .iter()
            .map(|v| v[dim])
            .collect::<Vec<f32>>()
            .into_iter()
            .reduce(f32::min);

        match min {
            Some(val) => Ok(val),
            None => Err("No vertices found"), // Handle the case where vertices is empty
        }
    }

    pub fn vert_max(&self, dim: usize) -> Result<f32, &'static str> {
        if dim > 2 {
            return Err("Dimension must be 0, 1, or 2");
        }

        let min = self
            .vertices
            .iter()
            .map(|v| v[dim])
            .collect::<Vec<f32>>()
            .into_iter()
            .reduce(f32::max);

        match min {
            Some(val) => Ok(val),
            None => Err("No vertices found"), // Handle the case where vertices is empty
        }
    }

    /// Computes the maximum z-distance to the vertices of another.
    /// This is defined as the lowest vertex in the subject to the highest
    /// vertex in the other.
    /// This is used to determine if any part of the other is visible along
    /// the projection direction, in which case the result is positive
    pub fn z_max(&self, other: &Face, proj: &Vector3<f32>) -> f32 {
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
    /// Determines if all vertices of a Face are in front of the plane
    /// of another Face.
    pub fn is_in_front_of(&self, face: &Face) -> bool {
        let origin = face.vertices[0]; // choose point in plane of face
        for point in &self.vertices {
            let vector = point - origin;
            if vector.dot(&face.normal) > 0.05 {
                // if point is not above the plane
                return false;
            }
        }
        true
    }

    /// Transforms a Face in place using a `nalgebra` matrix transformation.
    pub fn transform(&mut self, model_view: &Matrix4<f32>) {
        for point in &mut self.vertices {
            // Iterate mutably
            let vertex4 = Vector4::new(point.x, point.y, point.z, 1.0);
            let projected_vertex = model_view * vertex4;
            point.x = projected_vertex.x;
            point.y = projected_vertex.y;
            point.z = projected_vertex.z;
        }
        self.set_midpoint();
        self.set_normal();
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

    fn from_model(model: &Model) -> Shape {
        let mesh = &model.mesh;

        let vertices = mesh
            .positions
            .chunks_exact(3)
            .map(|v| Point3::new(v[0], v[1], v[2]))
            .collect::<Vec<_>>();

        let mut shape = Shape::new();
        shape.num_vertices = vertices.len();
        shape.vertices = vertices;

        let face_arities = if mesh.face_arities.is_empty() {
            vec![3; mesh.indices.len() / 3]
        } else {
            mesh.face_arities.clone()
        };

        let mut next_face = 0;
        for arity in face_arities {
            let end = next_face + arity as usize;
            let face_indices = &mesh.indices[next_face..end];

            let face_vertices: Vec<_> = face_indices
                .iter()
                .map(|&i| shape.vertices[i as usize])
                .collect();
            shape.add_face(Face::new(face_vertices));

            next_face = end;
        }
        shape
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

    pub fn transform(&mut self, transform: &Matrix4<f32>) {
        for face in &mut self.faces {
            // Iterate mutably
            face.transform(transform); // Call the in-place project method
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Geom {
    pub shapes: Vec<Shape>,
    pub num_shapes: usize,
}

impl Geom {
    pub fn from_file(filename: &str) -> Result<Self, String> {
        // Log current directory only in debug builds
        #[cfg(debug_assertions)]
        match std::env::current_dir() {
            Ok(path) => println!("Current directory: {}", path.display()),
            Err(e) => eprintln!("Error getting current directory: {}", e),
        }

        let path = Path::new(filename);
        let resolved_filename = if path.is_absolute() {
            filename.to_string()
        } else {
            std::env::current_dir()
                .map(|p| p.join(path).display().to_string())
                .map_err(|e| format!("Could not resolve path: {}", e))?
        };

        let (models, _) = tobj::load_obj(&resolved_filename, &tobj::LoadOptions::default())
            .map_err(|e| format!("Failed to load OBJ file '{}': {}", filename, e))?;

        if models.is_empty() {
            return Err("No models found in OBJ file".to_string());
        }

        let shapes: Vec<Shape> = models.iter().map(Shape::from_model).collect();

        Ok(Self {
            num_shapes: shapes.len(),
            shapes,
        })
    }

    pub fn transform(&mut self, transform: &Matrix4<f32>) {
        for shape in &mut self.shapes {
            shape.transform(transform);
        }
    }
}
