use geo::Area;
use geo_types::Coord;
use geo_types::LineString;
use geo_types::Polygon;
use nalgebra::Matrix4;
use nalgebra::Point3;
use nalgebra::Vector3;
use nalgebra::Vector4;
use std::cell::Ref;
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
        match &shape.faces[0] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, 5.0);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
        match &shape.faces[4] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, 2.5);
                assert_eq!(data.exterior[4].z, 5.0);
                assert_eq!(data.num_vertices, 6);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }

        let geom = Geom::from_file("./examples/data/hex.obj").unwrap();
        assert_eq!(geom.num_shapes, 1);
        assert_eq!(geom.shapes[0].num_faces, 8);
        assert_eq!(geom.shapes[0].num_vertices, 12);
        match &geom.shapes[0].faces[0] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, 5.0);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
        match &geom.shapes[0].faces[4] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, 2.5);
                assert_eq!(data.exterior[4].z, 5.0);
                assert_eq!(data.num_vertices, 6);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
    }

    #[test]
    fn load_multiple_geom() {
        let geom = Geom::from_file("./examples/data/multiple.obj").unwrap();
        assert_eq!(geom.num_shapes, 2);
        assert_eq!(geom.shapes[0].num_faces, 8);
        assert_eq!(geom.shapes[0].num_vertices, 12);
        match &geom.shapes[0].faces[4] {
            Face::Simple(data) => {
                assert_eq!(data.num_vertices, 6);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }

        assert_eq!(geom.shapes[1].num_faces, 8);
        assert_eq!(geom.shapes[1].num_vertices, 12);
        match &geom.shapes[1].faces[4] {
            Face::Simple(data) => {
                assert_eq!(data.num_vertices, 6);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
    }

    #[test]
    fn polygon_clip() {
        let shape = &Geom::from_file("./examples/data/hex2.obj").unwrap().shapes[0];

        let face1 = &shape.faces[4];
        let face2 = &shape.faces[7];

        let mut exterior = Vec::new();
        match face1 {
            Face::Simple(data) => {
                for vertex in &data.exterior {
                    exterior.push(Coord {
                        x: vertex.x,
                        y: vertex.y,
                    });
                }
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
        exterior.reverse();
        let subject = Polygon::new(LineString(exterior), vec![]);

        let mut exterior = Vec::new();
        match face2 {
            Face::Simple(data) => {
                for vertex in &data.exterior {
                    exterior.push(Coord {
                        x: vertex.x,
                        y: vertex.y,
                    });
                }
            }
            Face::Complex { .. } => {
                panic!();
            }
        }

        let clip = Polygon::new(LineString(exterior), vec![]);

        let result = subject.intersection(&clip, 100000.0);

        assert!(!result.0.is_empty());
    }
}

trait Point3Extensions {
    fn transform(&mut self, model_view: &Matrix4<f32>);
    fn to_xy(&self) -> Coord<f32>;
}

impl Point3Extensions for Point3<f32> {
    /// Transforms a Point3 type to another coordinate system.
    fn transform(&mut self, model_view: &Matrix4<f32>) {
        let vertex4 = Vector4::new(self.x, self.y, self.z, 1.0);
        let projected_vertex = model_view * vertex4;
        self.x = projected_vertex.x;
        self.y = projected_vertex.y;
        self.z = projected_vertex.z;
    }
    fn to_xy(&self) -> Coord<f32> {
        Coord {
            x: self.x,
            y: self.y,
        }
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

/// Represents a closed line of exterior points of a polygon 3D.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceData {
    pub exterior: Vec<Point3<f32>>, // List of exterior vertices
    pub normal: Vector3<f32>,       // Normal vector of the facet
    pub midpoint: Point3<f32>,      // Midpoint
    pub num_vertices: usize,        // Number of vertices
    pub area: Option<f32>,          // Unsigned area
    pub parent_id: Option<usize>,   // An optional parent shape id number
}

impl FaceData {
    pub fn new(vertices: Vec<Point3<f32>>, parent_id: Option<usize>) -> Self {
        let vertices = vertices.clone();
        let num_vertices = vertices.len();

        let mut face = Self {
            exterior: vertices,
            num_vertices,
            normal: Vector3::zeros(),
            midpoint: Point3::origin(),
            area: None, // compute as needed
            parent_id,
        };

        face.set_normal();
        face.set_midpoint();

        face
    }

    /// Compute the normal vector for the face.
    fn set_normal(&mut self) {
        let vertices = &self.exterior;

        // Using the first three vertices to compute the normal.
        let v1 = &vertices[0];
        let v2 = &vertices[1];
        let v3 = &vertices[2];

        // Compute edge vectors
        let u = v2 - v1;
        let v = v3 - v2;

        // Compute cross product u × v
        let mut normal = u.cross(&v);

        normal.normalize_mut(); // normalise to unit vector

        #[cfg(debug_assertions)]
        {
            assert!(
                u.dot(&normal) < 0.01,
                "u: {:?}, v: {:?}, dot: {}, verts: {:?}",
                u,
                v,
                u.dot(&normal),
                vertices
            );
            assert!(
                v.dot(&normal) < 0.01,
                "u: {:?}, v: {:?}, dot: {}, verts: {:?}",
                u,
                v,
                v.dot(&normal),
                vertices
            );
        }

        self.normal = normal;
    }

    /// Compute the midpoint of the facet.
    fn set_midpoint(&mut self) {
        let vertices = &self.exterior;
        let len = vertices.len() as f32;
        // let mut mid = vertices.iter().copied();
        let mut sum: Point3<f32> = vertices
            .iter()
            .fold(Point3::origin(), |acc, point| acc + point.coords);

        sum /= len;

        self.midpoint = sum;
    }

    /// Computes the plane containing the face.
    /// The components of the normal are a, b, and c, and the offset is d,
    /// such that ax + by + cz + d = 0
    pub fn plane(&self) -> Plane {
        Plane {
            normal: self.normal,
            offset: -self.normal.dot(&self.exterior[0].coords),
        }
    }

    /// Computes the z-distance from one facet to another.
    /// This is defined as the dot product of the position vector between
    ///     their centroids and a given projection vector.
    #[allow(dead_code)]
    fn z_distance(&self, other: &FaceData, proj: &Vector3<f32>) -> f32 {
        let vec = &other.midpoint - &self.midpoint;
        vec.dot(&proj)
    }

    /// Returns the minimum value of the vertices in a `FaceData` along the
    /// specified dimension.
    pub fn vert_min(&self, dim: usize) -> Result<f32, &'static str> {
        if dim > 2 {
            return Err("Dimension must be 0, 1, or 2");
        }

        let min = self
            .exterior
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

    /// Returns the maximum value of the vertices in a `FaceData` along the
    /// specified dimension.
    pub fn vert_max(&self, dim: usize) -> Result<f32, &'static str> {
        if dim > 2 {
            return Err("Dimension must be 0, 1, or 2");
        }

        let min = self
            .exterior
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
    pub fn z_max(&self, other: &FaceData, proj: &Vector3<f32>) -> f32 {
        let lowest = self
            .exterior
            .iter()
            .map(|v| v.coords.dot(&proj))
            .collect::<Vec<f32>>()
            .into_iter()
            .reduce(f32::min)
            .unwrap();

        let highest = other
            .exterior
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
    pub fn is_in_front_of(&self, face: &FaceData) -> bool {
        let origin = face.exterior[0]; // choose point in plane of face
        for point in &self.exterior {
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
        for point in &mut self.exterior {
            point.transform(model_view);
        }
        self.set_midpoint();
        self.set_normal();
    }
}

/// An enum for 2 different types of polygon in 3D.
/// `Face::Simple` represents a polygon with only exterior vertices.
/// `Face::Complex` represents a polygon that may also contain interior vertices (holes).

#[derive(Debug, Clone, PartialEq)]
pub enum Face {
    Simple(FaceData),
    Complex {
        data: FaceData,
        interiors: Vec<Vec<Point3<f32>>>,
    },
}

impl Face {
    pub fn new_simple(exterior: Vec<Point3<f32>>, parent_id: Option<usize>) -> Self {
        Face::Simple(FaceData::new(exterior, parent_id))
    }

    pub fn new_complex(
        exterior: Vec<Point3<f32>>,
        interiors: Vec<Vec<Point3<f32>>>,
        parent_id: Option<usize>,
    ) -> Self {
        Face::Complex {
            data: FaceData::new(exterior, parent_id),
            interiors,
        }
    }

    /// Transform a `Face` to another coordinate system.
    pub fn transform(&mut self, model_view: &Matrix4<f32>) {
        match self {
            Face::Simple(data) => data.transform(model_view),
            Face::Complex { data, interiors } => {
                data.transform(model_view);
                for interior in interiors {
                    for point in interior {
                        point.transform(model_view);
                    }
                }
            }
        }
    }

    pub fn midpoint(&self) -> Point3<f32> {
        match self {
            Face::Simple(data) => data.midpoint,
            Face::Complex { data, .. } => data.midpoint,
        }
    }

    pub fn to_polygon(&self) -> Polygon<f32> {
        match self {
            Face::Simple(data) => {
                let mut exterior = Vec::new();
                for vertex in &data.exterior {
                    exterior.push(vertex.to_xy());
                }
                // exterior.reverse();
                Polygon::new(LineString(exterior), vec![])
            }
            Face::Complex { data, interiors } => {
                let mut exterior = Vec::new();
                for vertex in &data.exterior {
                    exterior.push(vertex.to_xy());
                }
                // exterior.reverse();
                let mut holes = Vec::new();
                for interior in interiors {
                    let mut hole = Vec::new();
                    for vertex in interior {
                        hole.push(vertex.to_xy());
                    }
                    holes.push(LineString(hole));
                }
                Polygon::new(LineString(exterior), holes)
            }
        }
    }

    pub fn plane(&self) -> Plane {
        match self {
            Face::Simple(data) => data.plane(),
            Face::Complex { data, .. } => data.plane(),
        }
    }

    /// Setter for the area of a `Face`.
    pub fn set_area(&mut self, area: f32) {
        match self {
            Face::Simple(data) => data.area = Some(area),
            Face::Complex { data, .. } => data.area = Some(area),
        }
    }

    /// Creates a `Face` struct from a `Polygon`
    #[allow(dead_code)]
    fn from_polygon(polygon: &Polygon<f32>) -> Self {
        // do the exterior
        let mut exterior = Vec::new();
        for coord in polygon
            .exterior()
            .0
            .iter()
            .take(polygon.exterior().0.len() - 1)
        {
            exterior.push(Point3::new(coord.x, coord.y, 0.0));
        }

        if polygon.interiors().is_empty() {
            let mut face = Face::new_simple(exterior, None);
            match face {
                Face::Simple(ref mut data) => data.area = Some(polygon.unsigned_area()),
                Face::Complex { .. } => {}
            }
            face
        } else {
            let mut interiors = Vec::new();
            for interior in polygon.interiors() {
                let mut vertices = Vec::new();
                for coord in interior.0.iter().take(interior.0.len() - 1) {
                    vertices.push(Point3::new(coord.x, coord.y, 0.0));
                }
                interiors.push(vertices);
            }
            let mut face = Face::new_complex(exterior, interiors, None);
            match face {
                Face::Simple { .. } => {}
                Face::Complex { ref mut data, .. } => data.area = Some(polygon.unsigned_area()),
            }
            face
        }
    }

    pub fn data(&self) -> &FaceData {
        match self {
            Face::Simple(data) => data,
            Face::Complex { data, .. } => data,
        }
    }

    pub fn data_mut(&mut self) -> &mut FaceData {
        match self {
            Face::Simple(data) => data,
            Face::Complex { data, .. } => data,
        }
    }
}

/// A refractive index
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefrIndex {
    pub real: f32,
    pub imag: f32,
}

impl RefrIndex {
    pub fn new(real: f32, imag: f32) -> Self {
        Self { real, imag }
    }
}

/// Represents a 3D surface mesh.
#[derive(Debug, PartialEq)]
pub struct Shape {
    pub vertices: Vec<Point3<f32>>, // List of all vertices in the mesh
    pub num_vertices: usize,        // Number of vertices in the mesh
    pub faces: Vec<Face>,           // List of all facets in the mesh
    pub num_faces: usize,           // Number of facets in the mesh
    pub refr_index: RefrIndex,      // Refractive index of this shape
    pub id: Option<usize>,          // an id number
    pub parent_id: Option<usize>,   // An optional parent shape index, which encompasses this one
}

impl Shape {
    pub fn new(id: Option<usize>, parent_id: Option<usize>) -> Self {
        Self {
            vertices: Vec::new(),
            num_vertices: 0,
            faces: Vec::new(),
            num_faces: 0,
            refr_index: RefrIndex {
                real: 1.31, // some default settings, temporary
                imag: 0.0,
            },
            id,
            parent_id,
        }
    }

    fn from_model(model: Model, id: Option<usize>) -> Shape {
        let mesh = &model.mesh;

        let vertices = mesh
            .positions
            .chunks_exact(3)
            .map(|v| Point3::new(v[0], v[1], v[2]))
            .collect::<Vec<_>>();

        let mut shape = Shape::new(id, None);
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
            shape.add_face(Face::new_simple(face_vertices, id));

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

#[derive(Debug, PartialEq)]
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

        let shapes: Vec<Shape> = models
            .into_iter()
            .enumerate()
            .map(|(i, model)| Shape::from_model(model, Some(i)))
            .collect();

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

/// An intersection in the context of a `Geom`.
#[derive(Debug, PartialEq)]
pub struct Intersection {
    pub face: Face,
    pub shape_id: usize,
}

impl Intersection {
    /// A new `Intersection`.
    pub fn new(face: Face, shape_id: usize) -> Self {
        Self { face, shape_id }
    }
}
