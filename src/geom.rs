//! Geometric representation and operations for electromagnetic scattering simulations.
//!
//! This module provides the core geometric foundation for the GOAD (Geometric Optics
//! with Aperture Diffraction) method. It handles 3D shape representation, geometric
//! transformations, containment relationships, and mesh operations required for
//! electromagnetic beam propagation through complex particle geometries.
//!
//! The geometric representation supports:
//! - Complex polyhedral particles with arbitrary shapes
//! - Nested geometries with inclusions and coatings
//! - Material properties attached to geometric regions  
//! - Coordinate transformations for beam alignment
//! - Mesh validation and geometric quality checks
//!
//! # Key Components
//!
//! - [`Shape`]: Individual closed 3D surfaces with material properties
//! - [`Geom`]: Collections of shapes with containment relationships
//! - [`Face`]: Polygonal surface elements for beam-surface interactions
//! - Geometric transformations and scaling operations
//! - OBJ file loading and mesh processing

use crate::containment::{ContainmentGraph, AABB};
use crate::orientation::*;
use crate::settings;
use anyhow::Result;
use geo::{Area, TriangulateEarcut};
use geo_types::{Coord, LineString, Polygon};
use nalgebra::{self as na, Complex, Isometry3, Matrix4, Point3, Vector3, Vector4};
use pyo3::prelude::*;
use std::path::Path;
use tobj::{self, Model};

#[cfg(test)]
mod tests {

    use super::*;
    use geo_clipper::Clipper;
    use geo_types::{Coord, LineString, Polygon};

    #[test]
    fn earcut_xy() {
        let mut geom = Geom::from_file("./examples/data/plane_xy.obj").unwrap();

        let face = geom.shapes[0].faces.remove(0);
        assert_eq!(face.data().exterior.len(), 4);
        assert_eq!(face.data().exterior[0], Point3::new(1.0, 1.0, 0.0));

        let triangles = Face::earcut(&face);

        assert_eq!(triangles.len(), 2);
        assert_eq!(triangles[0].data().normal, face.data().normal);
    }

    #[test]
    fn earcut_zy() {
        let mut geom = Geom::from_file("./examples/data/plane_yz.obj").unwrap();

        let face = geom.shapes[0].faces.remove(0);
        assert_eq!(face.data().exterior.len(), 4);
        assert_eq!(face.data().exterior[0], Point3::new(0.0, 1.0, 1.0));

        let triangles = Face::earcut(&face);

        assert_eq!(triangles.len(), 2);
        assert_eq!(triangles[0].data().normal, face.data().normal);
    }

    #[test]
    fn rescale_hex() {
        let mut geom = Geom::from_file("./examples/data/hex2.obj").unwrap();
        let x_dim = geom.shapes[0].aabb.as_ref().unwrap().max.x
            - geom.shapes[0].aabb.as_ref().unwrap().min.x;

        geom.shapes[0].rescale(0.5);
        let rescaled_x_dim = geom.shapes[0].aabb.as_ref().unwrap().max.x
            - geom.shapes[0].aabb.as_ref().unwrap().min.x;

        assert!((rescaled_x_dim / x_dim - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_com() {
        let geom = Geom::from_file("./examples/data/hex.obj").unwrap();
        let com = geom.centre_of_mass();
        println!("{:?}", com);
        assert!(com.coords.norm() < 1e-6);
        assert!(geom.is_centered());

        let geom = Geom::from_file("./examples/data/multiple2.obj").unwrap();
        let com = geom.centre_of_mass();
        println!("{:?}", com);
        assert!(com.coords.norm() < 1e-6);
        assert!(geom.is_centered());

        let geom = Geom::from_file("./examples/data/multiple3.obj").unwrap();
        let com = geom.centre_of_mass();
        println!("{:?}", com);
        assert!(com.coords.norm() - 5.0 < 1e-6);
        assert!(com.y.abs() - 5.0 < 1e-6);
        assert!(com.z.abs() < 1e-6);
        assert!(!geom.is_centered());

        let mut recentred = geom.clone();
        recentred.recentre();
        let com = recentred.centre_of_mass();
        println!("{:?}", com);
        assert!(com.coords.norm() < 1e-6);
        assert!(recentred.is_centered());
    }

    #[test]
    fn load_hex_shape() {
        let shape = &Geom::from_file("./examples/data/hex.obj").unwrap().shapes[0];
        assert_eq!(shape.num_faces, 8);
        assert_eq!(shape.num_vertices, 12);
        match &shape.faces[0] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, -0.0);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
        match &shape.faces[4] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, -4.330127);
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
                assert_eq!(data.exterior[0].x, -0.0);
            }
            Face::Complex { .. } => {
                panic!();
            }
        }
        match &geom.shapes[0].faces[4] {
            Face::Simple(data) => {
                assert_eq!(data.exterior[0].x, -4.330127);
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

    #[test]
    fn shape_within() {
        let geom = &Geom::from_file("./examples/data/cubes.obj").unwrap();

        assert_eq!(geom.num_shapes, 6);
        assert!(geom.shapes[1].is_within(&geom, Some(0)));
        assert!(!geom.shapes[2].is_within(&geom, Some(1)));
        assert!(!geom.shapes[1].is_within(&geom, Some(2)));
        assert!(geom.shapes[3].is_within(&geom, Some(0)));
        assert!(geom.shapes[3].is_within(&geom, Some(1)));
        assert!(!geom.shapes[4].is_within(&geom, Some(0)));
        assert!(geom.shapes[5].is_within(&geom, Some(3)));
    }
}

trait Coord3Extensions {
    fn projected_z(&self, plane: &Plane) -> f32;
}
impl Coord3Extensions for Coord<f32> {
    /// Returns the z-coordinate of a `Coord` projected onto a plane in 3D
    fn projected_z(&self, plane: &Plane) -> f32 {
        -(plane.normal.x * self.x + plane.normal.y * self.y + plane.offset) / plane.normal.z
    }
}

pub trait PolygonExtensions {
    fn project(&self, plane: &Plane) -> Result<Face>;
}

impl PolygonExtensions for Polygon<f32> {
    /// Projects the xy coordinates of a polygon onto a plane in 3D
    ///  the last vertex, which is a duplicate of the first
    fn project(&self, plane: &Plane) -> Result<Face> {
        let area = self.unsigned_area() / plane.normal.z.abs();

        // condition to enforce that all normals point outwards,
        // assuming the initial planes were correctly oriented
        let reverse = if plane.normal.z < 0.0 { true } else { false };

        let project_coords = |coords: &Vec<Coord<f32>>| -> Vec<Point3<f32>> {
            coords
                .iter()
                .take(coords.len() - 1)
                .map(|coord| Point3::new(coord.x, coord.y, coord.projected_z(plane)))
                .collect()
        };

        let mut exterior = project_coords(&self.exterior().0);
        if reverse {
            exterior.reverse()
        }

        if self.interiors().is_empty() {
            let mut face = Face::new_simple(exterior, None, None)?;
            face.set_area(area);
            Ok(face)
        } else {
            let mut interiors: Vec<_> = self
                .interiors()
                .iter()
                .rev()
                .map(|interior| project_coords(&interior.0))
                .collect();
            if reverse {
                interiors.iter_mut().for_each(|interior| interior.reverse());
            }
            let mut face = Face::new_complex(exterior, interiors, None)?;
            face.set_area(area);
            Ok(face)
        }
    }
}

trait Point3Extensions {
    fn transform(&mut self, model_view: &Matrix4<f32>) -> Result<()>;
    fn to_xy(&self) -> Coord<f32>;
}

impl Point3Extensions for Point3<f32> {
    /// Transforms a Point3 type to another coordinate system.
    fn transform(&mut self, model_view: &Matrix4<f32>) -> Result<()> {
        let vertex4 = Vector4::new(self.x, self.y, self.z, 1.0);
        let projected_vertex = model_view * vertex4;
        self.x = projected_vertex.x;
        self.y = projected_vertex.y;
        self.z = projected_vertex.z;

        Ok(())
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
    pub exterior: Vec<Point3<f32>>,           // List of exterior vertices
    pub exterior_indices: Option<Vec<usize>>, // List of exterior vertex indices
    pub normal: Vector3<f32>,                 // Normal vector of the facet
    pub midpoint: Point3<f32>,                // Midpoint
    pub num_vertices: usize,                  // Number of vertices
    pub area: Option<f32>,                    // Unsigned area
    pub shape_id: Option<usize>,              // An optional parent shape id number
}

impl FaceData {
    pub fn new(
        vertices: Vec<Point3<f32>>,
        shape_id: Option<usize>,
        indices: Option<Vec<usize>>,
    ) -> Result<Self> {
        let vertices = vertices.clone();
        let num_vertices = vertices.len();

        let mut face = Self {
            exterior: vertices,
            exterior_indices: indices,
            num_vertices,
            normal: Vector3::zeros(),
            midpoint: Point3::origin(),
            area: None, // compute as needed
            shape_id,
        };

        face.set_midpoint();
        face.set_normal()?; // midpoint should be set first

        Ok(face)
    }

    /// Computes the normal vector for the face.
    /// 
    /// **Context**: Face normals are required for lighting calculations,
    /// ray-surface intersection tests, and determining face orientation
    /// in electromagnetic field computations.
    /// 
    /// **How it Works**: Finds two vertices with sufficient separation,
    /// uses the face midpoint as a third point, then computes the cross
    /// product to get the normal vector. Normalizes the result and verifies
    /// orthogonality to the face plane.
    fn set_normal(&mut self) -> Result<()> {
        let vertices = &self.exterior;

        if vertices.len() < 2 {
            return Err(anyhow::anyhow!(
                "Not enough vertices to compute the normal."
            ));
        }

        // Find a pair of vertices with a distance greater than the threshold
        let mut v1 = None;
        let mut v2 = None;
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                if (vertices[j] - vertices[i]).magnitude() > settings::VEC_LENGTH_THRESHOLD {
                    v1 = Some(&vertices[i]);
                    v2 = Some(&vertices[j]);
                    break;
                }
            }
            if v1.is_some() && v2.is_some() {
                break;
            }
        }

        // Return an error if no suitable pair is found
        let v1 = v1.ok_or_else(|| {
            anyhow::anyhow!("No vertex pair found with a distance greater than the threshold.")
        })?;
        let v2 = v2.ok_or_else(|| {
            anyhow::anyhow!("No vertex pair found with a distance greater than the threshold.")
        })?;

        let v3 = self.midpoint;

        // Compute edge vectors
        let u = v2 - v1;
        let v = v3 - v1;

        // Compute the cross product
        let mut normal = u.cross(&v);

        if normal.magnitude() == 0.0 {
            return Err(anyhow::anyhow!(
                "Degenerate face detected; the cross product is zero. u: {u}, v: {v}"
            ));
        }

        normal.normalize_mut();

        // Verify the normal
        if u.dot(&normal).abs() < 0.01 && v.dot(&normal).abs() < 0.01 {
            self.normal = normal;
            Ok(())
        } else {
            Err(anyhow::anyhow!(
                "Normal could not be computed correctly. u: {u}, v: {v}, face: {:?}",
                self
            ))
        }
    }

    /// Computes the geometric center of the face.
    /// 
    /// **Context**: Face midpoints serve as reference points for normal
    /// vector calculations and provide representative positions for
    /// face-based operations in the electromagnetic solver.
    /// 
    /// **How it Works**: Averages the coordinates of all exterior vertices
    /// to find the centroid of the face polygon.
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

    /// Computes the plane equation containing this face.
    /// 
    /// **Context**: Many geometric algorithms require the plane equation
    /// (ax + by + cz + d = 0) for intersection tests, distance calculations,
    /// and coordinate transformations.
    /// 
    /// **How it Works**: Uses the face normal as the plane normal (a,b,c)
    /// and computes the offset d from the dot product of the normal with
    /// any point on the face. The components of the normal are a, b, and c,
    /// and the offset is d, such that ax + by + cz + d = 0.
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
    pub fn vert_min(&self, dim: usize) -> Result<f32> {
        if dim > 2 {
            return Err(anyhow::anyhow!("Dimension must be 0, 1, or 2"));
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
            None => Err(anyhow::anyhow!("No vertices found")), // Handle the case where vertices is empty
        }
    }

    /// Returns the maximum value of the vertices in a `FaceData` along the
    /// specified dimension.
    pub fn vert_max(&self, dim: usize) -> Result<f32> {
        if dim > 2 {
            return Err(anyhow::anyhow!("Dimension must be 0, 1, or 2"));
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
            None => Err(anyhow::anyhow!("No vertices found")), // Handle the case where vertices is empty
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
    /// Tests if all vertices of this face lie in front of another face's plane.
    /// 
    /// **Context**: Visibility and occlusion calculations in ray tracing require
    /// determining the spatial relationship between faces. This test helps
    /// establish which faces can potentially obstruct others.
    /// 
    /// **How it Works**: Computes the signed distance from each vertex to the
    /// other face's plane using dot products. Returns true only if all vertices
    /// are on the positive side of the plane (in front).
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

    /// Applies a transformation matrix to the face geometry.
    /// 
    /// **Context**: Coordinate system transformations are required for aligning
    /// geometry with beam propagation directions and for various geometric
    /// operations in the electromagnetic solver.
    /// 
    /// **How it Works**: Transforms all exterior vertices using the 4x4 matrix,
    /// then recomputes the face midpoint and normal vector to maintain geometric
    /// consistency after the transformation.
    pub fn transform(&mut self, model_view: &Matrix4<f32>) -> Result<()> {
        for point in &mut self.exterior {
            point.transform(model_view)?;
        }
        self.set_midpoint();
        self.set_normal()
    }

    /// Checks if the face polygon intersects itself.
    /// 
    /// **Context**: Self-intersecting faces can cause problems in electromagnetic
    /// field calculations and ray tracing algorithms. Detecting these geometric
    /// anomalies helps ensure mesh quality.
    /// 
    /// **How it Works**: Tests all pairs of non-adjacent edges for 3D intersection
    /// using parametric line segment intersection tests. Returns true if any
    /// edge pair intersects within their segments.
    pub fn self_intersects(&self) -> bool {
        let vertices = &self.exterior;
        let n = vertices.len();

        // Need at least 4 vertices for self-intersection
        if n < 4 {
            return false;
        }

        // Check each pair of non-adjacent edges for intersection
        for i in 0..n {
            let i_next = (i + 1) % n;
            let edge1_start = &vertices[i];
            let edge1_end = &vertices[i_next];

            // Start j at i+2 to avoid adjacent edges
            for j in (i + 2)..n {
                // Skip if this would create adjacent edges
                if j == n - 1 && i == 0 {
                    continue;
                }

                let j_next = (j + 1) % n;
                // Also skip if edges would be adjacent
                if j_next == i {
                    continue;
                }

                let edge2_start = &vertices[j];
                let edge2_end = &vertices[j_next];

                // Check if these two edges intersect in 3D
                if Self::segments_intersect_3d(edge1_start, edge1_end, edge2_start, edge2_end) {
                    return true;
                }
            }
        }

        false
    }

    // Helper function to check if two line segments intersect in 3D
    fn segments_intersect_3d(
        p1: &Point3<f32>,
        p2: &Point3<f32>,
        p3: &Point3<f32>,
        p4: &Point3<f32>,
    ) -> bool {
        // Convert points to vectors for calculations
        let p13 = p3 - p1;
        let p43 = p3 - p4;
        let p21 = p1 - p2;

        // Check if lines are parallel
        let d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
        let d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
        let d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
        let d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
        let d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

        let denom = d2121 * d4343 - d4321 * d4321;

        // If denominator is close to 0, lines are parallel
        if denom.abs() < settings::COLINEAR_THRESHOLD {
            return false;
        }

        let numer = d1343 * d4321 - d1321 * d4343;

        let mua = numer / denom;
        let mub = (d1343 + d4321 * mua) / d4343;

        // Check if intersection point is within both line segments
        if mua >= 0.0 && mua <= 1.0 && mub >= 0.0 && mub <= 1.0 {
            // Calculate intersection point
            let pa = Point3::new(
                p1.x + mua * (p2.x - p1.x),
                p1.y + mua * (p2.y - p1.y),
                p1.z + mua * (p2.z - p1.z),
            );

            let pb = Point3::new(
                p3.x + mub * (p4.x - p3.x),
                p3.y + mub * (p4.y - p3.y),
                p3.z + mub * (p4.z - p3.z),
            );

            // Check if intersection points are close enough
            let dist = (pa - pb).norm();
            return dist < settings::VEC_LENGTH_THRESHOLD;
        }

        false
    }

    /// Determines if the face represents a convex polygon.
    /// 
    /// **Context**: Some algorithms work more efficiently with convex polygons,
    /// and triangulation methods may behave differently for convex vs concave
    /// faces. This provides a geometric classification of the face.
    /// 
    /// **How it Works**: Projects the face onto its best-fitting 2D plane
    /// based on the largest component of the face normal, then checks if all 
    /// cross products of consecutive edges have the same sign, which indicates
    /// convexity.
    pub fn is_convex(&self) -> bool {
        let vertices = &self.exterior;
        let n = vertices.len();

        // Need at least 3 vertices for a valid polygon
        if n < 3 {
            return true; // Technically, a line or point is trivially convex
        }

        // Find the best projection plane based on the face normal
        let abs_normal = Vector3::new(
            self.normal.x.abs(),
            self.normal.y.abs(),
            self.normal.z.abs(),
        );

        // We'll project the face onto the plane where the normal has the largest component
        let (i1, i2) = if abs_normal.x >= abs_normal.y && abs_normal.x >= abs_normal.z {
            // Project onto YZ plane
            (1, 2) // Use Y and Z coordinates
        } else if abs_normal.y >= abs_normal.x && abs_normal.y >= abs_normal.z {
            // Project onto XZ plane
            (0, 2) // Use X and Z coordinates
        } else {
            // Project onto XY plane
            (0, 1) // Use X and Y coordinates
        };

        // Check convexity using the cross product method
        // A polygon is convex if all cross products of consecutive edges have the same sign
        let mut sign = 0; // 0 = uninitialized, 1 = positive, -1 = negative

        for i in 0..n {
            let p1 = &vertices[i];
            let p2 = &vertices[(i + 1) % n];
            let p3 = &vertices[(i + 2) % n];

            // Form 2D vectors for two consecutive edges
            let v1 = [p2[i1] - p1[i1], p2[i2] - p1[i2]];
            let v2 = [p3[i1] - p2[i1], p3[i2] - p2[i2]];

            // Compute the 2D cross product
            let cross = v1[0] * v2[1] - v1[1] * v2[0];

            // If cross product is close to zero, these points are collinear
            if cross.abs() < settings::COLINEAR_THRESHOLD {
                continue;
            }

            // Initialize sign with first non-zero cross product
            if sign == 0 {
                sign = if cross > 0.0 { 1 } else { -1 };
            } else if (cross > 0.0 && sign < 0) || (cross < 0.0 && sign > 0) {
                // If sign changes, the polygon is not convex
                return false;
            }
        }

        // If we reach here, the polygon is convex
        true
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
    pub fn new_simple(
        exterior: Vec<Point3<f32>>,
        parent_id: Option<usize>,
        indices: Option<Vec<usize>>,
    ) -> Result<Self> {
        Ok(Face::Simple(FaceData::new(exterior, parent_id, indices)?))
    }

    pub fn new_complex(
        exterior: Vec<Point3<f32>>,
        interiors: Vec<Vec<Point3<f32>>>,
        parent_id: Option<usize>,
    ) -> Result<Self> {
        Ok(Face::Complex {
            data: FaceData::new(exterior, parent_id, None)?,
            interiors,
        })
    }

    /// Transform a `Face` to another coordinate system.
    pub fn transform(&mut self, model_view: &Matrix4<f32>) -> Result<()> {
        match self {
            Face::Simple(data) => data.transform(model_view),
            Face::Complex { data, interiors } => {
                data.transform(model_view)?;

                for interior in interiors {
                    for point in interior {
                        point.transform(model_view)?;
                    }
                }
                Ok(())
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

    // /// Creates a `Face` struct from a `Polygon`
    // fn from_polygon(polygon: &Polygon<f32>) -> Result<Self> {
    //     // do the exterior
    //     let mut exterior = Vec::new();
    //     for coord in polygon
    //         .exterior()
    //         .0
    //         .iter()
    //         .take(polygon.exterior().0.len() - 1)
    //     {
    //         exterior.push(Point3::new(coord.x, coord.y, 0.0));
    //     }

    //     if polygon.interiors().is_empty() {
    //         let mut face = Face::new_simple(exterior, None)?;
    //         if let Face::Simple(ref mut data) = face {
    //             data.area = Some(polygon.unsigned_area());
    //         }
    //         Ok(face)
    //     } else {
    //         let mut interiors = Vec::new();
    //         for interior in polygon.interiors() {
    //             let mut vertices = Vec::new();
    //             for coord in interior.0.iter().take(interior.0.len() - 1) {
    //                 vertices.push(Point3::new(coord.x, coord.y, 0.0));
    //             }
    //             interiors.push(vertices);
    //         }
    //         let mut face = Face::new_complex(exterior, interiors, None)?;
    //         if let Face::Complex { ref mut data, .. } = face {
    //             data.area = Some(polygon.unsigned_area());
    //         }
    //         Ok(face)
    //     }
    // }

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

    /// Triangulates a face using the ear clipping algorithm.
    /// 
    /// **Context**: Many algorithms in computational geometry work more efficiently
    /// with triangular faces rather than arbitrary polygons. Triangulation breaks
    /// complex faces into simpler triangular elements while preserving the surface.
    /// 
    /// **How it Works**: Projects the face onto a 2D plane aligned with its normal,
    /// applies the earcut triangulation algorithm, then projects the resulting
    /// triangles back to 3D space. Returns a vector of triangular faces that
    /// collectively represent the original face.
    /// 
    /// # Example
    /// ```rust
    /// #[test]
    /// fn earcut_xy() {
    ///     let mut geom = Geom::from_file("./examples/data/plane_xy.obj").unwrap();
    ///     let face = geom.shapes[0].faces.remove(0);
    ///     let triangles = Face::earcut(&face);  // Triangulation call
    ///     assert_eq!(triangles.len(), 2);
    ///     assert_eq!(triangles[0].data().normal, face.data().normal);
    /// }
    /// ```
    pub fn earcut(face: &Face) -> Vec<Face> {
        let mut face = face.clone();
        // use nalgebra to get transform to xy plane
        let model = Isometry3::new(Vector3::zeros(), na::zero()); // do some sort of projection - set to nothing
        let origin = Point3::origin(); // camera location

        let target = Point3::new(
            face.data().normal.x,
            face.data().normal.y,
            face.data().normal.z,
        ); // projection direction, defines negative z-axis in new coords

        let up: Vector3<f32> =
            if face.data().normal.cross(&Vector3::y()).norm() < settings::COLINEAR_THRESHOLD {
                Vector3::x()
            } else {
                Vector3::y()
            };

        let view = Isometry3::look_at_rh(&origin, &target, &up);

        let transform = (view * model).to_homogeneous(); // transform to clipping system
        let itransform = transform.try_inverse().unwrap(); // inverse transform

        face.transform(&transform).unwrap();
        let poly = face.to_polygon();
        let triangles = poly.earcut_triangles();
        let outputs = triangles
            .iter()
            .filter_map(|tri| {
                let poly = tri.to_polygon();

                let mut face = match poly.project(&face.plane()) {
                    Ok(face) => face,
                    Err(_) => return None,
                };
                face.data_mut().exterior.reverse();

                if let Err(_) = face.transform(&itransform) {
                    return None;
                }

                Some(face)
            })
            .collect();

        outputs
    }
}

/// A closed 3D surface representing a single optical scatterer.
/// 
/// **Context**: Light scattering simulations require geometric representations
/// of particles. Each particle must be modeled as a closed surface with well-defined 
/// interior and exterior regions, along with its optical properties for [`crate::beam::Beam`]
/// propagation and [`crate::clip::Clipping`] operations.
/// 
/// **How it Works**: A [`Shape`] stores geometric and material information for
/// electromagnetic field calculations. The mesh contains [`nalgebra::Point3<f32>`] vertices 
/// defining 3D positions, [`Face`] collections describing surface polygons (which should 
/// be planar - error checking for non-planar faces will be added in future versions), 
/// and a [`nalgebra::Complex<f32>`] refractive index characterizing how light interacts 
/// with the material. Each [`Face`] maintains its own copy of vertex positions (rather 
/// than shared references) to enable independent geometric operations like 
/// [`crate::distortion::distort`] - though this requires careful synchronization when vertices change.
/// 
/// The [`crate::containment::AABB`] enables spatial queries for determining containment 
/// relationships between nested shapes, necessary for hierarchical geometries managed 
/// by [`crate::containment::ContainmentGraph`].
/// 
/// **Initialization**: [`Shape`] objects are created by loading Wavefront OBJ files through
/// [`Geom::from_file()`] which calls [`Shape::from_model()`], or directly via the
/// [`pyo3`] Python bindings for programmatic geometry generation.
#[pyclass]
#[derive(Debug, Clone, PartialEq)]
pub struct Shape {
    pub vertices: Vec<Point3<f32>>, // List of all vertices in the mesh
    pub num_vertices: usize,        // Number of vertices in the mesh
    pub faces: Vec<Face>,           // List of all [`Face`] elements in the mesh
    pub num_faces: usize,           // Number of facets in the mesh
    pub refr_index: Complex<f32>,   // Refractive index of this shape
    pub id: Option<usize>,          // an id number for [`crate::containment::ContainmentGraph`]
    pub parent_id: Option<usize>,   // An optional parent shape index, which encompasses this one
    pub aabb: Option<AABB>,         // axis-aligned bounding box for [`crate::containment`] operations
}

impl Shape {
    pub fn new(id: Option<usize>, parent_id: Option<usize>) -> Self {
        Self {
            vertices: Vec::new(),
            num_vertices: 0,
            faces: Vec::new(),
            num_faces: 0,
            refr_index: Complex { re: 1.31, im: 0.0 },
            id,
            parent_id,
            aabb: None,
        }
    }

    /// Converts a [`tobj::Model`] into a [`Shape`].
    /// 
    /// **Context**: The [`tobj`] library provides raw mesh data (vertices, indices, face counts)
    /// but this needs to be restructured into the [`Shape`] format required for electromagnetic
    /// calculations. The conversion must handle arbitrary polygon faces, establish proper 
    /// vertex-face relationships, and set up the geometric properties needed for subsequent
    /// operations.
    /// 
    /// **How it Works**: Extracts vertex positions from the model's position array, then
    /// processes face definitions using the face_arities to group indices into individual
    /// faces. Each face becomes a [`Face::Simple`] containing its vertex data and optional
    /// index references. The axis-aligned bounding box is computed from all vertices
    /// to enable spatial queries.
    fn from_model(model: Model, id: Option<usize>) -> Result<Shape> {
        let mesh = &model.mesh;

        let vertices = mesh
            .positions
            .chunks_exact(3)
            .map(|v| Point3::new(v[0] as f32, v[1] as f32, v[2] as f32))
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

            // Convert face indices to usize
            let usize_indices: Vec<usize> = face_indices.iter().map(|&i| i as usize).collect();

            let face_vertices: Vec<_> = usize_indices.iter().map(|&i| shape.vertices[i]).collect();
            shape.add_face(Face::new_simple(face_vertices, id, Some(usize_indices))?);

            next_face = end;
        }

        shape.set_aabb();

        Ok(shape)
    }

    pub fn set_aabb(&mut self) {
        let (min, max) = self.vertices.iter().fold(
            ([f32::INFINITY; 3], [-f32::INFINITY; 3]),
            |(min_acc, max_acc), v| {
                (
                    [
                        min_acc[0].min(v[0]),
                        min_acc[1].min(v[1]),
                        min_acc[2].min(v[2]),
                    ],
                    [
                        max_acc[0].max(v[0]),
                        max_acc[1].max(v[1]),
                        max_acc[2].max(v[2]),
                    ],
                )
            },
        );

        let min = Point3::from(min);
        let max = Point3::from(max);

        self.aabb = Some(AABB { min, max });
    }

    /// Uniformly scales all vertices and faces by a given factor.
    /// 
    /// **Context**: Shape-level scaling is needed during distortion recovery
    /// and for maintaining consistent geometry sizes after various operations.
    /// This operates on individual shapes rather than the entire geometry.
    /// 
    /// **How it Works**: Multiplies all vertex coordinates and face vertex
    /// coordinates by the scale factor, then recomputes the axis-aligned
    /// bounding box to reflect the new dimensions.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // From distortion recovery
    /// // let rescale_fac = max_dim / new_max_dim;
    /// // shape.rescale(rescale_fac);  // Rescale after distortion to maintain size
    /// ```
    pub fn rescale(&mut self, scale: f32) {
        for vertex in &mut self.vertices {
            vertex.coords *= scale;
        }
        for face in self.faces.iter_mut() {
            for vertex in face.data_mut().exterior.iter_mut() {
                vertex.coords *= scale;
            }
        }
        self.set_aabb(); // recompute axis-aligned bounding box
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

    pub fn transform(&mut self, transform: &Matrix4<f32>) -> Result<()> {
        for face in &mut self.faces {
            // Iterate mutably
            face.transform(transform)?; // Call the in-place project method
        }
        Ok(())
    }

    /// Determines if the axis-aligned bounding box of this shape contains
    /// that of another.
    /// Tests if this shape's bounding box completely contains another's.
    /// 
    /// **Context**: Determining hierarchical relationships between shapes
    /// requires testing spatial containment. This is used during containment
    /// graph construction to establish parent-child relationships.
    /// 
    /// **How it Works**: Compares axis-aligned bounding boxes, checking that
    /// the other shape's bounding box lies entirely within this shape's
    /// bounding box across all three dimensions.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // Iterate over distinct pairs of shapes
    /// // for (id_a, a) in &shapes_with_ids {
    /// //     for (id_b, b) in &shapes_with_ids {
    /// //         if id_a != id_b && a.contains(b) {  // Containment testing
    /// //             containment_graph.set_parent(*id_b, *id_a);
    /// //         }
    /// //     }
    /// // }
    /// ```
    pub fn contains(&self, other: &Shape) -> bool {
        match (&self.aabb, &other.aabb) {
            (Some(a), Some(b)) => (0..3).all(|i| b.min[i] > a.min[i] && a.max[i] > b.max[i]),
            (_, _) => false,
        }
    }

    /// determines if a shape in a geometry is inside another. Returns `true`
    /// if the two shapes have the same id.
    /// Checks if this shape is contained within another shape in the hierarchy.
    /// 
    /// **Context**: Ray tracing and clipping operations need to determine
    /// whether faces belong to internal or external boundaries. This requires
    /// understanding the containment relationships between shapes.
    /// 
    /// **How it Works**: Traverses up the containment graph starting from this
    /// shape's ID, checking each parent until either finding the target shape ID
    /// or reaching the top level. Returns true if the target shape is found
    /// in the ancestry chain.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // From clipping logic
    /// // for shape in self.geom.shapes.iter() {
    /// //     if internal && !shape.is_within(&self.geom, clip_shape_id) {
    /// //         continue;  // Skip shapes not within the clipping shape
    /// //     }
    /// //     // Process faces...
    /// // }
    /// ```
    pub fn is_within(&self, geom: &Geom, other_id: Option<usize>) -> bool {
        if other_id.is_none() {
            return false;
        } else if other_id.unwrap() == self.id.unwrap() {
            return true;
        }

        // traverse up the parents:
        let mut current = self.id; // get current shape id
        while current.is_some() {
            // while current shape id exists
            let parent_id = geom.containment_graph.get_parent(current.unwrap()); // try to get parent id
            if parent_id == other_id {
                // if parent id matches
                return true; // other must contain this shape
            }

            current = parent_id; // else, move up and try again
        }
        false
    }
}

/// Python bindings for the `Shape` struct.
#[pymethods]
impl Shape {
    #[new]
    fn py_new(
        vertices: Vec<(f32, f32, f32)>,
        face_indices: Vec<Vec<usize>>,
        id: usize,
        refr_index_re: f32,
        refr_index_im: f32,
    ) -> PyResult<Self> {
        let vertices = vertices
            .into_iter()
            .map(|(x, y, z)| Point3::new(x, y, z))
            .collect::<Vec<_>>();

        let mut shape = Shape::new(Some(id), None);
        shape.num_vertices = vertices.len();
        shape.vertices = vertices;

        const BODGE_SHAPE_ID: usize = 0;

        for indices in face_indices {
            let face_vertices: Vec<_> = indices.into_iter().map(|i| shape.vertices[i]).collect();
            shape.add_face(
                Face::new_simple(face_vertices, Some(BODGE_SHAPE_ID), None)
                    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?,
            );
        }

        shape.set_aabb();
        shape.refr_index = Complex {
            re: refr_index_re,
            im: refr_index_im,
        };

        Ok(shape)
    }
}

/// A geometric model for light scattering simulations.
/// 
/// **Context**: Light scattering scenarios involve multi-particle systems where
/// particles may contain hierarchical relationships, with some particles existing 
/// inside others. Each [`crate::problem::Problem`] requires working with an ensemble of shapes 
/// that maintain these spatial relationships.
/// 
/// **How it Works**: A [`Geom`] manages a collection of closed-surface [`Shape`] objects and their
/// spatial relationships through a [`crate::containment::ContainmentGraph`]. This graph tracks
/// parent-child relationships based on bounding box containment, enabling handling
/// of nested geometries during ray tracing and field calculations. The containment
/// information determines refractive indices at shape boundaries and manages
/// optical interfaces with [`crate::fresnel`] calculations.
/// 
/// Each [`crate::problem::Problem`] operates on exactly one [`Geom`] instance, but that geometry can
/// represent hierarchical structures. All shapes must be closed surfaces to ensure 
/// well-defined interior/exterior boundaries for electromagnetic field calculations.
#[pyclass]
#[derive(Debug, Clone, PartialEq)]
pub struct Geom {
    pub shapes: Vec<Shape>,                     // Vector of [`Shape`] objects
    pub containment_graph: ContainmentGraph,    // Hierarchical relationships between shapes
    pub num_shapes: usize,                      // Number of shapes in the collection
}

impl Geom {
    /// Loads a geometric model from a Wavefront OBJ file.
    /// 
    /// **Context**: Geometry for light scattering simulations typically originates
    /// from external sources - CAD models, 3D scanners, or mesh generation software.
    /// The standard OBJ format provides a widely-supported way to import these geometries,
    /// but the raw mesh data needs to be processed into the internal representation
    /// required for electromagnetic field calculations.
    /// 
    /// **How it Works**: Parses the OBJ file using the [`tobj`] library to extract vertices
    /// and face definitions. Each mesh object in the file becomes a separate [`Shape`] with
    /// automatically assigned IDs. Face polygons are converted to the internal [`Face`]
    /// representation, and axis-aligned bounding boxes are computed for each shape.
    /// A [`crate::containment::ContainmentGraph`] is constructed by testing bounding box relationships between
    /// all shape pairs, establishing parent-child hierarchies for nested geometries.
    /// 
    /// # Example
    /// ```rust,no_run
    /// use goad::geom::Geom;
    /// use goad::problem::Problem;
    /// let mut geom = Geom::from_file("./examples/data/hex2.obj").unwrap();
    /// geom.shapes[0].refr_index.re = 1.5;
    /// geom.shapes[0].refr_index.im = 0.0001;
    /// let mut problem = Problem::new(Some(geom), None);
    /// ```
    pub fn from_file(filename: &str) -> Result<Self> {
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
                .map_err(|e| anyhow::anyhow!("Could not resolve path: {}", e))?
        };

        let (models, _) = tobj::load_obj(&resolved_filename, &tobj::LoadOptions::default())
            .map_err(|e| anyhow::anyhow!("Failed to load OBJ file '{}': {}", filename, e))?;

        if models.is_empty() {
            return Err(anyhow::anyhow!("No models found in OBJ file"));
        }

        let shapes = Self::shapes_from_models(models)?;

        // Create containment graph
        let mut containment_graph = ContainmentGraph::new(shapes.len());

        // Ensure all shapes have valid IDs upfront
        let shapes_with_ids: Vec<_> = shapes
            .iter()
            .filter_map(|shape| {
                shape
                    .id
                    .map(|id| (id, shape))
                    .or_else(|| panic!("Shape cannot be added to containment graph without an id"))
            })
            .collect();

        // Iterate over distinct pairs of shapes
        for (id_a, a) in &shapes_with_ids {
            for (id_b, b) in &shapes_with_ids {
                if id_a != id_b && a.contains(b) {
                    containment_graph.set_parent(*id_b, *id_a);
                }
            }
        }

        Ok(Self {
            num_shapes: shapes.len(),
            shapes,
            containment_graph,
        })
    }

    fn shapes_from_models(models: Vec<Model>) -> Result<Vec<Shape>> {
        models
            .into_iter()
            .enumerate()
            .map(|(i, model)| Shape::from_model(model, Some(i)))
            .collect()
    }

    /// Applies a 4x4 transformation matrix to all geometry.
    /// 
    /// **Context**: Computing how [`crate::beam::Beam`] segments intersect with faces requires rotating 
    /// beam cross-section polygons and faces in the geometry into the plane normal 
    /// to the beam propagation direction. This coordinate transformation is essential 
    /// for ray-tracing calculations and [`crate::field`] propagation algorithms.
    /// 
    /// **How it Works**: Applies the transformation matrix to every [`Face`] in every [`Shape`].
    /// Each face handles its own transformation, updating vertex positions, midpoints,
    /// and normals according to the [`nalgebra::Matrix4<f32>`]. This enables arbitrary rotations, scaling,
    /// and translations of the entire geometric model.
    pub fn transform(&mut self, transform: &Matrix4<f32>) -> Result<()> {
        for shape in &mut self.shapes {
            shape.transform(transform)?;
        }
        Ok(())
    }

    /// Computes the center of mass of all vertices across all shapes.
    /// 
    /// **Context**: Many geometric operations require a reference point for
    /// the geometry, particularly rotations which should occur around the
    /// geometric center rather than an arbitrary origin. This function provides
    /// that reference point for centering operations.
    /// 
    /// **How it Works**: Averages the positions of all vertices in all shapes
    /// to find the centroid. This is not a true center of mass calculation
    /// (which would require mass/density information) but rather a geometric
    /// centroid suitable for coordinate system operations.
    pub fn centre_of_mass(&self) -> Point3<f32> {
        let mut centre = Point3::origin();

        for shape in &self.shapes {
            centre += calculate_center_of_mass(&shape.vertices).coords;
        }

        centre / self.num_shapes as f32
    }

    /// Returns the refractive outside a shape
    /// Returns the refractive index outside a given shape.
    /// 
    /// **Context**: Fresnel reflection and refraction calculations at shape
    /// boundaries require knowing the refractive indices on both sides of
    /// the interface. For nested geometries, the "outside" medium may be
    /// another shape rather than the background medium.
    /// 
    /// **How it Works**: Uses the containment graph to find the parent shape
    /// of the given shape ID. If a parent exists, returns its refractive index;
    /// otherwise returns the provided medium refractive index for the outermost
    /// background.
    pub fn n_out(&self, shape_id: usize, medium_refr_index: Complex<f32>) -> Complex<f32> {
        self.containment_graph
            .get_parent(shape_id)
            .map_or(medium_refr_index, |parent_id| {
                self.shapes[parent_id].refr_index
            })
    }

    /// Returns the axis-aligned bounding box encompassing all shapes.
    /// 
    /// **Context**: Beam initialization requires knowing the spatial extent
    /// of the geometry to create appropriate illumination areas. The bounding
    /// box provides the minimum and maximum coordinates needed for this setup.
    /// 
    /// **How it Works**: Iterates through all shape bounding boxes to find
    /// the global minimum and maximum coordinates across all three dimensions.
    /// Returns (min_point, max_point) defining the overall geometric extent.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // fn basic_initial_beam(geom: &Geom, wavelength: f32, medium_refractive_index: Complex<f32>) -> Beam {
    /// //     const FAC: f32 = 1.1; // scale factor to stretch beam to cover geometry
    /// //     let bounds = geom.bounds();
    /// //     let (min, max) = (bounds.0.map(|v| v * FAC), bounds.1.map(|v| v * FAC));
    /// //     // Create beam cross-section to fully illuminate the geometry
    /// //     // ...
    /// // }
    /// ```
    pub fn bounds(&self) -> (Point3<f32>, Point3<f32>) {
        let (min, max) = self.shapes.iter().fold(
            ([f32::INFINITY; 3], [-f32::INFINITY; 3]),
            |(min_acc, max_acc), shape| {
                let aabb = shape.aabb.as_ref().unwrap();
                (
                    [
                        min_acc[0].min(aabb.min[0]),
                        min_acc[1].min(aabb.min[1]),
                        min_acc[2].min(aabb.min[2]),
                    ],
                    [
                        max_acc[0].max(aabb.max[0]),
                        max_acc[1].max(aabb.max[1]),
                        max_acc[2].max(aabb.max[2]),
                    ],
                )
            },
        );

        (Point3::from(min), Point3::from(max))
    }

    /// Rescales the geometry so that the largest dimension is 1. Returns the
    /// scaling factor.
    /// Rescales the geometry so that the largest dimension is 1.
    /// 
    /// **Context**: Electromagnetic field calculations can suffer from numerical
    /// instability when geometry dimensions vary widely. Normalizing to unit scale
    /// provides better numerical conditioning for the solver algorithms.
    /// 
    /// **How it Works**: Computes the bounding box of all shapes, finds the
    /// maximum dimension, then scales all vertices and faces by 1/max_dimension.
    /// Returns the scaling factor applied, which can be used to scale results
    /// back to original units if needed.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // self.geom.recentre();
    /// // self.settings.scale = self.geom.rescale();  // Returns scaling factor
    /// ```
    pub fn rescale(&mut self) -> f32 {
        let bounds = self.bounds();
        let max_dim = bounds.1.iter().fold(0.0, |acc: f32, &x| acc.max(x));
        let scale = 1.0 / max_dim;

        for shape in self.shapes.iter_mut() {
            shape.rescale(scale);
        }

        scale
    }

    /// Translates the geometry so its center of mass is at the origin.
    /// 
    /// **Context**: Rotations and other transformations typically need to occur
    /// around the geometric center rather than an arbitrary coordinate system
    /// origin. This function ensures the geometry is properly centered before
    /// such operations.
    /// 
    /// **How it Works**: Computes the center of mass, then subtracts this offset
    /// from all vertex positions in both the shape vertex arrays and the face
    /// vertex data. This dual update maintains consistency between the two
    /// vertex representations.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // From problem initialization workflow
    /// // if let Some(distortion) = self.settings.distortion {
    /// //     self.geom.distort(distortion, self.settings.seed);
    /// // }
    /// // self.geom.recentre();  // Center before rescaling
    /// // self.settings.scale = self.geom.rescale();
    /// ```
    pub fn recentre(&mut self) {
        let com = self.centre_of_mass();

        for shape in self.shapes.iter_mut() {
            for vertex in shape.vertices.iter_mut() {
                vertex.coords -= com.coords;
            }

            for face in shape.faces.iter_mut() {
                match face {
                    Face::Simple(data) => {
                        for vertex in data.exterior.iter_mut() {
                            vertex.coords -= com.coords;
                        }
                    }
                    Face::Complex { data, interiors } => {
                        for vertex in data.exterior.iter_mut() {
                            vertex.coords -= com.coords;
                        }

                        for interior in interiors.iter_mut() {
                            for vertex in interior.iter_mut() {
                                vertex.coords -= com.coords;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Checks if the geometry's center of mass is at the origin.
    /// 
    /// **Context**: Some operations, particularly rotations, require the
    /// geometry to be centered at the origin. This function provides a
    /// precondition check to ensure transformations will behave correctly.
    /// 
    /// **How it Works**: Computes the center of mass and checks if its
    /// magnitude is below a tolerance threshold (1e-6).
    /// 
    /// # Example
    /// ```rust,no_run
    /// // pub fn euler_rotate(&mut self, euler: &Euler, convention: EulerConvention) -> Result<()> {
    /// //     if !self.is_centered() {
    /// //         return Err(anyhow::anyhow!(
    /// //             "Geometry must be centred before rotation can be applied. HINT: Try geom.recentre()"
    /// //         ));
    /// //     }
    /// //     // ... proceed with rotation
    /// // }
    /// ```
    pub fn is_centered(&self) -> bool {
        self.centre_of_mass().coords.norm() < 1e-6
    }

    /// Rotates the geometry by the Euler angles alpha, beta, and gamma (in degrees)
    /// Uses Mishchenko's Euler rotation matrix convention.
    /// Rotates the geometry using Euler angles.
    /// 
    /// **Context**: Light scattering calculations often require analyzing
    /// particles at different orientations. Euler rotations provide a
    /// systematic way to orient the geometry for scattering analysis.
    /// 
    /// **How it Works**: Converts the Euler angles to a rotation matrix
    /// using the specified convention, then applies this rotation to all
    /// vertices, face midpoints, and normals. Requires the geometry to be
    /// centered at the origin first to ensure rotation occurs around the
    /// geometric center.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // pub fn orient(&mut self, euler: &orientation::Euler) {
    /// //     if let Err(error) = self
    /// //         .geom
    /// //         .euler_rotate(euler, self.settings.orientation.euler_convention)
    /// //     {
    /// //         panic!("Error rotating geometry: {}", error);
    /// //     }
    /// // }
    /// ```
    pub fn euler_rotate(&mut self, euler: &Euler, convention: EulerConvention) -> Result<()> {
        if !self.is_centered() {
            return Err(anyhow::anyhow!(
                "Geometry must be centred before rotation can be applied. HINT: Try geom.recentre()"
            ));
        }

        let rotation = euler.rotation_matrix(convention);

        for shape in self.shapes.iter_mut() {
            for vertex in shape.vertices.iter_mut() {
                vertex.coords = rotation * vertex.coords;
            }

            for face in shape.faces.iter_mut() {
                match face {
                    Face::Simple(data) => {
                        data.midpoint = rotation * data.midpoint;
                        data.normal = rotation * data.normal;
                        for vertex in data.exterior.iter_mut() {
                            vertex.coords = rotation * vertex.coords;
                        }
                    }

                    Face::Complex { data, interiors } => {
                        for vertex in data.exterior.iter_mut() {
                            vertex.coords = rotation * vertex.coords;
                        }

                        for interior in interiors.iter_mut() {
                            for vertex in interior.iter_mut() {
                                vertex.coords = rotation * vertex.coords;
                            }
                        }
                    }
                }
            }
            shape.set_aabb();
        }

        Ok(())
    }

    /// Exports the geometry to a Wavefront OBJ file.
    /// 
    /// **Context**: After geometric processing and transformations, it's often
    /// necessary to export the modified geometry for visualization, verification,
    /// or use in other software. The OBJ format provides a standard interchange
    /// format for 3D geometry.
    /// 
    /// **How it Works**: Writes vertices, normals, and face definitions for all
    /// shapes in standard OBJ format. Each shape becomes a separate group with
    /// metadata comments for shape ID and refractive index. Does not support
    /// complex faces with interior holes.
    /// 
    /// # Example
    /// ```rust,no_run
    /// // From distortion example
    /// // geom.distort(0.5, None);
    /// // geom.recentre();
    /// // let euler = Euler::new(0.0, 30.0, 0.0);
    /// // let result = geom.euler_rotate(&euler, goad::orientation::EulerConvention::XYX);
    /// 
    /// // write the distorted object to a file
    /// // geom.write_obj("hex_distorted.obj").unwrap();
    /// ```
    pub fn write_obj<P: AsRef<Path>>(&self, filename: P) -> Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};

        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "# OBJ file generated by GOAD")?;
        writeln!(writer, "# Total shapes: {}", self.num_shapes)?;

        for shape_idx in 0..self.num_shapes {
            let vertex_offset = 1; // OBJ indices start at 1
            let normal_offset = 1;
            let shape = &self.shapes[shape_idx];

            writeln!(writer, "g shape_{}", shape_idx)?;
            writeln!(writer, "# Shape ID: {:?}", shape.id)?;
            writeln!(
                writer,
                "# Refractive index: {} + {}i",
                shape.refr_index.re, shape.refr_index.im
            )?;

            // Write vertices
            for vertex in &shape.vertices {
                writeln!(writer, "v {} {} {}", vertex.x, vertex.y, vertex.z)?;
            }

            // Write normals from faces
            let mut normals = Vec::new();
            for face in &shape.faces {
                let normal = match face {
                    Face::Simple(data) => &data.normal,
                    Face::Complex { .. } => {
                        panic!("Complex faces with interior holes are not supported in OBJ export")
                    }
                };

                normals.push(normal);
                writeln!(writer, "vn {} {} {}", normal.x, normal.y, normal.z)?;
            }

            // Write faces (using vertex indices and normal indices)
            for (face_idx, face) in shape.faces.iter().enumerate() {
                match face {
                    Face::Simple(data) => {
                        write!(writer, "f")?;

                        // If we have explicit indices in the face data, use those
                        if let Some(indices) = &data.exterior_indices {
                            for &idx in indices {
                                write!(
                                    writer,
                                    " {}//{}",
                                    idx + vertex_offset,
                                    face_idx + normal_offset
                                )?;
                            }
                        } else {
                            // Otherwise find vertex indices by matching vertices
                            for vertex in &data.exterior {
                                // Find the index of this vertex in the shape's vertices
                                if let Some(idx) = shape.vertices.iter().position(|v| v == vertex) {
                                    write!(
                                        writer,
                                        " {}//{}",
                                        idx + vertex_offset,
                                        face_idx + normal_offset
                                    )?;
                                } else {
                                    return Err(anyhow::anyhow!(
                                        "Vertex not found in shape vertices"
                                    ));
                                }
                            }
                        }
                        writeln!(writer)?;
                    }
                    Face::Complex { .. } => {
                        panic!("Complex faces with interior holes are not supported in OBJ export");
                    }
                }
            }
        }

        Ok(())
    }

    /// Applies non-uniform scaling with different factors for each dimension.
    /// 
    /// **Context**: Some simulations require anisotropic deformation of particles
    /// to study aspect ratio effects or simulate non-spherical geometries.
    /// This enables independent scaling along X, Y, and Z axes.
    /// 
    /// **How it Works**: Multiplies vertex coordinates by the corresponding
    /// scale factors, updates face midpoints and normals accordingly, then
    /// renormalizes the normals and recomputes bounding boxes. Requires a
    /// 3-element scale vector.
    pub fn vector_scale(&mut self, scale: &Vec<f32>) {
        if scale.len() != 3 {
            panic!("Scale vector must have length 3");
        }

        let scale_vec = Vector3::new(scale[0], scale[1], scale[2]);

        for shape in self.shapes.iter_mut() {
            for vertex in shape.vertices.iter_mut() {
                vertex.coords.component_mul_assign(&scale_vec);
            }

            for face in shape.faces.iter_mut() {
                match face {
                    Face::Simple(data) => {
                        data.midpoint.coords.component_mul_assign(&scale_vec);
                        data.normal.component_mul_assign(&scale_vec);
                        data.normal.normalize_mut(); // Re-normalize after scaling

                        for vertex in data.exterior.iter_mut() {
                            vertex.coords.component_mul_assign(&scale_vec);
                        }
                    }

                    Face::Complex { data, interiors } => {
                        data.midpoint.coords.component_mul_assign(&scale_vec);
                        data.normal.component_mul_assign(&scale_vec);
                        data.normal.normalize_mut(); // Re-normalize after scaling

                        for vertex in data.exterior.iter_mut() {
                            vertex.coords.component_mul_assign(&scale_vec);
                        }

                        for interior in interiors.iter_mut() {
                            for vertex in interior.iter_mut() {
                                vertex.coords.component_mul_assign(&scale_vec);
                            }
                        }
                    }
                }
            }
            shape.set_aabb();
        }
    }
}

/// Python bindings for the `Geom` struct.
#[pymethods]
impl Geom {
    #[new]
    fn py_new(shapes: Vec<Shape>) -> Self {
        let num_shapes = shapes.len();
        let mut containment_graph = ContainmentGraph::new(num_shapes);

        // Ensure all shapes have valid IDs upfront
        let shapes_with_ids: Vec<_> = shapes
            .iter()
            .filter_map(|shape| {
                shape
                    .id
                    .map(|id| (id, shape))
                    .or_else(|| panic!("Shape cannot be added to containment graph without an id"))
            })
            .collect();

        // Iterate over distinct pairs of shapes
        for (id_a, a) in &shapes_with_ids {
            for (id_b, b) in &shapes_with_ids {
                if id_a != id_b && a.contains(b) {
                    containment_graph.set_parent(*id_b, *id_a);
                }
            }
        }

        Self {
            shapes,
            containment_graph,
            num_shapes,
        }
    }

    /// Getter for the vertices of the first shape
    #[getter]
    fn get_first_shape_vertices(&self) -> Vec<(f32, f32, f32)> {
        self.shapes[0]
            .vertices
            .iter()
            .map(|v| (v.x, v.y, v.z))
            .collect()
    }
}

/// Computes the geometric centroid of a vertex set.
/// 
/// **Context**: Various geometric operations require finding a representative
/// center point for a collection of vertices. While not a true center of mass
/// (which would require density information), this centroid is sufficient
/// for coordinate system operations.
/// 
/// **How it Works**: Averages the coordinates of all vertices to find the
/// geometric center point of the vertex cloud.
pub fn calculate_center_of_mass(verts: &[Point3<f32>]) -> Point3<f32> {
    Point3::from(
        verts
            .iter()
            .map(|vert| vert.coords)
            .fold(Vector3::zeros(), |acc, coords| acc + coords)
            / verts.len() as f32,
    )
}

/// Translates vertices relative to a reference point.
/// 
/// **Context**: Coordinate system transformations often require expressing
/// vertex positions relative to a reference point rather than the global
/// origin. This is commonly used for centering operations.
/// 
/// **How it Works**: Subtracts the reference point coordinates from each
/// vertex position, returning the translated positions as vectors.
pub fn translate(verts: &[Point3<f32>], center_of_mass: &Point3<f32>) -> Vec<Vector3<f32>> {
    verts
        .iter()
        .map(|point| point.coords - center_of_mass.coords)
        .collect()
}
