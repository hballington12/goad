use super::geom::{Face, Geom, Plane};
use super::settings;
use crate::geom::PolygonExtensions;
use anyhow::Result;
use geo::{Area, BooleanOps, Simplify};

use nalgebra::{self as na, Isometry3, Matrix4, Point3, Vector3};
use std::cmp::Ordering;
use std::fmt;

#[cfg(test)]
mod tests {

    use geo::polygon;
    use geo::{CoordsIter, MultiPolygon, Simplify};

    use super::*;
    const AREA_THRESHOLD: f32 = 0.01;

    /// Helper to verify vertices match expected values in cyclic order.
    /// Finds the first expected vertex in `actual`, then checks that subsequent
    /// vertices match in order (wrapping around). Also checks reverse order
    /// in case winding direction differs.
    fn assert_vertices_match_cyclic(
        actual: &[Point3<f32>],
        expected: &[(f32, f32)],
        expected_z: f32,
        tolerance: f32,
    ) {
        assert_eq!(
            actual.len(),
            expected.len(),
            "Vertex count mismatch: got {}, expected {}",
            actual.len(),
            expected.len()
        );

        // Find starting index in actual that matches expected[0]
        let start_idx = actual.iter().position(|v| {
            (v.x - expected[0].0).abs() < tolerance && (v.y - expected[0].1).abs() < tolerance
        });

        let start_idx = match start_idx {
            Some(idx) => idx,
            None => {
                // Try finding any expected vertex as start
                let mut found_start = None;
                for (ei, (ex, ey)) in expected.iter().enumerate() {
                    if let Some(ai) = actual
                        .iter()
                        .position(|v| (v.x - ex).abs() < tolerance && (v.y - ey).abs() < tolerance)
                    {
                        found_start = Some((ai, ei));
                        break;
                    }
                }
                match found_start {
                    Some((_ai, ei)) => {
                        // Rotate expected to start from ei
                        let rotated: Vec<_> = expected
                            .iter()
                            .cycle()
                            .skip(ei)
                            .take(expected.len())
                            .cloned()
                            .collect();
                        return assert_vertices_match_cyclic(actual, &rotated, expected_z, tolerance);
                    }
                    None => panic!(
                        "Could not find any expected vertex in actual vertices.\nExpected: {:?}\nActual: {:?}",
                        expected,
                        actual.iter().map(|v| (v.x, v.y)).collect::<Vec<_>>()
                    ),
                }
            }
        };

        // Try forward order
        let forward_match = expected.iter().enumerate().all(|(i, (ex, ey))| {
            let v = &actual[(start_idx + i) % actual.len()];
            (v.x - ex).abs() < tolerance && (v.y - ey).abs() < tolerance
        });

        // Try reverse order (different winding)
        let reverse_match = expected.iter().enumerate().all(|(i, (ex, ey))| {
            let idx = (start_idx + actual.len() - i) % actual.len();
            let v = &actual[idx];
            (v.x - ex).abs() < tolerance && (v.y - ey).abs() < tolerance
        });

        assert!(
            forward_match || reverse_match,
            "Vertices do not match in cyclic order (forward or reverse).\nExpected: {:?}\nActual: {:?}",
            expected,
            actual.iter().map(|v| (v.x, v.y)).collect::<Vec<_>>()
        );

        // Check z coordinates
        for v in actual {
            assert!(
                (v.z - expected_z).abs() < tolerance,
                "Vertex z={} does not match expected z={}",
                v.z,
                expected_z
            );
        }
    }

    #[test]
    #[should_panic]
    fn concave_clip() {
        let geoms = Geom::load("./examples/data/concave1.obj").unwrap();
        let mut geom = geoms[0].clone();

        let clip_index = 4; // the index of the face to be used as the clip
        let projection = Vector3::new(-0.3, 0.0, -1.0);
        let mut clip = geom.shapes[0].faces.remove(clip_index); // choose a face be the clip

        // start function `do_clip` here:
        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        let _ = clipping.clip(AREA_THRESHOLD);
        let _ = clipping.clip(AREA_THRESHOLD); // cannot redo clipping
    }

    #[test]
    fn remove_duplicate_vertices() {
        // Define a MultiPolygon
        let multipolygon = MultiPolygon(vec![polygon![
            (x: 0.0, y: 0.0),
            (x: 5.0, y: 0.0),
            (x: 5.0, y: 5.0),
            (x: 5.0, y: 4.99),
            (x: 0.0, y: 5.0),
            (x: 0.0, y: 0.0),
        ]]);

        println!("Original MultiPolygon: {:?}", multipolygon);

        let cleaned = Simplify::simplify(&multipolygon, 0.01);

        // Print the cleaned polygon
        println!("Cleaned MultiPolygon: {:?}", cleaned);

        // Assert that the number of vertices in the cleaned exterior is 5
        let cleaned_exterior = &cleaned.0[0].exterior();
        assert_eq!(cleaned_exterior.coords_count(), 5);
    }

    // =========================================================================
    // Regression tests for clipping behavior (for i_overlay migration)
    // Based on old macroquad visualization examples
    // =========================================================================

    /// Test based on projection-debug.rs example
    /// hex.obj, face 4, projection (0,0,-1)
    #[test]
    fn projection_debug_hex() {
        let geoms = Geom::load("./examples/data/hex.obj").unwrap();
        let mut geom = geoms[0].clone();

        assert_eq!(geom.shapes.len(), 1);
        assert_eq!(geom.shapes[0].faces.len(), 8);

        let mut clip = geom.shapes[0].faces.remove(4);
        let projection = Vector3::new(0.0, 0.0, -1.0);

        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip(AREA_THRESHOLD).unwrap();

        let stats = clipping.stats.as_ref().unwrap();

        // Expected values from old implementation
        assert!(
            (stats.clipping_area - 64.951904).abs() < 0.01,
            "clipping_area: expected 64.951904, got {}",
            stats.clipping_area
        );
        assert!(
            (stats.intersection_area - 64.951904).abs() < 0.01,
            "intersection_area: expected 64.951904, got {}",
            stats.intersection_area
        );
        assert!(
            stats.remaining_area.abs() < 0.01,
            "remaining_area: expected 0, got {}",
            stats.remaining_area
        );
        assert!(
            (stats.total_consvtn - 1.0).abs() < 0.001,
            "total_consvtn: expected 1.0, got {}",
            stats.total_consvtn
        );

        // Check counts
        assert_eq!(clipping.intersections.len(), 1, "Expected 1 intersection");
        assert_eq!(clipping.remaining.len(), 0, "Expected 0 remaining");

        // Check intersection geometry
        let intsn = &clipping.intersections[0];
        assert_eq!(
            intsn.data().num_vertices,
            6,
            "Intersection should have 6 vertices"
        );

        let midpoint = intsn.data().midpoint;
        assert!(
            (midpoint.x - 0.0).abs() < 0.01,
            "midpoint.x: expected 0.0, got {}",
            midpoint.x
        );
        assert!(
            (midpoint.y - 0.0).abs() < 0.01,
            "midpoint.y: expected 0.0, got {}",
            midpoint.y
        );
        assert!(
            (midpoint.z - (-5.0)).abs() < 0.01,
            "midpoint.z: expected -5.0, got {}",
            midpoint.z
        );

        let normal = intsn.data().normal;
        assert!(
            (normal.z - (-1.0)).abs() < 0.01,
            "normal.z: expected -1.0, got {}",
            normal.z
        );

        // Check vertex positions (hexagonal shape at z=-5)
        let verts = &intsn.data().exterior;
        let expected_verts = [
            (0.0, -5.0),
            (-4.330127, -2.5),
            (-4.330127, 2.5),
            (0.0, 5.0),
            (4.330127, 2.5),
            (4.330127, -2.5),
        ];
        assert_vertices_match_cyclic(verts, &expected_verts, -5.0, 0.01);
    }

    /// Test based on projection1.rs example
    /// concave1.obj, face 4, projection (-0.3, 0, -1)
    #[test]
    fn projection1_concave() {
        let geoms = Geom::load("./examples/data/concave1.obj").unwrap();
        let mut geom = geoms[0].clone();

        assert_eq!(geom.shapes.len(), 1);
        assert_eq!(geom.shapes[0].faces.len(), 8);

        let mut clip = geom.shapes[0].faces.remove(4);
        let projection = Vector3::new(-0.3, 0.0, -1.0);

        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip(AREA_THRESHOLD).unwrap();

        let stats = clipping.stats.as_ref().unwrap();

        // Expected values from old implementation
        assert!(
            (stats.clipping_area - 43.419937).abs() < 0.01,
            "clipping_area: expected 43.419937, got {}",
            stats.clipping_area
        );
        assert!(
            (stats.intersection_area - 43.41993).abs() < 0.01,
            "intersection_area: expected 43.41993, got {}",
            stats.intersection_area
        );
        assert!(
            (stats.total_consvtn - 1.0).abs() < 0.001,
            "total_consvtn: expected 1.0, got {}",
            stats.total_consvtn
        );

        // Check counts
        assert_eq!(clipping.intersections.len(), 4, "Expected 4 intersections");
        assert_eq!(clipping.remaining.len(), 0, "Expected 0 remaining");

        // Check vertex counts for each intersection
        assert_eq!(clipping.intersections[0].data().num_vertices, 4);
        assert_eq!(clipping.intersections[1].data().num_vertices, 4);
        assert_eq!(clipping.intersections[2].data().num_vertices, 4);
        assert_eq!(clipping.intersections[3].data().num_vertices, 7);

        // Check first intersection midpoint
        let mid0 = clipping.intersections[0].data().midpoint;
        assert!(
            (mid0.x - 1.166282).abs() < 0.01,
            "intersection[0] midpoint.x"
        );
        assert!(
            (mid0.y - 1.495947).abs() < 0.01,
            "intersection[0] midpoint.y"
        );
        assert!(
            (mid0.z - 0.033455).abs() < 0.01,
            "intersection[0] midpoint.z"
        );

        // Check last intersection (7 vertices) midpoint
        let mid3 = clipping.intersections[3].data().midpoint;
        assert!(
            (mid3.x - (-2.315247)).abs() < 0.01,
            "intersection[3] midpoint.x"
        );
        assert!(
            (mid3.y - (-0.006494)).abs() < 0.01,
            "intersection[3] midpoint.y"
        );
        assert!(
            (mid3.z - (-4.682372)).abs() < 0.01,
            "intersection[3] midpoint.z"
        );
    }

    /// Test based on projection2.rs example
    /// cube_inside_ico.obj, face 5 from shape 0, projection (-0.2, 0, -1)
    #[test]
    fn projection2_cube_inside_ico() {
        let geoms = Geom::load("./examples/data/cube_inside_ico.obj").unwrap();
        let mut geom = geoms[0].clone();

        assert_eq!(geom.shapes.len(), 2);
        assert_eq!(geom.shapes[0].faces.len(), 6);
        assert_eq!(geom.shapes[1].faces.len(), 20);

        let mut clip = geom.shapes[0].faces.remove(5);
        let projection = Vector3::new(-0.2, 0.0, -1.0);

        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip(AREA_THRESHOLD).unwrap();

        let stats = clipping.stats.as_ref().unwrap();

        // Expected values from old implementation
        assert!(
            (stats.clipping_area - 3.385843).abs() < 0.01,
            "clipping_area: expected 3.385843, got {}",
            stats.clipping_area
        );
        assert!(
            (stats.intersection_area - 3.385844).abs() < 0.01,
            "intersection_area: expected 3.385844, got {}",
            stats.intersection_area
        );
        assert!(
            (stats.total_consvtn - 1.0).abs() < 0.001,
            "total_consvtn: expected 1.0, got {}",
            stats.total_consvtn
        );

        // Check counts
        assert_eq!(clipping.intersections.len(), 3, "Expected 3 intersections");
        assert_eq!(clipping.remaining.len(), 0, "Expected 0 remaining");

        // Check vertex counts
        assert_eq!(clipping.intersections[0].data().num_vertices, 4);
        assert_eq!(clipping.intersections[1].data().num_vertices, 4);
        assert_eq!(clipping.intersections[2].data().num_vertices, 4);

        // Check intersection midpoints
        let mid0 = clipping.intersections[0].data().midpoint;
        assert!(
            (mid0.x - 0.264427).abs() < 0.01,
            "intersection[0] midpoint.x"
        );
        assert!(
            (mid0.y - (-0.852768)).abs() < 0.01,
            "intersection[0] midpoint.y"
        );
        assert!(
            (mid0.z - (-0.508341)).abs() < 0.01,
            "intersection[0] midpoint.z"
        );

        let mid2 = clipping.intersections[2].data().midpoint;
        assert!(
            (mid2.x - (-0.231690)).abs() < 0.01,
            "intersection[2] midpoint.x"
        );
        assert!((mid2.y - 0.0).abs() < 0.01, "intersection[2] midpoint.y");
        assert!(
            (mid2.z - (-1.158448)).abs() < 0.01,
            "intersection[2] midpoint.z"
        );
    }

    /// Test based on projection_multi.rs example
    /// multiple.obj, shape[0].face[5], projection (-1, 0, 0)
    #[test]
    fn projection_multi() {
        let geoms = Geom::load("./examples/data/multiple.obj").unwrap();
        let mut geom = geoms[0].clone();

        assert_eq!(geom.shapes.len(), 2);
        assert_eq!(geom.shapes[0].faces.len(), 8);
        assert_eq!(geom.shapes[1].faces.len(), 8);

        let mut clip = geom.shapes[0].faces.remove(5);
        let projection = Vector3::new(-1.0, 0.0, 0.0);

        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip(AREA_THRESHOLD).unwrap();

        let stats = clipping.stats.as_ref().unwrap();

        // Expected values from old implementation
        assert!(
            (stats.clipping_area - 49.999989).abs() < 0.01,
            "clipping_area: expected 49.999989, got {}",
            stats.clipping_area
        );
        assert!(
            (stats.intersection_area - 41.416775).abs() < 0.1,
            "intersection_area: expected 41.416775, got {}",
            stats.intersection_area
        );
        assert!(
            (stats.remaining_area - 8.583215).abs() < 0.1,
            "remaining_area: expected 8.583215, got {}",
            stats.remaining_area
        );
        assert!(
            (stats.total_consvtn - 1.0).abs() < 0.001,
            "total_consvtn: expected 1.0, got {}",
            stats.total_consvtn
        );

        // Check counts
        assert_eq!(clipping.intersections.len(), 3, "Expected 3 intersections");
        assert_eq!(clipping.remaining.len(), 2, "Expected 2 remaining");

        // Check vertex counts (allowing +1 for collinear points that may be preserved)
        assert!(clipping.intersections[0].data().num_vertices >= 4);
        assert!(clipping.intersections[1].data().num_vertices >= 7);
        assert!(clipping.intersections[2].data().num_vertices >= 4);
        assert!(clipping.remaining[0].data().num_vertices >= 3);
        assert!(clipping.remaining[1].data().num_vertices >= 3);

        // Check intersection midpoints
        let mid0 = clipping.intersections[0].data().midpoint;
        assert!(
            (mid0.x - (-9.606176)).abs() < 0.01,
            "intersection[0] midpoint.x"
        );
        assert!(
            (mid0.y - (-0.560584)).abs() < 0.01,
            "intersection[0] midpoint.y"
        );
        assert!(
            (mid0.z - 1.742397).abs() < 0.01,
            "intersection[0] midpoint.z"
        );

        // Check that expected remaining midpoints exist (order-independent)
        let expected_remaining_midpoints = [
            (-4.330127, 1.951482, -2.436753),
            (-4.330127, -5.016, -0.710), // approximate second remaining
        ];
        for (ex, ey, ez) in &expected_remaining_midpoints {
            let found = clipping.remaining.iter().any(|rem| {
                let mid = rem.data().midpoint;
                (mid.x - ex).abs() < 0.1 && (mid.y - ey).abs() < 0.1 && (mid.z - ez).abs() < 0.1
            });
            assert!(
                found,
                "Expected remaining midpoint ({}, {}, {}) not found",
                ex, ey, ez
            );
        }
    }

    /// Test based on clip_test.rs example
    /// clip_test.obj, shape[1].face[1], projection (1, 1, 0)
    #[test]
    fn clip_test_two_shapes() {
        let geoms = Geom::load("./examples/data/clip_test.obj").unwrap();
        let mut geom = geoms[0].clone();

        assert_eq!(geom.shapes.len(), 2);
        assert_eq!(geom.shapes[0].faces.len(), 6);
        assert_eq!(geom.shapes[1].faces.len(), 6);

        let mut clip = geom.shapes[1].faces.remove(1);
        let projection = Vector3::new(1.0, 1.0, 0.0);

        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip(AREA_THRESHOLD).unwrap();

        let stats = clipping.stats.as_ref().unwrap();

        // Expected values from old implementation
        assert!(
            (stats.clipping_area - 2.828427).abs() < 0.01,
            "clipping_area: expected 2.828427, got {}",
            stats.clipping_area
        );
        assert!(
            (stats.intersection_area - 2.828427).abs() < 0.01,
            "intersection_area: expected 2.828427, got {}",
            stats.intersection_area
        );
        assert!(
            (stats.total_consvtn - 1.0).abs() < 0.001,
            "total_consvtn: expected 1.0, got {}",
            stats.total_consvtn
        );

        // Check counts
        assert_eq!(clipping.intersections.len(), 3, "Expected 3 intersections");
        assert_eq!(clipping.remaining.len(), 0, "Expected 0 remaining");

        // Check vertex counts (allowing +1 for collinear points)
        assert!(clipping.intersections[0].data().num_vertices >= 4);
        assert!(clipping.intersections[1].data().num_vertices >= 4);
        assert!(clipping.intersections[2].data().num_vertices >= 4);

        // Check intersection midpoints exist (order-independent)
        let expected_midpoints = [(2.0515, 3.0, 0.0), (3.0, 2.909945, 0.0)];
        for (ex, ey, ez) in &expected_midpoints {
            let found = clipping.intersections.iter().any(|intsn| {
                let mid = intsn.data().midpoint;
                (mid.x - ex).abs() < 0.1 && (mid.y - ey).abs() < 0.1 && (mid.z - ez).abs() < 0.1
            });
            assert!(
                found,
                "Expected intersection midpoint ({}, {}, {}) not found",
                ex, ey, ez
            );
        }

        // Check that expected vertices exist in some intersection (order-independent)
        let expected_vertex = (2.922890, 3.0);
        let found = clipping.intersections.iter().any(|intsn| {
            intsn.data().exterior.iter().any(|v| {
                (v.x - expected_vertex.0).abs() < 0.01 && (v.y - expected_vertex.1).abs() < 0.01
            })
        });
        assert!(
            found,
            "Expected vertex ({}, {}) not found",
            expected_vertex.0, expected_vertex.1
        );
    }

    /// Test based on remainder.rs example
    /// multiple.obj with custom rectangular clip at z=10, projection (0, 0, -1)
    #[test]
    fn remainder_custom_clip() {
        let geoms = Geom::load("./examples/data/multiple.obj").unwrap();
        let mut geom = geoms[0].clone();

        assert_eq!(geom.shapes.len(), 2);

        let projection = Vector3::new(0.0, 0.0, -1.0);

        // Create custom clip rectangle (vertices reversed for correct normal)
        let mut clip_vertices = vec![
            Point3::new(-19.0, 3.0, 10.0),
            Point3::new(-19.0, -3.0, 10.0),
            Point3::new(10.0, -3.0, 10.0),
            Point3::new(10.0, 3.0, 10.0),
        ];
        clip_vertices.reverse();
        let mut clip = Face::new_simple(clip_vertices, None, None).unwrap();

        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip(AREA_THRESHOLD).unwrap();

        let stats = clipping.stats.as_ref().unwrap();

        // Expected values from old implementation
        assert!(
            (stats.clipping_area - 174.0).abs() < 0.1,
            "clipping_area: expected 174.0, got {}",
            stats.clipping_area
        );
        assert!(
            (stats.intersection_area - 107.045547).abs() < 0.5,
            "intersection_area: expected 107.045547, got {}",
            stats.intersection_area
        );
        assert!(
            (stats.remaining_area - 66.954453).abs() < 0.5,
            "remaining_area: expected 66.954453, got {}",
            stats.remaining_area
        );
        assert!(
            (stats.total_consvtn - 1.0).abs() < 0.001,
            "total_consvtn: expected 1.0, got {}",
            stats.total_consvtn
        );

        // Check counts
        assert_eq!(clipping.intersections.len(), 4, "Expected 4 intersections");
        assert_eq!(clipping.remaining.len(), 4, "Expected 4 remaining");

        // Check vertex counts (allowing +1 for collinear points)
        for intsn in &clipping.intersections {
            assert!(intsn.data().num_vertices >= 4);
        }
        for rem in &clipping.remaining {
            assert!(rem.data().num_vertices >= 3);
        }

        // Check intersection midpoints exist (order-independent)
        let expected_intersection_midpoints =
            [(-11.720262, 0.0, 7.021183), (2.165063, 0.0, 3.921344)];
        for (ex, ey, ez) in &expected_intersection_midpoints {
            let found = clipping.intersections.iter().any(|intsn| {
                let mid = intsn.data().midpoint;
                (mid.x - ex).abs() < 0.1 && (mid.y - ey).abs() < 0.1 && (mid.z - ez).abs() < 0.1
            });
            assert!(
                found,
                "Expected intersection midpoint ({}, {}, {}) not found",
                ex, ey, ez
            );
        }

        // Check remaining midpoints (all at z=10)
        for rem in &clipping.remaining {
            assert!(
                (rem.data().midpoint.z - 10.0).abs() < 0.01,
                "Remaining faces should be at z=10"
            );
        }

        // Check expected remaining midpoint exists (order-independent)
        let expected_remaining = (7.165063, 0.0, 10.0);
        let found = clipping.remaining.iter().any(|rem| {
            let mid = rem.data().midpoint;
            (mid.x - expected_remaining.0).abs() < 0.1
                && (mid.y - expected_remaining.1).abs() < 0.1
                && (mid.z - expected_remaining.2).abs() < 0.1
        });
        assert!(
            found,
            "Expected remaining midpoint ({}, {}, {}) not found",
            expected_remaining.0, expected_remaining.1, expected_remaining.2
        );
    }
}
trait Point3Extensions {
    fn ray_cast_z(&self, plane: &Plane) -> f32;
}

impl Point3Extensions for Point3<f32> {
    /// Returns the ray-cast distance along the -z axis from a point to its intersection with a plane in 3D
    fn ray_cast_z(&self, plane: &Plane) -> f32 {
        -(plane.normal.x * self.x + plane.normal.y * self.y + plane.offset) / plane.normal.z
            - self.z
    }
}

/// Statistics for a `Clipping` object.
#[derive(Debug, PartialEq, Clone, Default)] // Added Default derive
pub struct Stats {
    pub clipping_area: f32,     // the total input clipping area
    pub intersection_area: f32, // the total intersection area
    pub remaining_area: f32,    // the total remaining area
    pub consvtn: f32,           // the ratio of intersection to clipping area
    pub total_consvtn: f32,     // the ratio of (intersection + remaining) to clipping area
    pub area_loss: f32,         // the total area loss
}

impl Stats {
    pub fn new(clip: &Face, intersection: &Vec<Face>, remaining: &Vec<Face>) -> Self {
        let clipping_area = clip.to_polygon().unsigned_area();
        let intersection_area = intersection
            .iter()
            .fold(0.0, |acc, i| acc + i.to_polygon().unsigned_area());
        let remaining_area = remaining
            .iter()
            .fold(0.0, |acc, i| acc + i.to_polygon().unsigned_area());

        let consvtn = if clipping_area == 0.0 {
            0.0 // Avoid division by zero
        } else {
            intersection_area / clipping_area
        };

        let total_consvtn = if clipping_area == 0.0 {
            0.0 // Avoid division by zero
        } else {
            (intersection_area + remaining_area) / clipping_area
        };

        let area_loss = clipping_area - intersection_area - remaining_area;

        Self {
            clipping_area,
            intersection_area,
            remaining_area,
            consvtn,
            total_consvtn,
            area_loss,
        }
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Clipping Area: {}\nIntersection Area: {}\nRemaining Area: {}\nConservation (Intersection/Clipping): {}\nTotal Conservation ((Intersection + Remaining)/Clipping): {}",
            self.clipping_area,
            self.intersection_area,
            self.remaining_area,
            self.consvtn,
            self.total_consvtn
        )
    }
}

/// A clipping object
#[derive(Debug, PartialEq)]
pub struct Clipping<'a> {
    pub geom: &'a mut Geom,       // a geometry holding subjects to clip against
    pub clip: &'a mut Face,       // a clipping face
    pub proj: &'a Vector3<f32>,   // a projection vector
    pub intersections: Vec<Face>, // a list of intersection faces
    pub remaining: Vec<Face>,     // a list of remaining clips
    transform: Matrix4<f32>,      // a transform matrix to the clipping system
    itransform: Matrix4<f32>,     // a transform matrix from the clipping system
    is_done: bool,                // whether or not the clipping has been computed
    pub stats: Option<Stats>,     // statistics about the clipping result
}

impl<'a> Clipping<'a> {
    /// A new clipping object.
    /// If `clip` exists inside `geom`, it is ommitted from the subjects.
    pub fn new(geom: &'a mut Geom, clip: &'a mut Face, proj: &'a Vector3<f32>) -> Self {
        let mut clipping = Self {
            geom,
            clip,
            proj,
            intersections: Vec::new(),
            remaining: Vec::new(),
            transform: Matrix4::zeros(),
            itransform: Matrix4::zeros(),
            is_done: false,
            stats: None,
        };
        clipping.set_transform();

        clipping
    }

    /// Sets the forward and inverse transform for the clipping
    fn set_transform(&mut self) {
        let model = Isometry3::new(Vector3::zeros(), na::zero()); // do some sort of projection - set to nothing
        let origin = Point3::origin(); // camera location
        let target = Point3::new(self.proj.x, self.proj.y, self.proj.z); // projection direction, defines negative z-axis in new coords

        let up: Vector3<f32> =
            if self.proj.cross(&Vector3::y()).norm() < settings::COLINEAR_THRESHOLD {
                Vector3::x()
            } else {
                Vector3::y()
            };

        let view = Isometry3::look_at_rh(&origin, &target, &up);

        self.transform = (view * model).to_homogeneous(); // transform to clipping system
        self.itransform = self.transform.try_inverse().unwrap(); // inverse transform
    }

    pub fn init_clip(&mut self) -> Result<(&Face, Vec<&Face>)> {
        if self.is_done {
            panic!("Method clip() called, but the clipping was already done previously.");
        }

        self.geom.transform(&self.transform)?; // transform to clipping coordinate system
        self.clip.transform(&self.transform)?;

        let mut subjects = Vec::new();

        let clip_shape_id = self.clip.data().shape_id;
        let internal = if self.clip.data().normal.z > 0.0 {
            true
        } else {
            false
        };

        // create a mapping where each element links a subject to its shape and
        // face in the geometry
        for shape in self.geom.shapes.iter() {
            if internal && !shape.is_within(&self.geom, clip_shape_id) {
                continue;
            }

            for face in shape.faces.iter() {
                if face == self.clip {
                    // don't include the clip in the subjects
                    continue;
                }
                subjects.push(face);
            }
        }

        Ok((self.clip, subjects))
    }

    pub fn finalise_clip(
        &mut self,
        mut intersection: Vec<Face>,
        mut remaining: Vec<Face>,
    ) -> Result<()> {
        // transform back to original coordinate system
        self.geom.transform(&self.itransform)?;
        intersection
            .iter_mut()
            .try_for_each(|x| x.transform(&self.itransform))?;
        remaining
            .iter_mut()
            .try_for_each(|face| face.transform(&self.itransform))?;
        self.clip.transform(&self.itransform)?;

        // append the remapped intersections to the struct
        self.intersections.extend(intersection);
        self.remaining.extend(remaining);
        self.is_done = true;
        Ok(())
    }

    /// Performs the clip on a `Clipping` object.
    pub fn clip(&mut self, area_threshold: f32) -> Result<()> {
        if self.is_done {
            panic!("Method clip() called, but the clipping was already done previously.");
        }

        let (clip, mut subjects) = self.init_clip()?;

        // compute remapped intersections, converting to Intersection structs
        let (intersection, remaining) = clip_faces(&clip, &mut subjects, area_threshold)?;

        // compute statistics in clipping system
        self.set_stats(&intersection, &remaining);

        self.finalise_clip(intersection, remaining)?;

        Ok(())
    }

    fn set_stats(&mut self, intersection: &Vec<Face>, remaining: &Vec<Face>) {
        self.stats = Some(Stats::new(self.clip, intersection, remaining));
    }
}

/// Determines if a subject face can possibly intersect the clip face.
/// Returns false if:
/// - The subject is entirely behind the clip (in z-coordinate)
/// - The 2D bounding boxes (x,y) don't overlap (with tolerance)
fn can_subject_clip(subject: &Face, clip_in: &Face) -> bool {
    use crate::settings::constants::BBOX_TOLERANCE;

    // Check z-coordinate: subject must not be entirely behind clip
    let z_ok = match (subject.data().vert_min(2), clip_in.data().vert_max(2)) {
        (Ok(subj_min), Ok(clip_max)) => subj_min <= clip_max + BBOX_TOLERANCE,
        _ => return false,
    };
    if !z_ok {
        return false;
    }

    // Check 2D bounding box overlap (x and y dimensions)
    let subj_data = subject.data();
    let clip_data = clip_in.data();

    let (subj_x_min, subj_x_max) = match (subj_data.vert_min(0), subj_data.vert_max(0)) {
        (Ok(min), Ok(max)) => (min, max),
        _ => return true, // Can't determine, assume possible overlap
    };
    let (subj_y_min, subj_y_max) = match (subj_data.vert_min(1), subj_data.vert_max(1)) {
        (Ok(min), Ok(max)) => (min, max),
        _ => return true,
    };
    let (clip_x_min, clip_x_max) = match (clip_data.vert_min(0), clip_data.vert_max(0)) {
        (Ok(min), Ok(max)) => (min, max),
        _ => return true,
    };
    let (clip_y_min, clip_y_max) = match (clip_data.vert_min(1), clip_data.vert_max(1)) {
        (Ok(min), Ok(max)) => (min, max),
        _ => return true,
    };

    // Check if bounding boxes overlap (with tolerance for floating-point precision)
    let x_overlap =
        subj_x_min <= clip_x_max + BBOX_TOLERANCE && subj_x_max >= clip_x_min - BBOX_TOLERANCE;
    let y_overlap =
        subj_y_min <= clip_y_max + BBOX_TOLERANCE && subj_y_max >= clip_y_min - BBOX_TOLERANCE;

    x_overlap && y_overlap
}

/// Clips the `clip_in` against the `subjects_in`, in the current coordinate system.
pub fn clip_faces<'a>(
    clip_in: &Face,
    subjects_in: &Vec<&'a Face>,
    area_threshold: f32,
) -> Result<(Vec<Face>, Vec<Face>)> {
    if subjects_in.is_empty() {
        return Ok((Vec::new(), vec![clip_in.clone()]));
    }

    let clip_polygon = clip_in.to_polygon();
    let mut intersections = Vec::new();
    let mut remaining_clips = vec![clip_polygon];

    // Sort subjects by their Z-coordinate midpoint, descending.
    let sorted_subjects = {
        let mut subjects = subjects_in.clone();
        subjects.sort_by(|a, b| {
            b.midpoint()
                .z
                .partial_cmp(&a.midpoint().z)
                .unwrap_or(Ordering::Equal)
        });
        subjects
    };

    for subject in sorted_subjects
        .iter()
        .filter(|subj| can_subject_clip(subj, clip_in))
    {
        let subject_poly = subject.to_polygon();
        let mut next_clips = Vec::new();

        for clip in &remaining_clips {
            let mut intersection = Simplify::simplify(
                &subject_poly.intersection(clip),
                settings::VERTEX_MERGE_DISTANCE,
            );
            let mut difference = Simplify::simplify(
                &clip.difference(&subject_poly),
                settings::VERTEX_MERGE_DISTANCE,
            );

            // Retain only meaningful intersections and differences.
            intersection
                .0
                .retain(|f| f.unsigned_area() > area_threshold);
            difference.0.retain(|f| f.unsigned_area() > area_threshold);

            for poly in intersection.0.into_iter() {
                // try to project the polygon onto the subject plane
                let mut face = match poly.project(&subject.plane()) {
                    Ok(face) => face,
                    Err(_) => continue, // skip face if poly project failed
                };

                // cast a ray to determine if the intersection was in front
                if face.data().midpoint.ray_cast_z(&clip_in.plane())
                    > settings::RAYCAST_MINIMUM_DISTANCE
                {
                    face.data_mut().shape_id = subject.data().shape_id;
                    intersections.push(face);
                } else {
                    difference.0.push(poly);
                }
            }
            next_clips.extend(difference.0);
        }

        remaining_clips = next_clips;
        if remaining_clips.is_empty() {
            break;
        }
    }

    let remaining: Vec<_> = remaining_clips
        .into_iter()
        .map(|poly| poly.project(&clip_in.plane()))
        .collect::<Result<Vec<_>>>()?;

    Ok((intersections, remaining))
}
