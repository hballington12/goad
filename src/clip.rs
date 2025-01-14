use super::geom::{Face, Plane};
use geo_clipper::Clipper;
use geo_types::Polygon;
use macroquad::prelude::*;
use nalgebra::Point3;
use std::cmp::Ordering;

trait PolygonExtensions {
    fn project(&self, plane: &Plane) -> Face;
}

impl PolygonExtensions for Polygon<f32> {
    /// Projects the xy coordinates of a polygon onto a plane in 3D
    /// Ignores the last vertex, which is a duplicate of the first
    fn project(&self, plane: &Plane) -> Face {
        let mut vertices = Vec::new();

        #[cfg(debug_assertions)]
        {
            assert!(
                !self.interiors().len() != 0,
                "hold detected, please fix: {:?}",
                self.interiors()
            );

            assert!(
                plane.normal.z != 0.0,
                "plane normal cannot have 0 z-component: {:?}",
                plane.normal
            );
        }

        for coord in self.exterior().0.iter().take(self.exterior().0.len() - 1) {
            // compute z intersection via z = -(ax + by + d)/c
            let z = -(plane.normal.x * coord.x + plane.normal.y * coord.y + plane.offset)
                / plane.normal.z;

            vertices.push(Point3::new(coord.x, coord.y, z));
        }

        Face::new(vertices)
    }
}

/// Clips the `clip_in` against the `subjects_in`, in the current coordinate system.
pub fn clip_faces(clip_in: &Face, subjects_in: &[Face]) -> Vec<Face> {
    if subjects_in.is_empty() {
        return Vec::new();
    }

    // Filter and sort subjects by maximum vertex z-coordinate.
    let mut subjects: Vec<&Face> = subjects_in
        .iter()
        .filter(|subj| subj.vert_max(2) < clip_in.vert_max(2))
        .collect();

    subjects.sort_by(|&b, &a| {
        a.midpoint
            .z // may wish to interchange vert_max(2) <-> midpoint.z at some point
            .partial_cmp(&b.midpoint.z)
            .unwrap_or(Ordering::Equal)
    });

    let clip = clip_in.polygon();
    let mut intersections = Vec::new();
    let mut remaining_clips = vec![clip];

    for subject_in in subjects {
        // Iterate by reference to avoid cloning
        let subject = subject_in.polygon();
        let mut next_clips = Vec::new();

        for clip in &remaining_clips {
            // Iterate by reference
            let intersection = subject.intersection(clip, 100000.0);
            let difference = clip.difference(&subject, 100000.0);

            intersections.extend(
                intersection
                    .into_iter()
                    .map(|intsn| intsn.project(&subject_in.plane())),
            );
            next_clips.extend(difference);
        }

        remaining_clips = next_clips;
        if remaining_clips.is_empty() {
            break;
        }
    }

    intersections
}
