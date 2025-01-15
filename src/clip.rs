use super::geom::{Face, Geom, Plane};
use geo::Area;
use geo_clipper::Clipper;
use geo_types::Polygon;
use macroquad::prelude::*;
use nalgebra::{self as na, Isometry3, Point3, Vector3};
use std::cmp::Ordering;

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn concave_clip() {
        let mut geom = Geom::from_file("./examples/data/concave1.obj").unwrap();

        let clip_index = 4; // the index of the face to be used as the clip
        let projection = Vector3::new(-0.3, 0.0, -1.0);
        let mut clip = geom.shapes[0].faces.remove(clip_index); // choose a face be the clip

        // start function `do_clip` here:
        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip();

        panic!();
    }
}

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

/// A clipping object
#[derive(Debug, PartialEq)]
pub struct Clipping<'a> {
    pub geom: &'a mut Geom,       // a geometry holding subjects to clip against
    pub clip: &'a mut Face,       // a clipping face
    pub proj: &'a Vector3<f32>,   // a projection vector
    pub intersections: Vec<Face>, // a list of intersection faces
}

impl<'a> Clipping<'a> {
    /// A new cllipping object.
    pub fn new(geom: &'a mut Geom, clip: &'a mut Face, proj: &'a Vector3<f32>) -> Self {
        Self {
            geom: geom,
            clip: clip,
            proj: proj,
            intersections: Vec::new(),
        }
    }

    /// Performs the clip on a `Clipping` object.
    pub fn clip(&mut self) {
        let model = Isometry3::new(Vector3::zeros(), na::zero()); // do some sort of projection - set to nothing
        let origin = Point3::origin(); // camera location
        let target = Point3::new(self.proj.x, self.proj.y, self.proj.z); // projection direction, defines negative z-axis in new coords

        let up: Vector3<f32> = if self.proj.cross(&Vector3::y()).norm() < 0.01 {
            Vector3::x()
        } else {
            Vector3::y()
        };

        let view = Isometry3::look_at_rh(&origin, &target, &up);

        let transform = (view * model).to_homogeneous(); // transform to clipping system
        let itransform = transform.try_inverse().unwrap(); // inverse transform

        self.geom.transform(&transform); // transform to clipping coordinate system
        self.clip.transform(&transform);

        let subjects: Vec<&Face> = self
            .geom
            .shapes
            .iter()
            .flat_map(|f| f.faces.iter())
            .collect();

        // compute remapped intersections
        let mut remapped: Vec<Face> = clip_faces(&self.clip, &subjects);

        // transform back to original coordinate system
        self.geom.transform(&itransform);
        remapped
            .iter_mut()
            .for_each(|face| face.transform(&itransform));
        self.clip.transform(&itransform);

        // append the remapped intersections to the struct
        self.intersections.extend(remapped);
    }
}

const CLIP_TOLERANCE: f32 = 1e6; // Named constant for tolerance
const AREA_THRESHOLD: f32 = 1e-2; // Named constant for minimum intersection area
/// Clips the `clip_in` against the `subjects_in`, in the current coordinate system.
pub fn clip_faces(clip_in: &Face, subjects_in: &Vec<&Face>) -> Vec<Face> {
    if subjects_in.is_empty() {
        return Vec::new();
    }

    // Use a single pass to filter and sort by midpoint.z
    let mut subjects = subjects_in
        .iter()
        .filter(|subj| subj.midpoint.z < clip_in.midpoint.z)
        .copied()
        .collect::<Vec<_>>();
    subjects.sort_by(|a, b| {
        b.midpoint
            .z
            .partial_cmp(&a.midpoint.z)
            .unwrap_or(Ordering::Equal)
    });

    let clip = clip_in.polygon();
    let mut intersections = Vec::new();
    let mut remaining_clips = vec![clip];

    for subject in subjects {
        // subject is now &Face
        let subject_poly = subject.polygon(); // Or better: if Face has a &Polygon field, use that.
        let mut next_clips = Vec::new();

        for clip in &remaining_clips {
            let mut intersection = subject_poly.intersection(clip, CLIP_TOLERANCE);
            let difference = clip.difference(&subject_poly, CLIP_TOLERANCE);

            // filter out "bad" intersections
            intersection
                .0
                .retain(|f| f.unsigned_area() > AREA_THRESHOLD);

            intersections.extend(
                intersection
                    .into_iter()
                    .map(|intsn| intsn.project(&subject.plane())),
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
