use super::geom::{Face, Geom, Plane};
use super::settings;
use crate::geom::PolygonExtensions;
use anyhow::Result;
use geo::{Area, Simplify};
use geo_clipper::Clipper;
#[cfg(feature = "macroquad")]
use macroquad::prelude::*;
use nalgebra::{self as na, Isometry3, Matrix4, Point3, Vector3};
use std::cmp::Ordering;
use std::fmt;

#[cfg(test)]
mod tests {

    use geo::polygon;
    use geo::{CoordsIter, MultiPolygon, Simplify};

    use super::*;
    const AREA_THRESHOLD: f32 = 0.01;

    #[test]
    #[should_panic]
    fn concave_clip() {
        let mut geom = Geom::from_file("./examples/data/concave1.obj").unwrap();

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

        let cleaned = Simplify::simplify(&multipolygon, &0.01);

        // Print the cleaned polygon
        println!("Cleaned MultiPolygon: {:?}", cleaned);

        // Assert that the number of vertices in the cleaned exterior is 5
        let cleaned_exterior = &cleaned.0[0].exterior();
        assert_eq!(cleaned_exterior.coords_count(), 5);
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

    for subject in sorted_subjects.iter().filter(|subj| {
        match (subj.data().vert_min(2), clip_in.data().vert_max(2)) {
            (Ok(subj_min), Ok(clip_max)) => subj_min <= clip_max,
            _ => false,
        }
    }) {
        let subject_poly = subject.to_polygon();
        let mut next_clips = Vec::new();

        for clip in &remaining_clips {
            let mut intersection = Simplify::simplify(
                &subject_poly.intersection(clip, settings::CLIP_TOLERANCE),
                &settings::VERTEX_MERGE_DISTANCE,
            );
            let mut difference = Simplify::simplify(
                &clip.difference(&subject_poly, settings::CLIP_TOLERANCE),
                &settings::VERTEX_MERGE_DISTANCE,
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
