use super::config;
use super::geom::{Face, Geom, Intersection, Plane};
use geo::{Area, Point};
use geo_clipper::Clipper;
use geo_types::Coord;
use geo_types::Polygon;
use macroquad::prelude::*;
use nalgebra::{self as na, Isometry3, Matrix4, Point3, Vector3};
use std::cmp::Ordering;
use std::fmt;
use std::iter::repeat;

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    #[should_panic]
    fn concave_clip() {
        let mut geom = Geom::from_file("./examples/data/concave1.obj").unwrap();

        let clip_index = 4; // the index of the face to be used as the clip
        let projection = Vector3::new(-0.3, 0.0, -1.0);
        let mut clip = geom.shapes[0].faces.remove(clip_index); // choose a face be the clip

        // start function `do_clip` here:
        let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
        clipping.clip();
        clipping.clip(); // try to clip again, which should panic
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
trait Coord3Extensions {
    fn projected_z(&self, plane: &Plane) -> f32;
}
impl Coord3Extensions for Coord<f32> {
    /// Returns the z-coordinate of a `Coord` projected onto a plane in 3D
    fn projected_z(&self, plane: &Plane) -> f32 {
        -(plane.normal.x * self.x + plane.normal.y * self.y + plane.offset) / plane.normal.z
    }
}

trait PolygonExtensions {
    fn project(&self, plane: &Plane) -> Face;
}

impl PolygonExtensions for Polygon<f32> {
    /// Projects the xy coordinates of a polygon onto a plane in 3D
    ///  the last vertex, which is a duplicate of the first
    fn project(&self, plane: &Plane) -> Face {
        let area = self.unsigned_area() / plane.normal.z;

        let project_coords = |coords: &Vec<Coord<f32>>| -> Vec<Point3<f32>> {
            coords
                .iter()
                .take(coords.len() - 1)
                .map(|coord| Point3::new(coord.x, coord.y, coord.projected_z(plane)))
                .collect()
        };

        let exterior = project_coords(&self.exterior().0);

        if self.interiors().is_empty() {
            let mut face = Face::new_simple(exterior, None);
            face.set_area(area);
            face
        } else {
            let interiors = self
                .interiors()
                .iter()
                .map(|interior| project_coords(&interior.0))
                .collect();
            let mut face = Face::new_complex(exterior, interiors, None);
            face.set_area(area);
            face
        }
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
}

impl Stats {
    pub fn new(clip: &Face, intersection: &Vec<Intersection>, remaining: &Vec<Face>) -> Self {
        let clipping_area = clip.to_polygon().unsigned_area();
        let intersection_area = intersection
            .iter()
            .fold(0.0, |acc, i| acc + i.face.to_polygon().unsigned_area());
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
        Self {
            clipping_area,
            intersection_area,
            remaining_area,
            consvtn,
            total_consvtn,
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
    pub geom: &'a mut Geom,     // a geometry holding subjects to clip against
    pub clip: &'a mut Face,     // a clipping face
    pub proj: &'a Vector3<f32>, // a projection vector
    pub intersections: Vec<Intersection>, // a list of intersection faces
    pub source_mapping: Vec<(usize, usize)>, // a list of references to the source of each intersection
    pub remaining: Vec<Face>,                // a list of remaining clips
    transform: Matrix4<f32>,                 // a transform matrix to the clipping system
    itransform: Matrix4<f32>,                // a transform matrix from the clipping system
    is_done: bool,                           // whether or not the clipping has been computed
    pub stats: Option<Stats>,                // statistics about the clipping result
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
            source_mapping: Vec::new(),
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

        let up: Vector3<f32> = if self.proj.cross(&Vector3::y()).norm() < 0.01 {
            Vector3::x()
        } else {
            Vector3::y()
        };

        let view = Isometry3::look_at_rh(&origin, &target, &up);

        self.transform = (view * model).to_homogeneous(); // transform to clipping system
        self.itransform = self.transform.try_inverse().unwrap(); // inverse transform
    }

    pub fn init_clip(&mut self) -> (&Face, Vec<&Face>, Vec<(usize, usize)>) {
        if self.is_done {
            panic!("Method clip() called, but the clipping was already done previously.");
        }

        self.geom.transform(&self.transform); // transform to clipping coordinate system
        self.clip.transform(&self.transform);

        let mut subjects = Vec::new();
        let mut mapping = Vec::new();

        for (i, shape) in self.geom.shapes.iter().enumerate() {
            for (j, face) in shape.faces.iter().enumerate() {
                if face == self.clip {
                    // don't include the clip in the subjects
                    continue;
                }
                subjects.push(face);
                mapping.push((i, j));
            }
        }

        (self.clip, subjects, mapping)
    }

    pub fn finalise_clip(
        &mut self,
        mut intersection: Vec<Intersection>,
        mut remaining: Vec<Face>,
        source_mapping: Vec<(usize, usize)>,
    ) {
        // transform back to original coordinate system
        self.geom.transform(&self.itransform);
        intersection
            .iter_mut()
            .for_each(|x| x.face.transform(&self.itransform));
        remaining
            .iter_mut()
            .for_each(|face| face.transform(&self.itransform));
        self.clip.transform(&self.itransform);

        // append the remapped intersections to the struct
        self.intersections.extend(intersection);
        self.remaining.extend(remaining);
        self.source_mapping = source_mapping;
        self.is_done = true;
    }

    /// Performs the clip on a `Clipping` object.
    pub fn clip(&mut self) {
        if self.is_done {
            panic!("Method clip() called, but the clipping was already done previously.");
        }

        let (clip, mut subjects, mapping) = self.init_clip();

        // compute remapped intersections, converting to Intersection structs
        let (intersection, remaining, sources) = clip_faces(&clip, &mut subjects);
        // get mapping back to geometry
        let source_mapping: Vec<_> = sources.iter().map(|&i| mapping[i]).collect();
        let intersection: Vec<Intersection> = intersection
            .into_iter()
            .enumerate()
            .map(|f| Intersection::new(f.1, mapping[f.0].0))
            .collect();

        // compute statistics in clipping system
        self.set_stats(&intersection, &remaining);

        self.finalise_clip(intersection, remaining, source_mapping);
    }

    fn set_stats(&mut self, intersection: &Vec<Intersection>, remaining: &Vec<Face>) {
        self.stats = Some(Stats::new(self.clip, intersection, remaining));
    }
}

const AREA_THRESHOLD: f32 = 1e-2;

/// Clips the `clip_in` against the `subjects_in`, in the current coordinate system.
pub fn clip_faces<'a>(
    clip_in: &Face,
    subjects_in: &mut Vec<&'a Face>,
) -> (Vec<Face>, Vec<Face>, Vec<usize>) {
    if subjects_in.is_empty() {
        return (Vec::new(), vec![clip_in.clone()], Vec::new());
    }

    let clip_polygon = clip_in.to_polygon();
    let mut intersections = Vec::new();
    let mut remaining_clips = vec![clip_polygon];
    let mut intsn_sources = Vec::new();

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

    for (i, subject) in sorted_subjects
        .into_iter()
        .enumerate()
        .filter(|(_, subj)| subj.data().vert_min(2) <= clip_in.data().vert_max(2))
    {
        let subject_poly = subject.to_polygon();
        let mut next_clips = Vec::new();

        for clip in &remaining_clips {
            let mut intersection = subject_poly.intersection(clip, config::CLIP_TOLERANCE);
            let mut difference = clip.difference(&subject_poly, config::CLIP_TOLERANCE);

            if intersection.0.iter().any(|poly| {
                let face = poly.project(&subject.plane());
                face.data().midpoint.ray_cast_z(&clip_in.plane()) < 0.0
            }) {
                // Include intersections that are unphysical back into the difference.
                difference.0.extend(intersection.0);
            } else {
                // Retain only meaningful intersections and differences.
                intersection
                    .0
                    .retain(|f| f.unsigned_area() > AREA_THRESHOLD);
                difference.0.retain(|f| f.unsigned_area() > AREA_THRESHOLD);

                intsn_sources.extend(std::iter::repeat(i).take(intersection.0.len()));
                intersections.extend(
                    intersection
                        .0
                        .into_iter()
                        .map(|poly| poly.project(&subject.plane())),
                );
            }

            next_clips.extend(difference.0);
        }

        remaining_clips = next_clips;
        if remaining_clips.is_empty() {
            break;
        }
    }

    let remaining = remaining_clips
        .into_iter()
        .map(|poly| poly.project(&clip_in.plane()))
        .collect();

    (intersections, remaining, intsn_sources)
}
