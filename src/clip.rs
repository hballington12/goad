use super::geom::{Face, Geom, Plane};
use super::settings;
use crate::geom::PolygonExtensions;
use anyhow::Result;
use geo::{Area, Simplify};
use geo_clipper::Clipper;
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
    /// Computes ray casting distance for occlusion testing in clipping operations.
    /// 
    /// **Context**: During clipping, intersections must be tested to determine
    /// whether they occur in front of or behind the clipping plane. This
    /// ray casting operation enables proper depth ordering.
    /// 
    /// **How it Works**: Calculates the z-distance from a point to its projection
    /// onto a plane along the negative z-axis direction used in clipping coordinates.
    fn ray_cast_z(&self, plane: &Plane) -> f32 {
        -(plane.normal.x * self.x + plane.normal.y * self.y + plane.offset) / plane.normal.z
            - self.z
    }
}

/// Area conservation statistics for geometric clipping operations.
/// 
/// **Context**: Geometric clipping can introduce numerical errors that result
/// in area loss or gain. Tracking area conservation is essential for power
/// conservation verification in electromagnetic simulations where geometric
/// cross-sections determine beam power.
/// 
/// **How it Works**: Computes area statistics by comparing the original clipping
/// face area with the sum of intersection and remaining areas after clipping.
/// Conservation ratios indicate the fraction of area properly accounted for
/// versus numerical losses from polygon operations.
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
    /// Computes area conservation statistics for a clipping operation.
    /// 
    /// **Context**: Geometric clipping algorithms can introduce numerical errors
    /// that cause area loss or gain. These statistics quantify conservation
    /// to verify that beam power is properly tracked through the simulation.
    /// 
    /// **How it Works**: Sums areas of intersection and remaining faces,
    /// computes conservation ratios relative to the original clipping area,
    /// and calculates total area loss from numerical operations.
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

/// Geometric clipping engine for beam-surface intersection calculations.
/// 
/// **Context**: Electromagnetic beam propagation requires determining which portions
/// of beam cross-sections intersect with particle surfaces and which portions
/// continue propagating. This geometric clipping operation is fundamental to
/// the geometric optics approach, where beams are discretized into polygonal
/// cross-sections that interact with particle faces.
/// 
/// **How it Works**: Transforms both the beam cross-section (clip) and particle
/// geometry into a coordinate system aligned with the beam propagation direction.
/// Uses 2D polygon clipping algorithms to compute intersections between the beam
/// and each relevant particle face, sorted by depth. Tracks area conservation
/// and handles complex polygon operations including holes and non-convex shapes.
/// Returns intersection faces (portions interacting with surfaces) and remaining
/// faces (portions continuing propagation).
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
    /// Creates a clipping operation for beam-surface intersection analysis.
    /// 
    /// **Context**: Each beam propagation step requires determining intersections
    /// with particle surfaces along the propagation direction. The clipping
    /// operation must be set up with appropriate coordinate transformations
    /// for numerical stability.
    /// 
    /// **How it Works**: Establishes coordinate transformation matrices to align
    /// the clipping face with the xy-plane and propagation direction with the
    /// negative z-axis. This reduces the 3D clipping problem to efficient 2D
    /// polygon operations.
    /// 
    /// # Example
    /// ```rust
    /// let mut clipping = Clipping::new(geom, &mut self.face, &self.prop);
    /// clipping.clip(area_threshold)?;
    /// ```
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

    /// Establishes coordinate transformation for optimal clipping geometry.
    /// 
    /// **Context**: 2D polygon clipping algorithms are most efficient and
    /// numerically stable when the clipping face lies in the xy-plane.
    /// The transformation aligns the coordinate system with the beam
    /// propagation direction.
    /// 
    /// **How it Works**: Creates a look-at transformation that places the
    /// clipping face in the xy-plane with the propagation direction along
    /// the negative z-axis. Handles colinear cases by choosing appropriate
    /// up vectors.
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

    /// Initializes clipping operation by transforming geometry and selecting subjects.
    /// 
    /// **Context**: Not all faces in the geometry are relevant for clipping
    /// against a given beam. Subject selection depends on containment relationships
    /// and the beam's interaction with internal versus external surfaces.
    /// 
    /// **How it Works**: Transforms geometry to clipping coordinates, determines
    /// whether the beam is internal or external based on face normal direction,
    /// and selects appropriate subject faces based on containment rules.
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

    /// Transforms clipping results back to world coordinates and finalizes operation.
    /// 
    /// **Context**: Clipping operations are performed in the transformed coordinate
    /// system for numerical efficiency. Results must be transformed back to
    /// world coordinates for use in beam propagation.
    /// 
    /// **How it Works**: Applies inverse transformation to geometry, intersection
    /// faces, remaining faces, and the clipping face itself. Stores results
    /// in the clipping structure and marks the operation as complete.
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

    /// Executes the geometric clipping operation to find beam-surface intersections.
    /// 
    /// **Context**: This is the main clipping operation that determines which portions
    /// of a beam cross-section intersect with particle surfaces and which portions
    /// remain for continued propagation. Area thresholding prevents processing
    /// of tiny polygon fragments that contribute negligible power.
    /// 
    /// **How it Works**: Transforms geometry to clipping coordinates, identifies
    /// relevant subject faces based on containment rules, performs depth-sorted
    /// polygon clipping operations, computes area statistics, and transforms
    /// results back to world coordinates.
    /// 
    /// # Example
    /// ```rust
    /// let mut clipping = Clipping::new(&mut geom, &mut self.face, &self.prop);
    /// clipping.clip(area_threshold)?;
    /// let (intersections, remainders) = (
    ///     clipping.intersections.into_iter().collect(),
    ///     clipping.remaining.into_iter().collect(),
    /// );
    /// ```
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

    /// Computes and stores area conservation statistics for the clipping operation.
    /// 
    /// **Context**: Area conservation statistics are computed in the clipping
    /// coordinate system before transformation back to world coordinates
    /// to avoid transformation-induced numerical errors.
    /// 
    /// **How it Works**: Creates a Stats object from the clipping face and
    /// computed intersection and remaining faces.
    fn set_stats(&mut self, intersection: &Vec<Face>, remaining: &Vec<Face>) {
        self.stats = Some(Stats::new(self.clip, intersection, remaining));
    }
}

/// Performs 2D polygon clipping between beam cross-section and particle faces.
/// 
/// **Context**: Once transformed to the clipping coordinate system, the 3D
/// beam-surface intersection problem reduces to 2D polygon clipping operations.
/// Multiple particle faces may intersect the beam, requiring depth-sorted
/// processing to handle occlusion correctly.
/// 
/// **How it Works**: Sorts subject faces by depth, iteratively clips the remaining
/// beam area against each face using Boolean polygon operations (intersection
/// and difference). Applies area thresholding to filter insignificant fragments.
/// Uses ray casting to distinguish intersections in front of the beam from
/// those behind, ensuring correct occlusion handling.
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
