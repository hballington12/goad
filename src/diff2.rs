//! Reimplemented near-to-far field aperture diffraction.

use anyhow::Result;
use nalgebra::{Complex, Matrix3, Matrix4, Point3, Vector3};

use crate::beam::Beam;
use crate::bins::SolidAngleBin;
use crate::field::Ampl;
#[cfg(debug_assertions)]
use crate::geom::Face;
use crate::geom::{self};
use crate::settings::constants::COLINEAR_THRESHOLD;
#[cfg(debug_assertions)]
use crate::settings::constants::PLANARITY_TOLERANCE;

/// Incident beam parameters for the prerotation matrix.
/// These define the reference frame of the original illumination source.
pub struct IncidentBeam {
    /// Perpendicular polarisation direction of incident beam
    pub e_perp: Vector3<f32>,
    /// Propagation direction of incident beam
    pub prop: Vector3<f32>,
}

impl Default for IncidentBeam {
    /// Default incident beam: e_perp along +y, prop along -z
    fn default() -> Self {
        Self {
            e_perp: Vector3::new(0.0, 1.0, 0.0),
            prop: Vector3::new(0.0, 0.0, -1.0),
        }
    }
}

// =============================================================================
// Debug assertion helpers
// =============================================================================

/// Panic if the face vertices are not coplanar.
#[cfg(debug_assertions)]
fn assert_face_planar(verts: &[Point3<f32>]) {
    if verts.len() < 3 {
        panic!("Face must have at least 3 vertices, got {}", verts.len());
    }

    let v0: Vector3<f32> = verts[1] - verts[0];
    let v1: Vector3<f32> = verts[2] - verts[0];
    let normal = v0.cross(&v1).normalize();

    for (i, vert) in verts.iter().enumerate().skip(3) {
        let v = vert - verts[0];
        let dist = v.dot(&normal).abs();
        if dist > PLANARITY_TOLERANCE {
            panic!(
                "Face is non-planar: vertex {} deviates by {} from plane (tolerance: {})",
                i, dist, PLANARITY_TOLERANCE
            );
        }
    }
}

/// Panic if the face is not a Simple variant.
#[cfg(debug_assertions)]
fn assert_face_simple(face: &Face) {
    match face {
        Face::Simple(_) => {}
        Face::Complex { .. } => {
            panic!("Face must be Simple variant, got Complex");
        }
    }
}

/// Panic if e_perp is not perpendicular to prop.
#[cfg(debug_assertions)]
fn assert_e_perp_perpendicular_to_prop(e_perp: Vector3<f32>, prop: Vector3<f32>) {
    let dot = e_perp.dot(&prop).abs();
    if dot > COLINEAR_THRESHOLD {
        panic!(
            "e_perp is not perpendicular to prop (|e_perp · prop| = {} > {})",
            dot, COLINEAR_THRESHOLD
        );
    }
}

/// Panic if e_perp is not perpendicular to the face normal.
#[cfg(debug_assertions)]
fn assert_e_perp_perpendicular_to_normal(e_perp: Vector3<f32>, normal: Vector3<f32>) {
    let dot = e_perp.dot(&normal).abs();
    if dot > COLINEAR_THRESHOLD {
        panic!(
            "e_perp is not perpendicular to face normal (|e_perp · normal| = {} > {})",
            dot, COLINEAR_THRESHOLD
        );
    }
}

// =============================================================================
// Transform helper functions
// =============================================================================

/// Build a 4x4 translation matrix that moves the given point to the origin.
fn translation_to_origin(center_of_mass: &Point3<f32>) -> Matrix4<f32> {
    Matrix4::new_translation(&(-center_of_mass.coords))
}

/// Build a 4x4 rotation matrix that rotates around z-axis to put e_perp along +y.
/// Takes the xy-plane beam (face already in xy plane).
/// Since e_perp ⊥ prop, this also ensures prop lies in the xz plane.
/// Using e_perp instead of prop avoids degeneracy at normal incidence.
#[allow(dead_code)]
fn rotation_e_perp_to_y(beam: &Beam) -> Matrix4<f32> {
    let e_perp = beam.field.e_perp();
    // Rotate so that e_perp aligns with +y axis
    // e_perp is in xy plane at angle atan2(x, y) from +y axis
    // Rotate by that angle (positive) to bring it to +y
    let angle = e_perp.x.atan2(e_perp.y);
    let (sin_angle, cos_angle) = angle.sin_cos();

    let rot3 = Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    );
    rot3.to_homogeneous()
}

/// Build a 4x4 rotation matrix that rotates the aperture into the xy plane.
/// Takes the centered beam (face already at origin).
fn rotation_to_xy_plane(beam: &Beam) -> Matrix4<f32> {
    #[cfg(debug_assertions)]
    {
        let midpoint = beam.face.midpoint();
        let dist = midpoint.coords.norm();
        assert!(
            dist < PLANARITY_TOLERANCE,
            "Beam face midpoint must be at origin, got distance {} from origin",
            dist
        );
    }

    let verts: Vec<Vector3<f32>> = beam.face.data().exterior.iter().map(|p| p.coords).collect();
    let rot3 = get_rotation_matrix(&verts);
    rot3.to_homogeneous()
}

// =============================================================================
// Fraunhofer integral pre-computation
// =============================================================================

use std::f32::consts::PI;

/// Pre-computed edge data for the Fraunhofer integral.
/// Each edge of the aperture polygon has associated slopes and adjusted values.
struct EdgeData {
    x: Vec<f32>,     // x coordinates of vertices
    y: Vec<f32>,     // y coordinates of vertices
    dx: Vec<f32>,    // delta x for each edge
    dy: Vec<f32>,    // delta y for each edge
    m_adj: Vec<f32>, // adjusted slope (clamped for numerical stability)
    n_adj: Vec<f32>, // adjusted inverse slope
}

impl EdgeData {
    /// Compute edge data from aperture vertices (already in xy plane).
    fn from_vertices(verts: &[Point3<f32>]) -> Self {
        let nv = verts.len();

        let x: Vec<f32> = verts.iter().map(|v| v.x).collect();
        let y: Vec<f32> = verts.iter().map(|v| v.y).collect();

        let mut dx_vec = Vec::with_capacity(nv);
        let mut dy_vec = Vec::with_capacity(nv);
        let mut m = Vec::with_capacity(nv);
        let mut n = Vec::with_capacity(nv);
        let mut m_adj = Vec::with_capacity(nv);
        let mut n_adj = Vec::with_capacity(nv);

        for j in 0..nv {
            let next_j = (j + 1) % nv;
            let mut dx = x[next_j] - x[j];
            let mut dy = y[next_j] - y[j];

            // Calculate slope, handling near-zero dx
            let mj = if dx.abs() < crate::settings::DIFF_DMIN {
                if dy.signum() == dx.signum() {
                    1e6
                } else {
                    -1e6
                }
            } else {
                dy / dx
            };
            m.push(mj);

            // Calculate inverse slope
            let nj = if mj.abs() < 1e-6 {
                if mj.signum() > 0.0 {
                    1e6
                } else {
                    -1e6
                }
            } else {
                1.0 / mj
            };
            n.push(nj);

            // Adjust dx/dy for numerical stability
            dx = if dx.abs() < crate::settings::DIFF_DMIN {
                crate::settings::DIFF_DMIN * dx.signum()
            } else {
                dx
            };
            dy = if dy.abs() < crate::settings::DIFF_DMIN {
                crate::settings::DIFF_DMIN * dy.signum()
            } else {
                dy
            };
            dx_vec.push(dx);
            dy_vec.push(dy);

            // Pre-calculate adjusted m and n
            let (adj_mj, adj_nj) = adjust_mj_nj(mj, nj);
            m_adj.push(adj_mj);
            n_adj.push(adj_nj);
        }

        Self {
            x,
            y,
            dx: dx_vec,
            dy: dy_vec,
            m_adj,
            n_adj,
        }
    }
}

/// Adjust m and n values for numerical stability.
#[inline]
fn adjust_mj_nj(mj: f32, nj: f32) -> (f32, f32) {
    if mj.abs() > 1e6 || nj.abs() < 1e-6 {
        (1e6, 1e-6)
    } else if nj.abs() > 1e6 || mj.abs() < 1e-6 {
        (1e-6, 1e6)
    } else {
        (mj, nj)
    }
}

/// Calculate the field of view cosine threshold for filtering.
fn calculate_fov_cosine(verts: &[Point3<f32>], wavenumber: f32, fov_factor: Option<f32>) -> f32 {
    let aperture_dimension = verts.iter().map(|v| v.coords.norm()).fold(0.0, f32::max);
    (fov_factor.unwrap_or(1.0) * 2.0 * 2.0 * PI / (wavenumber * aperture_dimension)).cos()
}

// =============================================================================
// Fraunhofer integral helper functions
// =============================================================================

use crate::settings;
use nalgebra::Matrix2;

#[inline]
fn calculate_kxx_kyy(kinc: &[f32; 2], k: &Vector3<f32>, wavenumber: f32) -> (f32, f32) {
    let kxx = kinc[0] - wavenumber * k.x;
    let kyy = kinc[1] - wavenumber * k.y;

    let kxx = if kxx.abs() < settings::KXY_EPSILON {
        settings::KXY_EPSILON
    } else {
        kxx
    };
    let kyy = if kyy.abs() < settings::KXY_EPSILON {
        settings::KXY_EPSILON
    } else {
        kyy
    };

    (kxx, kyy)
}

#[inline]
fn calculate_deltas(kxx: f32, kyy: f32, xj: f32, yj: f32, mj: f32, nj: f32) -> (f32, f32, f32) {
    let delta = kxx * xj + kyy * yj;
    let delta1 = kyy * mj + kxx;
    let delta2 = kxx * nj + kyy;
    (delta, delta1, delta2)
}

#[inline]
fn calculate_omegas(dx: f32, dy: f32, delta1: f32, delta2: f32) -> (f32, f32) {
    let omega1 = dx * delta1;
    let omega2 = dy * delta2;
    (omega1, omega2)
}

#[inline]
fn calculate_alpha_beta(delta1: f32, delta2: f32, kxx: f32, kyy: f32) -> (f32, f32) {
    let alpha = 1.0 / (2.0 * kyy * delta1);
    let beta = 1.0 / (2.0 * kxx * delta2);
    (alpha, beta)
}

#[inline]
fn calculate_summand(
    bvsk: f32,
    delta: f32,
    omega1: f32,
    omega2: f32,
    alpha: f32,
    beta: f32,
    inv_denom: Complex<f32>,
) -> Complex<f32> {
    let (sin_delta, cos_delta) = delta.sin_cos();
    let (sin_delta_omega1, cos_delta_omega1) = (delta + omega1).sin_cos();
    let (sin_delta_omega2, cos_delta_omega2) = (delta + omega2).sin_cos();

    let sumim = alpha * (cos_delta - cos_delta_omega1) - beta * (cos_delta - cos_delta_omega2);
    let sumre = -alpha * (sin_delta - sin_delta_omega1) + beta * (sin_delta - sin_delta_omega2);

    let exp_factor = Complex::cis(bvsk);

    exp_factor * Complex::new(sumre, sumim) * inv_denom
}

// =============================================================================
// Karczewski polarisation matrix
// =============================================================================

/// Compute the Karczewski polarisation matrix for aperture diffraction.
///
/// # Arguments
/// * `prop` - Propagation vector of the beam in aperture system
/// * `k` - Rotated observation direction (bin vector) in aperture system
///
/// # Returns
/// * Karczewski polarisation matrix (2x2)
/// * `m` - Vector perpendicular to the scattering plane
pub fn karczewski(prop: &Vector3<f32>, k: &Vector3<f32>) -> (Matrix2<f32>, Vector3<f32>) {
    let big_kx = prop.x;
    let big_ky = prop.y;
    let big_kz = prop.z;

    let one_minus_k2y2 = (1.0 - k.y.powi(2)).max(0.0);
    let sqrt_1_minus_k2y2 = one_minus_k2y2.sqrt();
    let sqrt_1_minus_k2y2 = if sqrt_1_minus_k2y2.abs() < settings::DIFF_EPSILON {
        settings::DIFF_EPSILON
    } else {
        sqrt_1_minus_k2y2
    };

    let m = Vector3::new(
        -k.x * k.y / sqrt_1_minus_k2y2,
        sqrt_1_minus_k2y2,
        -k.y * k.z / sqrt_1_minus_k2y2,
    );

    let frac = (one_minus_k2y2 / (1.0 - big_ky.powi(2))).sqrt();
    let frac = if frac.abs() < settings::DIFF_EPSILON {
        settings::DIFF_EPSILON
    } else {
        frac
    };

    let a1m = -big_kz * frac;
    let b2m = -k.z / frac;
    let a1e = b2m;
    let b2e = a1m;
    let b1m = -k.x * k.y / frac + big_kx * big_ky * frac;
    let a2e = -b1m;

    let a1em = 0.5 * (a1m + a1e);
    let a2em = 0.5 * a2e;
    let b1em = 0.5 * b1m;
    let b2em = 0.5 * (b2m + b2e);

    let diff_ampl = Matrix2::new(a1em, b1em, a2em, b2em);

    (diff_ampl, m)
}

// =============================================================================
// Rotation matrix computation
// =============================================================================

/// Compute a 3x3 rotation matrix for rotating vertices into the xy plane.
///
/// This follows the Fortran implementation: three sequential rotations
/// (about z, x, then y axes) to bring the aperture into the xy plane,
/// with an optional 180° flip about y if the face is pointing downward.
pub fn get_rotation_matrix(verts: &[Vector3<f32>]) -> Matrix3<f32> {
    let a1 = verts[0];
    let b1 = verts[1];

    // First rotation about z-axis
    let theta1 = a1.x.atan2(a1.y);
    let (sin1, cos1) = theta1.sin_cos();
    let rot1 = Matrix3::new(cos1, -sin1, 0.0, sin1, cos1, 0.0, 0.0, 0.0, 1.0);

    let a2 = rot1 * a1;
    let b2 = rot1 * b1;

    // Second rotation about x-axis
    let theta2 = -a2.z.atan2(a2.y);
    let (sin2, cos2) = theta2.sin_cos();
    let rot2 = Matrix3::new(1.0, 0.0, 0.0, 0.0, cos2, -sin2, 0.0, sin2, cos2);

    let a3 = rot2 * a2;
    let b3 = rot2 * b2;

    // Third rotation about y-axis
    let theta3 = b3.z.atan2(b3.x);
    let (sin3, cos3) = theta3.sin_cos();
    let rot3 = Matrix3::new(cos3, 0.0, sin3, 0.0, 1.0, 0.0, -sin3, 0.0, cos3);

    let a4 = rot3 * a3;
    let b4 = rot3 * b3;

    // Check if face is pointing down (cross product z-component > 0)
    // If so, flip 180° about y-axis
    if a4.x * b4.y - a4.y * b4.x > 0.0 {
        let rot4 = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
        rot4 * rot3 * rot2 * rot1
    } else {
        rot3 * rot2 * rot1
    }
}

// =============================================================================
// Main function
// =============================================================================

/// Perform near-to-far field aperture diffraction.
///
/// # Arguments
/// * `beam` - The beam to diffract
/// * `bins` - List of solid angle bins to compute diffraction for
/// * `incident` - Incident beam parameters (e_perp and prop of original illumination)
/// * `fov_factor` - Optional field of view factor for filtering
///
/// # Returns
/// Vector of (bin_index, amplitude_matrix) pairs
pub fn n2f_aperture_diffraction(
    beam: &Beam,
    bins: &[SolidAngleBin],
    incident: &IncidentBeam,
    fov_factor: Option<f32>,
) -> Result<Vec<(usize, Ampl)>> {
    let verts = &beam.face.data().exterior;

    // Apply small perturbation to prop to reduce numerical errors (matches original diff.rs)
    let prop = (beam.field.prop()
        + Vector3::new(
            settings::constants::PROP_PERTURBATION,
            settings::constants::PROP_PERTURBATION,
            settings::constants::PROP_PERTURBATION,
        ))
    .normalize();

    #[cfg(debug_assertions)]
    {
        let e_perp = beam.field.e_perp();
        let normal = beam.face.data().normal;
        if beam.field.prop().dot(&beam.face.data().normal) < 0.0 {
            log::warn!("prop should be pointing away from the face but the dot product with face normal is {}",
                beam.field.prop().dot(&beam.face.data().normal)
                );
        };

        assert_face_simple(&beam.face);
        assert_face_planar(verts);
        // assert_clockwise_winding(verts, prop);
        assert_e_perp_perpendicular_to_prop(e_perp, prop);
        assert_e_perp_perpendicular_to_normal(e_perp, normal);
    }

    // Step 1: Translate beam so face is centered at origin
    let center_of_mass = geom::calculate_center_of_mass(verts);
    let beam_centered = beam.transformed(&translation_to_origin(&center_of_mass))?;

    // Step 2: Rotate beam into xy plane (vertices are now centered at origin)
    let rot_to_xy = rotation_to_xy_plane(&beam_centered);
    let beam_xy = beam_centered.transformed(&rot_to_xy)?;

    // Step 3: Rotate around z-axis to put e_perp along +y (aperture system)
    // This also ensures prop lies in the xz plane since e_perp ⊥ prop
    let rot_e_perp_to_y = rotation_e_perp_to_y(&beam_xy);
    let beam_aperture = beam_xy.transformed(&rot_e_perp_to_y)?;

    // Combined rotation matrix (3x3) for rotating vectors from original to aperture frame
    let rot3_to_xy: Matrix3<f32> = rot_to_xy.fixed_view::<3, 3>(0, 0).into_owned();
    let rot3_e_perp_to_y: Matrix3<f32> = rot_e_perp_to_y.fixed_view::<3, 3>(0, 0).into_owned();
    let rot3 = rot3_e_perp_to_y * rot3_to_xy;

    // Check for amplitude sign flip (matches original diff.rs)
    let e_perp = beam.field.e_perp();
    let perp2 = rot3 * e_perp;
    let prop2 = rot3 * prop;
    let e_par2 = perp2.cross(&prop2).normalize();

    #[cfg(debug_assertions)]
    {
        let aperture_verts = &beam_aperture.face.data().exterior;
        let aperture_prop = beam_aperture.field.prop();
        let aperture_e_perp = beam_aperture.field.e_perp();
        let aperture_normal = beam_aperture.face.data().normal;

        assert_face_simple(&beam_aperture.face);
        assert_face_planar(aperture_verts);
        // assert_clockwise_winding(aperture_verts, aperture_prop);
        assert_e_perp_perpendicular_to_prop(aperture_e_perp, aperture_prop);
        assert_e_perp_perpendicular_to_normal(aperture_e_perp, aperture_normal);
    }

    // Get wavenumber from beam and scale amplitude matrix
    let wavenumber = beam_aperture.wavenumber();
    let mut ampl = beam_aperture.field.ampl();
    ampl *= Complex::new(wavenumber, 0.0);

    // Apply amplitude sign flip if needed (matches original diff.rs)
    if e_par2.z > COLINEAR_THRESHOLD {
        ampl = -ampl;
    }

    // Get aperture vertices and prop in aperture system
    let aperture_verts = &beam_aperture.face.data().exterior;
    let prop_aperture = beam_aperture.field.prop();

    // Pre-compute edge data for Fraunhofer integral
    let edge_data = EdgeData::from_vertices(aperture_verts);

    // Calculate field of view cosine for filtering
    let cos_fov = calculate_fov_cosine(aperture_verts, wavenumber, fov_factor);

    // Incident wave vector in aperture frame
    let kinc = prop_aperture * wavenumber;

    // Pre-calculate constant for Fraunhofer integral
    let inv_denom = Complex::new(wavenumber / (2.0 * PI), 0.0);

    // Output amplitude for each bin
    let mut ampl_cs = vec![Ampl::zeros(); bins.len()];

    // Phase offset vector: uses ORIGINAL center of mass (before transforms)
    // This is the displacement from aperture center to far-field reference
    let r_offset = -center_of_mass.coords;

    // Main loop over scattering bins
    for (index, bin) in bins.iter().enumerate() {
        // Calculate observation direction in the original (lab) frame
        // Uses inverted z-axis convention: theta=0 is forward (along -z)
        let k_obs = bin.unit_vector();

        // Rotate observation direction into aperture frame
        let k = rot3 * k_obs;

        // Phase calculation: path difference times wavenumber
        // Uses k_obs (original frame) with r_offset (original center of mass)
        let path_difference = k_obs.dot(&r_offset);
        let bvsk = path_difference * wavenumber;

        // Field of view filtering: skip bins outside the valid scattering cone
        if fov_factor.is_some() && k.dot(&prop_aperture) < cos_fov {
            continue;
        }

        // Get mutable reference to output amplitude for this bin
        let ampl_far_field = &mut ampl_cs[index];

        // Compute Karczewski polarisation matrix and scattering plane normal
        let (karczewski_matrix, karczewski_e_perp) = karczewski(&prop_aperture, &k);

        // Precompute sin/cos phi for rotation matrices
        let (sin_phi, cos_phi) = bin.phi.center.to_radians().sin_cos();

        // rot4: rotation from Karczewski scattering plane to aperture system scattering plane
        // hc is the vector perpendicular to the scattering plane, rotated into aperture system
        let scattering_e_perp = Vector3::new(-sin_phi, cos_phi, 0.0);
        let hc = rot3 * scattering_e_perp;
        let rot4 = crate::field::Field::rotation_matrix(karczewski_e_perp, hc, k);

        // prerotation: rotation of initial incidence reference frame
        // Uses the incident beam's e_perp and prop to define the reference frame
        let prerotation =
            crate::field::Field::rotation_matrix(incident.e_perp, scattering_e_perp, incident.prop)
                .transpose();

        // Compute amplitude: rot4 * karczewski * ampl * prerotation
        let ampl_temp = rot4.map(Complex::from)
            * karczewski_matrix.map(Complex::from)
            * ampl
            * prerotation.map(Complex::from);

        *ampl_far_field = ampl_temp;

        // Calculate Fraunhofer factor for this direction
        let mut fraunhofer_sum = Complex::new(0.0, 0.0);

        let (kxx, kyy) = calculate_kxx_kyy(
            &kinc
                .fixed_rows::<2>(0)
                .into_owned()
                .as_slice()
                .try_into()
                .unwrap(),
            &k,
            wavenumber,
        );

        // Loop over aperture edges
        let nv = edge_data.x.len();
        for j in 0..nv {
            let xj = edge_data.x[j];
            let yj = edge_data.y[j];
            let dx = edge_data.dx[j];
            let dy = edge_data.dy[j];
            let mj = edge_data.m_adj[j];
            let nj = edge_data.n_adj[j];

            let (delta, delta1, delta2) = calculate_deltas(kxx, kyy, xj, yj, mj, nj);
            let (omega1, omega2) = calculate_omegas(dx, dy, delta1, delta2);
            let (alpha, beta) = calculate_alpha_beta(delta1, delta2, kxx, kyy);

            // Skip invalid cases
            if alpha.is_infinite() || beta.is_infinite() || alpha.is_nan() || beta.is_nan() {
                continue;
            }

            let summand = calculate_summand(bvsk, delta, omega1, omega2, alpha, beta, inv_denom);

            // Final check
            if summand.is_nan() {
                continue;
            }

            fraunhofer_sum += summand;
        }

        *ampl_far_field *= fraunhofer_sum;
    }

    // Collect non-zero results
    let results: Vec<(usize, Ampl)> = ampl_cs
        .into_iter()
        .enumerate()
        .filter(|(_, a)| a.iter().any(|c| c.norm() > 0.0))
        .collect();

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::beam::{Beam, BeamVariant, DefaultBeamVariant};
    use crate::field::Field;
    use crate::geom::Face;
    use crate::settings::{default_e_perp, default_prop};
    use nalgebra::{Complex, Point3, Vector3};

    /// Create a simple triangular beam in the aperture system for testing.
    /// - Face: triangle in xy plane, centered at origin
    /// - Prop: along +z axis
    /// - e_perp: along +y axis
    fn make_test_beam_in_aperture_system() -> Beam {
        // Triangle in xy plane, clockwise when viewed from +z
        let verts = vec![
            Point3::new(-0.5, 0.866, 0.0),
            Point3::new(-0.5, -0.866, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];

        let face = Face::new_simple(verts, None, None).unwrap();

        // Prop in xz plane (not along +z to avoid collinear issues), e_perp along +y
        let prop = Vector3::new(0.5, 0.0, 0.866).normalize(); // ~30 degrees from +z
        let e_perp = prop.cross(&face.data().normal).normalize();
        let field = Field::new_identity(e_perp, prop).unwrap();

        Beam::new(
            face,
            Complex::new(1.0, 0.0), // refr_index
            0,                      // rec_count
            0,                      // tir_count
            field,
            BeamVariant::Default(DefaultBeamVariant::Refr),
            0.532, // wavelength
        )
    }

    #[test]
    fn test_null_beam_already_in_aperture_system() {
        let beam = make_test_beam_in_aperture_system();
        let bins = vec![]; // empty bins for now

        let result = n2f_aperture_diffraction(&beam, &bins, &IncidentBeam::default(), None);
        assert!(result.is_ok(), "Null test failed: {:?}", result.err());
    }

    /// Create a beam with normal incidence (prop along -z, matching default IncidentBeam).
    fn make_test_beam_normal_incidence() -> Beam {
        // Triangle in xy plane, clockwise when viewed from -z (i.e., from below)
        let verts = vec![
            Point3::new(-0.5, -0.866, 0.0),
            Point3::new(-0.5, 0.866, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];

        let face = Face::new_simple(verts, None, None).unwrap();

        // Prop along -z (forward scattering direction), e_perp along +y
        let prop = default_prop();
        let e_perp = default_e_perp();
        let field = Field::new_identity(e_perp, prop).unwrap();

        Beam::new(
            face,
            Complex::new(1.0, 0.0), // refr_index
            0,                      // rec_count
            0,                      // tir_count
            field,
            BeamVariant::Default(DefaultBeamVariant::Refr),
            0.532, // wavelength
        )
    }

    #[test]
    fn test_normal_incidence_beam() {
        let beam = make_test_beam_normal_incidence();
        let bins = vec![]; // empty bins for now

        let result = n2f_aperture_diffraction(&beam, &bins, &IncidentBeam::default(), None);
        assert!(
            result.is_ok(),
            "Normal incidence test failed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normal_incidence_forward_scattering_is_imaginary() {
        use crate::bins::{AngleBin, SolidAngleBin};

        let beam = make_test_beam_normal_incidence();

        // Forward scattering bin: use small but non-zero theta to avoid exact forward degeneracy
        let forward_bin = SolidAngleBin::new(
            AngleBin::from_center_width(0.001, 0.002),
            AngleBin::from_center_width(0.0, 1.0),
        );
        let bins = vec![forward_bin];

        let result = n2f_aperture_diffraction(&beam, &bins, &IncidentBeam::default(), None);
        assert!(
            result.is_ok(),
            "Forward scattering test failed: {:?}",
            result.err()
        );

        let results = result.unwrap();
        println!("Number of results: {}", results.len());
        for (idx, a) in &results {
            println!("Result {}: {:?}", idx, a);
        }
        assert_eq!(
            results.len(),
            1,
            "Expected exactly one result for forward bin"
        );

        let (index, ampl) = &results[0];
        assert_eq!(*index, 0, "Expected result for bin index 0");

        // Check that amplitude matrix is mostly imaginary
        // For each element, |imag| should be >> |real|
        for (i, elem) in ampl.iter().enumerate() {
            let ratio = if elem.im.abs() > 1e-10 {
                elem.re.abs() / elem.im.abs()
            } else {
                0.0 // If imaginary is ~0, ratio is 0
            };
            println!(
                "ampl[{}]: re = {:.6}, im = {:.6}, |re/im| = {:.6}",
                i, elem.re, elem.im, ratio
            );
            assert!(
                ratio < 0.1,
                "Element {} is not mostly imaginary: re = {}, im = {}, ratio = {}",
                i,
                elem.re,
                elem.im,
                ratio
            );
        }
    }

    #[test]
    fn test_cross_polarised_input_gives_cross_polarised_output() {
        use crate::bins::{AngleBin, SolidAngleBin};

        // Triangle in xy plane, clockwise when viewed from -z (along prop direction)
        let verts = vec![
            Point3::new(-0.5, -0.866, 0.0),
            Point3::new(-0.5, 0.866, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];

        let face = crate::geom::Face::new_simple(verts, None, None).unwrap();

        // Normal incidence
        let prop = default_prop();
        let e_perp = default_e_perp();

        // Cross-polarised amplitude: only off-diagonal elements
        // [[0, 1], [1, 0]] - swaps perp and par components
        let cross_pol_ampl = nalgebra::Matrix2::new(
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
        );

        let field = crate::field::Field::new(e_perp, prop, cross_pol_ampl, 0.0).unwrap();

        let beam = Beam::new(
            face,
            Complex::new(1.0, 0.0),
            0,
            0,
            field,
            BeamVariant::Default(DefaultBeamVariant::Refr),
            0.532,
        );

        // Near-forward bin (avoiding exact theta=0 numerical issues)
        let bins = vec![SolidAngleBin::new(
            AngleBin::from_center_width(0.01, 0.02),
            AngleBin::from_center_width(45.0, 1.0), // avoid phi=0 instability
        )];

        let results =
            n2f_aperture_diffraction(&beam, &bins, &IncidentBeam::default(), None).unwrap();

        assert_eq!(results.len(), 1, "Expected one result");
        let (_, ampl) = &results[0];

        println!("Cross-polarised input, forward scattering output:");
        println!("  ampl[0,0] = {:?}", ampl[(0, 0)]);
        println!("  ampl[0,1] = {:?}", ampl[(0, 1)]);
        println!("  ampl[1,0] = {:?}", ampl[(1, 0)]);
        println!("  ampl[1,1] = {:?}", ampl[(1, 1)]);

        // For cross-polarised input, output should also be cross-polarised
        // i.e., diagonal elements should be small, off-diagonal should dominate
        let diag_norm = ampl[(0, 0)].norm() + ampl[(1, 1)].norm();
        let off_diag_norm = ampl[(0, 1)].norm() + ampl[(1, 0)].norm();

        println!("  diag_norm = {:.6e}", diag_norm);
        println!("  off_diag_norm = {:.6e}", off_diag_norm);

        assert!(
            off_diag_norm > diag_norm * 10.0,
            "Cross-polarised input should give cross-polarised output: off_diag={:.3e}, diag={:.3e}",
            off_diag_norm,
            diag_norm
        );
    }

    #[test]
    fn test_max_amplitude_follows_prop_direction() {
        use crate::bins::{AngleBin, SolidAngleBin};

        // Test various incidence angles
        let test_angles = [
            (0.5_f32, 0.0_f32), // near-normal incidence
            (10.0, 0.0),        // 10 degrees in xz plane
            (20.0, 0.0),        // 20 degrees in xz plane
            (10.0, 90.0),       // 10 degrees in yz plane
            (15.0, 45.0),       // 15 degrees at phi=45
        ];

        for (inc_theta, inc_phi) in test_angles {
            println!("\n=== Incidence: theta={}, phi={} ===", inc_theta, inc_phi);

            // Create prop direction from incidence angles
            // Using same convention as bins: theta from -z, phi from +x
            let theta_rad = inc_theta.to_radians();
            let phi_rad = inc_phi.to_radians();
            let prop = Vector3::new(
                theta_rad.sin() * phi_rad.cos(),
                theta_rad.sin() * phi_rad.sin(),
                -theta_rad.cos(),
            );
            println!("prop = {:?}", prop);

            // e_perp must be perpendicular to both prop and face normal
            let normal = Vector3::new(0.0, 0.0, 1.0);
            let e_perp = prop.cross(&normal).normalize();

            // Triangle in xy plane, clockwise when viewed from -z (along prop direction)
            let verts = vec![
                Point3::new(-0.5, -0.866, 0.0),
                Point3::new(-0.5, 0.866, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ];

            let face = crate::geom::Face::new_simple(verts, None, None).unwrap();
            let field = crate::field::Field::new_identity(e_perp, prop).unwrap();

            let beam = Beam::new(
                face,
                Complex::new(1.0, 0.0),
                0,
                0,
                field,
                BeamVariant::Default(DefaultBeamVariant::Refr),
                0.532,
            );

            // Create bins covering a range of theta and phi
            let mut bins = Vec::new();
            for theta_idx in 0..30 {
                let theta = 0.5 + theta_idx as f32 * 1.0; // 0.5 to 29.5 degrees
                for phi_idx in 0..8 {
                    let phi = phi_idx as f32 * 45.0; // 0, 45, 90, ..., 315
                    bins.push(SolidAngleBin::new(
                        AngleBin::from_center_width(theta, 1.0),
                        AngleBin::from_center_width(phi, 45.0),
                    ));
                }
            }

            let results =
                n2f_aperture_diffraction(&beam, &bins, &IncidentBeam::default(), None).unwrap();

            // Find bin with maximum amplitude
            let mut max_ampl = 0.0_f32;
            let mut max_bin_idx = 0;
            for (idx, ampl) in &results {
                let norm = ampl.norm();
                if norm > max_ampl {
                    max_ampl = norm;
                    max_bin_idx = *idx;
                }
            }

            let max_bin = &bins[max_bin_idx];
            println!(
                "Max amplitude at theta={:.1}, phi={:.1} (norm={:.3e})",
                max_bin.theta.center, max_bin.phi.center, max_ampl
            );

            // Check that the maximum is near the propagation direction
            // Allow some tolerance since bins are discrete
            let theta_diff = (max_bin.theta.center - inc_theta).abs();
            let phi_diff = if inc_theta < 1.0 {
                0.0 // phi is undefined at theta=0
            } else {
                let d = (max_bin.phi.center - inc_phi).abs();
                d.min(360.0 - d) // handle wrap-around
            };

            println!(
                "Difference from expected: theta_diff={:.1}, phi_diff={:.1}",
                theta_diff, phi_diff
            );

            // The maximum should be within a few bins of the propagation direction
            assert!(
                theta_diff < 5.0,
                "Max amplitude theta ({:.1}) too far from incidence theta ({:.1})",
                max_bin.theta.center,
                inc_theta
            );
            if inc_theta >= 1.0 {
                assert!(
                    phi_diff < 50.0,
                    "Max amplitude phi ({:.1}) too far from incidence phi ({:.1})",
                    max_bin.phi.center,
                    inc_phi
                );
            }
        }
    }

    #[test]
    fn test_prerotation_preserves_physical_field() {
        use nalgebra::Matrix2;

        // Test that prerotation correctly rotates amplitude into scattering plane
        // while preserving the physical E field.
        //
        // For unpolarised light: E_perp = E_par = 1
        // Physical E = E_perp * e_perp + E_par * e_par
        //
        // After rotation to scattering plane, the physical E should be unchanged.

        let phi = 45.0_f32;
        let (sin_phi, cos_phi) = phi.to_radians().sin_cos();

        // Scattering plane basis vectors (at phi=45)
        let e_perp_scat = Vector3::new(-sin_phi, cos_phi, 0.0); // perpendicular to scattering plane
        let prop = Vector3::new(0.0, 0.0, -1.0);
        let e_par_scat = e_perp_scat.cross(&prop).normalize(); // parallel to scattering plane

        println!(
            "Scattering plane (phi={}): e_perp={:?}, e_par={:?}",
            phi, e_perp_scat, e_par_scat
        );

        // Test different incident e_perp directions
        let test_cases = [
            ("e_perp along +y", Vector3::new(0.0, 1.0, 0.0).normalize()),
            ("e_perp along +x", Vector3::new(1.0, 0.0, 0.0).normalize()),
            (
                "e_perp at 30deg",
                Vector3::new(
                    30.0_f32.to_radians().sin(),
                    30.0_f32.to_radians().cos(),
                    0.0,
                ),
            ),
            ("e_perp at -45deg", Vector3::new(1.0, -1.0, 0.0).normalize()),
        ];

        for (name, e_perp_inc) in test_cases {
            println!("\n=== {} ===", name);
            let e_par_inc = e_perp_inc.cross(&prop).normalize();
            println!("Incident: e_perp={:?}, e_par={:?}", e_perp_inc, e_par_inc);

            // Unpolarised input: E_perp = E_par = 1
            let e_perp_component = 1.0_f32;
            let e_par_component = 1.0_f32;

            // Physical E field in incident frame
            let e_physical = e_perp_inc * e_perp_component + e_par_inc * e_par_component;
            println!("Physical E field: {:?}", e_physical);

            // Compute prerotation matrix
            let prerotation =
                crate::field::Field::rotation_matrix(e_perp_inc, e_perp_scat, prop).transpose();
            println!("Prerotation matrix: {:?}", prerotation);

            // Apply prerotation to amplitude (identity for unpolarised)
            let ampl_inc = Matrix2::<f32>::identity();
            let ampl_scat = prerotation * ampl_inc;
            println!("Amplitude after prerotation: {:?}", ampl_scat);

            // The rotated amplitude should give the same physical E field
            // when we use the scattering plane basis vectors
            let e_perp_out =
                ampl_scat[(0, 0)] * e_perp_component + ampl_scat[(0, 1)] * e_par_component;
            let e_par_out =
                ampl_scat[(1, 0)] * e_perp_component + ampl_scat[(1, 1)] * e_par_component;
            let e_physical_out = e_perp_scat * e_perp_out + e_par_scat * e_par_out;
            println!("Physical E field after rotation: {:?}", e_physical_out);

            // Check that physical E field is preserved
            let diff = (e_physical - e_physical_out).norm();
            println!("Difference: {:.6e}", diff);

            assert!(
                diff < 1e-4,
                "{}: Physical E field not preserved! diff={:.6e}",
                name,
                diff
            );
        }
    }

    #[test]
    fn test_debug_phi_zero_instability() {
        use crate::bins::{AngleBin, SolidAngleBin};
        use nalgebra::Matrix2;

        // For normal incidence beam in aperture system:
        // prop_aperture = (0, 0, 1) approximately (with small perturbation)
        // We'll use a simplified prop for debugging
        let prop_aperture = Vector3::new(0.0, 0.0, 1.0_f32);

        // Compare phi=0 vs phi=45
        for phi in [0.0_f32, 45.0] {
            println!("\n=== phi = {} ===", phi);

            let bin = SolidAngleBin::new(
                AngleBin::from_center_width(0.01, 0.02),
                AngleBin::from_center_width(phi, 1.0),
            );

            let k_obs = bin.unit_vector();
            println!("k_obs = {:?}", k_obs);

            // For aperture system, rot3 is identity (beam already in aperture system)
            let rot3 = Matrix3::<f32>::identity();
            let k = rot3 * k_obs;
            println!("k (in aperture) = {:?}", k);

            // Karczewski matrix
            let (karczewski_matrix, m) = karczewski(&prop_aperture, &k);
            println!("karczewski_matrix = {:?}", karczewski_matrix);
            println!("m = {:?}", m);

            // Simulate what happens in the main function
            let (sin_phi, cos_phi) = phi.to_radians().sin_cos();
            println!("sin_phi = {}, cos_phi = {}", sin_phi, cos_phi);

            // hc vector (rot3 is identity here)
            let hc = rot3 * Vector3::new(sin_phi, -cos_phi, 0.0);
            println!("hc = {:?}", hc);

            // evo2
            let evo2 = k.cross(&m);
            println!("evo2 = {:?}", evo2);

            // rot4
            let rot4 = Matrix2::new(hc.dot(&m), -hc.dot(&evo2), hc.dot(&evo2), hc.dot(&m));
            println!("rot4 = {:?}", rot4);

            // prerotation
            let e_perp_in = Vector3::<f32>::y();
            let e_perp_out = Vector3::new(-sin_phi, cos_phi, 0.0);
            let prop = -Vector3::<f32>::z();
            let prerotation =
                crate::field::Field::rotation_matrix(e_perp_in, e_perp_out, prop).transpose();
            println!("prerotation = {:?}", prerotation);

            // Combined rotation effect on identity amplitude
            let ampl = Matrix2::<f32>::identity();
            let combined = rot4 * karczewski_matrix * ampl * prerotation;
            println!(
                "combined (rot4 * karczewski * I * prerotation) = {:?}",
                combined
            );
        }
    }
}
