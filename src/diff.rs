use anyhow::{anyhow, Result};
use clap::ValueEnum;
use nalgebra::{Complex, Matrix2, Matrix3, Point3, Vector3};
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

use crate::beam::Beam;
use crate::bins::{get_n_linear_search, get_n_simple, BinningScheme, Scheme, SolidAngleBin};
use crate::field::{Ampl, Field};
use crate::{geom, settings};

/// Tolerance for planarity check - maximum allowed deviation from the aperture plane
const PLANARITY_TOLERANCE: f32 = 1e-4;

/// Check if a set of vertices are coplanar.
/// Returns Ok(()) if planar, Err with details if not.
pub fn check_planarity(verts: &[Point3<f32>]) -> Result<()> {
    if verts.len() < 3 {
        return Err(anyhow!(
            "Aperture must have at least 3 vertices, got {}",
            verts.len()
        ));
    }

    if verts.len() == 3 {
        // Three points are always coplanar
        return Ok(());
    }

    // Compute normal from first three vertices
    let v0 = verts[1] - verts[0];
    let v1 = verts[2] - verts[0];
    let normal = v0.cross(&v1);

    if normal.norm() < settings::VEC_LENGTH_THRESHOLD {
        return Err(anyhow!(
            "First three vertices are collinear, cannot determine plane"
        ));
    }

    let normal = normal.normalize();

    // Check all remaining vertices lie in the plane
    for (i, vert) in verts.iter().enumerate().skip(3) {
        let v = vert - verts[0];
        let distance = v.dot(&normal).abs();

        if distance > PLANARITY_TOLERANCE {
            return Err(anyhow!(
                "Aperture is non-planar: vertex {} deviates from plane by {:.6} (tolerance: {:.6})",
                i,
                distance,
                PLANARITY_TOLERANCE
            ));
        }
    }

    Ok(())
}

/// Enum representing different mapping methods from near to far field.
#[pyclass]
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, ValueEnum, Copy)]
pub enum Mapping {
    GeometricOptics,
    ApertureDiffraction,
}

#[pymethods]
impl Mapping {
    #[new]
    pub fn py_new(str: &str) -> PyResult<Self> {
        match str.to_lowercase().as_str() {
            "go" => Ok(Mapping::GeometricOptics),
            "ad" => Ok(Mapping::ApertureDiffraction),
            _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "'{}' is not a valid Mapping method. Valid options are: 'go' (Geometric Optics), 'ad' (Aperture Diffraction)",
                str
            ))),
        }
    }
}

/// Map a beam to the far-field using geometric optics. Assumes delta theta and delta phi are provided if the binning scheme is Simple. Returns a single-element vector containing the bin index and amplitude matrix.
pub fn n2f_go(binning: &BinningScheme, bins: &[SolidAngleBin], beam: &Beam) -> Vec<(usize, Ampl)> {
    // Use the precomputed theta and phi spacings if using Simple binning
    let (delta_theta, delta_phi) = match binning.scheme {
        Scheme::Simple {
            num_theta: _,
            num_phi: _,
            delta_theta,
            delta_phi,
        } => (Some(delta_theta), Some(delta_phi)),
        Scheme::Interval { .. } => (None, None),
        Scheme::Custom { .. } => (None, None),
    };
    // Get beam scattering angles
    let (theta, phi) = beam.get_scattering_angles();

    // Map scattering angles to corresponding bin
    let Some(n) = (match &binning.scheme {
        Scheme::Simple {
            num_theta, num_phi, ..
        } => {
            // Safe to unwrap because we know the scheme is Simple
            get_n_simple(
                *num_theta,
                *num_phi,
                delta_theta.unwrap(),
                delta_phi.unwrap(),
                theta,
                phi,
            )
        }
        Scheme::Interval { .. } => get_n_linear_search(bins, theta, phi),
        Scheme::Custom { .. } => get_n_linear_search(bins, theta, phi),
    }) else {
        return vec![];
    };

    // Get amplitude rotation matrices
    let (rotation, prerotation) = get_mapping_rotations(beam, phi);

    // Calculate the reference phase correction
    let phase_correction = get_reference_phase(beam);

    // Compute solid angle
    let solid_angle = &bins[n].solid_angle();

    // Compute scaling factor (sqrt <- amplitude, not intensity)
    let scale_factor = beam.csa().sqrt() // account for beam cross-sectional area
        / solid_angle.sqrt() // account for Jacobian: Cartesian to spherical
        * 5.34464802915 // bodge empirical factor (probably slight underestimate)
        / beam.wavelength;
    // account for scaled wavelength

    // Compute far-field amplitude matrix
    let ampl = rotation // rotation from beam plane to scattering plane
        * beam.field.ampl() // outgoing beam amplitude matrix
        * prerotation // pre-rotation of the initial incidence
        * Complex::new(scale_factor, 0.0) // amplitude scaling factor
        * phase_correction; // reference phase correction

    vec![(n, ampl)]
}

/// Mapping from near to far field using aperture diffraction theory.
///
/// # Errors
/// Returns an empty vector if the aperture is non-planar.
pub fn n2f_aperture_diffraction(
    verts: &[Point3<f32>],
    mut ampl: Ampl,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
    bins: &[SolidAngleBin],
    wavenumber: f32,
    fov_factor: Option<f32>,
) -> Vec<Matrix2<Complex<f32>>> {
    // Translate to aperture system, rotate, and transform propagation and auxiliary vectors.
    let (center_of_mass, relative_vertices, rot3, prop2) =
        match init_diff(verts, &mut ampl, prop, vk7, wavenumber) {
            Ok(result) => result,
            Err(e) => {
                eprintln!("n2f_aperture_diffraction: {}", e);
                return vec![Matrix2::<Complex<f32>>::zeros(); bins.len()];
            }
        };

    // --- Optimizations Start ---
    // Pre-calculate transformed vertices and related quantities outside the main loop
    let nv = relative_vertices.len();
    let mut v1_data = Vec::with_capacity(nv * 3);
    let mut transformed_vertices_vec = Vec::with_capacity(nv);
    for vertex in &relative_vertices {
        let transformed_vertex = rot3 * vertex;
        transformed_vertices_vec.push(transformed_vertex);
        v1_data.push(transformed_vertex.x);
        v1_data.push(transformed_vertex.y);
        v1_data.push(transformed_vertex.z);
    }

    // Get estimated field of view cosine
    let aperture_dimension = transformed_vertices_vec
        .iter()
        .map(|v| v.norm())
        .fold(0.0, f32::max);

    // fov_factor * lambda / 2 * r
    let cos_fov =
        (fov_factor.unwrap_or(1.0) * 2.0 * 2.0 * PI / (wavenumber * aperture_dimension)).cos();

    // Extract x and y coordinates from v1_data (which is stored as [x0, y0, z0, x1, y1, z1, ...])
    let x: Vec<f32> = (0..nv).map(|i| v1_data[i * 3]).collect();
    let y: Vec<f32> = (0..nv).map(|i| v1_data[i * 3 + 1]).collect();

    // Pre-calculate dx, dy, m, n, m_adj, n_adj
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

        let mj = if dx.abs() < settings::DIFF_DMIN {
            if dy.signum() == dx.signum() {
                1e6
            } else {
                -1e6
            }
        } else {
            dy / dx
        };
        m.push(mj);

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

        // Adjust dx/dy based on DIFF_DMIN *after* calculating m
        dx = if dx.abs() < settings::DIFF_DMIN {
            settings::DIFF_DMIN * dx.signum()
        } else {
            dx
        };
        dy = if dy.abs() < settings::DIFF_DMIN {
            settings::DIFF_DMIN * dy.signum()
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
    // --- Optimizations End ---

    // Define the output variables.
    let mut ampl_cs = vec![Matrix2::<Complex<f32>>::default(); bins.len()];
    let kinc = prop2 * wavenumber;

    // Pre-calculate constants dependent only on wavenumber
    // let radius = settings::RADIUS * 2.0 * PI / wavenumber;
    let inv_denom = Complex::new(wavenumber / (2.0 * PI), 0.0); // Pre-calculate 1.0 / (2*PI/wavenumber)

    // Iterate over the flattened combinations
    for (index, bin) in bins.iter().enumerate() {
        // Compute sin and cos values for current theta and phi bin centers
        let (sin_theta, cos_theta) = bin.theta_bin.center.to_radians().sin_cos();
        let (sin_phi, cos_phi) = bin.phi_bin.center.to_radians().sin_cos();

        // Calculate observation direction in original frame
        let k_obs = Vector3::new(sin_theta * cos_phi, sin_theta * sin_phi, -cos_theta);

        // Rotate observation direction to aperture frame
        let k = rot3 * k_obs;

        // Phase calculation: relative to far-field reference distance
        // The phase reference is at the aperture center, so we compute
        // the additional phase due to the center of mass displacement
        let r_offset = -center_of_mass.coords;
        let path_difference = k_obs.dot(&r_offset);
        let bvsk = path_difference * wavenumber;

        // Apply filtering based on field of view if specified
        if fov_factor.is_some() && k.dot(&prop2) < cos_fov {
            continue;
        }

        let ampl_far_field = &mut ampl_cs[index];

        let (karczewski, rot4, prerotation) = get_rotations(rot3, prop2, sin_phi, cos_phi, k);

        let exp_factor = Complex::cis(bvsk); // Use cis for complex exponential

        // Debug hook for forward scattering (theta ~ 0)
        if sin_theta.abs() < 0.01 {
            // debug prints here
        }

        // in direct forwards, ampl_temp is diagonal for external diffraction
        let ampl_temp = rot4.map(Complex::from)
            * karczewski.map(Complex::from)
            * ampl
            // * exp_factor
            * prerotation.map(Complex::from);

        if sin_theta.abs() < 0.01 {
            println!("ampl_temp (before Fraunhofer):\n{}", ampl_temp);
            println!("--- End forward scattering debug ---");
        }

        *ampl_far_field = ampl_temp;

        // Calculate fraunhofer factor for this direction
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

        for j in 0..nv {
            // Use pre-calculated values
            let mj = m_adj[j];
            let nj = n_adj[j];
            let xj = x[j];
            let yj = y[j];
            let dx = dx_vec[j];
            let dy = dy_vec[j];

            // mj, nj are already adjusted, dx/dy already handled DIFF_DMIN

            let (delta, delta1, delta2) = calculate_deltas(kxx, kyy, xj, yj, mj, nj);
            let (omega1, omega2) = calculate_omegas(dx, dy, delta1, delta2);
            let (alpha, beta) = calculate_alpha_beta(delta1, delta2, kxx, kyy);

            // Initial checks for frequent cases before calculate_summand()
            if alpha.is_infinite() || beta.is_infinite() || alpha.is_nan() || beta.is_nan() {
                continue;
            }

            let summand = calculate_summand(delta, omega1, omega2, alpha, beta, inv_denom); // Pass inv_denom

            // Final check just to be sure
            if summand.is_nan() {
                continue;
            }

            fraunhofer_sum += summand;
        }

        *ampl_far_field *= fraunhofer_sum;
    }
    ampl_cs
}

// Other functions remain unchanged
#[inline]
fn get_rotations(
    rot3: Matrix3<f32>,
    prop2: Vector3<f32>,
    sin_phi: f32,
    cos_phi: f32,
    k: Vector3<f32>,
) -> (Matrix2<f32>, Matrix2<f32>, Matrix2<f32>) {
    let (karczewski, m) = karczewski(&prop2, &k); // compute karczweski polarisation matrix

    let hc = rot3 * Vector3::new(sin_phi, -cos_phi, 0.0); // rotate the vector perpendicular to the scattering plane into the aperture system
    let evo2 = k.cross(&m); // compute the perpendicular component from the product of scattering direction and the parallel vector
    let rot4 = Matrix2::new(hc.dot(&m), -hc.dot(&evo2), hc.dot(&evo2), hc.dot(&m)); // compute the rotation matrix

    let prerotation = Field::rotation_matrix(
        Vector3::x(),
        Vector3::new(-sin_phi, cos_phi, 0.0),
        -Vector3::z(),
    )
    .transpose();
    (karczewski, rot4, prerotation)
}

/// Tolerance for prop direction check - how close to parallel with aperture is allowed
const PROP_DIRECTION_TOLERANCE: f32 = 0.1;

/// Initialize diffraction calculation by transforming aperture to standard frame.
/// Returns Result containing (center_of_mass, relative_vertices, rot3, prop2) where:
/// - rot3: combined rotation matrix that maps lab frame -> aperture frame
/// - prop2: propagation direction in aperture frame (should have +z, lie in xz plane)
///
/// # Errors
/// Returns an error if:
/// - The aperture vertices are non-planar
/// - The propagation direction is nearly parallel to the aperture (not pointing into it)
pub fn init_diff(
    verts: &[Point3<f32>],
    ampl: &mut Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
    wavenumber: f32,
) -> Result<(Point3<f32>, Vec<Vector3<f32>>, Matrix3<f32>, Vector3<f32>)> {
    // Check that the aperture is planar
    check_planarity(verts)?;

    // Compute aperture normal to validate prop direction
    let v0: Vector3<f32> = verts[1] - verts[0];
    let v1: Vector3<f32> = verts[2] - verts[0];
    let normal = v0.cross(&v1).normalize();

    // Check that prop has a significant component into or out of the aperture
    // (i.e., prop is not nearly parallel to the aperture plane)
    let prop_dot_normal = prop.dot(&normal).abs();
    if prop_dot_normal < PROP_DIRECTION_TOLERANCE {
        return Err(anyhow!(
            "Propagation direction is nearly parallel to aperture plane (|prop · normal| = {:.4} < {})",
            prop_dot_normal,
            PROP_DIRECTION_TOLERANCE
        ));
    }

    // Check that e_perp (vk7) is perpendicular to prop
    let perp_dot_prop = vk7.dot(&prop).abs();
    if perp_dot_prop > PROP_DIRECTION_TOLERANCE {
        return Err(anyhow!(
            "e_perp is not perpendicular to propagation direction (|e_perp · prop| = {:.4} > {})",
            perp_dot_prop,
            PROP_DIRECTION_TOLERANCE
        ));
    }

    // Check that e_perp (vk7) is perpendicular to face normal
    let perp_dot_normal = vk7.dot(&normal).abs();
    if perp_dot_normal > PROP_DIRECTION_TOLERANCE {
        return Err(anyhow!(
            "e_perp is not perpendicular to face normal (|e_perp · normal| = {:.4} > {})",
            perp_dot_normal,
            PROP_DIRECTION_TOLERANCE
        ));
    }

    let prop = (prop
        + Vector3::new(
            settings::PROP_PERTURBATION,
            settings::PROP_PERTURBATION,
            settings::PROP_PERTURBATION,
        ))
    .normalize();

    *ampl *= Complex::new(wavenumber, 0.0);

    let center_of_mass = geom::calculate_center_of_mass(verts);

    let relative_vertices = geom::translate(verts, &center_of_mass);

    // Use the new rotation function that accounts for prop direction
    // Step 1: Rotate aperture into xy plane (normal along ±z)
    let rot1 = get_rotation_matrix_for_aperture(&relative_vertices, prop);
    let prop1 = rot1 * prop;
    let perp1 = rot1 * vk7;

    // Step 2: Rotate around z-axis to put e_perp along ±y axis
    // Since e_perp is perpendicular to both prop and face normal, and face normal is now along z,
    // e_perp lies in the xy plane. Rotating around z aligns it with y.
    let rot2 = calculate_rotation_to_align_perp_with_y(perp1);
    let rot_combined = rot2 * rot1;

    let prop2 = rot2 * prop1;
    let perp2 = rot2 * perp1;

    println!(
        "init_diff vk7 (lab e_perp): ({:.3},{:.3},{:.3})",
        vk7.x, vk7.y, vk7.z
    );
    let e_par2 = perp2.cross(&prop2).normalize();
    println!(
        "init_diff prop2 (aperture): ({:.3},{:.3},{:.3})",
        prop2.x, prop2.y, prop2.z
    );
    println!(
        "init_diff e_par2 (aperture): ({:.3},{:.3},{:.3})",
        e_par2.x, e_par2.y, e_par2.z
    );
    println!(
        "init_diff perp2 (aperture e_perp): ({:.3},{:.3},{:.3})",
        perp2.x, perp2.y, perp2.z
    );

    if e_par2.z > settings::COLINEAR_THRESHOLD {
        *ampl = -*ampl;
    }
    Ok((center_of_mass, relative_vertices, rot_combined, prop2))
}

/// Compute a rotation matrix that transforms the aperture plane into the xy plane,
/// ensuring the propagation direction ends up with a +z component.
///
/// The approach:
/// 1. Compute the aperture normal from vertices using cross product
/// 2. Determine target based on prop direction (prop should end up with +z)
/// 3. Compute a rotation that aligns the normal appropriately
///
/// Uses Rodrigues' rotation formula to find the rotation matrix.
pub fn get_rotation_matrix_for_aperture(
    verts: &Vec<Vector3<f32>>,
    prop: Vector3<f32>,
) -> Matrix3<f32> {
    // Compute aperture normal from first three vertices
    let v0 = verts[1] - verts[0];
    let v1 = verts[2] - verts[0];
    let normal = v0.cross(&v1).normalize();

    // We want prop to end up pointing towards +z after rotation.
    // Convention: both normal and prop point outward from the face.
    // So prop · normal > 0 means they point in the same hemisphere.
    // To make prop point towards +z, we align normal with +z when prop · normal > 0.
    let target = if prop.dot(&normal) > 0.0 {
        // prop and normal point in same direction (both outward) - align normal with +z
        Vector3::new(0.0, 0.0, 1.0)
    } else {
        // prop points opposite to normal - align normal with -z
        Vector3::new(0.0, 0.0, -1.0)
    };

    // If normal is already aligned with target, handle specially
    let dot = normal.dot(&target);

    if dot > 1.0 - settings::COLINEAR_THRESHOLD {
        // Normal already aligned with target, no rotation needed
        return Matrix3::identity();
    }

    if dot < -1.0 + settings::COLINEAR_THRESHOLD {
        // Normal points opposite to target, rotate 180° around x-axis
        return Matrix3::new(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0);
    }

    // General case: use Rodrigues' rotation formula
    // Rotation axis is the cross product of normal and target
    let axis = normal.cross(&target).normalize();

    // Rotation angle from dot product
    let cos_angle = dot;
    let sin_angle = (1.0 - cos_angle * cos_angle).sqrt();

    // Rodrigues' formula: R = I + sin(θ)K + (1-cos(θ))K²
    // where K is the skew-symmetric cross-product matrix of the axis
    let k = Matrix3::new(
        0.0, -axis.z, axis.y, axis.z, 0.0, -axis.x, -axis.y, axis.x, 0.0,
    );

    let k_squared = k * k;

    Matrix3::identity() + k * sin_angle + k_squared * (1.0 - cos_angle)
}

/// Legacy wrapper - use get_rotation_matrix_for_aperture instead
pub fn get_rotation_matrix2(verts: &Vec<Vector3<f32>>) -> Matrix3<f32> {
    // Default behavior: assume prop points in -z direction
    get_rotation_matrix_for_aperture(verts, Vector3::new(0.0, 0.0, -1.0))
}

#[inline]
pub fn karczewski(prop2: &Vector3<f32>, bvk: &Vector3<f32>) -> (Matrix2<f32>, Vector3<f32>) {
    let big_kx = prop2.x;
    let big_ky = prop2.y;
    let big_kz = prop2.z;

    let one_minus_k2y2 = (1.0 - bvk.y.powi(2)).max(0.0);
    let sqrt_1_minus_k2y2 = one_minus_k2y2.sqrt();
    let sqrt_1_minus_k2y2 = if sqrt_1_minus_k2y2.abs() < settings::DIFF_EPSILON {
        settings::DIFF_EPSILON
    } else {
        sqrt_1_minus_k2y2
    };

    let m = Vector3::new(
        -bvk.x * bvk.y / sqrt_1_minus_k2y2,
        sqrt_1_minus_k2y2,
        -bvk.y * bvk.z / sqrt_1_minus_k2y2,
    );

    let frac = (one_minus_k2y2 / (1.0 - big_ky.powi(2))).sqrt();
    let frac = if frac.abs() < settings::DIFF_EPSILON {
        settings::DIFF_EPSILON
    } else {
        frac
    };

    let a1m = -big_kz * frac;
    let b2m = -bvk.z / frac;
    let a1e = b2m;
    let b2e = a1m;
    let b1m = -bvk.x * bvk.y / frac + big_kx * big_ky * frac;
    let a2e = -b1m;

    let a1em = 0.5 * (a1m + a1e);
    let a2em = 0.5 * a2e;
    let b1em = 0.5 * b1m;
    let b2em = 0.5 * (b2m + b2e);

    let diff_ampl = Matrix2::new(a1em, b1em, a2em, b2em);

    (diff_ampl, m)
}

#[inline]
pub fn calculate_rotation_matrix(prop1: Vector3<f32>) -> Matrix3<f32> {
    let angle = -prop1.y.atan2(prop1.x);
    let (sin_angle, cos_angle) = angle.sin_cos();

    Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    )
}

/// Calculate a rotation matrix around the z-axis to align e_perp with +y axis.
///
/// Precondition: e_perp must be perpendicular to both prop and face normal.
/// After rot1 puts the face normal along z, e_perp will lie in the xy plane (z ≈ 0).
/// This function rotates around z to put e_perp along +y.
#[inline]
pub fn calculate_rotation_to_align_perp_with_y(perp: Vector3<f32>) -> Matrix3<f32> {
    // perp is in xy plane (z ≈ 0), we want to rotate it to align with +y
    // Current angle of perp in xy plane: atan2(y, x)
    // To put perp along +y, we rotate by (π/2 - current_angle)
    let current_angle = perp.y.atan2(perp.x);
    let target_angle = std::f32::consts::FRAC_PI_2; // π/2 = +y direction
    let rotation_angle = target_angle - current_angle;

    let (sin_angle, cos_angle) = rotation_angle.sin_cos();

    Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    )
}

#[inline]
pub fn adjust_mj_nj(mj: f32, nj: f32) -> (f32, f32) {
    if mj.abs() > 1e6 || nj.abs() < 1e-6 {
        (1e6, 1e-6)
    } else if nj.abs() > 1e6 || mj.abs() < 1e-6 {
        (1e-6, 1e6)
    } else {
        (mj, nj)
    }
}

#[inline]
pub fn calculate_kxx_kyy(kinc: &[f32; 2], k: &Vector3<f32>, wavenumber: f32) -> (f32, f32) {
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
pub fn calculate_deltas(kxx: f32, kyy: f32, xj: f32, yj: f32, mj: f32, nj: f32) -> (f32, f32, f32) {
    let delta = kxx * xj + kyy * yj;
    let delta1 = kyy * mj + kxx;
    let delta2 = kxx * nj + kyy;
    (delta, delta1, delta2)
}

#[inline]
pub fn calculate_omegas(dx: f32, dy: f32, delta1: f32, delta2: f32) -> (f32, f32) {
    let omega1 = dx * delta1;
    let omega2 = dy * delta2;
    (omega1, omega2)
}

#[inline]
pub fn calculate_alpha_beta(delta1: f32, delta2: f32, kxx: f32, kyy: f32) -> (f32, f32) {
    let alpha = 1.0 / (2.0 * kyy * delta1);
    let beta = 1.0 / (2.0 * kxx * delta2);
    (alpha, beta)
}

#[inline]
pub fn calculate_summand(
    // bvsk: f32,
    delta: f32,
    omega1: f32,
    omega2: f32,
    alpha: f32,
    beta: f32,
    inv_denom: Complex<f32>, // Accept pre-calculated inverse denominator
) -> Complex<f32> {
    let (sin_delta, cos_delta) = delta.sin_cos();
    let (sin_delta_omega1, cos_delta_omega1) = (delta + omega1).sin_cos();
    let (sin_delta_omega2, cos_delta_omega2) = (delta + omega2).sin_cos();

    let sumim = alpha * (cos_delta - cos_delta_omega1) - beta * (cos_delta - cos_delta_omega2);
    let sumre = -alpha * (sin_delta - sin_delta_omega1) + beta * (sin_delta - sin_delta_omega2);

    // let exp_factor = Complex::cis(bvsk); // Use cis for complex exponential

    // exp_factor * Complex::new(sumre, sumim) * inv_denom // Multiply by inverse denominator
    Complex::new(sumre, sumim) * inv_denom // Multiply by inverse denominator
}

/// Returns the reference phase correction for accounting for how far the beam must travel to reach a point on the scattering sphere in the far-field.
fn get_reference_phase(beam: &Beam) -> Complex<f32> {
    let exp_factor = {
        let position = beam.face.data().midpoint.coords;
        let correction = -beam.field.prop().dot(&position) * beam.wavenumber();
        Complex::cis(correction)
    };
    exp_factor
}

/// Returns the prerotation and rotation matrices for rotating a beam into the scattering plane based on phi angle.
fn get_mapping_rotations(beam: &Beam, phi: f32) -> (Matrix2<Complex<f32>>, Matrix2<Complex<f32>>) {
    let (sin_phi, cos_phi) = phi.to_radians().sin_cos();
    let hc = Vector3::new(sin_phi, -cos_phi, 0.0);
    let rotation = Field::rotation_matrix(beam.field.e_perp(), hc, beam.field.prop());

    let prerotation = Field::rotation_matrix(
        Vector3::x(),
        Vector3::new(-sin_phi, cos_phi, 0.0),
        -Vector3::z(),
    )
    .transpose();
    (rotation.map(Complex::from), prerotation.map(Complex::from))
}
