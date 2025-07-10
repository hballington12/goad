use nalgebra::{Complex, Matrix2, Matrix3, Point3, Vector3};
use ndarray::Array2;
use std::f32::consts::PI;

use crate::field::Field;
use crate::{geom, settings};

// Diffraction. Face must be convex.
pub fn diffraction(
    verts: &[Point3<f32>],
    mut ampl: Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
    theta_phi_combinations: &[(f32, f32)],
    wavenumber: f32,
    fov_factor: Option<f32>,
) -> Vec<Matrix2<Complex<f32>>> {
    // Translate to aperture system, rotate, and transform propagation and auxiliary vectors.
    let (center_of_mass, relative_vertices, rot3, prop2) =
        init_diff(verts, &mut ampl, prop, vk7, wavenumber);

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

    let v1 = Array2::from_shape_vec((nv, 3), v1_data).expect("Shape error creating v1");

    let x: Vec<f32> = v1.column(0).iter().cloned().collect();
    let y: Vec<f32> = v1.column(1).iter().cloned().collect();

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
    let mut ampl_cs = vec![Matrix2::<Complex<f32>>::default(); theta_phi_combinations.len()];
    let kinc = prop2 * wavenumber;

    // Pre-calculate constants dependent only on wavenumber
    // let radius = settings::RADIUS * 2.0 * PI / wavenumber;
    let inv_denom = Complex::new(wavenumber / (2.0 * PI), 0.0); // Pre-calculate 1.0 / (2*PI/wavenumber)

    // Iterate over the flattened combinations
    for (index, (theta, phi)) in theta_phi_combinations.iter().enumerate() {
        // Compute sin and cos values for current theta and phi
        let (sin_theta, cos_theta) = theta.to_radians().sin_cos();
        let (sin_phi, cos_phi) = phi.to_radians().sin_cos();

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

        let ampl_temp = rot4.map(Complex::from)
            * karczewski.map(Complex::from)
            * ampl
            * prerotation.map(Complex::from);

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

            if alpha.is_infinite() || beta.is_infinite() || alpha.is_nan() || beta.is_nan() {
                continue;
            }

            let summand = calculate_summand(bvsk, delta, omega1, omega2, alpha, beta, inv_denom); // Pass inv_denom

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
    let (karczewski, m) = karczewski(&prop2, &k);

    let hc = rot3 * Vector3::new(sin_phi, -cos_phi, 0.0);
    let evo2 = k.cross(&m);
    let rot4 = Matrix2::new(hc.dot(&m), -hc.dot(&evo2), hc.dot(&evo2), hc.dot(&m));

    let prerotation = Field::rotation_matrix(
        Vector3::x(),
        Vector3::new(-sin_phi, cos_phi, 0.0),
        -Vector3::z(),
    )
    .transpose();
    (karczewski, rot4, prerotation)
}

fn init_diff(
    verts: &[Point3<f32>],
    ampl: &mut Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
    wavenumber: f32,
) -> (Point3<f32>, Vec<Vector3<f32>>, Matrix3<f32>, Vector3<f32>) {
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

    let rot1 = get_rotation_matrix2(&relative_vertices);
    let prop1 = rot1 * prop;
    let perp1 = rot1 * vk7;
    let rot2 = calculate_rotation_matrix(prop1);
    let rot3 = rot2 * rot1;

    let prop2 = rot2 * prop1;
    let perp2 = rot2 * perp1;
    let e_par2 = perp2.cross(&prop2).normalize();

    if e_par2.z > settings::COLINEAR_THRESHOLD {
        *ampl = -*ampl;
    }
    (center_of_mass, relative_vertices, rot3, prop2)
}

pub fn get_rotation_matrix2(verts: &Vec<Vector3<f32>>) -> Matrix3<f32> {
    let a1 = verts[0];
    let b1 = verts[1];

    let theta1 = if a1.y.abs() > settings::COLINEAR_THRESHOLD {
        (a1[0] / a1[1]).atan()
    } else {
        PI / 4.0
    };

    let rot1 = Matrix3::new(
        theta1.cos(),
        -theta1.sin(),
        0.0,
        theta1.sin(),
        theta1.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    );

    let a2 = rot1 * a1;
    let b2 = rot1 * b1;

    let theta2 = if a2.y.abs() > settings::COLINEAR_THRESHOLD {
        -(a2[2] / a2[1]).atan()
    } else {
        -PI / 4.0
    };

    let rot2 = Matrix3::new(
        1.0,
        0.0,
        0.0,
        0.0,
        theta2.cos(),
        -theta2.sin(),
        0.0,
        theta2.sin(),
        theta2.cos(),
    );

    let a3 = rot2 * a2;
    let b3 = rot2 * b2;

    let theta3 = if b3.x.abs() > settings::COLINEAR_THRESHOLD {
        (b3[2] / b3[0]).atan()
    } else {
        PI / 4.0
    };

    let rot3 = Matrix3::new(
        theta3.cos(),
        0.0,
        theta3.sin(),
        0.0,
        1.0,
        0.0,
        -theta3.sin(),
        0.0,
        theta3.cos(),
    );

    let a4 = rot3 * a3;
    let b4 = rot3 * b3;

    let rot = if a4[0] * b4[1] - a4[1] * b4[0] > 0.0 {
        let rot4 = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
        rot4 * rot3 * rot2 * rot1
    } else {
        rot3 * rot2 * rot1
    };

    rot
}

#[inline]
pub fn karczewski(prop2: &Vector3<f32>, bvk: &Vector3<f32>) -> (Matrix2<f32>, Vector3<f32>) {
    let big_kx = prop2.x;
    let big_ky = prop2.y;
    let big_kz = prop2.z;

    let sqrt_1_minus_k2y2 = (1.0 - bvk.y.powi(2)).sqrt();
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

    let frac = ((1.0 - bvk.y.powi(2)) / (1.0 - big_ky.powi(2))).sqrt();
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
    bvsk: f32,
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

    let exp_factor = Complex::cis(bvsk); // Use cis for complex exponential

    exp_factor * Complex::new(sumre, sumim) * inv_denom // Multiply by inverse denominator
}
