use nalgebra::{Complex, Matrix2, Matrix3, Point3, Vector3};
use ndarray::Array2;
use std::f32::consts::PI;

use crate::field::Field;
use crate::{config, geom};

// Diffraction. Face must be convex.
pub fn diffraction(
    verts: &[Point3<f32>],
    mut ampl: Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
    theta_phi_combinations: &[(f32, f32)],
) -> Vec<Matrix2<Complex<f32>>> {
    // Translate to aperture system, rotate, and transform propagation and auxiliary vectors.
    let (center_of_mass, relative_vertices, rot3, prop2) = init_diff(verts, &mut ampl, prop, vk7);

    // Define the output variables.
    let mut ampl_cs = vec![Matrix2::<Complex<f32>>::default(); theta_phi_combinations.len()];
    let mut area_facs2 = vec![Complex::<f32>::default(); theta_phi_combinations.len()];

    // Iterate over the flattened combinations
    for (index, (theta, phi)) in theta_phi_combinations.iter().enumerate() {
        // Compute sin and cos values for current theta and phi
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        // Calculate pos (xfar, yfar, zfar) for the current (theta, phi)
        let r_sin_theta = config::RADIUS * sin_theta;
        let rotated_pos = rot3
            * (Vector3::new(
                r_sin_theta * cos_phi,
                r_sin_theta * sin_phi,
                -config::RADIUS * cos_theta,
            ) - center_of_mass.coords);

        // Calculate distance to bins and bin unit vectors
        let bvs = rotated_pos.norm();
        let k = rotated_pos / bvs;
        let bvsk = bvs * config::WAVENO;
        let ampl_far_field = &mut ampl_cs[index];

        let (polarisation, rot4, prerotation) = get_rotations(rot3, prop2, sin_phi, cos_phi, k);

        let ampl_temp = rot4.map(Complex::from)
            * polarisation.map(Complex::from)
            * ampl
            * prerotation.map(Complex::from);

        ampl_far_field[(0, 0)] = ampl_temp[(0, 0)];
        ampl_far_field[(1, 0)] = ampl_temp[(1, 0)];
        ampl_far_field[(0, 1)] = ampl_temp[(0, 1)];
        ampl_far_field[(1, 1)] = ampl_temp[(1, 1)];

        let nv = relative_vertices.len();
        let mut v1 = Array2::<f32>::zeros((nv, 3));
        for (i, vertex) in relative_vertices.iter().enumerate() {
            let transformed_vertex = rot3 * vertex;
            v1[[i, 0]] = transformed_vertex.x;
            v1[[i, 1]] = transformed_vertex.y;
            v1[[i, 2]] = transformed_vertex.z;
        }

        let kinc = prop2 * config::WAVENO;
        let x: Vec<f32> = v1.column(0).iter().cloned().collect();
        let y: Vec<f32> = v1.column(1).iter().cloned().collect();
        let m: Vec<f32> = (0..nv)
            .map(|j| {
                if j == nv - 1 {
                    (y[0] - y[j]) / (x[0] - x[j])
                } else {
                    (y[j + 1] - y[j]) / (x[j + 1] - x[j])
                }
            })
            .collect();
        let n: Vec<f32> = m.iter().map(|&mj| 1.0 / mj).collect();

        let fraunhofer = &mut area_facs2[index];
        let mut fraunhofer_sum = Complex::new(0.0, 0.0); // Example, modify as needed

        for (j, vertex) in v1.rows().into_iter().enumerate() {
            let mj = m[j];
            let nj = n[j];
            let xj = vertex[0];
            let yj = vertex[1];

            // Adjust mj and nj
            let (mj, nj) = adjust_mj_nj(mj, nj);

            let (xj_plus1, yj_plus1) = if j == nv - 1 {
                (x[0], y[0])
            } else {
                (x[j + 1], y[j + 1])
            };

            let dx = xj_plus1 - xj;
            let dy = yj_plus1 - yj;

            let dx = if dx.abs() < config::DIFF_EPSILON {
                config::DIFF_EPSILON
            } else {
                dx
            };
            let dy = if dy.abs() < config::DIFF_EPSILON {
                config::DIFF_EPSILON
            } else {
                dy
            };

            let (kxx, kyy) = calculate_kxx_kyy(&kinc.fixed_rows::<2>(0).into(), &rotated_pos, bvsk);
            let (delta, delta1, delta2) = calculate_deltas(kxx, kyy, xj, yj, mj, nj);
            let (omega1, omega2) = calculate_omegas(dx, dy, delta1, delta2);
            let (alpha, beta) = calculate_alpha_beta(delta1, delta2, kxx, kyy);

            if alpha.is_infinite() || beta.is_infinite() {
                continue;
            }

            let summand = calculate_summand(bvsk, delta, omega1, omega2, alpha, beta);

            fraunhofer_sum += summand;
        }

        *fraunhofer = fraunhofer_sum;

        ampl_far_field[(0, 0)] *= *fraunhofer;
        ampl_far_field[(1, 0)] *= *fraunhofer;
        ampl_far_field[(0, 1)] *= *fraunhofer;
        ampl_far_field[(1, 1)] *= *fraunhofer;
    }
    ampl_cs
}

fn get_rotations(
    rot3: Matrix3<f32>,
    prop2: Vector3<f32>,
    sin_phi: f32,
    cos_phi: f32,
    k: Vector3<f32>,
) -> (Matrix2<f32>, Matrix2<f32>, Matrix2<f32>) {
    // Compute Karczewski polarisation matrix for each element
    let (polarisation, m) = karczewski(&prop2, &k);

    // Vector perpendicular to the scattering plane in the aperture system
    let hc = rot3 * Vector3::new(sin_phi, -cos_phi, 0.0);
    let evo2 = k.cross(&m);
    let rot4 = Matrix2::new(hc.dot(&m), -hc.dot(&evo2), hc.dot(&evo2), hc.dot(&m));

    let prerotation = Field::rotation_matrix(
        -Vector3::x(),
        Vector3::new(sin_phi, -cos_phi, 0.0),
        Vector3::z(),
    );
    (polarisation, rot4, prerotation)
}

fn init_diff(
    verts: &[Point3<f32>],
    ampl: &mut Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
) -> (Point3<f32>, Vec<Vector3<f32>>, Matrix3<f32>, Vector3<f32>) {
    // -1. Apply a small perturbation to the propagation vector to reduce numerical errors
    let prop = (prop
        + Vector3::new(
            config::DIFF_EPSILON,
            config::DIFF_EPSILON,
            config::DIFF_EPSILON,
        ))
    .normalize();

    // 0. Account for 1/waveno factor in Bohren & Huffman eq 3.12
    *ampl *= Complex::new(config::WAVENO, 0.0);

    // 1. Compute the center of mass
    let center_of_mass = geom::calculate_center_of_mass(verts);

    // 2. Transform vertices relative to the center of mass
    let relative_vertices = geom::transform_to_center_of_mass(verts, &center_of_mass);

    // 3. Compute rotation matrices
    let rot1 = get_rotation_matrix2(&relative_vertices);
    let prop1 = rot1 * prop;
    let perp1 = rot1 * vk7;
    let rot2 = calculate_rotation_matrix(prop1);
    let rot3 = rot2 * rot1;

    // 4. Transform propagation and auxiliary vectors
    let prop2 = rot2 * prop1;
    let perp2 = rot2 * perp1;
    let e_par2 = perp2.cross(&prop2).normalize();

    // 5. Update amplitude based on anti-parallel condition
    if e_par2.z > config::COLINEAR_THRESHOLD {
        *ampl = -*ampl;
    }
    (center_of_mass, relative_vertices, rot3, prop2)
}

pub fn get_rotation_matrix2(verts: &Vec<Vector3<f32>>) -> Matrix3<f32> {
    let a1 = verts[0];
    let b1 = verts[1];

    let theta1 = if a1.y.abs() > config::COLINEAR_THRESHOLD {
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

    let theta2 = if a2.y.abs() > config::COLINEAR_THRESHOLD {
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

    let theta3 = if b3.x.abs() > config::COLINEAR_THRESHOLD {
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

pub fn karczewski(
    prop2: &Vector3<f32>, // Outgoing propagation direction in aperture system
    bvk: &Vector3<f32>,   // bin unit vector
) -> (Matrix2<f32>, Vector3<f32>) {
    // Compute the diff ampl matrix for polarisation of far-field diffraction at a far-field bin

    // Propagation direction in aperture system
    let big_kx = prop2.x;
    let big_ky = prop2.y;
    let big_kz = prop2.z;

    // Perpendicular field direction
    let sqrt_1_minus_k2y2 = (1.0 - bvk.y.powi(2)).sqrt();
    let m = Vector3::new(
        -bvk.x * bvk.y / sqrt_1_minus_k2y2,
        sqrt_1_minus_k2y2,
        -bvk.y * bvk.z / sqrt_1_minus_k2y2,
    );

    // Pre-calculate factor
    let frac = ((1.0 - bvk.y.powi(2)) / (1.0 - big_ky.powi(2))).sqrt();

    // KW coefficients
    let a1m = -big_kz * frac;
    let b2m = -bvk.z / frac;
    let a1e = b2m;
    let b2e = a1m;
    let b1m = -bvk.x * bvk.y / frac + big_kx * big_ky * frac;
    let a2e = -b1m;

    // Combined (e-m theory) KW coefficients
    let a1em = 0.5 * (a1m + a1e);
    let a2em = 0.5 * a2e;
    let b1em = 0.5 * b1m;
    let b2em = 0.5 * (b2m + b2e);

    // Fill the diff_ampl matrix (e-m theory)
    let diff_ampl = Matrix2::new(a1em, b1em, a2em, b2em);

    // Return the outputs
    (diff_ampl, m)
}

pub fn calculate_rotation_matrix(prop1: Vector3<f32>) -> Matrix3<f32> {
    let angle = -prop1.y.atan2(prop1.x);
    let (sin_angle, cos_angle) = angle.sin_cos();

    Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    )
}

pub fn adjust_mj_nj(mj: f32, nj: f32) -> (f32, f32) {
    if mj.abs() > 1e6 || nj.abs() < 1e-6 {
        (1e6, 1e6)
    } else if nj.abs() > f32::MAX || mj.abs() < 1e-6 {
        (1e6, 1e6)
    } else {
        (mj, nj)
    }
}

pub fn calculate_bvsk(rotated_pos: &Vector3<f32>) -> f32 {
    config::WAVENO * (rotated_pos.x.powi(2) + rotated_pos.y.powi(2) + rotated_pos.z.powi(2)).sqrt()
}

pub fn calculate_kxx_kyy(kinc: &[f32; 2], rotated_pos: &Vector3<f32>, bvsk: f32) -> (f32, f32) {
    let kxx = kinc[0] - config::WAVENO.powi(2) * rotated_pos.x / bvsk;
    let kyy = kinc[1] - config::WAVENO.powi(2) * rotated_pos.y / bvsk;

    // numerical fix for kxx and kyy
    let kxx = if kxx.abs() < config::DIFF_EPSILON {
        config::DIFF_EPSILON
    } else {
        kxx
    };
    let kyy = if kyy.abs() < config::DIFF_EPSILON {
        config::DIFF_EPSILON
    } else {
        kyy
    };

    (kxx, kyy)
}

pub fn calculate_deltas(kxx: f32, kyy: f32, xj: f32, yj: f32, mj: f32, nj: f32) -> (f32, f32, f32) {
    let delta = kxx * xj + kyy * yj;
    let delta1 = kyy * mj + kxx;
    let delta2 = kxx * nj + kyy;
    (delta, delta1, delta2)
}

pub fn calculate_omegas(dx: f32, dy: f32, delta1: f32, delta2: f32) -> (f32, f32) {
    let omega1 = dx * delta1;
    let omega2 = dy * delta2;
    (omega1, omega2)
}

pub fn calculate_alpha_beta(delta1: f32, delta2: f32, kxx: f32, kyy: f32) -> (f32, f32) {
    let alpha = 1.0 / (2.0 * kyy * delta1);
    let beta = 1.0 / (2.0 * kxx * delta2);
    (alpha, beta)
}

pub fn calculate_summand(
    bvsk: f32,
    delta: f32,
    omega1: f32,
    omega2: f32,
    alpha: f32,
    beta: f32,
) -> Complex<f32> {
    let sumim = alpha * (delta.cos() - (delta + omega1).cos())
        - beta * (delta.cos() - (delta + omega2).cos());
    let sumre = -alpha * (delta.sin() - (delta + omega1).sin())
        + beta * (delta.sin() - (delta + omega2).sin());

    Complex::new((bvsk).cos(), (bvsk).sin()) * Complex::new(sumre, sumim)
        / Complex::new(config::WAVELENGTH, 0.0)
}
