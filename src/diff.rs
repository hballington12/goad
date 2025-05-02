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
    let v1 = Array2::from_shape_vec((nv, 3), v1_data).expect("Shape error creating v1");

    let x: Vec<f32> = v1.column(0).iter().cloned().collect();
    let y: Vec<f32> = v1.column(1).iter().cloned().collect();
    let m: Vec<f32> = (0..nv)
        .map(|j| {
            let next_j = (j + 1) % nv;
            let dx = x[next_j] - x[j];
            let dy = y[next_j] - y[j];
            if dx.abs() < settings::DIFF_DMIN {
                if dy.signum() == dx.signum() {
                    1e6
                } else {
                    -1e6
                }
            } else {
                dy / dx
            }
        })
        .collect();
    let n: Vec<f32> = m
        .iter()
        .map(|&mj| {
            if mj.abs() < 1e-6 {
                if mj.signum() > 0.0 {
                    1e6
                } else {
                    -1e6
                }
            } else {
                1.0 / mj
            }
        })
        .collect();
    // --- Optimizations End ---

    // Define the output variables.
    let mut ampl_cs = vec![Matrix2::<Complex<f32>>::default(); theta_phi_combinations.len()];
    let mut area_facs2 = vec![Complex::<f32>::default(); theta_phi_combinations.len()];
    let kinc = prop2 * wavenumber;

    // Iterate over the flattened combinations
    for (index, (theta, phi)) in theta_phi_combinations.iter().enumerate() {
        // Compute sin and cos values for current theta and phi
        let (sin_theta, cos_theta) = theta.to_radians().sin_cos();
        let (sin_phi, cos_phi) = phi.to_radians().sin_cos();

        // Calculate pos (xfar, yfar, zfar) for the current (theta, phi)
        let radius = settings::RADIUS * 2.0 * PI / wavenumber;
        let r_sin_theta = radius * sin_theta;
        let rotated_pos = rot3
            * (Vector3::new(
                r_sin_theta * cos_phi,
                r_sin_theta * sin_phi,
                -radius * cos_theta,
            ) - center_of_mass.coords);

        // Calculate distance to bins and bin unit vectors
        let bvs = rotated_pos.norm();
        let k = rotated_pos / bvs;
        let bvsk = bvs * wavenumber;
        let ampl_far_field = &mut ampl_cs[index];

        let (karczewski, rot4, prerotation) = get_rotations(rot3, prop2, sin_phi, cos_phi, k);

        let ampl_temp = rot4.map(Complex::from)
            * karczewski.map(Complex::from)
            * ampl
            * prerotation.map(Complex::from);

        *ampl_far_field = ampl_temp;

        let fraunhofer = &mut area_facs2[index];
        let mut fraunhofer_sum = Complex::new(0.0, 0.0);

        let (kxx, kyy) = calculate_kxx_kyy(
            &kinc
                .fixed_rows::<2>(0)
                .into_owned()
                .as_slice()
                .try_into()
                .unwrap(),
            &rotated_pos,
            bvsk,
            wavenumber,
        );

        for j in 0..nv {
            let mj = m[j];
            let nj = n[j];
            let xj = x[j];
            let yj = y[j];

            let (mj, nj) = adjust_mj_nj(mj, nj);

            let next_j = (j + 1) % nv;
            let xj_plus1 = x[next_j];
            let yj_plus1 = y[next_j];

            let dx = xj_plus1 - xj;
            let dy = yj_plus1 - yj;

            let dx = if dx.abs() < settings::DIFF_DMIN {
                settings::DIFF_DMIN * dx.signum()
            } else {
                dx
            };
            let dy = if dy.abs() < settings::DIFF_DMIN {
                settings::DIFF_DMIN * dy.signum()
            } else {
                dy
            };

            let (delta, delta1, delta2) = calculate_deltas(kxx, kyy, xj, yj, mj, nj);
            let (omega1, omega2) = calculate_omegas(dx, dy, delta1, delta2);
            let (alpha, beta) = calculate_alpha_beta(delta1, delta2, kxx, kyy);

            if alpha.is_infinite() || beta.is_infinite() || alpha.is_nan() || beta.is_nan() {
                continue;
            }

            let summand = calculate_summand(bvsk, delta, omega1, omega2, alpha, beta, wavenumber);

            fraunhofer_sum += summand;
        }

        *fraunhofer = fraunhofer_sum;

        *ampl_far_field *= *fraunhofer;
    }
    ampl_cs
}

// Other functions remain unchanged
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

pub fn calculate_rotation_matrix(prop1: Vector3<f32>) -> Matrix3<f32> {
    let angle = -prop1.y.atan2(prop1.x);
    let (sin_angle, cos_angle) = angle.sin_cos();

    Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    )
}

pub fn adjust_mj_nj(mj: f32, nj: f32) -> (f32, f32) {
    if mj.abs() > 1e6 || nj.abs() < 1e-6 {
        (1e6, 1e-6)
    } else if nj.abs() > 1e6 || mj.abs() < 1e-6 {
        (1e-6, 1e6)
    } else {
        (mj, nj)
    }
}

pub fn calculate_bvsk(rotated_pos: &Vector3<f32>, wavenumber: f32) -> f32 {
    wavenumber * (rotated_pos.x.powi(2) + rotated_pos.y.powi(2) + rotated_pos.z.powi(2)).sqrt()
}

pub fn calculate_kxx_kyy(
    kinc: &[f32; 2],
    rotated_pos: &Vector3<f32>,
    bvsk: f32,
    wavenumber: f32,
) -> (f32, f32) {
    let kxx = kinc[0] - wavenumber.powi(2) * rotated_pos.x / bvsk;
    let kyy = kinc[1] - wavenumber.powi(2) * rotated_pos.y / bvsk;

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
    wavenumber: f32,
) -> Complex<f32> {
    let (sin_delta, cos_delta) = delta.sin_cos();
    let (sin_delta_omega1, cos_delta_omega1) = (delta + omega1).sin_cos();
    let (sin_delta_omega2, cos_delta_omega2) = (delta + omega2).sin_cos();

    let sumim = alpha * (cos_delta - cos_delta_omega1) - beta * (cos_delta - cos_delta_omega2);
    let sumre = -alpha * (sin_delta - sin_delta_omega1) + beta * (sin_delta - sin_delta_omega2);

    let exp_factor = Complex::new(bvsk.cos(), bvsk.sin());
    let denom = Complex::new(2.0 * PI / wavenumber, 0.0);

    exp_factor * Complex::new(sumre, sumim) / denom
}
