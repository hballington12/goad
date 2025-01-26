use macroquad::prelude::*;
use miniquad::ElapsedQuery;
use nalgebra::{Complex, ComplexField, Matrix2, Matrix3, Point3, Vector, Vector3};
use ndarray::{Array1, Array2, Array3, Array4, Axis};
use pbt::field::Field;
use pbt::problem::Problem;
use pbt::{
    beam::Beam,
    config,
    geom::{self, Face},
};
use std::f32::consts::PI;
use std::f32::MAX;
use std::fs::File;
use std::io::BufWriter;
use std::io::{self, Write};

fn main() {
    let mut geom = geom::Geom::from_file("./examples/data/hex2.obj").unwrap();

    let projection = Vector3::new(0.0, -1.0, 0.0).normalize();
    // let e_perp = Vector3::z(); // choose e_perp along z-axis for now

    let lower_left = vec![-10.0, -2.0];
    let upper_right = vec![10.0, 2.0];
    let clip_vertices = vec![
        Point3::new(lower_left[0], 10.0, upper_right[1]),
        Point3::new(lower_left[0], 10.0, lower_left[1]),
        Point3::new(upper_right[0], 10.0, lower_left[1]),
        Point3::new(upper_right[0], 10.0, upper_right[1]),
    ];
    let mut clip = Face::new_simple(clip_vertices, None).unwrap();
    clip.data_mut().area =
        Some((upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]));
    geom.shapes[0].refr_index.re = 1.5;
    geom.shapes[0].refr_index.im = 0.00001;

    // let mut problem = Problem::new(
    //     geom,
    //     Beam::new_initial(clip, projection, Complex::new(1.00, 0.0), e_perp).unwrap(),
    // );

    // problem.solve_near();

    // pull rectangular face and print vertices
    let face = geom.shapes[0].faces[4].clone();

    let m11: Complex<f32> = Complex::new(0.5, 0.25);
    let m12: Complex<f32> = Complex::new(0.25, -0.45);
    let m21: Complex<f32> = Complex::new(0.85, 0.2);
    let m22: Complex<f32> = Complex::new(-0.5, 0.5);
    let ampl = Matrix2::new(m11, m12, m21, m22);
    let prop: Vector3<f32> = Vector3::new(0.5, 0.3, -0.2).normalize();
    // let prop: Vector3<f32> = Vector3::new(0.0, 0.0, 1.0).normalize();
    let vk7: Vector3<f32> = Vector3::new(1.0, 0.0, 0.0);
    let vk7 = vk7.cross(&prop).normalize();
    let verts = face.data().exterior.clone();

    let theta_phi_combinations = generate_theta_phi_combinations();
    let ampls = diffraction(&verts, ampl, prop, vk7, &theta_phi_combinations);
    writeup(&theta_phi_combinations, &ampls);
}

/// Generate theta and phi combinations
fn generate_theta_phi_combinations() -> Vec<(f32, f32)> {
    let thetas = Array1::linspace(0.2, std::f32::consts::PI, 50).insert_axis(ndarray::Axis(1)); // Reshape to (50, 1)
    let phis = Array1::linspace(0.0, 2.0 * std::f32::consts::PI, 60).insert_axis(ndarray::Axis(0)); // Reshape to (1, 60)

    // Flatten the combinations of theta and phi into a 1D array of tuples
    thetas
        .iter()
        .flat_map(|&theta| phis.iter().map(move |&phi| (theta, phi)))
        .collect()
}

/// Diffraction. face in must be convex!
fn diffraction(
    verts: &[Point3<f32>],
    mut ampl: Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
    theta_phi_combinations: &[(f32, f32)],
) -> Vec<Matrix2<Complex<f32>>> {
    // Example theta and phi (replace later with actual values)
    // let thetas = Array2::from_shape_vec((2, 1), vec![0.1, 0.2]).unwrap(); // Theta
    // let thetas = Array1::linspace(0.2, std::f32::consts::PI, 50).insert_axis(ndarray::Axis(1)); // Reshape to (50, 1)
    // let phis = Array1::linspace(0.0, 2.0 * std::f32::consts::PI, 60).insert_axis(ndarray::Axis(0)); // Reshape to (1, 60)

    // Flatten the combinations of theta and phi into a 1D array of tuples
    // let theta_phi_combinations: Vec<(f32, f32)> = thetas
    //     .iter()
    //     .flat_map(|&theta| phis.iter().map(move |&phi| (theta, phi)))
    //     .collect();

    // 1. Compute the center of mass
    let center_of_mass = calculate_center_of_mass(verts);

    // 2. Transform vertices relative to the center of mass
    let relative_vertices = transform_to_center_of_mass(verts, &center_of_mass);

    // 3. Compute rotation matrices
    let rot = get_rotation_matrix2(&relative_vertices);
    let prop1 = rot * prop;
    let perp1 = rot * vk7;
    let rot2 = calculate_rotation_matrix(prop1);

    // 4. Transform propagation and auxiliary vectors
    let prop2 = rot2 * prop1;
    let perp2 = rot2 * perp1;
    let e_par2 = perp2.cross(&prop2).normalize();

    // 5. Update amplitude based on anti-parallel condition
    if e_par2.z > config::COLINEAR_THRESHOLD {
        ampl = -ampl;
    }

    let rot3 = rot2 * rot;
    let incidence = Vector3::new(0.0, 0.0, -1.0); // TODO: generalise
    let incidence2 = rot3 * incidence;

    // Constants
    const RADIUS: f32 = 1e4;

    // Define a 1D array with length theta_phi_combinations.len() for Complex<f32>
    let mut ampl_cs = vec![Matrix2::<Complex<f32>>::default(); theta_phi_combinations.len()];
    let mut area_facs2 = vec![Complex::<f32>::default(); theta_phi_combinations.len()];

    // Iterate over the flattened combinations
    for (index, (theta, phi)) in theta_phi_combinations.iter().enumerate() {
        // Compute sin and cos values for current theta and phi
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        // Calculate xfar, yfar, zfar for the current (theta, phi)
        let r_sin_theta = RADIUS * sin_theta;
        let xfar = r_sin_theta * cos_phi;
        let yfar = r_sin_theta * sin_phi;
        let zfar = RADIUS * cos_theta;

        // Translate far-field bins to the aperture system
        let x1 = xfar - center_of_mass[0];
        let y1 = yfar - center_of_mass[1];
        let z1 = zfar - center_of_mass[2];

        // Rotate the bins using rot3 matrix (assuming rot3 is a 3x3 matrix)
        let pos = Vector3::new(x1, y1, z1);
        let rotated_pos = rot3 * pos; // Use matrix-vector multiplication

        // Use the calculated x3, y3, z3 for further processing in the loop
        let amplc = &mut ampl_cs[index];
        *amplc = Matrix2::identity(); // Example, modify as needed

        // Call karczewski for each element
        let (diff_ampl, m, k) =
            karczewski(&prop2, rotated_pos.x, rotated_pos.y, rotated_pos.z, RADIUS);

        let hc = if incidence2.dot(&k).abs() < 1.0 - config::COLINEAR_THRESHOLD {
            incidence2.cross(&k).normalize()
        } else {
            print!("warn");
            incidence2.cross(&k).normalize()
        };

        let evo2 = k.cross(&m);

        // Initialize a 2x2 matrix for rot4
        let rot4 = Matrix2::new(hc.dot(&m), -hc.dot(&evo2), hc.dot(&evo2), hc.dot(&m));

        let temp_vec3 = Vector3::new(cos_phi, sin_phi, 0.0);
        let temp_rot1 = Field::rotation_matrix(
            Vector3::x(),
            Vector3::new(-temp_vec3[1], temp_vec3[0], 0.0),
            -Vector3::z(),
        );

        let temp_rot1 = temp_rot1.transpose();

        let rot4_complex = rot4.map(|x| Complex::new(x, 0.0));
        let diff_ampl_complex = diff_ampl.map(|x| Complex::new(x, 0.0));
        let temp_rot1_complex = temp_rot1.map(|x| Complex::new(x, 0.0));

        let ampl_temp2 = rot4_complex * diff_ampl_complex * ampl * temp_rot1_complex;

        amplc[(0, 0)] = ampl_temp2[(0, 0)];
        amplc[(1, 0)] = ampl_temp2[(1, 0)];
        amplc[(0, 1)] = ampl_temp2[(0, 1)];
        amplc[(1, 1)] = ampl_temp2[(1, 1)];

        let nv = relative_vertices.len();

        let mut v1 = Array2::<f32>::zeros((relative_vertices.len(), relative_vertices[0].len()));

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

        let area_fac = &mut area_facs2[index];
        let mut area_fac_sum = Complex::new(0.0, 0.0); // Example, modify as needed

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

            let bvsk = calculate_bvsk(&rotated_pos);
            let (kxx, kyy) = calculate_kxx_kyy(&kinc.fixed_rows::<2>(0).into(), &rotated_pos, bvsk);
            let (delta, delta1, delta2) = calculate_deltas(kxx, kyy, xj, yj, mj, nj);
            let (omega1, omega2) = calculate_omegas(dx, dy, delta1, delta2);
            let (alpha, beta) = calculate_alpha_beta(delta1, delta2, kxx, kyy);
            let summand = calculate_summand(bvsk, delta, omega1, omega2, alpha, beta);

            area_fac_sum += summand;
        }

        *area_fac = area_fac_sum;

        amplc[(0, 0)] *= *area_fac;
        amplc[(1, 0)] *= *area_fac;
        amplc[(0, 1)] *= *area_fac;
        amplc[(1, 1)] *= *area_fac;
    }
    ampl_cs
}
fn writeup(theta_phi_combinations: &[(f32, f32)], ampl_cs: &Vec<Matrix2<Complex<f32>>>) {
    println!("done.");
    let mut s11 = Array1::<f32>::zeros(theta_phi_combinations.len());

    for (index, amplc) in ampl_cs.iter().enumerate() {
        s11[index] = (Complex::new(0.5, 0.0)
            * (amplc[(0, 0)] * amplc[(0, 0)].conj()
                + amplc[(0, 1)] * amplc[(0, 1)].conj()
                + amplc[(1, 0)] * amplc[(1, 0)].conj()
                + amplc[(1, 1)] * amplc[(1, 1)].conj()))
        .real()
    }

    // Open a file for writing
    let file = File::create("output.txt").unwrap();
    let mut writer = BufWriter::new(file);

    // Write header
    // writeln!(writer, "theta,phi,s11");

    // Iterate over the array and write data to the file
    for (index, val) in s11.iter().enumerate() {
        let (theta, phi) = theta_phi_combinations[index];
        writeln!(
            writer,
            "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            theta * 180.0 / PI, // Theta at index i
            phi * 180.0 / PI,   // Phi at index j
            val,                // s11 value
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
        );
    }
}

fn get_rotation_matrix2(verts: &Vec<Vector3<f32>>) -> Matrix3<f32> {
    let a1 = verts[0];
    let b1 = verts[1];

    println!("a1: {}", a1);
    println!("b1: {}", b1);

    let theta1 = if a1.y.abs() > config::COLINEAR_THRESHOLD {
        (a1[0] / a1[1]).atan()
    } else {
        PI / 4.0
    };

    println!("theta1: {}", theta1);

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

    println!("a2: {}", a2);
    println!("b2: {}", b2);

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

    println!("a3: {}", a3);
    println!("b3: {}", b3);

    let theta3 = if b3.x.abs() > config::COLINEAR_THRESHOLD {
        (b3[2] / b3[0]).atan()
    } else {
        PI / 4.0
    };

    println!("theta3: {}", theta3);

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

    println!("a4: {}", a4);
    println!("b4: {}", b4);

    let rot = if a4[0] * b4[1] - a4[1] * b4[0] > 0.0 {
        let rot4 = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
        rot4 * rot3 * rot2 * rot1
    } else {
        rot3 * rot2 * rot1
    };

    rot
}

fn karczewski(
    prop2: &Vector3<f32>, // Outgoing propagation direction in aperture system
    x3: f32,              // x coordinate of far-field bin
    y3: f32,              // y coordinate of far-field bin
    z3: f32,              // z coordinate of far-field bin
    r: f32,               // Distance to bin vector
) -> (Matrix2<f32>, Vector3<f32>, Vector3<f32>) {
    // Compute the diff ampl matrix for polarisation of far-field diffraction at a far-field bin
    // TODO: Add debug mode with numerical checks (as mentioned in the original MATLAB code)

    // Far-field bin distance
    let bin_vec_size = r;

    // Propagation vector components for each bin vector in the aperture system
    let mut k = Vector3::new(x3 / bin_vec_size, y3 / bin_vec_size, z3 / bin_vec_size);

    // Ensure k.y is within bounds
    if k.y.abs() > 0.999_999 {
        k.y = 0.999_999_f32.copysign(k.y);
    }

    // Propagation direction in aperture system
    let big_kx = prop2.x;
    let big_ky = prop2.y;
    let big_kz = prop2.z;

    // Perpendicular field direction
    let sqrt_1_minus_k2y2 = (1.0 - k.y.powi(2)).sqrt();
    let m = Vector3::new(
        -k.x * k.y / sqrt_1_minus_k2y2,
        sqrt_1_minus_k2y2,
        -k.y * k.z / sqrt_1_minus_k2y2,
    );

    // Pre-calculate factor
    let frac = ((1.0 - k.y.powi(2)) / (1.0 - big_ky.powi(2))).sqrt();

    // KW coefficients
    let a1m = -big_kz * frac;
    let b2m = -k.z / frac;
    let a1e = b2m;
    let b2e = a1m;
    let b1m = -k.x * k.y / frac + big_kx * big_ky * frac;
    let a2e = -b1m;

    // Combined (e-m theory) KW coefficients
    let a1em = 0.5 * (a1m + a1e);
    let a2em = 0.5 * a2e;
    let b1em = 0.5 * b1m;
    let b2em = 0.5 * (b2m + b2e);

    // Fill the diff_ampl matrix (e-m theory)
    let diff_ampl = Matrix2::new(a1em, b1em, a2em, b2em);

    // Return the outputs
    (diff_ampl, m, k)
}

fn calculate_center_of_mass(verts: &[Point3<f32>]) -> Point3<f32> {
    Point3::from(
        verts
            .iter()
            .map(|vert| vert.coords)
            .fold(Vector3::zeros(), |acc, coords| acc + coords)
            / verts.len() as f32,
    )
}

fn transform_to_center_of_mass(
    verts: &[Point3<f32>],
    center_of_mass: &Point3<f32>,
) -> Vec<Vector3<f32>> {
    verts
        .iter()
        .map(|point| point.coords - center_of_mass.coords)
        .collect()
}

fn calculate_rotation_matrix(prop1: Vector3<f32>) -> Matrix3<f32> {
    let angle = -prop1.y.atan2(prop1.x);
    let (sin_angle, cos_angle) = angle.sin_cos();

    Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    )
}

fn adjust_mj_nj(mj: f32, nj: f32) -> (f32, f32) {
    if mj.abs() > 1e6 || nj.abs() < 1e-6 {
        (1e6, 1e6)
    } else if nj.abs() > f32::MAX || mj.abs() < 1e-6 {
        (1e6, 1e6)
    } else {
        (mj, nj)
    }
}

fn calculate_bvsk(rotated_pos: &Vector3<f32>) -> f32 {
    config::WAVENO * (rotated_pos.x.powi(2) + rotated_pos.y.powi(2) + rotated_pos.z.powi(2)).sqrt()
}

fn calculate_kxx_kyy(kinc: &[f32; 2], rotated_pos: &Vector3<f32>, bvsk: f32) -> (f32, f32) {
    let kxx = kinc[0] - config::WAVENO.powi(2) * rotated_pos.x / bvsk;
    let kyy = kinc[1] - config::WAVENO.powi(2) * rotated_pos.y / bvsk;
    (kxx, kyy)
}

fn calculate_deltas(kxx: f32, kyy: f32, xj: f32, yj: f32, mj: f32, nj: f32) -> (f32, f32, f32) {
    let delta = kxx * xj + kyy * yj;
    let delta1 = kyy * mj + kxx;
    let delta2 = kxx * nj + kyy;
    (delta, delta1, delta2)
}

fn calculate_omegas(dx: f32, dy: f32, delta1: f32, delta2: f32) -> (f32, f32) {
    let omega1 = dx * delta1;
    let omega2 = dy * delta2;
    (omega1, omega2)
}

fn calculate_alpha_beta(delta1: f32, delta2: f32, kxx: f32, kyy: f32) -> (f32, f32) {
    let alpha = 1.0 / (2.0 * kyy * delta1);
    let beta = 1.0 / (2.0 * kxx * delta2);
    (alpha, beta)
}

fn calculate_summand(
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
