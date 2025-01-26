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

    diffraction(&verts, ampl, prop, vk7);
}

fn diffraction(
    verts: &[Point3<f32>],
    mut ampl: Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
) {
    let num_verts = verts.len();

    // Sum all vertex coordinates and compute the center of mass
    // Compute the center of mass directly without intermediate storage
    let center_of_mass = Point3::from(
        verts
            .iter()
            .map(|vert| vert.coords)
            .fold(Vector3::zeros(), |acc, coords| acc + coords)
            / num_verts as f32,
    );

    // Transform vertices relative to the center of mass
    let v20: Vec<Vector3<f32>> = verts
        .iter()
        .map(|point| point.coords - center_of_mass.coords)
        .collect();

    // Rotation matrix from custom function
    let rot = get_rotation_matrix2(&v20);

    // Transform propagation and auxiliary vectors
    let prop1 = rot * prop;
    let perp1 = rot * vk7;

    let angle = -prop1.y.atan2(prop1.x);
    let cos_angle = angle.cos();
    let sin_angle = angle.sin();

    // Create the rotation matrix directly
    let rot2 = Matrix3::new(
        cos_angle, -sin_angle, 0.0, sin_angle, cos_angle, 0.0, 0.0, 0.0, 1.0,
    );

    let prop2 = rot2 * prop1;
    let perp2 = rot2 * perp1;
    let e_par2 = perp2.cross(&prop2).normalize();

    // Check anti-parallel condition and modify amplitude in-place
    if e_par2.z > config::COLINEAR_THRESHOLD {
        ampl = -ampl;
    }

    let incidence = Vector3::new(0.0, 0.0, -1.0); // !todo: generalise

    let incidence2 = rot2 * rot * incidence;

    // Constants
    let r = 1e4;

    // Example theta and phi (replace later with actual values)
    // let thetas = Array2::from_shape_vec((2, 1), vec![0.1, 0.2]).unwrap(); // Theta
    // let phis = Array2::from_shape_vec((1, 2), vec![0.1, 0.4]).unwrap(); // Phi
    let thetas = Array1::linspace(0.2, std::f32::consts::PI, 50).insert_axis(ndarray::Axis(1)); // Reshape to (100, 1)

    // Phi: 100 steps from 0 to 2pi
    let phis = Array1::linspace(0.0, 2.0 * std::f32::consts::PI, 50).insert_axis(ndarray::Axis(0)); // Reshape to (1, 100)

    // Define a 4D array with shape (thetas.len(), phis.len(), 2, 2) for Complex<f32>
    let mut amplCs = Array2::<Matrix2<Complex<f32>>>::default((thetas.len(), phis.len()));

    let mut area_facs2 = Array2::<Complex<f32>>::zeros((thetas.len(), phis.len()));

    // Compute xfar, yfar, zfar
    let sin_theta = thetas.mapv(f32::sin);
    let cos_theta = thetas.mapv(f32::cos);
    let sin_phi = phis.mapv(f32::sin);
    let cos_phi = phis.mapv(f32::cos);

    let xfar = r * &sin_theta * &cos_phi;
    let yfar = r * &sin_theta * &sin_phi;
    let zfar = r * &cos_theta;

    // Translate far-field bins to the aperture system
    let x1 = &xfar - center_of_mass[0];
    let y1 = &yfar - center_of_mass[1];
    let z1 = &zfar - center_of_mass[2];

    // Compute r1 (distance from the center of mass)
    let r1 = (&x1.mapv(|v| v.powi(2)) + &y1.mapv(|v| v.powi(2)) + &z1.mapv(|v| v.powi(2)))
        .mapv(f32::sqrt);

    // TODO: numerical bodge for rot4 matrix here
    // rotate bins

    let rot3 = rot2 * rot;

    let x3 = rot3[(0, 0)] * &x1 + rot3[(0, 1)] * &y1 + rot3[(0, 2)] * &z1;
    let y3 = rot3[(1, 0)] * &x1 + rot3[(1, 1)] * &y1 + rot3[(1, 2)] * &z1;
    let z3 = rot3[(2, 0)] * &x1 + rot3[(2, 1)] * &y1 + rot3[(2, 2)] * &z1;

    for ((i, j), amplc) in amplCs.indexed_iter_mut() {
        let x3_val = &x3[(i, j)];
        let y3_val = &y3[(i, j)];
        let z3_val = &z3[(i, j)];
        let phi_val = &phis[(0, i)];
        let area_fac2 = &mut area_facs2[(i, j)];

        // Call karczewski for each element
        let (diff_ampl, m, k) = karczewski(&prop2, *x3_val, *y3_val, *z3_val, r);

        let hc = if incidence2.dot(&k).abs() < 0.999 {
            incidence2.cross(&k).normalize()
        } else {
            print!("warn: need to implement this");
            incidence2.cross(&k).normalize()
        };

        let evo2 = k.cross(&m);

        // Initialize a 2x2 matrix for rot4
        let mut rot4 = Matrix2::zeros();

        // Populate rot4 using the dot products
        rot4[(0, 0)] = hc.dot(&m); // rot4(1,1) = dot_product(hc, m)
        rot4[(1, 1)] = hc.dot(&m); // rot4(2,2) = dot_product(hc, m)
        rot4[(0, 1)] = -hc.dot(&evo2); // rot4(1,2) = -dot_product(hc, evo2)
        rot4[(1, 0)] = hc.dot(&evo2); // rot4(2,1) = +dot_product(hc, evo2)

        let temp_vec3 = Vector3::new(phi_val.cos(), phi_val.sin(), 0.0);

        // TODO this has normalisation, which is not present in Fortran version,
        // but it appears like normalisation should be present, so probably
        // uncomment this back in later...
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

        let rot = rot2 * rot;

        let nv = v20.len();

        let mut v1 = Array2::<f32>::zeros((v20.len(), v20[0].len()));

        for (i, vertex) in v20.iter().enumerate() {
            let transformed_vertex = rot * vertex;
            v1[[i, 0]] = transformed_vertex.x;
            v1[[i, 1]] = transformed_vertex.y;
            v1[[i, 2]] = transformed_vertex.z;
        }

        let kinc = prop2 * config::WAVENO;

        let mut x = vec![0.0; nv];
        let mut y = vec![0.0; nv];
        for i in 0..nv {
            x[i] = v1[[i, 0]];
            y[i] = v1[[i, 1]];
        }

        let mut m = vec![0.0; nv];
        let mut n = vec![0.0; nv];
        for j in 0..nv {
            if j == nv - 1 {
                // Special case for the last vertex
                m[j] = (y[0] - y[j]) / (x[0] - x[j]);
            } else {
                // Compute gradient for other vertices
                m[j] = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
            }
            // Compute inverse gradient
            n[j] = 1.0 / m[j];
        }

        for (j, vertex) in v1.rows().into_iter().enumerate() {
            let mut mj = m[j];
            let mut nj = n[j];
            let mut xj = x[j];
            let mut yj = y[j];

            let mj = if mj.abs() > f32::MAX { 1e6 } else { mj };
            let nj = if mj.abs() > f32::MAX { 1e6 } else { nj };
            let mj = if nj.abs() < 1e-9 { 1e6 } else { mj };
            let nj = if mj.abs() < 1e-9 { 1e6 } else { nj };

            let (xj_plus1, yj_plus1) = if j == nv - 1 {
                (x[0], y[0])
            } else {
                (x[j + 1], y[j + 1])
            };

            let dx = xj_plus1 - xj;
            let dy = yj_plus1 - yj;

            let bvsk = config::WAVENO * (x3_val.powi(2) + y3_val.powi(2) + z3_val.powi(2)).sqrt();

            let kxx = kinc[0] - config::WAVENO.powi(2) * x3_val / bvsk;
            let kyy = kinc[1] - config::WAVENO.powi(2) * y3_val / bvsk;

            let delta = kxx * xj + kyy * yj;
            let delta1 = kyy * mj + kxx;
            let delta2 = kxx * nj + kyy;
            let omega1 = dx * delta1;
            let omega2 = dy * delta2;

            let alpha = 1.0 / (2.0 * kyy * delta1);
            let beta = 1.0 / (2.0 * kxx * delta2);

            let sumim = alpha * (delta.cos() - (delta + omega1).cos())
                - beta * (delta.cos() - (delta + omega2).cos());
            let sumre = -alpha * (delta.sin() - (delta + omega1).sin())
                + beta * (delta.sin() - (delta + omega2).sin());

            let summand = Complex::new((bvsk).cos(), (bvsk).sin()) * Complex::new(sumre, sumim)
                / Complex::new(config::WAVELENGTH, 0.0);

            *area_fac2 += summand;
        }

        amplc[(0, 0)] *= *area_fac2;
        amplc[(1, 0)] *= *area_fac2;
        amplc[(0, 1)] *= *area_fac2;
        amplc[(1, 1)] *= *area_fac2;
    }

    println!("done.");
    let mut s11 = Array2::<f32>::zeros((thetas.len(), phis.len()));

    for ((i, j), amplc) in amplCs.indexed_iter_mut() {
        s11[(i, j)] = (Complex::new(0.5, 0.0)
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
    for ((i, j), val) in s11.indexed_iter_mut() {
        writeln!(
            writer,
            "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            thetas[[i, 0]] * 180.0 / PI, // Theta at index i
            phis[[0, j]] * 180.0 / PI,   // Phi at index j
            val,                         // s11 value
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
