use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Matrix3, Point3, Vector, Vector3};
use pbt::problem::Problem;
use pbt::{
    beam::Beam,
    config,
    geom::{self, Face},
};
use std::f32::consts::PI;
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
    let face = geom.shapes[0].faces[0].clone();

    let m11: Complex<f32> = Complex::new(0.5, 0.25);
    let m12: Complex<f32> = Complex::new(0.25, -0.45);
    let m21: Complex<f32> = Complex::new(0.85, 0.2);
    let m22: Complex<f32> = Complex::new(-0.5, 0.5);
    let ampl = Matrix2::new(m11, m12, m21, m22);
    let prop: Vector3<f32> = Vector3::new(0.5, 0.3, -0.2).normalize();
    let vk7: Vector3<f32> = Vector3::new(0.0, 0.0, 0.1);
    let vk7 = vk7.cross(&prop).normalize();
    let verts = face.data().exterior.clone();

    diffraction(verts, ampl, prop, vk7);
}

fn diffraction(
    verts: Vec<Point3<f32>>,
    ampl: Matrix2<Complex<f32>>,
    prop: Vector3<f32>,
    vk7: Vector3<f32>,
) {
    println!("face vertices: {:?}", verts);
    println!(" ampl: {:?}", ampl);
    println!(" prop: {:?}", prop);
    println!(" vk7: {:?}", vk7);
    let num_verts = 4;

    // Sum all vertex coordinates
    let sum: Vector3<f32> = verts
        .iter()
        .fold(Vector3::zeros(), |acc, vert| acc + vert.coords);

    // Compute center of mass
    let center_of_mass = Point3::from(sum / num_verts as f32);

    println!("com: {}", center_of_mass);

    let v20: Vec<_> = verts
        .iter()
        .map(|point| point.coords - center_of_mass.coords)
        .collect();

    println!("v20: {:?}", v20);

    let rot = get_rotation_matrix2(v20);
    println!("rot: {}", rot);

    let prop1 = rot * prop;
    let perp1 = rot * vk7;
    println!("prop1: {}", prop1);
    println!("per1: {}", perp1);

    let angle = -prop1.y.atan2(prop1.x);
    println!("angle is: {}", angle);

    let rot2 = Matrix3::new(
        angle.cos(),
        -angle.sin(),
        0.0,
        angle.sin(),
        angle.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    );

    println!("second rotation matrix rot2 {}", rot2);

    let prop2 = rot2 * prop1;
    let perp2 = rot2 * perp1;

    println!("prop2: {}", prop2);
    println!("per2: {}", perp2);

    let e_par2 = perp2.cross(&prop2).normalize();
    println!("epar2 {}", e_par2);

    let anti_parallel = if e_par2.z > config::COLINEAR_THRESHOLD {
        true
    } else {
        false
    };
    println!("antiparallel: {}", anti_parallel);

    let ampl = if anti_parallel {
        println!("reversed ampl: {}", -ampl);
        -ampl
    } else {
        println!("non-reversed ampl: {}", ampl);
        ampl
    };

    let incidence = Vector3::new(0.0, 0.0, -1.0); // !todo: generalise

    let incidence2 = rot2 * rot * incidence;
    println!("incidence2: {}", incidence2);

    // make some far field bins
    let r = 1e6;
    let help1s: Vec<f32> = vec![0.1, 0.2]; // replace with theta later
    let help2s: Vec<f32> = vec![0.1, 0.4]; // replace with phi later
    let mut xfar: [[f32; 2]; 2] = [[0.0; 2]; 2];
    let mut yfar: [[f32; 2]; 2] = [[0.0; 2]; 2];
    let mut zfar: [[f32; 2]; 2] = [[0.0; 2]; 2];
    for (i, help1) in help1s.iter().enumerate() {
        for (j, help2) in help2s.iter().enumerate() {
            xfar[i][j] = r * help1.sin() * help2.cos();
            yfar[i][j] = r * help1.sin() * help2.sin();
            zfar[i][j] = r * help1.cos();
        }
    }
    println!("xfar: {:?}", xfar);
    println!("yfar: {:?}", yfar);
    println!("zfar: {:?}", zfar);

    // translate far-field bins to aperture system
    let mut x1 = xfar.clone();
    let mut y1 = yfar.clone();
    let mut z1 = zfar.clone();
    for row in x1.iter_mut() {
        for elem in row.iter_mut() {
            *elem -= center_of_mass.x;
        }
    }
    for row in y1.iter_mut() {
        for elem in row.iter_mut() {
            *elem -= center_of_mass.y;
        }
    }
    for row in z1.iter_mut() {
        for elem in row.iter_mut() {
            *elem -= center_of_mass.z;
        }
    }

    println!("x1: {:?}", x1);
    println!("y1: {:?}", y1);
    println!("z1: {:?}", z1);

    let mut r1: [[f32; 2]; 2] = [[0.0; 2]; 2];
    for (i, row) in x1.iter().enumerate() {
        for (j, elem) in row.iter().enumerate() {
            r1[i][j] = (x1[i][j].powi(2) + y1[i][j].powi(2) + z1[i][j].powi(2)).sqrt();
        }
    }
    println!("r1: {:?}", r1);

    todo!()
}

fn get_rotation_matrix2(verts: Vec<Vector3<f32>>) -> Matrix3<f32> {
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
