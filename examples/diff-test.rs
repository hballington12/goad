use anyhow::Result;
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
use pbt::{bins, diff, output};
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
    let face = geom.shapes[0].faces[0].clone();
    println!("face vertices: {:?}", face.data().exterior);

    let m11: Complex<f32> = Complex::new(0.5, 0.25);
    let m12: Complex<f32> = Complex::new(0.25, -0.45);
    let m21: Complex<f32> = Complex::new(0.85, 0.2);
    let m22: Complex<f32> = Complex::new(-0.5, 0.5);
    let ampl = Matrix2::new(m11, m12, m21, m22);
    let prop: Vector3<f32> = Vector3::new(0.5, -0.3, -0.2).normalize();
    // let prop: Vector3<f32> = Vector3::new(0.0, 0.0, 1.0).normalize();
    let vk7: Vector3<f32> = Vector3::new(0.0, 0.0, 1.0);
    let vk7 = vk7.cross(&prop).normalize();
    let verts = face.data().exterior.clone();

    let theta_phi_combinations = bins::generate_theta_phi_combinations();
    let ampls = diff::diffraction(&verts, ampl, prop, vk7, &theta_phi_combinations);
    let _ = output::writeup(&theta_phi_combinations, &ampls);
}
