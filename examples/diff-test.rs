use goad::geom::{self};
use goad::{bins, diff, output};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Vector3};

fn main() {
    let geom = geom::Geom::from_file("./examples/data/hex.obj").unwrap();

    // pull rectangular face and print vertices
    let face = geom.shapes[0].faces[4].clone();
    println!("face vertices: {:?}", face.data().exterior);

    let m11: Complex<f32> = Complex::new(0.5, 0.25);
    let m12: Complex<f32> = Complex::new(0.25, -0.45);
    let m21: Complex<f32> = Complex::new(0.85, 0.2);
    let m22: Complex<f32> = Complex::new(-0.5, 0.5);
    let ampl = Matrix2::new(m11, m12, m21, m22);
    let prop: Vector3<f32> = Vector3::new(0.0, 0.5, -0.9).normalize();
    let vk7: Vector3<f32> = Vector3::new(1.0, 0.0, 0.0);
    let vk7 = vk7.cross(&prop).normalize();
    let verts = face.data().exterior.clone();
    let theta_phi_combinations = bins::simple_bins(180, 180);
    let ampl_far_field =
        diff::diffraction(&verts, ampl, prop, vk7, &theta_phi_combinations, 1.0, None);
    let mueller = output::ampl_to_mueller(&theta_phi_combinations, &ampl_far_field);
    let cwd = std::env::current_dir().expect("Failed to get current directory");
    let _ = output::write_mueller(&theta_phi_combinations, &mueller, "", &cwd);
}
