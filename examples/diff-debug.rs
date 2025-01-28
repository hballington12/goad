use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};
use pbt::geom::{self};
use pbt::{bins, diff, output};

fn main() {
    let geom = geom::Geom::from_file("./examples/data/hex.obj").unwrap();

    // pull rectangular face and print vertices
    let face = geom.shapes[0].faces[4].clone();
    println!("face vertices: {:?}", face.data().exterior);

    let m11: Complex<f32> = Complex::new(0.02172605, 0.010047015);
    let m12: Complex<f32> = Complex::new(-0.022072379, -0.010206628);
    let m21: Complex<f32> = Complex::new(-0.009356209, -0.004326305);
    let m22: Complex<f32> = Complex::new(0.052926008, 0.024474261);
    let ampl = Matrix2::new(m11, m12, m21, m22);
    let prop: Vector3<f32> = Vector3::new(0.33373073, 0.7969924, -0.5034153).normalize();
    let vk7: Vector3<f32> = Vector3::new(0.8347903, -0.0017971284, 0.5505651);
    let verts = vec![
        Point3::new(-1.0991076, 2.5884433, -6.0657525),
        Point3::new(-1.1482885, 3.0610435, -5.934424),
        Point3::new(-1.1696503, 3.266319, -5.8773804),
        Point3::new(-0.6184374, 2.8900104, -5.64386),
    ];
    let theta_phi_combinations = bins::generate_theta_phi_combinations();
    let ampl_far_field = diff::diffraction(&verts, ampl, prop, vk7, &theta_phi_combinations);
    let _ = output::writeup(&theta_phi_combinations, &ampl_far_field);
}
