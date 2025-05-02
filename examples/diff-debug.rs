use goad::geom::{self};
use goad::{bins, diff, output};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};

fn main() {
    let geom = geom::Geom::from_file("./examples/data/hex.obj").unwrap();

    // pull rectangular face and print vertices
    let face = geom.shapes[0].faces[4].clone();
    println!("face vertices: {:?}", face.data().exterior);

    let m11: Complex<f32> = Complex::new(-0.16685463, -0.013710017);
    let m12: Complex<f32> = Complex::new(0.958626, 0.078767836);
    let m21: Complex<f32> = Complex::new(-0.91943324, -0.07554753);
    let m22: Complex<f32> = Complex::new(-0.16977832, -0.0139502855);
    let ampl = Matrix2::new(m11, m12, m21, m22);
    let prop: Vector3<f32> = Vector3::new(-0.6380149, 0.3386723, 0.6915476).normalize();
    let vk7: Vector3<f32> = Vector3::new(0.38914767, 0.91677344, -0.08994864);
    let verts = vec![
        Point3::new(-17.16684, 3.791843, 5.771179),
        Point3::new(-15.896614, 6.7843, 5.7711763),
        Point3::new(-15.896614, 6.7843, 4.7025113),
        Point3::new(-15.896613, 6.7843, 3.7215571),
        Point3::new(-16.952244, 4.2973986, 5.4249063),
    ];

    let theta_phi_combinations = bins::simple_bins(100, 100);
    let ampl_far_field =
        diff::diffraction(&verts, ampl, prop, vk7, &theta_phi_combinations, 1.0, None);
    let mueller = output::ampl_to_mueller(&theta_phi_combinations, &ampl_far_field);
    let cwd = std::env::current_dir().expect("Failed to get current directory");
    let _ = output::write_mueller(&theta_phi_combinations, &mueller, "", &cwd);
}
