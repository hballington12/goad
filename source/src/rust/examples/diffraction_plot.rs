//! Example: Plot diffraction amplitude vs theta angle
//!
//! Run with: cargo run --example diffraction_plot
//! Or with non-normal incidence: cargo run --example diffraction_plot -- --oblique

use goad::beam::{Beam, BeamVariant, DefaultBeamVariant};
use goad::bins::{AngleBin, SolidAngleBin};
use goad::diff;
use goad::diff2;
use goad::field::Field;
use goad::geom::Face;
use nalgebra::{Complex, Point3, Vector3};
use plotters::backend::BitMapBackend;
use plotters::prelude::*;
use std::f32::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let oblique = std::env::args().any(|a| a == "--oblique");

    // Create a triangular beam
    // Vertices are clockwise when viewed along prop direction
    let verts = vec![
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(-0.5, 0.866, 0.0),
        Point3::new(-0.5, -0.866, 0.0),
    ];

    let face = Face::new_simple(verts, None, None)?;

    let (prop, e_perp, output_file) = if oblique {
        // Non-normal incidence: prop at 30 degrees from -z in xz plane
        let theta_inc = 30.0_f32.to_radians();
        let prop = Vector3::new(theta_inc.sin(), 0.0, -theta_inc.cos());
        let e_perp = Vector3::new(0.0, 1.0, 0.0);
        println!("Using OBLIQUE incidence: prop = {:?}", prop);
        (prop, e_perp, "diffraction_oblique.png")
    } else {
        // Normal incidence: prop along -z
        let prop = Vector3::new(0.0, 0.0, -1.0);
        let e_perp = Vector3::new(0.0, 1.0, 0.0);
        println!("Using NORMAL incidence: prop = {:?}", prop);
        (prop, e_perp, "diffraction_amplitude.png")
    };

    let field = Field::new_identity(e_perp, prop)?;

    let beam = Beam::new(
        face,
        Complex::new(1.0, 0.0),
        0,
        0,
        field,
        BeamVariant::Default(DefaultBeamVariant::Refr),
        0.532,
    );

    // Create bins spanning theta from 0.1 to 30 degrees, phi = 0
    let num_theta = 100;
    let theta_min = 0.1_f32;
    let theta_max = 30.0_f32;
    let d_theta = (theta_max - theta_min) / num_theta as f32;

    let bins: Vec<SolidAngleBin> = (0..num_theta)
        .map(|i| {
            let theta_center = theta_min + (i as f32 + 0.5) * d_theta;
            SolidAngleBin::new(
                AngleBin::from_center_width(theta_center, d_theta),
                AngleBin::from_center_width(0.0, 1.0),
            )
        })
        .collect();

    // Compute diffraction using diff2 (new implementation)
    println!("Computing diffraction for {} theta bins...", bins.len());
    let incident = diff2::IncidentBeam::default();
    let results_new = diff2::n2f_aperture_diffraction(&beam, &bins, &incident, None)?;

    // Extract amplitude magnitudes from new implementation (use S11 element)
    let mut data_new: Vec<(f32, f32)> = Vec::new();
    for (idx, ampl) in &results_new {
        let theta = bins[*idx].theta.center;
        let magnitude = ampl[(0, 0)].norm(); // S11 element
        data_new.push((theta, magnitude));
    }
    data_new.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    println!("diff2: Got {} non-zero results", data_new.len());

    // Compute diffraction using original diff implementation
    let wavelength = 0.532_f32;
    let wavenumber = 2.0 * PI / wavelength;
    let verts = beam.face.data().exterior.clone();
    let ampl_matrix = beam.field.ampl();
    let prop = beam.field.prop();
    let e_perp = beam.field.e_perp();

    let results_old =
        diff::n2f_aperture_diffraction(&verts, ampl_matrix, prop, e_perp, &bins, wavenumber, None);

    // Extract amplitude magnitudes from old implementation
    let mut data_old: Vec<(f32, f32)> = Vec::new();
    for (idx, ampl) in results_old.iter().enumerate() {
        if ampl.iter().any(|c| c.norm() > 0.0) {
            let theta = bins[idx].theta.center;
            let magnitude = ampl[(0, 0)].norm(); // S11 element
            data_old.push((theta, magnitude));
        }
    }
    data_old.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    println!("diff (original): Got {} non-zero results", data_old.len());

    // Find max for scaling
    let max_new = data_new.iter().map(|(_, m)| *m).fold(0.0_f32, f32::max);
    let max_old = data_old.iter().map(|(_, m)| *m).fold(0.0_f32, f32::max);
    let max_mag = max_new.max(max_old);
    println!("Max amplitude (new): {:.6e}", max_new);
    println!("Max amplitude (old): {:.6e}", max_old);

    // Plot
    let root = BitMapBackend::new(output_file, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Diffraction Amplitude vs Theta", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(theta_min..theta_max, 0.0_f32..max_mag * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Theta (degrees)")
        .y_desc("|S11| Amplitude")
        .draw()?;

    // Plot new implementation (diff2) in blue
    chart
        .draw_series(LineSeries::new(data_new.iter().cloned(), &BLUE))?
        .label("diff2 (new)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Plot old implementation (diff) in red
    chart
        .draw_series(LineSeries::new(data_old.iter().cloned(), &RED))?
        .label("diff (original)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Draw legend
    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    println!("Plot saved to {}", output_file);

    Ok(())
}
