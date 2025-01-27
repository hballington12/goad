use std::{f32::consts::PI, fs::File, io::BufWriter};

use anyhow::Result;
use nalgebra::{Complex, ComplexField, Matrix2};
use ndarray::Array2;
use std::io::Write;

pub fn writeup(
    theta_phi_combinations: &[(f32, f32)],
    ampl_cs: &Vec<Matrix2<Complex<f32>>>,
) -> Result<()> {
    println!("done.");
    let mut mueller = Array2::<f32>::zeros((theta_phi_combinations.len(), 16));

    for (index, amplc) in ampl_cs.iter().enumerate() {
        mueller[[index, 0]] = (Complex::new(0.5, 0.0)
            * (amplc[(0, 0)] * amplc[(0, 0)].conj()
                + amplc[(0, 1)] * amplc[(0, 1)].conj()
                + amplc[(1, 0)] * amplc[(1, 0)].conj()
                + amplc[(1, 1)] * amplc[(1, 1)].conj()))
        .real();

        mueller[[index, 1]] = (Complex::new(0.5, 0.0)
            * (amplc[(0, 0)] * amplc[(0, 0)].conj() - amplc[(0, 1)] * amplc[(0, 1)].conj()
                + amplc[(1, 0)] * amplc[(1, 0)].conj()
                - amplc[(1, 1)] * amplc[(1, 1)].conj()))
        .real();
        // s22
        mueller[[index, 5]] = (Complex::new(0.5, 0.0)
            * (amplc[(0, 0)] * amplc[(0, 0)].conj()
                - amplc[(0, 1)] * amplc[(0, 1)].conj()
                - amplc[(1, 0)] * amplc[(1, 0)].conj()
                + amplc[(1, 1)] * amplc[(1, 1)].conj()))
        .real();
    }

    // Open a file for writing
    let file = File::create("output.txt").unwrap();
    let mut writer = BufWriter::new(file);

    // Write header

    // Iterate over the array and write data to the file
    for (index, row) in mueller.outer_iter().enumerate() {
        let (theta, phi) = theta_phi_combinations[index];
        writeln!(
            writer,
            "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            theta * 180.0 / PI, // Theta at index i
            phi * 180.0 / PI,   // Phi at index j
            row[0],             // s11 value
            row[1],
            0.0,
            0.0,
            0.0,
            row[5],
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
        )?;
    }
    Ok(())
}
