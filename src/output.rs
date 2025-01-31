use std::{f64::consts::PI, fs::File, io::BufWriter};

use anyhow::Result;
use nalgebra::{Complex, Matrix2};
use ndarray::Array2;
use std::io::Write;

pub fn writeup(bins: &[(f64, f64)], mueller: &Array2<f64>) -> Result<()> {
    // Open a file for writing
    let file = File::create("mueller_scatgrid").unwrap();
    let mut writer = BufWriter::new(file);

    // Iterate over the array and write data to the file
    for (index, row) in mueller.outer_iter().enumerate() {
        let (theta, phi) = bins[index];
        write!(writer, "{} {} ", theta * 180.0 / PI, phi * 180.0 / PI)?;
        for value in row.iter() {
            write!(writer, "{} ", value)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

pub fn ampl_to_mueller(
    theta_phi_combinations: &[(f64, f64)],
    ampl_cs: &[Matrix2<Complex<f64>>],
) -> Array2<f64> {
    let mut mueller = Array2::<f64>::zeros((theta_phi_combinations.len(), 16));

    for (index, amplc) in ampl_cs.iter().enumerate() {
        mueller[[index, 0]] = (Complex::new(0.5, 0.0)
            * (amplc[(0, 0)] * amplc[(0, 0)].conj()
                + amplc[(0, 1)] * amplc[(0, 1)].conj()
                + amplc[(1, 0)] * amplc[(1, 0)].conj()
                + amplc[(1, 1)] * amplc[(1, 1)].conj()))
        .re;

        mueller[[index, 1]] = (Complex::new(0.5, 0.0)
            * (amplc[(0, 0)] * amplc[(0, 0)].conj() - amplc[(0, 1)] * amplc[(0, 1)].conj()
                + amplc[(1, 0)] * amplc[(1, 0)].conj()
                - amplc[(1, 1)] * amplc[(1, 1)].conj()))
        .re;
        mueller[[index, 2]] =
            (amplc[(0, 0)] * amplc[(0, 1)].conj() + amplc[(1, 1)] * amplc[(1, 0)].conj()).re;
        mueller[[index, 3]] =
            (amplc[(0, 0)] * amplc[(0, 1)].conj() - amplc[(1, 1)] * amplc[(1, 0)].conj()).im;
        mueller[[index, 4]] = 0.5
            * (amplc[(0, 0)] * amplc[(0, 0)].conj() + amplc[(0, 1)] * amplc[(0, 1)].conj()
                - amplc[(1, 0)] * amplc[(1, 0)].conj()
                - amplc[(1, 1)] * amplc[(1, 1)].conj())
            .re;
        mueller[[index, 5]] = 0.5
            * (amplc[(0, 0)] * amplc[(0, 0)].conj()
                - amplc[(0, 1)] * amplc[(0, 1)].conj()
                - amplc[(1, 0)] * amplc[(1, 0)].conj()
                + amplc[(1, 1)] * amplc[(1, 1)].conj())
            .re;
        mueller[[index, 6]] =
            (amplc[(0, 0)] * amplc[(0, 1)].conj() - amplc[(1, 1)] * amplc[(1, 0)].conj()).re;
        mueller[[index, 7]] =
            (amplc[(0, 0)] * amplc[(0, 1)].conj() + amplc[(1, 1)] * amplc[(1, 0)].conj()).im;
        mueller[[index, 8]] =
            (amplc[(0, 0)] * amplc[(1, 0)].conj() + amplc[(1, 1)] * amplc[(0, 1)].conj()).re;
        mueller[[index, 9]] =
            (amplc[(0, 0)] * amplc[(1, 0)].conj() - amplc[(1, 1)] * amplc[(0, 1)].conj()).re;
        mueller[[index, 10]] =
            (amplc[(0, 0)] * amplc[(1, 1)].conj() + amplc[(0, 1)] * amplc[(1, 0)].conj()).re;
        mueller[[index, 11]] =
            (amplc[(0, 0)] * amplc[(1, 1)].conj() + amplc[(1, 0)] * amplc[(0, 1)].conj()).im;
        mueller[[index, 12]] =
            (amplc[(1, 0)] * amplc[(0, 0)].conj() + amplc[(1, 1)] * amplc[(0, 1)].conj()).im;
        mueller[[index, 13]] =
            (amplc[(1, 0)] * amplc[(0, 0)].conj() - amplc[(1, 1)] * amplc[(0, 1)].conj()).im;
        mueller[[index, 14]] =
            (amplc[(1, 1)] * amplc[(0, 0)].conj() - amplc[(0, 1)] * amplc[(1, 0)].conj()).im;
        mueller[[index, 15]] =
            (amplc[(1, 1)] * amplc[(0, 0)].conj() - amplc[(0, 1)] * amplc[(1, 0)].conj()).re;
    }
    mueller
}
