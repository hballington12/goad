use std::{fs::File, io::BufWriter};

use anyhow::Result;
use nalgebra::{Complex, Matrix2};
// use ndarray::array;
use ndarray::{s, Array1, Array2, Axis};
// use ndarray_stats::QuantileExt;
use anyhow::anyhow;
use itertools::Itertools;
use std::io::Write;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn test_integrate_trapezoidal() {
        // Test the integral of sin(x) from 0 to pi
        let x = Array1::linspace(0.0, PI, 1000);
        let y = x.mapv(f32::sin);
        let result = integrate_trapezoidal(&x, &y);
        assert!((result - 2.0).abs() < 1e-4, "result: {}", result);

        // Test the integral of xsin(x) from 0 to pi
        let x = Array1::linspace(0.0, PI, 1000);
        let y = &x * &x.mapv(f32::sin);
        let result = integrate_trapezoidal(&x, &y);
        assert!((result - PI).abs() < 1e-4, "result: {}", result);
    }

    #[test]
    fn test_unique() {
        let mut arr = vec![1.0, 1.0, 1.2, 2.0, 2.0, 3.0];
        arr.sort_by(|a, b| a.partial_cmp(b).expect("NaN encountered"));
        arr.dedup();
        let expected = vec![1.0, 1.2, 2.0, 3.0];
        assert_eq!(arr, expected);

        let mut arr = vec![1.0, 1.2, 1.0, 3.0, 2.01, 2.0];
        arr.sort_by(|a, b| a.partial_cmp(b).expect("NaN encountered"));
        arr.dedup();
        let expected = vec![1.0, 1.2, 2.0, 2.01, 3.0];
        assert_eq!(arr, expected);
    }

    #[test]
    fn test_unique_grid() {
        let input = vec![(1.0, 2.0), (1.0, 2.0), (1.0, 3.0), (2.0, 3.0), (2.0, 3.0)];
        let expected = vec![(1.0, 2.0), (1.0, 3.0), (2.0, 2.0), (2.0, 3.0)];
        let result = unique_grid(&input);
        assert_eq!(result, expected);

        let input = vec![(1.0, 1.0), (1.0, 1.0), (2.0, 2.0), (2.0, 2.0)];
        let expected = vec![(1.0, 1.0), (1.0, 2.0), (2.0, 1.0), (2.0, 2.0)];
        let result = unique_grid(&input);
        assert_eq!(result, expected);
    }
}

fn integrate_trapezoidal(x: &Array1<f32>, y: &Array1<f32>) -> f32 {
    let dx = &x.slice(s![1..]) - &x.slice(s![..-1]);
    let avg_y = (&y.slice(s![1..]) + &y.slice(s![..-1])) / 2.0;
    (dx * avg_y).sum()
}

#[allow(dead_code)]
/// Given a slice of 2-element tuples, return a Vec of unique tuples.
fn unique_grid(tuple_slice: &[(f32, f32)]) -> Vec<(f32, f32)> {
    // separate the grid into two arrays
    let (mut thetas, mut phis): (Vec<f32>, Vec<f32>) = tuple_slice.iter().cloned().unzip();
    // sort each array
    thetas.sort_by(|a, b| a.partial_cmp(b).expect("NaN encountered"));
    phis.sort_by(|a, b| a.partial_cmp(b).expect("NaN encountered"));
    // deduplicate each array
    thetas.dedup();
    phis.dedup();
    // create a new Vec of tuples
    let mut unique_grid = Vec::new();
    for theta in thetas.iter() {
        for phi in phis.iter() {
            unique_grid.push((*theta, *phi));
        }
    }
    unique_grid
}

/// Write the Mueller matrix to a file against the theta and phi bins
pub fn writeup(bins: &[(f32, f32)], mueller: &Array2<f32>) -> Result<()> {
    // Open a file for writing
    let file = File::create("mueller_scatgrid").unwrap();
    let mut writer = BufWriter::new(file);

    // Iterate over the array and write data to the file
    for (index, row) in mueller.outer_iter().enumerate() {
        let (theta, phi) = bins[index];
        write!(writer, "{} {} ", theta, phi)?;
        for value in row.iter() {
            write!(writer, "{} ", value)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

pub fn ampl_to_mueller(
    theta_phi_combinations: &[(f32, f32)],
    ampl_cs: &[Matrix2<Complex<f32>>],
) -> Array2<f32> {
    let mut mueller = Array2::<f32>::zeros((theta_phi_combinations.len(), 16));

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

/// Integrate over phi (second bin of the tuple) to get the 1D Mueller matrix
/// Uses the trapezoidal rule
/// Returns a tuple of the theta bins and the 1D Mueller matrix
/// NOTE: Assumes phi is ordered
pub fn try_mueller_to_1d(
    bins: &[(f32, f32)],
    mueller: &Array2<f32>,
) -> Result<(Vec<f32>, Array2<f32>)> {
    // Check that the mueller matrix and bins are the same length
    if mueller.len_of(Axis(0)) != bins.len() {
        return Err(anyhow!(
            "Mueller matrix and bins must have the same length. Got {} and {}",
            mueller.len_of(Axis(0)),
            bins.len()
        ));
    }

    // Sort the bins and mueller by theta
    let mut bins = bins.to_owned();
    bins.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut mueller = mueller.to_owned();

    // Create a new sorted mueller matrix
    let mut sorted_mueller = mueller.clone();
    for (i, bin) in bins.iter().enumerate() {
        let index = bins.iter().position(|b| b == bin).unwrap();
        sorted_mueller
            .slice_mut(s![i, ..])
            .assign(&mueller.slice(s![index, ..]));
    }

    // zip the bins and mueller matrix
    // let mueller_vec = mueller
    //     .outer_iter()
    //     .map(|row| row.to_vec())
    //     .collect::<Vec<_>>();

    // zip the bins and mueller matrix
    let combined: Vec<_> = bins
        .iter()
        .zip(mueller.outer_iter())
        .map(|(bin, row)| (*bin, row.to_owned()))
        .collect();

    // group the combined Vec by theta
    let grouped: Vec<Vec<_>> = combined
        .into_iter()
        .chunk_by(|((key, _), _)| *key)
        .into_iter()
        .map(|(_, group)| group.map(|x| x).collect())
        .collect();

    println!("{:?}", grouped);

    // get number of thetas
    let num_thetas = grouped.len();

    // integrate over phi
    let mut thetas = Vec::new();
    let mut mueller_1d = Array2::<f32>::zeros((num_thetas, 16));
    // loop over vectors at each theta

    for (i, muellers) in grouped.iter().enumerate() {
        let thetas_phi: Vec<_> = muellers
            .iter()
            .map(|((theta, phi), _)| (*theta, *phi))
            .collect();
        let mueller_phi: Vec<_> = muellers.iter().map(|(_, mueller)| mueller).collect();
        let mut mueller_1d_row = Array1::<f32>::zeros(16);
        for j in 0..16 {
            let y = Array1::from(mueller_phi.iter().map(|row| row[j]).collect::<Vec<_>>());
            let x = Array1::from(thetas_phi.iter().map(|(_, phi)| *phi).collect::<Vec<_>>());
            mueller_1d_row[j] = integrate_trapezoidal(&x, &y); // integrate over phi
        }
        thetas.push(thetas_phi[0].0);
        mueller_1d
            .slice_mut(s![i, ..])
            .assign(&mueller_1d_row.slice(s![..]));
    }

    println!("");
    println!("{:?}", mueller_1d);

    !todo!()
}

fn split_at_values<T: PartialEq + Clone>(vec: &[T]) -> Vec<Vec<T>> {
    vec.iter()
        .chunk_by(|&x| x)
        .into_iter()
        .map(|(_, group)| group.cloned().collect())
        .collect()
}
