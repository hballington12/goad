use crate::result::MuellerMatrix;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::{fs::File, io::BufWriter};

use anyhow::Result;
use ndarray::{s, Array1};

use crate::bins::SolidAngleBin;
use crate::result::{Mueller, Results};

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn test_integrate_trapezoidal() {
        // Test the integral of sin(x) from 0 to pi
        let x = Array1::linspace(0.0, PI, 1000);
        let y = x.mapv(f32::sin);
        let result = integrate_trapezoidal(&x, &y, |_, y| y);
        assert!((result - 2.0).abs() < 1e-4, "result: {}", result);

        // Test the integral of xsin(x) from 0 to pi
        let x = Array1::linspace(0.0, PI, 1000);
        let y = &x * &x.mapv(f32::sin);
        let result = integrate_trapezoidal(&x, &y, |_, y| y);
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
}

pub fn integrate_trapezoidal<F>(x: &Array1<f32>, y: &Array1<f32>, transform: F) -> f32
where
    F: Fn(f32, f32) -> f32,
{
    let dx = &x.slice(s![1..]) - &x.slice(s![..-1]);
    let xl = x.slice(s![..-1]);
    let xr = x.slice(s![1..]);
    let yl = y.slice(s![..-1]);
    let yr = y.slice(s![1..]);

    let mut sum = 0.0;
    for i in 0..dx.len() {
        let left = transform(xl[i], yl[i]);
        let right = transform(xr[i], yr[i]);
        sum += dx[i] * (left + right) / 2.0;
    }

    sum
}

// /// Write the 1d Mueller matrix to a file against the theta and phi bins
// pub fn write_mueller_1d(
//     bins: &[f32],
//     mueller_1d: &[Mueller],
//     suffix: &str,
//     output_dir: &Path,
// ) -> Result<()> {
//     let file_name = format!("mueller_scatgrid_1d{}", suffix);
//     let path = output_path(Some(output_dir), &file_name)?;

//     let file = File::create(&path)?;
//     let mut writer = BufWriter::new(file);

//     // Iterate over the array and write data to the file
//     for (index, mueller) in mueller_1d.iter().enumerate() {
//         let theta = bins[index];
//         write!(writer, "{} ", theta)?;
//         for element in mueller.to_vec().into_iter() {
//             write!(writer, "{} ", element)?;
//         }
//         writeln!(writer)?;
//     }

//     Ok(())
// }

// /// Write the Mueller matrix to a file against the theta and phi bins
// pub fn write_mueller(
//     bins: &[SolidAngleBin],
//     muellers: &[&Mueller],
//     suffix: &str,
//     output_dir: &Path,
// ) -> Result<()> {
//     let file_name_total = "mueller_scatgrid";
//     let path_total = output_path(Some(output_dir), &file_name_total)?;
//     let file_total = File::create(&path_total)?;
//     let mut writer_total = BufWriter::new(file_total);

//     // Iterate over the array and write data to the file
//     for (index, mueller) in muellers.iter().enumerate() {
//         let bin = bins[index];
//         write!(writer, "{} {} ", bin.theta_bin.center, bin.phi_bin.center)?;
//         for element in mueller.to_vec().into_iter() {
//             write!(writer, "{} ", element)?;
//         }
//         writeln!(writer)?;
//     }

//     Ok(())
// }

/// Write the Mueller matrix to a file against the theta and phi bins
pub fn write_result(result: &Results, output_dir: &Path) -> Result<()> {
    let file_name = format!("results.dat");
    let path = output_path(Some(output_dir), &file_name)?;

    let file = File::create(&path)?;
    let mut writer = BufWriter::new(file);

    // Write the results to a file
    writeln!(writer, "# GOAD Simulation Results")?;
    writeln!(writer, "# ======================")?;

    // Write parameters section
    writeln!(writer, "\n# Parameters")?;
    writeln!(writer, "# ----------")?;
    if let Some(scat) = result.params.scat_cross {
        writeln!(writer, "Scattering Cross Section: {:.6}", scat)?;
    }
    if let Some(ext) = result.params.ext_cross {
        writeln!(writer, "Extinction Cross Section: {:.6}", ext)?;
    }
    if let Some(albedo) = result.params.albedo {
        writeln!(writer, "Single Scattering Albedo: {:.6}", albedo)?;
    }
    if let Some(asym) = result.params.asymettry {
        writeln!(writer, "Asymmetry Parameter: {:.6}", asym)?;
    }

    // Write powers section
    writeln!(writer, "\n# Power Distribution")?;
    writeln!(writer, "# ----------------")?;
    writeln!(writer, "Input Power:           {:.6}", result.powers.input)?;
    writeln!(writer, "Output Power:          {:.6}", result.powers.output)?;
    writeln!(
        writer,
        "Absorbed Power:        {:.6}",
        result.powers.absorbed
    )?;
    writeln!(
        writer,
        "Truncated Reflections: {:.6}",
        result.powers.trnc_ref
    )?;
    writeln!(
        writer,
        "Truncated Recursions:  {:.6}",
        result.powers.trnc_rec
    )?;
    writeln!(
        writer,
        "Truncated Clip Error:  {:.6}",
        result.powers.clip_err
    )?;
    writeln!(
        writer,
        "Truncated Energy:      {:.6}",
        result.powers.trnc_energy
    )?;
    writeln!(
        writer,
        "Truncated Area:        {:.6}",
        result.powers.trnc_area
    )?;
    writeln!(
        writer,
        "Truncated Cutoff:      {:.6}",
        result.powers.trnc_cop
    )?;
    writeln!(
        writer,
        "External Diffraction:  {:.6}",
        result.powers.ext_diff
    )?;
    writeln!(
        writer,
        "Missing Power:         {:.6}",
        result.powers.missing()
    )?;

    // Write ratios
    writeln!(writer, "\n# Power Ratios")?;
    writeln!(writer, "# ------------")?;
    let output_ratio = result.powers.output / result.powers.input;
    let absorbed_ratio = result.powers.absorbed / result.powers.input;
    let total_ratio = (result.powers.output + result.powers.absorbed) / result.powers.input;
    writeln!(writer, "Scattered/Input Ratio: {:.6}", output_ratio)?;
    writeln!(writer, "Absorbed/Input Ratio:  {:.6}", absorbed_ratio)?;
    writeln!(writer, "Total/Input Ratio:     {:.6}", total_ratio)?;

    // Write binning information
    writeln!(writer, "\n# Simulation Information")?;
    writeln!(writer, "# ---------------------")?;
    writeln!(writer, "Number of bins: {}", result.bins().len())?;
    if let Some(bins_1d) = &result.bins_1d {
        writeln!(writer, "Number of 1D bins: {}", bins_1d.len())?;
    }

    println!("Results written to: {}", path.display());
    Ok(())
}

// Helper function to construct the output path and ensure the directory exists
pub fn output_path(output_dir: Option<&Path>, file_name: &str) -> Result<PathBuf> {
    match output_dir {
        Some(dir) => {
            fs::create_dir_all(dir)?;
            Ok(dir.join(file_name))
        }
        None => Ok(PathBuf::from(file_name)),
    }
}
