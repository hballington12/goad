use crate::result::MuellerMatrix;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::{fs::File, io::BufWriter};

use anyhow::Result;

use crate::bins::SolidAngleBin;
use crate::result::{Mueller, Results};
use crate::settings::{OutputConfig, Settings};

/// Trait for writing output data to files
pub trait OutputWriter {
    /// Write data to a file in the output directory
    fn write(&self, output_dir: &Path) -> Result<()>;

    /// Get the filename this writer uses
    fn filename(&self) -> String;

    /// Check if this output is enabled in the configuration
    fn is_enabled(&self, config: &OutputConfig) -> bool;
}

/// Manager for coordinating all output operations
pub struct OutputManager<'a> {
    pub settings: &'a Settings,
    pub results: &'a Results,
}

impl<'a> OutputManager<'a> {
    pub fn new(settings: &'a Settings, results: &'a Results) -> Self {
        Self { settings, results }
    }

    /// Write all enabled outputs based on configuration
    pub fn write_all(&self) -> Result<()> {
        let output_dir = &self.settings.directory;
        fs::create_dir_all(output_dir)?;

        // Create all possible output writers
        let writers: Vec<Box<dyn OutputWriter>> = vec![
            Box::new(ResultsSummaryWriter::new(self.results)),
            Box::new(SettingsJsonWriter::new(self.settings)),
            Box::new(PowersJsonWriter::new(&self.results.powers)),
            Box::new(ParamsJsonWriter::new(&self.results.params)),
        ];

        // Write enabled outputs
        for writer in writers {
            if writer.is_enabled(&self.settings.output) {
                writer.write(output_dir)?;
            }
        }

        // Handle Mueller matrix outputs separately (they have custom logic)
        self.write_mueller_matrices()?;

        Ok(())
    }

    fn write_mueller_matrices(&self) -> Result<()> {
        if !self.settings.output.mueller_2d && !self.settings.output.mueller_1d {
            return Ok(());
        }

        let output_dir = &self.settings.directory;
        let config = &self.settings.output.mueller_components;

        // Write 2D Mueller matrices
        if self.settings.output.mueller_2d {
            if config.total {
                let muellers = &self
                    .results
                    .field_2d
                    .iter()
                    .map(|f| f.mueller_total)
                    .collect::<Vec<_>>();
                write_mueller(&self.results.bins(), muellers, "", output_dir)?;
            }
            if config.beam {
                let muellers = &self
                    .results
                    .field_2d
                    .iter()
                    .map(|f| f.mueller_beam)
                    .collect::<Vec<_>>();
                write_mueller(&self.results.bins(), muellers, "_beam", output_dir)?;
            }
            if config.external {
                let muellers = &self
                    .results
                    .field_2d
                    .iter()
                    .map(|f| f.mueller_ext)
                    .collect::<Vec<_>>();
                write_mueller(&self.results.bins(), muellers, "_ext", output_dir)?;
            }
        }

        // Write 1D Mueller matrices
        if self.settings.output.mueller_1d {
            if let Some(field_1d) = &self.results.field_1d {
                if config.total {
                    write_mueller_1d(
                        "",
                        field_1d,
                        &|r: &crate::result::ScattResult1D| r.mueller_total.clone(),
                        output_dir,
                    )?;
                }
                if config.beam {
                    write_mueller_1d(
                        "_beam",
                        field_1d,
                        &|r: &crate::result::ScattResult1D| r.mueller_beam.clone(),
                        output_dir,
                    )?;
                }
                if config.external {
                    write_mueller_1d(
                        "_ext",
                        field_1d,
                        &|r: &crate::result::ScattResult1D| r.mueller_ext.clone(),
                        output_dir,
                    )?;
                }
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {

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

/// Write the Mueller matrix to a file against the theta and phi bins
pub fn write_mueller(
    bins: &[SolidAngleBin],
    muellers: &[Mueller],
    suffix: &str,
    output_dir: &Path,
) -> Result<()> {
    let file_name_total = format!("mueller_scatgrid{}", suffix);
    let path_total = output_path(Some(output_dir), &file_name_total)?;
    let file_total = File::create(&path_total)?;
    let mut writer = BufWriter::new(file_total);

    // Iterate over the array and write data to the file
    for (index, mueller) in muellers.iter().enumerate() {
        let bin = bins[index];
        write!(writer, "{} {} ", bin.theta_bin.center, bin.phi_bin.center)?;
        for element in mueller.to_vec().into_iter() {
            write!(writer, "{} ", element)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Write the 1D Mueller matrix to a file
pub fn write_mueller_1d<F>(
    suffix: &str,
    field_1d: &[crate::result::ScattResult1D],
    mueller_getter: F,
    output_dir: &Path,
) -> Result<()>
where
    F: Fn(&crate::result::ScattResult1D) -> Mueller,
{
    let file_name = format!("mueller_scatgrid_1d{}", suffix);
    let path = output_path(Some(output_dir), &file_name)?;
    let file = File::create(&path)?;
    let mut writer = BufWriter::new(file);

    for result in field_1d {
        let bin = result.bin;
        write!(writer, "{} ", bin.center)?;

        let mueller = mueller_getter(result);
        for element in mueller.to_vec() {
            write!(writer, "{} ", element)?;
        }

        writeln!(writer)?;
    }

    Ok(())
}

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

    // Write parameters for Total component (backwards compatible)
    if let Some(scat) = result.params.scat_cross() {
        writeln!(writer, "Scattering Cross Section: {:.6}", scat)?;
    }
    if let Some(ext) = result.params.ext_cross() {
        writeln!(writer, "Extinction Cross Section: {:.6}", ext)?;
    }
    if let Some(albedo) = result.params.albedo() {
        writeln!(writer, "Single Scattering Albedo: {:.6}", albedo)?;
    }
    if let Some(asym) = result.params.asymmetry() {
        writeln!(writer, "Asymmetry Parameter: {:.6}", asym)?;
    }

    // Write component-specific parameters
    writeln!(writer, "\n# Component-Specific Parameters")?;
    writeln!(writer, "# ------------------------------")?;

    for component in [
        crate::result::GOComponent::Beam,
        crate::result::GOComponent::ExtDiff,
    ] {
        let comp_str = match component {
            crate::result::GOComponent::Beam => "Beam",
            crate::result::GOComponent::ExtDiff => "ExtDiff",
            _ => continue,
        };

        if let Some(scat) = result.params.scat_cross.get(&component) {
            writeln!(writer, "{} Scattering Cross Section: {:.6}", comp_str, scat)?;
        }
        if let Some(asym) = result.params.asymmetry.get(&component) {
            writeln!(writer, "{} Asymmetry Parameter: {:.6}", comp_str, asym)?;
        }
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
    // if let Some(bins_1d) = &result.bins_1d {
    //     writeln!(writer, "Number of 1D bins: {}", bins_1d.len())?;
    // }

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

// ========================================
// Individual Output Writer Implementations
// ========================================

/// Writer for the results.dat summary file (existing implementation)
pub struct ResultsSummaryWriter<'a> {
    results: &'a Results,
}

impl<'a> ResultsSummaryWriter<'a> {
    pub fn new(results: &'a Results) -> Self {
        Self { results }
    }
}

impl<'a> OutputWriter for ResultsSummaryWriter<'a> {
    fn write(&self, output_dir: &Path) -> Result<()> {
        write_result(self.results, output_dir)
    }

    fn filename(&self) -> String {
        "results.dat".to_string()
    }

    fn is_enabled(&self, config: &OutputConfig) -> bool {
        config.results_summary
    }
}

/// Writer for settings.json file
pub struct SettingsJsonWriter<'a> {
    settings: &'a Settings,
}

impl<'a> SettingsJsonWriter<'a> {
    pub fn new(settings: &'a Settings) -> Self {
        Self { settings }
    }
}

impl<'a> OutputWriter for SettingsJsonWriter<'a> {
    fn write(&self, output_dir: &Path) -> Result<()> {
        let path = output_path(Some(output_dir), &self.filename())?;
        let json = serde_json::to_string_pretty(self.settings)?;
        fs::write(path, json)?;
        Ok(())
    }

    fn filename(&self) -> String {
        "settings.json".to_string()
    }

    fn is_enabled(&self, config: &OutputConfig) -> bool {
        config.settings_json
    }
}

/// Writer for powers.json file
pub struct PowersJsonWriter<'a> {
    powers: &'a crate::powers::Powers,
}

impl<'a> PowersJsonWriter<'a> {
    pub fn new(powers: &'a crate::powers::Powers) -> Self {
        Self { powers }
    }
}

impl<'a> OutputWriter for PowersJsonWriter<'a> {
    fn write(&self, output_dir: &Path) -> Result<()> {
        let path = output_path(Some(output_dir), &self.filename())?;
        let json = serde_json::to_string_pretty(self.powers)?;
        fs::write(path, json)?;
        Ok(())
    }

    fn filename(&self) -> String {
        "powers.json".to_string()
    }

    fn is_enabled(&self, config: &OutputConfig) -> bool {
        config.powers_json
    }
}

/// Writer for params.json file
pub struct ParamsJsonWriter<'a> {
    params: &'a crate::params::Params,
}

impl<'a> ParamsJsonWriter<'a> {
    pub fn new(params: &'a crate::params::Params) -> Self {
        Self { params }
    }
}

impl<'a> OutputWriter for ParamsJsonWriter<'a> {
    fn write(&self, output_dir: &Path) -> Result<()> {
        let path = output_path(Some(output_dir), &self.filename())?;
        let json = serde_json::to_string_pretty(self.params)?;
        fs::write(path, json)?;
        Ok(())
    }

    fn filename(&self) -> String {
        "params.json".to_string()
    }

    fn is_enabled(&self, config: &OutputConfig) -> bool {
        config.params_json
    }
}
