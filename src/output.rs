use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::{fs::File, io::BufWriter};

use anyhow::Result;
use serde::Serialize;

use crate::bins::SolidAngleBin;
use crate::result::{Mueller, MuellerMatrix, Results};
use crate::settings::{OutputConfig, Settings};
use crate::zones::Zone;

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
            Box::new(SettingsJsonWriter::new(self.settings)),
            Box::new(ConsolidatedResultsWriter::new(self.results)),
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

        // Write Mueller matrices for each zone in its own directory
        for zone in self.results.zones.iter() {
            let zone_dir = output_dir.join(zone.display_name());
            fs::create_dir_all(&zone_dir)?;

            self.write_zone_mueller(zone, &zone_dir, config)?;
        }

        Ok(())
    }

    fn write_zone_mueller(
        &self,
        zone: &Zone,
        zone_dir: &Path,
        config: &crate::settings::MuellerComponentConfig,
    ) -> Result<()> {
        // Write 2D Mueller matrices
        if self.settings.output.mueller_2d && !zone.field_2d.is_empty() {
            let bins: Vec<_> = zone.bins.clone();

            if config.total {
                let muellers: Vec<_> = zone.field_2d.iter().map(|f| f.mueller_total).collect();
                write_mueller(&bins, &muellers, "", zone_dir)?;
            }
            if config.beam {
                let muellers: Vec<_> = zone.field_2d.iter().map(|f| f.mueller_beam).collect();
                write_mueller(&bins, &muellers, "_beam", zone_dir)?;
            }
            if config.external {
                let muellers: Vec<_> = zone.field_2d.iter().map(|f| f.mueller_ext).collect();
                write_mueller(&bins, &muellers, "_ext", zone_dir)?;
            }
        }

        // Write 1D Mueller matrices
        if self.settings.output.mueller_1d {
            if let Some(field_1d) = &zone.field_1d {
                if config.total {
                    write_mueller_1d(
                        "",
                        field_1d,
                        &|r: &crate::result::ScattResult1D| r.mueller_total.clone(),
                        zone_dir,
                    )?;
                }
                if config.beam {
                    write_mueller_1d(
                        "_beam",
                        field_1d,
                        &|r: &crate::result::ScattResult1D| r.mueller_beam.clone(),
                        zone_dir,
                    )?;
                }
                if config.external {
                    write_mueller_1d(
                        "_ext",
                        field_1d,
                        &|r: &crate::result::ScattResult1D| r.mueller_ext.clone(),
                        zone_dir,
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
        write!(writer, "{} {} ", bin.theta.center, bin.phi.center)?;
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

// ========================================
// Consolidated Results JSON Output
// ========================================

/// Serializable zone summary for results.json
#[derive(Serialize)]
struct ZoneOutput {
    label: String,
    zone_type: String,
    num_bins: usize,
    params: crate::params::Params,
    mueller_dir: String,
}

/// Serializable consolidated results
#[derive(Serialize)]
struct ConsolidatedResults {
    powers: crate::powers::Powers,
    zones: Vec<ZoneOutput>,
}

/// Writer for consolidated results.json file
pub struct ConsolidatedResultsWriter<'a> {
    results: &'a Results,
}

impl<'a> ConsolidatedResultsWriter<'a> {
    pub fn new(results: &'a Results) -> Self {
        Self { results }
    }
}

impl<'a> OutputWriter for ConsolidatedResultsWriter<'a> {
    fn write(&self, output_dir: &Path) -> Result<()> {
        let zones: Vec<ZoneOutput> = self
            .results
            .zones
            .iter()
            .map(|zone| {
                let label = zone.display_name();
                ZoneOutput {
                    zone_type: format!("{:?}", zone.zone_type),
                    num_bins: zone.bins.len(),
                    params: zone.params.clone(),
                    mueller_dir: label.clone(),
                    label,
                }
            })
            .collect();

        let consolidated = ConsolidatedResults {
            powers: self.results.powers.clone(),
            zones,
        };

        let path = output_path(Some(output_dir), &self.filename())?;
        let json = serde_json::to_string_pretty(&consolidated)?;
        fs::write(path, json)?;
        Ok(())
    }

    fn filename(&self) -> String {
        "results.json".to_string()
    }

    fn is_enabled(&self, _config: &OutputConfig) -> bool {
        // Always enabled for now - could add config option later
        true
    }
}
