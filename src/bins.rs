use ndarray::Array1;
use pyo3::prelude::*;
use serde::Deserialize;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_bins() {
        let values = vec![0.0, 1.0, 2.0];
        let spacings = vec![0.5, 0.5];
        let result = interval_spacings(&values, &spacings);
        let expected = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        assert_eq!(result, expected);
    }

    #[test]
    #[should_panic]
    fn test_interval_bins_bad_angle() {
        let values = vec![0.0, 1.0, 2.0];
        let spacings = vec![0.3, 0.5];
        interval_spacings(&values, &spacings);
    }

    #[test]
    fn test_simple_bins() {
        let num_theta = 3;
        let num_phi = 3;
        let result = simple_bins(num_theta, num_phi);
        let expected = vec![
            (0.0, 0.0),
            (0.0, 180.0),
            (0.0, 360.0),
            (90.0, 0.0),
            (90.0, 180.0),
            (90.0, 360.0),
            (180.0, 0.0),
            (180.0, 180.0),
            (180.0, 360.0),
        ];
        assert_eq!(result, expected);
    }
}

/// Angular binning schemes for far-field scattering calculations.
/// 
/// **Context**: Far-field scattering patterns require discretization of the
/// observation sphere into angular bins for numerical computation. Different
/// applications require different angular resolutions - uniform grids for
/// general analysis, fine spacing near forward/backward directions for
/// specific features, or custom patterns for experimental comparisons.
/// 
/// **How it Works**: Provides three discretization approaches: Simple creates
/// uniform grids in theta and phi, Interval allows variable spacing with
/// different resolutions in different angular regions, and Custom enables
/// arbitrary user-defined angular points.
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum Scheme {
    Simple {
        num_theta: usize,
        num_phi: usize,
    },
    Interval {
        thetas: Vec<f32>,
        theta_spacings: Vec<f32>,
        phis: Vec<f32>,
        phi_spacings: Vec<f32>,
    },
    Custom {
        bins: Vec<(f32, f32)>,
        file: Option<String>,
    },
}

/// Container for angular binning scheme configuration.
/// 
/// **Context**: Binning schemes need to be passed between different parts
/// of the simulation and exposed to Python interfaces for external control.
/// This wrapper provides the necessary serialization and Python bindings.
/// 
/// **How it Works**: Simple wrapper around the Scheme enum with Python
/// bindings for integration with external analysis tools.
#[pyclass]
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct BinningScheme {
    pub scheme: Scheme,
}

#[pymethods]
impl BinningScheme {
    /// Creates a custom binning scheme from a list of (theta, phi) pairs.
    /// 
    /// **Context**: Python interfaces require constructors for creating
    /// binning schemes from external data sources or analysis scripts.
    /// 
    /// **How it Works**: Wraps the provided bins in a Custom scheme variant
    /// for use in the simulation framework.
    #[new]
    fn py_new(bins: Vec<(f32, f32)>) -> Self {
        BinningScheme {
            scheme: Scheme::Custom { bins, file: None },
        }
    }
}

/// Creates variable-spacing angular grids with specified resolution intervals.
/// 
/// **Context**: Some scattering applications require fine angular resolution
/// in specific regions (e.g., near forward scattering) and coarser resolution
/// elsewhere. Variable spacing optimizes computational efficiency while
/// maintaining accuracy where needed.
/// 
/// **How it Works**: Takes split points and corresponding spacings to create
/// a variable-resolution grid. Validates that each interval is an integer
/// multiple of its spacing to ensure proper grid alignment.
pub fn interval_spacings(splits: &Vec<f32>, spacings: &Vec<f32>) -> Vec<f32> {
    let num_values = splits.len();
    let mut values = Vec::new();

    for i in 0..num_values - 1 {
        // Iterate over the splits

        // compute the number of values between the splits
        let jmax = ((splits[i + 1] - splits[i]) / spacings[i]).round() as usize;

        // validate that the split is close to an integer multiple of the spacing
        let remainder = (splits[i + 1] - splits[i]) % spacings[i];
        if remainder.abs() > 1e-3 && (spacings[i] - remainder).abs() > 1e-3 {
            panic!(
                "Invalid spacing: split at index {} (value: {}) to index {} (value: {}) is not an integer multiple of spacing {}. Computed remainder: {}",
                i,
                splits[i],
                i + 1,
                splits[i + 1],
                spacings[i],
                remainder
            );
        }

        for j in 0..=jmax {
            // Iterate over the number of values between the splits
            if i != num_values - 2 && j == jmax {
                // skip the last value unless it is the last split
                continue;
            } else {
                values.push(splits[i] + j as f32 * spacings[i]);
            }
        }
    }

    values
}

/// Generates 2D angular grid from variable-spacing theta and phi intervals.
/// 
/// **Context**: Variable spacing in both theta and phi dimensions allows
/// optimal angular resolution placement. Forward scattering features often
/// require fine theta resolution, while azimuthal symmetry may allow
/// coarser phi spacing.
/// 
/// **How it Works**: Creates variable-spacing grids for both theta and phi,
/// then generates all combinations to produce the full 2D angular grid.
pub fn interval_bins(
    theta_spacing: &Vec<f32>,
    theta_splits: &Vec<f32>,
    phi_spacing: &Vec<f32>,
    phi_splits: &Vec<f32>,
) -> Vec<(f32, f32)> {
    let thetas = interval_spacings(theta_splits, theta_spacing);
    let phis = interval_spacings(phi_splits, phi_spacing);

    let mut bins = Vec::new();
    for theta in thetas.iter() {
        for phi in phis.iter() {
            bins.push((*theta, *phi));
        }
    }

    bins
}

/// Creates uniform angular grid for general scattering analysis.
/// 
/// **Context**: Many scattering applications require uniform angular sampling
/// across the full observation sphere. This provides the simplest approach
/// for general-purpose analysis without specific angular features.
/// 
/// **How it Works**: Generates linearly-spaced theta values from 0° to 180°
/// and phi values from 0° to 360°, then creates all combinations for the
/// complete angular grid.
/// 
/// # Example
/// ```rust
/// let theta_phi_combinations = bins::simple_bins(180, 180);
/// ```
pub fn simple_bins(num_theta: usize, num_phi: usize) -> Vec<(f32, f32)> {
    let thetas = Array1::linspace(0.0, 180.0, num_theta).insert_axis(ndarray::Axis(1)); // Reshape to (50, 1)
    let phis = Array1::linspace(0.0, 360.0, num_phi).insert_axis(ndarray::Axis(0)); // Reshape to (1, 60)

    // Flatten the combinations of theta and phi into a 1D array of tuples
    thetas
        .iter()
        .flat_map(|&theta| phis.iter().map(move |&phi| (theta, phi)))
        .collect()
}

/// Generates angular bins according to the specified binning scheme.
/// 
/// **Context**: The simulation core needs a unified interface for generating
/// angular bins regardless of the chosen discretization scheme. This function
/// dispatches to appropriate generators based on scheme type.
/// 
/// **How it Works**: Pattern matches on the scheme type and calls the
/// corresponding bin generation function. For custom schemes, handles
/// file loading and parsing with fallback to default bins.
/// 
/// # Example
/// ```rust
/// let bins = generate_bins(&self.settings.bins.scheme);
/// ```
pub fn generate_bins(bin_type: &Scheme) -> Vec<(f32, f32)> {
    match bin_type {
        Scheme::Simple { num_theta, num_phi } => simple_bins(*num_theta, *num_phi),
        Scheme::Interval {
            thetas,
            theta_spacings,
            phis,
            phi_spacings,
        } => interval_bins(theta_spacings, thetas, phi_spacings, phis),
        Scheme::Custom { bins, file } => {
            println!("Loading custom bins from file: {:?}", file);
            if let Some(file) = file {
                let content = match std::fs::read_to_string(file) {
                    Ok(content) => content,
                    Err(e) => panic!("Could not read file '{}': {}", file, e),
                };

                // Parse the TOML file
                match toml::from_str::<CustomBins>(&content) {
                    Ok(custom_bins) => {
                        println!("Loaded {} custom bins from file", custom_bins.bins.len());
                        custom_bins.bins
                    }
                    Err(e) => {
                        eprintln!("Error parsing custom bins file: {}", e);
                        eprintln!("Falling back to default bins");
                        bins.to_vec()
                    }
                }
            } else {
                bins.to_vec()
            }
        }
    }
}

/// TOML file structure for custom angular bin definitions.
/// 
/// **Context**: External analysis workflows may require specific angular
/// sampling patterns that need to be loaded from configuration files.
/// 
/// **How it Works**: Defines the expected TOML file structure for custom
/// bin specifications with (theta, phi) angle pairs.
#[derive(Debug, Deserialize)]
struct CustomBins {
    bins: Vec<(f32, f32)>,
}
