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

/// Angular binning scheme for scattering calculations.
/// 
/// Defines how to discretize the scattering sphere into angular bins
/// for Mueller matrix and amplitude computations. Supports simple
/// regular grids, custom intervals, and arbitrary bin arrangements.
#[pyclass]
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct BinningScheme {
    pub scheme: Scheme,
}

#[pymethods]
impl BinningScheme {
    #[new]
    fn py_new(bins: Vec<(f32, f32)>) -> Self {
        BinningScheme {
            scheme: Scheme::Custom { bins, file: None },
        }
    }

    /// Create a simple binning scheme with uniform theta and phi spacing
    #[staticmethod]
    fn simple(num_theta: usize, num_phi: usize) -> PyResult<Self> {
        if num_theta == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "num_theta must be greater than 0"
            ));
        }
        if num_phi == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "num_phi must be greater than 0"
            ));
        }
        
        Ok(BinningScheme {
            scheme: Scheme::Simple { num_theta, num_phi },
        })
    }

    /// Create an interval binning scheme with variable spacing
    #[staticmethod]
    fn interval(
        thetas: Vec<f32>,
        theta_spacings: Vec<f32>,
        phis: Vec<f32>,
        phi_spacings: Vec<f32>,
    ) -> Self {
        BinningScheme {
            scheme: Scheme::Interval {
                thetas,
                theta_spacings,
                phis,
                phi_spacings,
            },
        }
    }

    /// Create a custom binning scheme with explicit bin positions
    #[staticmethod]
    fn custom(bins: Vec<(f32, f32)>) -> Self {
        BinningScheme {
            scheme: Scheme::Custom { bins, file: None },
        }
    }
}

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

/// Generate theta and phi combinations
pub fn simple_bins(num_theta: usize, num_phi: usize) -> Vec<(f32, f32)> {
    let thetas = Array1::linspace(0.0, 180.0, num_theta).insert_axis(ndarray::Axis(1)); // Reshape to (50, 1)
    let phis = Array1::linspace(0.0, 360.0, num_phi).insert_axis(ndarray::Axis(0)); // Reshape to (1, 60)

    // Flatten the combinations of theta and phi into a 1D array of tuples
    thetas
        .iter()
        .flat_map(|&theta| phis.iter().map(move |&phi| (theta, phi)))
        .collect()
}

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
            // println!("Loading custom bins from file: {:?}", file);
            if let Some(file) = file {
                let content = match std::fs::read_to_string(file) {
                    Ok(content) => content,
                    Err(e) => panic!("Could not read file '{}': {}", file, e),
                };

                // Parse the TOML file
                match toml::from_str::<CustomBins>(&content) {
                    Ok(custom_bins) => {
                        // println!("Loaded {} custom bins from file", custom_bins.bins.len());
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

// Define a struct to match the TOML structure
#[derive(Debug, Deserialize)]
struct CustomBins {
    bins: Vec<(f32, f32)>,
}
