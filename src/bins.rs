use pyo3::prelude::*;
use serde::{Deserialize, Deserializer, Serialize};

/// Represents a solid angle bin with theta and phi bins
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SolidAngleBin {
    pub theta_bin: AngleBin,
    pub phi_bin: AngleBin,
}

impl SolidAngleBin {
    /// Create a new bin from theta and phi bins
    pub fn new(theta_bin: AngleBin, phi_bin: AngleBin) -> Self {
        SolidAngleBin { theta_bin, phi_bin }
    }
    pub fn solid_angle(&self) -> f32 {
        2.0 * (self.theta_bin.center).to_radians().sin().abs()
            * (0.5 * self.theta_bin.width()).to_radians().sin()
            * self.phi_bin.width().to_radians()
    }
}

/// Represents an angular bin with edges and center. Fields: `min`, `max`, `center`
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize)]
pub struct AngleBin {
    pub min: f32,    // min edge
    pub max: f32,    // max edge
    pub center: f32, // center
}

impl AngleBin {
    /// Create a new bin from edges
    pub fn new(min: f32, max: f32) -> Self {
        AngleBin {
            min,
            max,
            center: (min + max) / 2.0,
        }
    }

    /// Create a bin from center and width
    pub fn from_center_width(center: f32, width: f32) -> Self {
        AngleBin {
            min: center - width / 2.0,
            max: center + width / 2.0,
            center,
        }
    }

    /// Get the width of the bin
    pub fn width(&self) -> f32 {
        self.max - self.min
    }
}

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
        // Check that we have the right number of bins
        assert_eq!(result.len(), 9);
        // Check first bin centers
        assert_eq!(result[0].phi_bin.center, 60.0);
        assert_eq!(result[0].phi_bin.center, 60.0);
        // Check bin edges for first theta bin
        assert_eq!(result[0].theta_bin.min, 0.0);
        assert_eq!(result[0].theta_bin.max, 60.0);
    }
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub enum Scheme {
    Simple {
        num_theta: usize,
        num_phi: usize,
        delta_theta: f32,
        delta_phi: f32,
    },
    Interval {
        thetas: Vec<f32>,
        theta_spacings: Vec<f32>,
        phis: Vec<f32>,
        phi_spacings: Vec<f32>,
    },
    Custom {
        bins: Vec<[[f32; 2]; 2]>, // Each bin is [[theta_min, theta_max], [phi_min, phi_max]]
        file: Option<String>,
    },
}

// Custom deserializer to handle missing delta_theta and delta_phi
impl<'de> Deserialize<'de> for Scheme {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct SimpleHelper {
            num_theta: usize,
            num_phi: usize,
            #[serde(default)]
            delta_theta: Option<f32>,
            #[serde(default)]
            delta_phi: Option<f32>,
        }

        #[derive(Deserialize)]
        struct IntervalHelper {
            thetas: Vec<f32>,
            theta_spacings: Vec<f32>,
            phis: Vec<f32>,
            phi_spacings: Vec<f32>,
        }

        #[derive(Deserialize)]
        struct CustomHelper {
            #[serde(default)]
            bins: Vec<[[f32; 2]; 2]>,
            file: Option<String>,
        }

        #[derive(Deserialize)]
        enum SchemeHelper {
            Simple(SimpleHelper),
            Interval(IntervalHelper),
            Custom(CustomHelper),
        }

        let helper = SchemeHelper::deserialize(deserializer)?;
        match helper {
            SchemeHelper::Simple(SimpleHelper {
                num_theta,
                num_phi,
                delta_theta,
                delta_phi,
            }) => {
                // Calculate deltas if not provided
                let delta_theta = delta_theta.unwrap_or(180.0 / num_theta as f32);
                let delta_phi = delta_phi.unwrap_or(360.0 / num_phi as f32);
                Ok(Scheme::Simple {
                    num_theta,
                    num_phi,
                    delta_theta,
                    delta_phi,
                })
            }
            SchemeHelper::Interval(IntervalHelper {
                thetas,
                theta_spacings,
                phis,
                phi_spacings,
            }) => Ok(Scheme::Interval {
                thetas,
                theta_spacings,
                phis,
                phi_spacings,
            }),
            SchemeHelper::Custom(CustomHelper { mut bins, file }) => {
                // If file is specified, load bins from file
                if let Some(ref filepath) = file {
                    #[derive(Deserialize)]
                    struct CustomBinsFile {
                        bins: Vec<[[f32; 2]; 2]>,
                    }

                    let content = std::fs::read_to_string(filepath).map_err(|e| {
                        serde::de::Error::custom(format!(
                            "Failed to read custom bins file '{}': {}",
                            filepath, e
                        ))
                    })?;

                    let file_data: CustomBinsFile = toml::from_str(&content).map_err(|e| {
                        serde::de::Error::custom(format!(
                            "Failed to parse custom bins file '{}': {}",
                            filepath, e
                        ))
                    })?;

                    bins = file_data.bins;
                }

                Ok(Scheme::Custom { bins, file })
            }
        }
    }
}

impl Scheme {
    pub fn new_simple(num_theta: usize, num_phi: usize) -> Self {
        let delta_theta = 180.0 / num_theta as f32;
        let delta_phi = 360.0 / num_phi as f32;
        Scheme::Simple {
            num_theta,
            num_phi,
            delta_theta,
            delta_phi,
        }
    }
}

/// Angular binning scheme for scattering calculations.
///
/// Defines how to discretize the scattering sphere into angular bins
/// for Mueller matrix and amplitude computations. Supports simple
/// regular grids, custom intervals, and arbitrary bin arrangements.
#[pyclass]
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
pub struct BinningScheme {
    pub scheme: Scheme,
}

#[pymethods]
impl BinningScheme {
    #[new]
    fn py_new(bins: Vec<[[f32; 2]; 2]>) -> Self {
        BinningScheme {
            scheme: Scheme::Custom { bins, file: None },
        }
    }

    /// Create a simple binning scheme with uniform theta and phi spacing
    #[staticmethod]
    fn simple(num_theta: usize, num_phi: usize) -> PyResult<Self> {
        if num_theta == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "num_theta must be greater than 0",
            ));
        }
        if num_phi == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "num_phi must be greater than 0",
            ));
        }

        Ok(BinningScheme {
            scheme: Scheme::new_simple(num_theta, num_phi),
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

    /// Create a custom binning scheme with explicit bin edges
    /// Each bin is specified as [[theta_min, theta_max], [phi_min, phi_max]]
    #[staticmethod]
    fn custom(bins: Vec<[[f32; 2]; 2]>) -> Self {
        BinningScheme {
            scheme: Scheme::Custom { bins, file: None },
        }
    }
}

pub fn interval_spacings(splits: &[f32], spacings: &[f32]) -> Vec<f32> {
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
            let val = splits[i] + j as f32 * spacings[i];

            // Iterate over the number of values between the splits
            if i != num_values - 2 && j == jmax {
                // skip the last value unless it is the last split
                continue;
            }

            values.push(val);
        }
    }

    values
}

pub fn interval_bins(
    theta_spacing: &Vec<f32>,
    theta_splits: &Vec<f32>,
    phi_spacing: &Vec<f32>,
    phi_splits: &Vec<f32>,
) -> Vec<SolidAngleBin> {
    // Get edge positions
    let theta_edges = interval_spacings(theta_splits, theta_spacing);
    let phi_edges = interval_spacings(phi_splits, phi_spacing);

    // Convert edges to bins
    let theta_bins: Vec<AngleBin> = theta_edges
        .windows(2)
        .map(|edges| AngleBin::new(edges[0], edges[1]))
        .collect();

    let phi_bins: Vec<AngleBin> = phi_edges
        .windows(2)
        .map(|edges| AngleBin::new(edges[0], edges[1]))
        .collect();

    let mut bins = Vec::new();
    for theta_bin in theta_bins.iter() {
        for phi_bin in phi_bins.iter() {
            bins.push(SolidAngleBin::new(*theta_bin, *phi_bin));
        }
    }

    bins
}

/// Generate theta and phi bin combinations
pub fn simple_bins(num_theta: usize, num_phi: usize) -> Vec<SolidAngleBin> {
    // Create theta bins
    let dtheta = 180.0 / (num_theta as f32);
    let theta_bins: Vec<AngleBin> = (0..num_theta)
        .map(|i| {
            let min = i as f32 * dtheta;
            let max = (i + 1) as f32 * dtheta;
            AngleBin::new(min, max)
        })
        .collect();

    // Create phi bins
    let dphi = 360.0 / (num_phi as f32);
    let phi_bins: Vec<AngleBin> = (0..num_phi)
        .map(|i| {
            let min = i as f32 * dphi;
            let max = (i + 1) as f32 * dphi;
            AngleBin::new(min, max)
        })
        .collect();

    let mut bins = Vec::new();
    for theta_bin in theta_bins.iter() {
        for phi_bin in phi_bins.iter() {
            bins.push(SolidAngleBin::new(*theta_bin, *phi_bin));
        }
    }

    bins
}

/// Generate custom bins from explicit edge specifications
/// Each bin is [[theta_min, theta_max], [phi_min, phi_max]]
pub fn custom_bins(bin_specs: &[[[f32; 2]; 2]]) -> Vec<SolidAngleBin> {
    bin_specs
        .iter()
        .map(|&[[theta_min, theta_max], [phi_min, phi_max]]| {
            let theta_bin = AngleBin::new(theta_min, theta_max);
            let phi_bin = AngleBin::new(phi_min, phi_max);
            SolidAngleBin::new(theta_bin, phi_bin)
        })
        .collect()
}

pub fn generate_bins(bin_type: &Scheme) -> Vec<SolidAngleBin> {
    match bin_type {
        Scheme::Simple {
            num_theta, num_phi, ..
        } => simple_bins(*num_theta, *num_phi),
        Scheme::Interval {
            thetas,
            theta_spacings,
            phis,
            phi_spacings,
        } => interval_bins(theta_spacings, thetas, phi_spacings, phis),
        Scheme::Custom { bins, .. } => custom_bins(bins),
    }
}

/// Gets the index of a theta-phi bin, assuming a `Simple` binning scheme, given an input theta and phi.
pub fn get_n_simple(
    num_theta: usize,
    num_phi: usize,
    delta_theta: f32,
    delta_phi: f32,
    theta: f32,
    phi: f32,
) -> Option<usize> {
    let n_theta = ((theta / delta_theta).floor() as usize).min(num_theta - 1);
    let n_phi = ((phi / delta_phi).floor() as usize).min(num_phi - 1);
    Some(n_theta * num_phi + n_phi)
}

/// Gets the index of a theta-phi bin by linearly searching through the bins until a match is found. Returns `None` if no match is found.
pub fn get_n_linear_search(bins: &[SolidAngleBin], theta: f32, phi: f32) -> Option<usize> {
    // Find the corresponding bin in the bins array
    let mut bin_idx = None;
    for (i, bin) in bins.iter().enumerate() {
        if theta >= bin.theta_bin.min
            && theta < bin.theta_bin.max
            && phi >= bin.phi_bin.min
            && phi < bin.phi_bin.max
        {
            bin_idx = Some(i);
            break;
        }
    }
    bin_idx
}
