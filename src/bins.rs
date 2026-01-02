use nalgebra::Vector3;
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use serde::{Deserialize, Deserializer, Serialize};

/// Represents a solid angle bin with theta and phi bins
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SolidAngleBin {
    pub theta: AngleBin,
    pub phi: AngleBin,
}

impl SolidAngleBin {
    /// Create a new bin from theta and phi bins
    pub fn new(theta_bin: AngleBin, phi_bin: AngleBin) -> Self {
        SolidAngleBin {
            theta: theta_bin,
            phi: phi_bin,
        }
    }

    pub fn solid_angle(&self) -> f32 {
        2.0 * (self.theta.center).to_radians().sin().abs()
            * (0.5 * self.theta.width()).to_radians().sin()
            * self.phi.width().to_radians()
    }

    /// Returns the unit observation vector for this bin's center direction.
    ///
    /// Uses the inverted z-axis convention where z → -z, converting from
    /// spherical (theta, phi) to Cartesian coordinates.
    pub fn unit_vector(&self) -> Vector3<f32> {
        let (sin_theta, cos_theta) = self.theta.center.to_radians().sin_cos();
        let (sin_phi, cos_phi) = self.phi.center.to_radians().sin_cos();
        Vector3::new(sin_theta * cos_phi, sin_theta * sin_phi, -cos_theta)
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
        assert_eq!(result[0].phi.center, 60.0);
        assert_eq!(result[0].phi.center, 60.0);
        // Check bin edges for first theta bin
        assert_eq!(result[0].theta.min, 0.0);
        assert_eq!(result[0].theta.max, 60.0);
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

    /// Returns the theta range (min, max) for this scheme.
    pub fn theta_range(&self) -> (f32, f32) {
        match self {
            Scheme::Simple { .. } => (0.0, 180.0),
            Scheme::Interval { thetas, .. } => {
                let min = thetas.first().copied().unwrap_or(0.0);
                let max = thetas.last().copied().unwrap_or(180.0);
                (min, max)
            }
            Scheme::Custom { bins, .. } => {
                if bins.is_empty() {
                    return (0.0, 0.0);
                }
                let min = bins.iter().map(|b| b[0][0]).fold(f32::INFINITY, f32::min);
                let max = bins
                    .iter()
                    .map(|b| b[0][1])
                    .fold(f32::NEG_INFINITY, f32::max);
                (min, max)
            }
        }
    }

    /// Generate the bins for this scheme.
    pub fn generate(&self) -> Vec<SolidAngleBin> {
        match self {
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
}

/// Angular binning scheme for scattering calculations.
///
/// Defines how to discretize the scattering sphere into angular bins
/// for Mueller matrix and amplitude computations. Supports simple
/// regular grids, custom intervals, and arbitrary bin arrangements.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
pub struct BinningScheme {
    pub scheme: Scheme,
}

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
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
        if num_theta <= 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "num_theta must be greater than 0",
            ));
        }
        if num_phi <= 0 {
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

    /// Returns a list of all theta bin centre values
    fn thetas(&self) -> Vec<f32> {
        let bins = match &self.scheme {
            Scheme::Simple { num_theta, .. } => simple_spacings(*num_theta, 180.0)
                .iter()
                .map(|&bin| bin.center)
                .collect(),
            Scheme::Interval {
                thetas,
                theta_spacings,
                ..
            } => edges_to_bins(interval_spacings(&thetas, &theta_spacings))
                .iter()
                .map(|bin| bin.center)
                .collect(),
            Scheme::Custom { bins, .. } => custom_bins(&bins)
                .iter()
                .map(|bin| bin.theta.center)
                .collect(),
        };
        bins
    }

    /// Returns a list of all phi bin centre values
    fn phis(&self) -> Vec<f32> {
        let bins = match &self.scheme {
            Scheme::Simple { num_phi, .. } => simple_spacings(*num_phi, 360.0)
                .iter()
                .map(|&bin| bin.center)
                .collect(),
            Scheme::Interval {
                phis, phi_spacings, ..
            } => edges_to_bins(interval_spacings(&phis, &phi_spacings))
                .iter()
                .map(|bin| bin.center)
                .collect(),
            Scheme::Custom { bins, .. } => custom_bins(&bins)
                .iter()
                .map(|bin| bin.phi.center)
                .collect(),
        };
        bins
    }

    /// Returns all 2D bins as a numpy array of shape (n_bins, 2) with columns [theta, phi]
    fn bins<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        use numpy::IntoPyArray;
        let solid_bins = self.scheme.generate();
        let flat: Vec<f32> = solid_bins
            .iter()
            .flat_map(|bin| vec![bin.theta.center, bin.phi.center])
            .collect();
        ndarray::Array2::from_shape_vec((solid_bins.len(), 2), flat)
            .unwrap()
            .into_pyarray(py)
    }

    /// Returns unique 1D theta bins as a numpy array
    fn bins_1d<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f32>> {
        use numpy::IntoPyArray;
        let thetas = self.thetas();
        // Get unique thetas (they may repeat for each phi)
        let mut unique_thetas: Vec<f32> = thetas.clone();
        unique_thetas.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        ndarray::Array1::from_vec(unique_thetas).into_pyarray(py)
    }

    /// Returns the number of bins
    fn num_bins(&self) -> usize {
        self.scheme.generate().len()
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
    let theta_bins = edges_to_bins(theta_edges);
    let phi_bins = edges_to_bins(phi_edges);

    let mut bins = Vec::new();
    for theta_bin in theta_bins.iter() {
        for phi_bin in phi_bins.iter() {
            bins.push(SolidAngleBin::new(*theta_bin, *phi_bin));
        }
    }

    bins
}

fn edges_to_bins(edges: Vec<f32>) -> Vec<AngleBin> {
    let bins: Vec<AngleBin> = edges
        .windows(2)
        .map(|edges| AngleBin::new(edges[0], edges[1]))
        .collect();
    bins
}

/// Helper function to generate evenly spaced angle bins
fn simple_spacings(num_bins: usize, limit: f32) -> Vec<AngleBin> {
    let dangle = limit / (num_bins as f32);
    (0..num_bins)
        .map(|i| {
            let min = i as f32 * dangle;
            let max = (i + 1) as f32 * dangle;
            AngleBin::new(min, max)
        })
        .collect()
}

/// Generate theta and phi bin combinations
pub fn simple_bins(num_theta: usize, num_phi: usize) -> Vec<SolidAngleBin> {
    let theta_bins = simple_spacings(num_theta, 180.0);
    let phi_bins = simple_spacings(num_phi, 360.0);

    // meshgrid
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

// pub fn generate_bins(bin_type: &Scheme) -> Vec<SolidAngleBin> {
//     match bin_type {
//         Scheme::Simple {
//             num_theta, num_phi, ..
//         } => simple_bins(*num_theta, *num_phi),
//         Scheme::Interval {
//             thetas,
//             theta_spacings,
//             phis,
//             phi_spacings,
//         } => interval_bins(theta_spacings, thetas, phi_spacings, phis),
//         Scheme::Custom { bins, .. } => custom_bins(bins),
//     }
// }

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
        if theta >= bin.theta.min
            && theta < bin.theta.max
            && phi >= bin.phi.min
            && phi < bin.phi.max
        {
            bin_idx = Some(i);
            break;
        }
    }
    bin_idx
}
