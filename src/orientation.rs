use nalgebra::Matrix3;
use serde::Serialize;
use std::{f32::consts::PI, str::FromStr};

use anyhow::Result;
use pyo3::prelude::*;
#[cfg(feature = "stub-gen")]
use pyo3_stub_gen::derive::*;
use rand::Rng;
use rand::SeedableRng;
use serde::Deserialize;

use crate::settings::DEFAULT_EULER_ORDER;

#[cfg_attr(feature = "stub-gen", gen_stub_pyclass_enum)]
#[pyclass(module = "goad._goad")]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum Scheme {
    /// Solve the problem by averaging over a uniform distribution of angles.
    /// Example: `uniform 100`
    Uniform { num_orients: usize },
    /// Solve the problem by averaging over a discrete set of angles (in degrees).
    /// Example: `discrete 0,0,0 20,30,40`
    Discrete { eulers: Vec<Euler> },
    /// Solve the problem using Sobol quasi-random sequence for faster convergence.
    /// Example: `sobol 100`
    Sobol { num_orients: usize },
    /// Solve the problem using Halton quasi-random sequence for faster convergence.
    /// Example: `halton 100`
    Halton { num_orients: usize },
}

/// Euler angle order for the discrete orientation scheme.
#[cfg_attr(feature = "stub-gen", gen_stub_pyclass_enum)]
#[pyclass(module = "goad._goad")]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Copy)]
pub enum EulerConvention {
    XZX,
    XYX,
    YXY,
    YZY,
    ZYZ,
    ZXZ,
    XZY,
    XYZ,
    YXZ,
    YZX,
    ZYX,
    ZXY,
}

#[pymethods]
impl EulerConvention {
    #[new]
    fn py_new(str: &str) -> PyResult<Self> {
        match str.to_lowercase().as_str() {
            "xzx" => Ok(EulerConvention::XZX),
            "xyx" => Ok(EulerConvention::XYX),
            "yxy" => Ok(EulerConvention::YXY),
            "yzy" => Ok(EulerConvention::YZY),
            "zyz" => Ok(EulerConvention::ZYZ),
            "zxz" => Ok(EulerConvention::ZXZ),
            "xzy" => Ok(EulerConvention::XZY),
            "xyz" => Ok(EulerConvention::XYZ),
            "yxz" => Ok(EulerConvention::YXZ),
            "yzx" => Ok(EulerConvention::YZX),
            "zyx" => Ok(EulerConvention::ZYX),
            "zxy" => Ok(EulerConvention::ZXY),
            _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "'{}' is not a valid Euler convention. Valid options are: 'XZX', 'XYX', 'YXY', 'YZY', 'ZYZ', 'ZXZ', 'XZY', 'XYZ', 'YXZ', 'YZX', 'ZYX', 'ZXY'",
                str
            ))),
        }
    }
}

#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Euler {
    #[pyo3(get, set)]
    pub alpha: f32,
    #[pyo3(get, set)]
    pub beta: f32,
    #[pyo3(get, set)]
    pub gamma: f32,
}

impl Euler {
    pub fn new(alpha: f32, beta: f32, gamma: f32) -> Self {
        Self { alpha, beta, gamma }
    }
    pub fn rotation_matrix(&self, convention: EulerConvention) -> Matrix3<f32> {
        let alpha = self.alpha.to_radians();
        let beta = self.beta.to_radians();
        let gamma = self.gamma.to_radians();

        let s1 = alpha.sin();
        let s2 = beta.sin();
        let s3 = gamma.sin();
        let c1 = alpha.cos();
        let c2 = beta.cos();
        let c3 = gamma.cos();

        match convention {
            EulerConvention::XZX => Matrix3::new(
                c2,
                -c3 * s2,
                s2 * s3,
                c1 * s2,
                c1 * c2 * c3 - s1 * s3,
                -c3 * s1 - c1 * c2 * s3,
                s1 * s2,
                c1 * s3 + c2 * c3 * s1,
                c1 * c3 - c2 * s1 * s3,
            ),
            EulerConvention::XYX => Matrix3::new(
                c2,
                s2 * s3,
                c3 * s2,
                s1 * s2,
                c1 * c3 - c2 * s1 * s3,
                -c1 * s3 - c2 * c3 * s1,
                -c1 * s2,
                c2 * s1 + c1 * c2 * s3,
                c1 * c2 * c3 - s1 * s3,
            ),
            EulerConvention::YXY => Matrix3::new(
                c1 * c3 - c2 * s1 * s3,
                s1 * s2,
                c1 * s3 + c2 * c3 * s1,
                s2 * s3,
                c2,
                -c3 * s2,
                -c3 * s1 - c1 * c2 * s3,
                c1 * s2,
                c1 * c2 * c3 - s1 * s3,
            ),
            EulerConvention::YZY => Matrix3::new(
                c1 * c2 * c3 - s1 * s3,
                -c1 * s2,
                c1 * s3 + c1 * c2 * s3,
                c3 * s2,
                c2,
                s2 * s3,
                -s1 * c2 * c3 - c1 * s3,
                s1 * s2,
                c1 * c3 - c2 * s1 * s3,
            ),
            EulerConvention::ZYZ => Matrix3::new(
                c1 * c2 * c3 - s1 * s3,
                -c1 * c2 * s3 - s1 * c3,
                c1 * s2,
                s1 * c2 * c3 + c1 * s3,
                -s1 * c2 * s3 + c1 * c3,
                s1 * s2,
                -s2 * c3,
                s2 * s3,
                c2,
            ),
            EulerConvention::ZXZ => Matrix3::new(
                c1 * c3 - s1 * s3 * c2,
                -c1 * s3 - s1 * c3 * c2,
                s1 * s2,
                s1 * c3 + c1 * s3 * c2,
                -s1 * s3 + c1 * c3 * c2,
                -c1 * s2,
                s3 * s2,
                c3 * s2,
                c2,
            ),
            EulerConvention::XZY => Matrix3::new(
                c2 * c3,
                -s2,
                c2 * s3,
                s1 * s3 + c1 * c3 * s2,
                c1 * c2,
                c1 * s2 * s3 - c3 * s1,
                c3 * s1 * s2 - c1 * s3,
                c2 * s1,
                c1 * c3 + s1 * s2 * s3,
            ),
            EulerConvention::XYZ => Matrix3::new(
                c2 * c3,
                -c2 * s3,
                s2,
                c1 * s3 + c3 * s1 * s2,
                c1 * c3 - s1 * s2 * s3,
                -c2 * s1,
                s1 * s3 - c1 * c3 * s2,
                c3 * s1 + c1 * s2 * s3,
                c1 * c2,
            ),
            EulerConvention::YXZ => Matrix3::new(
                c1 * c3 + s1 * s2 * s3,
                c3 * s1 * s2 - c1 * s3,
                c2 * s1,
                c2 * s3,
                c2 * c3,
                -s2,
                c1 * s2 * s3 - c3 * s1,
                s1 * s3 + c1 * c3 * s2,
                c1 * c2,
            ),
            EulerConvention::YZX => Matrix3::new(
                c1 * c2,
                c1 * s2 * s3 - c3 * s1,
                s1 * s3 + c1 * c3 * s2,
                s2,
                c2 * c3,
                -c2 * s3,
                -c2 * s1,
                c1 * s3 + c3 * s1 * s2,
                c1 * c3 - s1 * s2 * s3,
            ),
            EulerConvention::ZYX => Matrix3::new(
                c1 * c2,
                c1 * s2 * s3 - c3 * s1,
                s1 * s3 + c1 * c3 * s2,
                c2 * s1,
                c1 * c3 + s1 * s2 * s3,
                c3 * s1 * s2 - c1 * s3,
                -s2,
                c2 * s3,
                c2 * c3,
            ),
            EulerConvention::ZXY => Matrix3::new(
                c1 * c3 - s1 * s2 * s3,
                -c2 * s1,
                c1 * s3 + c3 * s1 * s2,
                c3 * s1 + c1 * s2 * s3,
                c1 * c2,
                s1 * s3 - c1 * c3 * s2,
                -c2 * s3,
                s2,
                c2 * c3,
            ),
        }
    }
}

impl FromStr for Euler {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 3 {
            return Err(anyhow::anyhow!("Invalid number of angles"));
        }
        let a = parts[0].trim().parse::<f32>()?;
        let b = parts[1].trim().parse::<f32>()?;
        let c = parts[2].trim().parse::<f32>()?;
        Ok(Euler {
            alpha: a,
            beta: b,
            gamma: c,
        })
    }
}

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl Euler {
    #[new]
    fn py_new(alpha: f32, beta: f32, gamma: f32) -> Self {
        Euler::new(alpha, beta, gamma)
    }

    fn __repr__(&self) -> String {
        format!(
            "Euler(alpha={}, beta={}, gamma={})",
            self.alpha, self.beta, self.gamma
        )
    }
}

#[cfg_attr(feature = "stub-gen", gen_stub_pyclass)]
#[pyclass(module = "goad._goad")]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Orientation {
    #[pyo3(get, set)]
    pub scheme: Scheme,
    #[pyo3(get, set)]
    pub euler_convention: EulerConvention,
}

#[cfg_attr(feature = "stub-gen", gen_stub_pymethods)]
#[pymethods]
impl Orientation {
    #[staticmethod]
    #[pyo3(name = "uniform", signature = (num_orients, euler_convention = None))]
    fn py_uniform(num_orients: usize, euler_convention: Option<EulerConvention>) -> Self {
        Orientation {
            scheme: Scheme::Uniform { num_orients },
            euler_convention: euler_convention.unwrap_or(DEFAULT_EULER_ORDER),
        }
    }
    #[staticmethod]
    #[pyo3(name = "discrete", signature = (eulers, euler_convention = None))]
    fn py_discrete(eulers: Vec<Euler>, euler_convention: Option<EulerConvention>) -> Self {
        Orientation {
            scheme: Scheme::Discrete { eulers },
            euler_convention: euler_convention.unwrap_or(DEFAULT_EULER_ORDER),
        }
    }

    #[staticmethod]
    #[pyo3(name = "sobol", signature = (num_orients, euler_convention = None))]
    fn py_sobol(num_orients: usize, euler_convention: Option<EulerConvention>) -> Self {
        Orientation {
            scheme: Scheme::Sobol { num_orients },
            euler_convention: euler_convention.unwrap_or(DEFAULT_EULER_ORDER),
        }
    }

    #[staticmethod]
    #[pyo3(name = "halton", signature = (num_orients, euler_convention = None))]
    fn py_halton(num_orients: usize, euler_convention: Option<EulerConvention>) -> Self {
        Orientation {
            scheme: Scheme::Halton { num_orients },
            euler_convention: euler_convention.unwrap_or(DEFAULT_EULER_ORDER),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Orientation(scheme={:?}, euler_convention={:?})",
            self.scheme, self.euler_convention
        )
    }
}

/// Compute the Halton sequence value for a given index and base.
/// Returns a value in [0, 1).
fn halton_sequence(index: u32, base: u32) -> f32 {
    let mut result = 0.0;
    let mut f = 1.0 / base as f32;
    let mut i = index;
    while i > 0 {
        result += f * (i % base) as f32;
        i /= base;
        f /= base as f32;
    }
    result
}

/// Orientation sampler that generates orientations on demand.
pub enum OrientationSampler {
    /// Generates uniform random orientations (never exhausts).
    Uniform {
        rng: rand::rngs::StdRng,
        seed: Option<u64>,
    },
    /// Iterates over a discrete set of orientations.
    Discrete { eulers: Vec<Euler>, index: usize },
    /// Generates orientations using Sobol quasi-random sequence (never exhausts).
    Sobol { seed: u32, index: u32 },
    /// Generates orientations using Halton quasi-random sequence (never exhausts).
    Halton { index: u32 },
}

impl OrientationSampler {
    /// Create a uniform random sampler.
    pub fn uniform(seed: Option<u64>) -> Self {
        let rng = if let Some(seed) = seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };
        Self::Uniform { rng, seed }
    }

    /// Create a discrete sampler from a list of Euler angles.
    pub fn discrete(eulers: Vec<Euler>) -> Self {
        Self::Discrete { eulers, index: 0 }
    }

    /// Create a Sobol quasi-random sampler.
    pub fn sobol(seed: Option<u64>) -> Self {
        Self::Sobol {
            seed: seed.unwrap_or(0) as u32,
            index: 0,
        }
    }

    /// Create a Halton quasi-random sampler.
    pub fn halton() -> Self {
        Self::Halton { index: 0 }
    }

    /// Get the next orientation. Returns None if exhausted (for discrete samplers).
    pub fn next(&mut self) -> Option<Euler> {
        match self {
            Self::Uniform { rng, .. } => {
                let alpha = rng.random_range(0.0..1.0) as f32 * 360.0;
                let beta = (1.0 - rng.random_range(0.0..1.0) as f32 * 2.0).acos() * 180.0 / PI;
                let gamma = rng.random_range(0.0..1.0) as f32 * 360.0;
                Some(Euler::new(alpha, beta, gamma))
            }
            Self::Discrete { eulers, index } => {
                if *index < eulers.len() {
                    let euler = eulers[*index].clone();
                    *index += 1;
                    Some(euler)
                } else {
                    None
                }
            }
            Self::Sobol { seed, index } => {
                let u1 = sobol_burley::sample(*index, 0, *seed);
                let u2 = sobol_burley::sample(*index, 1, *seed);
                let u3 = sobol_burley::sample(*index, 2, *seed);
                *index += 1;

                let alpha = u1 * 360.0;
                let beta = (1.0_f32 - u2 * 2.0).acos() * 180.0 / PI;
                let gamma = u3 * 360.0;
                Some(Euler::new(alpha, beta, gamma))
            }
            Self::Halton { index } => {
                // Use primes 2, 3, 5 for the three dimensions
                let u1 = halton_sequence(*index, 2);
                let u2 = halton_sequence(*index, 3);
                let u3 = halton_sequence(*index, 5);
                *index += 1;

                let alpha = u1 * 360.0;
                let beta = (1.0 - u2 * 2.0).acos() * 180.0 / PI;
                let gamma = u3 * 360.0;
                Some(Euler::new(alpha, beta, gamma))
            }
        }
    }

    /// Reset the sampler to its initial state (for reproducibility).
    pub fn reset(&mut self) {
        match self {
            Self::Uniform { rng, seed } => {
                *rng = if let Some(s) = seed {
                    rand::rngs::StdRng::seed_from_u64(*s)
                } else {
                    rand::rngs::StdRng::from_rng(&mut rand::rng())
                };
            }
            Self::Discrete { index, .. } => {
                *index = 0;
            }
            Self::Sobol { index, .. } | Self::Halton { index } => {
                *index = 0;
            }
        }
    }
}

/// Orientation scheme for problem averaging. Can either be a discrete list of angles
/// or a distribution.
#[derive(Debug, Clone, PartialEq)]
pub struct Orientations {
    pub num_orientations: usize,
    pub eulers: Vec<(f32, f32, f32)>,
}

impl Orientations {
    pub fn generate(scheme: &Scheme, seed: Option<u64>) -> Orientations {
        match &scheme {
            Scheme::Uniform {
                num_orients: num_orientations,
            } => Orientations::random_uniform(*num_orientations, seed),
            Scheme::Discrete { eulers } => {
                let alphas: Vec<f32> = eulers.iter().map(|e| e.alpha).collect();
                let betas: Vec<f32> = eulers.iter().map(|e| e.beta).collect();
                let gammas: Vec<f32> = eulers.iter().map(|e| e.gamma).collect();
                Orientations::new_discrete(alphas, betas, gammas).unwrap()
            }
            Scheme::Sobol { num_orients } => Orientations::sobol(*num_orients, seed),
            Scheme::Halton { num_orients } => Orientations::halton(*num_orients),
        }
    }

    /// Generate orientations using Sobol quasi-random sequence.
    pub fn sobol(num_orient: usize, seed: Option<u64>) -> Orientations {
        let seed = seed.unwrap_or(0) as u32;
        let eulers: Vec<(f32, f32, f32)> = (0..num_orient as u32)
            .map(|i| {
                let u1 = sobol_burley::sample(i, 0, seed);
                let u2 = sobol_burley::sample(i, 1, seed);
                let u3 = sobol_burley::sample(i, 2, seed);
                let alpha = u1 * 360.0;
                let beta = (1.0_f32 - u2 * 2.0).acos() * 180.0 / PI;
                let gamma = u3 * 360.0;
                (alpha, beta, gamma)
            })
            .collect();
        Orientations {
            num_orientations: num_orient,
            eulers,
        }
    }

    /// Generate orientations using Halton quasi-random sequence.
    pub fn halton(num_orient: usize) -> Orientations {
        let eulers: Vec<(f32, f32, f32)> = (0..num_orient as u32)
            .map(|i| {
                let u1 = halton_sequence(i, 2);
                let u2 = halton_sequence(i, 3);
                let u3 = halton_sequence(i, 5);
                let alpha = u1 * 360.0;
                let beta = (1.0 - u2 * 2.0).acos() * 180.0 / PI;
                let gamma = u3 * 360.0;
                (alpha, beta, gamma)
            })
            .collect();
        Orientations {
            num_orientations: num_orient,
            eulers,
        }
    }

    /// Creates a new orientation scheme with the given discrete angles.
    pub fn new_discrete(alphas: Vec<f32>, betas: Vec<f32>, gammas: Vec<f32>) -> Result<Self> {
        if alphas.is_empty() || betas.is_empty() || gammas.is_empty() {
            return Err(anyhow::anyhow!("Empty angle list"));
        }
        if alphas.len() != betas.len() || alphas.len() != gammas.len() {
            return Err(anyhow::anyhow!("Angle lists have different lengths"));
        }
        Ok(Self {
            num_orientations: alphas.len(),
            eulers: alphas
                .into_iter()
                .zip(betas.into_iter())
                .zip(gammas.into_iter())
                .map(|((alpha, beta), gamma)| (alpha, beta, gamma))
                .collect(),
        })
    }

    pub fn random_uniform(num_orient: usize, seed: Option<u64>) -> Orientations {
        let mut rng = if let Some(seed) = seed {
            rand::rngs::StdRng::seed_from_u64(seed)
        } else {
            rand::rngs::StdRng::from_rng(&mut rand::rng())
        };

        let alphas: Vec<f32> = (0..num_orient)
            .map(|_| rng.random_range(0.0..1.0) as f32 * 360.0)
            .collect();
        let betas: Vec<f32> = (0..num_orient)
            .map(|_| (1.0 - rng.random_range(0.0..1.0) as f32 * 2.0).acos() * 180.0 / PI)
            .collect();
        let gammas: Vec<f32> = (0..num_orient)
            .map(|_| rng.random_range(0.0..1.0) as f32 * 360.0)
            .collect();

        let orientations = Orientations::new_discrete(alphas, betas, gammas).unwrap();
        orientations
    }
}
