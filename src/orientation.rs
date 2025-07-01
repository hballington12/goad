use clap::Subcommand;
use nalgebra::Matrix3;
use std::{f32::consts::PI, str::FromStr};

use anyhow::Result;
use pyo3::prelude::*;
use rand::Rng;
use rand::SeedableRng;
use serde::Deserialize;

#[pyclass]
#[derive(Subcommand, Debug, Clone, Deserialize, PartialEq)]
pub enum Scheme {
    /// Solve the problem by averaging over a uniform distribution of angles.
    /// Example: `uniform 100`
    Uniform { num_orients: usize },
    /// Solve the problem by averaging over a discrete set of angles (in degrees).
    /// Example: `discrete 0,0,0 20,30,40`
    Discrete { eulers: Vec<Euler> },
}

/// Euler angle order for the discrete orientation scheme.
#[pyclass]
#[derive(Debug, Clone, Deserialize, PartialEq, Copy)]
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

#[pyclass]
#[derive(Debug, Clone, Deserialize, PartialEq)]
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

#[pymethods]
impl Euler {
    #[new]
    fn py_new(alpha: f32, beta: f32, gamma: f32) -> Self {
        Euler::new(alpha, beta, gamma)
    }

    fn __repr__(&self) -> String {
        format!("Euler(alpha={}, beta={}, gamma={})", self.alpha, self.beta, self.gamma)
    }
}

#[pyclass]
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Orientation {
    #[pyo3(get, set)]
    pub scheme: Scheme,
    #[pyo3(get, set)]
    pub euler_convention: EulerConvention,
}

#[pymethods]
impl Orientation {
    #[new]
    fn py_new(scheme: Scheme, euler_convention: Option<EulerConvention>) -> Self {
        Orientation {
            scheme,
            euler_convention: euler_convention.unwrap_or(EulerConvention::ZYZ),
        }
    }

    fn __repr__(&self) -> String {
        format!("Orientation(scheme={:?}, euler_convention={:?})", self.scheme, self.euler_convention)
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
