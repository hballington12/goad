use clap::Subcommand;
use std::{f32::consts::PI, str::FromStr};

use anyhow::Result;
use rand::Rng;
use rand::SeedableRng;
use serde::Deserialize;

#[derive(Subcommand, Debug, Clone, Deserialize, PartialEq)]
pub enum Scheme {
    /// Solve the problem by averaging over a uniform distribution of angles.
    /// Example: `uniform 100`
    Uniform { num_orients: usize },
    /// Solve the problem by averaging over a discrete set of angles (in degrees).
    /// Example: `discrete 0,0,0 20,30,40`
    Discrete { eulers: Vec<Euler> },
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Euler {
    alpha: f32,
    beta: f32,
    gamma: f32,
}

impl Euler {
    pub fn new(alpha: f32, beta: f32, gamma: f32) -> Self {
        Self { alpha, beta, gamma }
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

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct OrientationScheme {
    pub scheme: Scheme,
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
