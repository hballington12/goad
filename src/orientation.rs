use std::f32::consts::PI;

use anyhow::Result;
use rand::Rng;
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum OrientationScheme {
    /// Uniform distribution of angles.
    Uniform { num_orients: usize },
    /// Discrete list of angles in degrees.
    Discrete {
        alphas: Vec<f32>,
        betas: Vec<f32>,
        gammas: Vec<f32>,
    },
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Config {
    pub orient_scheme: OrientationScheme,
}

/// Orientation scheme for problem averaging. Can either be a discrete list of angles
/// or a distribution.
#[derive(Debug, Clone, PartialEq)]
pub struct Orientations {
    pub num_orientations: usize,
    pub eulers: Vec<(f32, f32, f32)>,
}

impl Orientations {
    pub fn generate(scheme: &OrientationScheme) -> Orientations {
        match &scheme {
            OrientationScheme::Uniform {
                num_orients: num_orientations,
            } => Orientations::random_uniform(*num_orientations),
            OrientationScheme::Discrete {
                alphas,
                betas,
                gammas,
            } => {
                let alphas = alphas.clone().iter().map(|x| x * PI / 180.0).collect();
                let betas = betas.clone().iter().map(|x| x * PI / 180.0).collect();
                let gammas = gammas.clone().iter().map(|x| x * PI / 180.0).collect();
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

    pub fn random_uniform(num_orient: usize) -> Orientations {
        let mut rng = rand::rng();
        let alphas: Vec<f32> = (0..num_orient)
            .map(|_| rng.random_range(0.0..1.0) as f32 * 2.0 * PI)
            .collect();
        let betas: Vec<f32> = (0..num_orient)
            .map(|_| (1.0 - rng.random_range(0.0..1.0) as f32 * 2.0).acos())
            .collect();
        let gammas: Vec<f32> = (0..num_orient)
            .map(|_| rng.random_range(0.0..1.0) as f32 * 2.0 * PI)
            .collect();

        let orientations = Orientations::new_discrete(alphas, betas, gammas).unwrap();
        orientations
    }
}
