//! Particle orientation handling and Euler angle operations.
//!
//! This module provides comprehensive support for 3D particle orientations using
//! Euler angles with multiple convention support. It handles both discrete orientation
//! sets for systematic studies and uniform random orientations for ensemble averaging
//! in electromagnetic scattering simulations.
//!
//! The orientation system provides:
//! - Multiple Euler angle conventions (proper and Tait-Bryan)
//! - Discrete orientation specifications
//! - Uniform random orientation generation on SO(3)
//! - 3D rotation matrix computation
//! - Command-line and file-based orientation input
//! - Proper statistical sampling for ensemble averaging
//!
//! # Key Components
//!
//! - [`Euler`]: Three-angle rotation representation
//! - [`EulerConvention`]: Standard rotation sequence specifications
//! - [`Scheme`]: Orientation generation methods (discrete/uniform)
//! - [`Orientations`]: Generated orientation sets for simulation
//! - Rotation matrix computation for all conventions

use clap::Subcommand;
use nalgebra::Matrix3;
use std::{f32::consts::PI, str::FromStr};

use anyhow::Result;
use rand::Rng;
use rand::SeedableRng;
use serde::Deserialize;

/// Particle orientation schemes for scattering pattern analysis.
/// 
/// **Context**: Scattering patterns depend strongly on particle orientation relative
/// to the incident beam. Different applications require different orientation approaches:
/// fixed orientations for specific studies, random orientations for ensemble averaging,
/// or systematic discrete sets for detailed pattern analysis.
/// 
/// **How it Works**: Provides two main schemes - Uniform generates random orientations
/// following uniform distributions on the rotation group, while Discrete uses
/// user-specified Euler angle sets for deterministic orientation control.
#[derive(Subcommand, Debug, Clone, Deserialize, PartialEq)]
pub enum Scheme {
    /// Solve the problem by averaging over a uniform distribution of angles.
    /// Example: `uniform 100`
    Uniform { num_orients: usize },
    /// Solve the problem by averaging over a discrete set of angles (in degrees).
    /// Example: `discrete 0,0,0 20,30,40`
    Discrete { eulers: Vec<Euler> },
}

/// Standard Euler angle rotation conventions for 3D orientation specification.
/// 
/// **Context**: Euler angles provide a three-parameter description of 3D rotations,
/// but the order of rotations affects the final orientation. Different fields
/// use different conventions, requiring support for multiple rotation sequences
/// to maintain compatibility with external data sources.
/// 
/// **How it Works**: Enumerates all standard Euler conventions with both proper
/// Euler angles (repeated axis) and Tait-Bryan angles (distinct axes).
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

/// Three-angle representation of 3D rotation using Euler angles.
/// 
/// **Context**: Euler angles provide an intuitive parameterization of 3D rotations
/// using three sequential rotations about coordinate axes. This representation
/// is widely used in crystallography, molecular dynamics, and engineering
/// applications for specifying particle orientations.
/// 
/// **How it Works**: Stores three rotation angles (alpha, beta, gamma) in degrees
/// that define sequential rotations according to the specified Euler convention.
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Euler {
    pub alpha: f32,
    pub beta: f32,
    pub gamma: f32,
}

impl Euler {
    /// Creates a new Euler angle specification.
    /// 
    /// **Context**: Euler angles are commonly specified from external sources
    /// or user input and need validation and storage in a structured format.
    /// 
    /// **How it Works**: Simple constructor storing the three angles.
    pub fn new(alpha: f32, beta: f32, gamma: f32) -> Self {
        Self { alpha, beta, gamma }
    }
    
    /// Computes the 3D rotation matrix from Euler angles according to the specified convention.
    /// 
    /// **Context**: Converting Euler angles to rotation matrices enables geometric
    /// transformations of particle coordinates and beam directions. Different
    /// Euler conventions require different matrix formulations.
    /// 
    /// **How it Works**: Implements the closed-form rotation matrix for each
    /// Euler convention by composing the three sequential rotations.
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

/// Complete orientation specification combining scheme and convention.
/// 
/// **Context**: Orientation analysis requires both the orientation scheme
/// (how orientations are generated) and the Euler convention (how angles
/// are interpreted), providing a complete specification for orientation
/// handling throughout the simulation.
/// 
/// **How it Works**: Pairs an orientation scheme with an Euler convention
/// to enable consistent orientation processing.
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Orientation {
    pub scheme: Scheme,
    pub euler_convention: EulerConvention,
}

/// Generated orientation set for multi-orientation simulations.
/// 
/// **Context**: Multi-orientation simulations require concrete sets of Euler
/// angles for each simulation run. This structure provides the generated
/// orientations in a consistent format regardless of the generation scheme.
/// 
/// **How it Works**: Stores the generated Euler angle tuples and count for
/// use in orientation averaging calculations.
#[derive(Debug, Clone, PartialEq)]
pub struct Orientations {
    pub num_orientations: usize,
    pub eulers: Vec<(f32, f32, f32)>,
}

impl Orientations {
    /// Generates orientation sets according to the specified scheme.
    /// 
    /// **Context**: Different orientation schemes require different generation
    /// algorithms. This factory method provides a unified interface for
    /// creating orientation sets from scheme specifications.
    /// 
    /// **How it Works**: Dispatches to appropriate generation methods based
    /// on scheme type, handling random number seeding for reproducibility.
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

    /// Creates an orientation set from explicit angle lists.
    /// 
    /// **Context**: Discrete orientation schemes require validation that
    /// all angle lists have consistent lengths and contain valid values.
    /// 
    /// **How it Works**: Validates input consistency and packages the
    /// angles into tuples for simulation use.
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

    /// Generates uniformly distributed random orientations on the rotation group.
    /// 
    /// **Context**: Uniform random orientation averaging requires proper sampling
    /// from the rotation group SO(3). Naive uniform sampling in Euler angles
    /// creates bias, requiring careful distribution design.
    /// 
    /// **How it Works**: Uses proper uniform distribution on SO(3) by sampling
    /// alpha and gamma uniformly on [0, 2π] and beta with distribution proportional
    /// to sin(beta) to ensure uniform coverage of the rotation group.
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
