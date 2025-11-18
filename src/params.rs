use crate::result::GOComponent;
use serde::Serialize;
use std::{
    collections::HashMap,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Debug, PartialEq, Clone, Serialize)]
pub struct Params {
    pub asymmetry: HashMap<GOComponent, f32>,
    pub scat_cross: HashMap<GOComponent, f32>,
    pub ext_cross: HashMap<GOComponent, f32>,
    pub albedo: HashMap<GOComponent, f32>,
}

impl Div<f32> for Params {
    type Output = Params;

    fn div(self, rhs: f32) -> Self::Output {
        let mut params = Params::new();
        if let Some(asymmetry) = self.asymmetry() {
            params.set_asymmetry(asymmetry / rhs);
        }
        if let Some(scat_cross) = self.scat_cross() {
            params.set_scat_cross(scat_cross / rhs);
        }
        if let Some(ext_cross) = self.ext_cross() {
            params.set_ext_cross(ext_cross / rhs);
        }
        if let Some(albedo) = self.albedo() {
            params.set_albedo(albedo / rhs);
        }
        params
    }
}

impl Mul for Params {
    type Output = Params;

    fn mul(self, other: Params) -> Self::Output {
        let mut params = Params::new();
        if let (Some(asymmetry), Some(other_asymmetry)) = (self.asymmetry(), other.asymmetry()) {
            params.set_asymmetry(asymmetry * other_asymmetry);
        }
        if let (Some(scatt), Some(other_scatt)) = (self.scat_cross(), other.scat_cross()) {
            params.set_scat_cross(scatt * other_scatt);
        }
        if let (Some(ext), Some(other_ext)) = (self.ext_cross(), other.ext_cross()) {
            params.set_ext_cross(ext * other_ext);
        }
        if let (Some(albedo), Some(other_albedo)) = (self.albedo(), other.albedo()) {
            params.set_albedo(albedo * other_albedo);
        }

        params
    }
}

impl Add for Params {
    type Output = Params;

    fn add(self, other: Params) -> Self::Output {
        let mut params = Params::new();
        if let (Some(asymmetry), Some(other_asymmetry)) = (self.asymmetry(), other.asymmetry()) {
            params.set_asymmetry(asymmetry + other_asymmetry);
        }
        if let (Some(scatt), Some(other_scatt)) = (self.scat_cross(), other.scat_cross()) {
            params.set_scat_cross(scatt + other_scatt);
        }
        if let (Some(ext), Some(other_ext)) = (self.ext_cross(), other.ext_cross()) {
            params.set_ext_cross(ext + other_ext);
        }
        if let (Some(albedo), Some(other_albedo)) = (self.albedo(), other.albedo()) {
            params.set_albedo(albedo + other_albedo);
        }

        params
    }
}

impl Sub for Params {
    type Output = Params;

    fn sub(self, other: Params) -> Self::Output {
        let mut params = Params::new();
        if let (Some(asymmetry), Some(other_asymmetry)) = (self.asymmetry(), other.asymmetry()) {
            params.set_asymmetry(asymmetry - other_asymmetry);
        }
        if let (Some(scatt), Some(other_scatt)) = (self.scat_cross(), other.scat_cross()) {
            params.set_scat_cross(scatt - other_scatt);
        }
        if let (Some(ext), Some(other_ext)) = (self.ext_cross(), other.ext_cross()) {
            params.set_ext_cross(ext - other_ext);
        }
        if let (Some(albedo), Some(other_albedo)) = (self.albedo(), other.albedo()) {
            params.set_albedo(albedo - other_albedo);
        }

        params
    }
}

impl Params {
    pub fn new() -> Self {
        Self {
            asymmetry: HashMap::new(),
            scat_cross: HashMap::new(),
            ext_cross: HashMap::new(),
            albedo: HashMap::new(),
        }
    }

    /// Get asymmetry for Total component (backwards compatibility)
    pub fn asymmetry(&self) -> Option<f32> {
        self.asymmetry.get(&GOComponent::Total).copied()
    }

    /// Set asymmetry for Total component (backwards compatibility)
    pub fn set_asymmetry(&mut self, value: f32) {
        self.asymmetry.insert(GOComponent::Total, value);
    }

    /// Get scat_cross for Total component (backwards compatibility)
    pub fn scat_cross(&self) -> Option<f32> {
        self.scat_cross.get(&GOComponent::Total).copied()
    }

    /// Set scat_cross for Total component (backwards compatibility)
    pub fn set_scat_cross(&mut self, value: f32) {
        self.scat_cross.insert(GOComponent::Total, value);
    }

    /// Get ext_cross for Total component (backwards compatibility)
    pub fn ext_cross(&self) -> Option<f32> {
        self.ext_cross.get(&GOComponent::Total).copied()
    }

    /// Set ext_cross for Total component (backwards compatibility)
    pub fn set_ext_cross(&mut self, value: f32) {
        self.ext_cross.insert(GOComponent::Total, value);
    }

    /// Get albedo for Total component (backwards compatibility)
    pub fn albedo(&self) -> Option<f32> {
        self.albedo.get(&GOComponent::Total).copied()
    }

    /// Set albedo for Total component (backwards compatibility)
    pub fn set_albedo(&mut self, value: f32) {
        self.albedo.insert(GOComponent::Total, value);
    }
}
