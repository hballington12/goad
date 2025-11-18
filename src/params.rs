use crate::result::GOComponent;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, PartialEq, Clone, Serialize)]
pub struct Params {
    pub asymmetry: HashMap<GOComponent, f32>,
    pub scat_cross: HashMap<GOComponent, f32>,
    pub ext_cross: HashMap<GOComponent, f32>,
    pub albedo: HashMap<GOComponent, f32>,
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
