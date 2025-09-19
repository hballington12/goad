use std::collections::HashMap;
use serde::Serialize;
use crate::result::GOComponent;

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
    
    /// Get scat_cross for Total component (backwards compatibility)
    pub fn scat_cross(&self) -> Option<f32> {
        self.scat_cross.get(&GOComponent::Total).copied()
    }
    
    /// Get ext_cross for Total component (backwards compatibility)
    pub fn ext_cross(&self) -> Option<f32> {
        self.ext_cross.get(&GOComponent::Total).copied()
    }
    
    /// Get albedo for Total component (backwards compatibility)
    pub fn albedo(&self) -> Option<f32> {
        self.albedo.get(&GOComponent::Total).copied()
    }
}
