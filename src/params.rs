use crate::result::GOComponent;
use serde::Serialize;
use std::{
    collections::HashMap,
    ops::{Add, Div, Sub},
};

#[derive(Debug, PartialEq, Clone, Serialize)]
pub struct Params {
    asymmetry: HashMap<GOComponent, f32>,
    scat_cross: HashMap<GOComponent, f32>,
    ext_cross: HashMap<GOComponent, f32>,
    albedo: HashMap<GOComponent, f32>,
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

impl Add for Params {
    type Output = Params;

    fn add(self, other: Params) -> Self::Output {
        let mut params = Params::new();
        match (self.scat_cross(), other.scat_cross()) {
            (Some(scat_cross), Some(other_scat_cross)) => {
                params.set_scat_cross(scat_cross + other_scat_cross);
                match (self.asymmetry(), other.asymmetry()) {
                    (Some(asymmetry), Some(other_asymmetry)) => {
                        params.set_asymmetry(
                            (asymmetry * scat_cross + other_asymmetry * other_scat_cross)
                                / (scat_cross + other_scat_cross),
                        );
                    }
                    (None, Some(other_asymmetry)) => {
                        params.set_asymmetry(other_asymmetry);
                    }
                    (Some(asymmetry), None) => {
                        params.set_asymmetry(asymmetry);
                    }
                    _ => {}
                }
            }
            (None, Some(other_scat_cross)) => {
                params.set_scat_cross(other_scat_cross);
                match (self.asymmetry(), other.asymmetry()) {
                    (None, Some(other_asymmetry)) => {
                        params.set_asymmetry(other_asymmetry);
                    }
                    _ => {}
                }
            }
            (Some(scat_cross), None) => {
                params.set_scat_cross(scat_cross);
                match (self.asymmetry(), other.asymmetry()) {
                    (Some(asymmetry), None) => {
                        params.set_asymmetry(asymmetry);
                    }
                    _ => {}
                }
            }
            _ => {}
        }

        match (self.ext_cross(), other.ext_cross()) {
            (Some(ext_cross), Some(other_ext_cross)) => {
                params.set_ext_cross(ext_cross + other_ext_cross);
                match (self.albedo(), other.albedo()) {
                    (Some(albedo), Some(other_albedo)) => {
                        params.set_albedo(
                            (albedo * ext_cross + other_albedo * other_ext_cross)
                                / (ext_cross + other_ext_cross),
                        );
                    }
                    (None, Some(other_albedo)) => {
                        params.set_albedo(other_albedo);
                    }
                    (Some(albedo), None) => {
                        params.set_albedo(albedo);
                    }
                    _ => {}
                }
            }
            (None, Some(other_ext_cross)) => {
                params.set_ext_cross(other_ext_cross);
                match (self.albedo(), other.albedo()) {
                    (None, Some(other_albedo)) => {
                        params.set_albedo(other_albedo);
                    }
                    _ => {}
                }
            }
            (Some(ext_cross), None) => {
                params.set_ext_cross(ext_cross);
                match (self.albedo(), other.albedo()) {
                    (Some(albedo), None) => {
                        params.set_albedo(albedo);
                    }
                    _ => {}
                }
            }
            _ => {}
        }

        params
    }
}

impl Sub for Params {
    type Output = Params;

    fn sub(self, other: Params) -> Self::Output {
        let mut params = Params::new();
        match (self.scat_cross(), other.scat_cross()) {
            (Some(scat_cross), Some(other_scat_cross)) => {
                params.set_scat_cross(scat_cross - other_scat_cross);
                match (self.asymmetry(), other.asymmetry()) {
                    (Some(asymmetry), Some(other_asymmetry)) => {
                        params.set_asymmetry(
                            (asymmetry * scat_cross - other_asymmetry * other_scat_cross)
                                / (scat_cross + other_scat_cross),
                        );
                    }
                    (None, Some(other_asymmetry)) => {
                        params.set_asymmetry(other_asymmetry);
                    }
                    (Some(asymmetry), None) => {
                        params.set_asymmetry(asymmetry);
                    }
                    _ => {}
                }
            }
            (None, Some(other_scat_cross)) => {
                params.set_scat_cross(other_scat_cross);
                match (self.asymmetry(), other.asymmetry()) {
                    (None, Some(other_asymmetry)) => {
                        params.set_asymmetry(other_asymmetry);
                    }
                    _ => {}
                }
            }
            (Some(scat_cross), None) => {
                params.set_scat_cross(scat_cross);
                match (self.asymmetry(), other.asymmetry()) {
                    (Some(asymmetry), None) => {
                        params.set_asymmetry(asymmetry);
                    }
                    _ => {}
                }
            }
            _ => {}
        }

        match (self.ext_cross(), other.ext_cross()) {
            (Some(ext_cross), Some(other_ext_cross)) => {
                params.set_ext_cross(ext_cross - other_ext_cross);
                match (self.albedo(), other.albedo()) {
                    (Some(albedo), Some(other_albedo)) => {
                        params.set_albedo(
                            (albedo * ext_cross - other_albedo * other_ext_cross)
                                / (ext_cross + other_ext_cross),
                        );
                    }
                    (None, Some(other_albedo)) => {
                        params.set_albedo(other_albedo);
                    }
                    (Some(albedo), None) => {
                        params.set_albedo(albedo);
                    }
                    _ => {}
                }
            }
            (None, Some(other_ext_cross)) => {
                params.set_ext_cross(other_ext_cross);
                match (self.albedo(), other.albedo()) {
                    (None, Some(other_albedo)) => {
                        params.set_albedo(other_albedo);
                    }
                    _ => {}
                }
            }
            (Some(ext_cross), None) => {
                params.set_ext_cross(ext_cross);
                match (self.albedo(), other.albedo()) {
                    (Some(albedo), None) => {
                        params.set_albedo(albedo);
                    }
                    _ => {}
                }
            }
            _ => {}
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
    pub fn asymmetry(&self, component: &GOComponent) -> Option<f32> {
        self.asymmetry.get(&component).copied()
    }

    /// Set asymmetry for Total component (backwards compatibility)
    pub fn set_asymmetry(&mut self, component: &GOComponent, value: f32) {
        self.asymmetry.insert(*component, value);
    }

    /// Get scat_cross for Total component (backwards compatibility)
    pub fn scat_cross(&self, component: &GOComponent) -> Option<f32> {
        self.scat_cross.get(&component).copied()
    }

    /// Set scat_cross for Total component (backwards compatibility)
    pub fn set_scat_cross(&mut self, component: &GOComponent, value: f32) {
        self.scat_cross.insert(*component, value);
    }

    /// Get ext_cross for Total component (backwards compatibility)
    pub fn ext_cross(&self, component: &GOComponent) -> Option<f32> {
        self.ext_cross.get(&component).copied()
    }

    /// Set ext_cross for Total component (backwards compatibility)
    pub fn set_ext_cross(&mut self, component: &GOComponent, value: f32) {
        self.ext_cross.insert(*component, value);
    }

    /// Get albedo for Total component (backwards compatibility)
    pub fn albedo(&self, component: &GOComponent) -> Option<f32> {
        self.albedo.get(&component).copied()
    }

    /// Set albedo for Total component (backwards compatibility)
    pub fn set_albedo(&mut self, component: &GOComponent, value: f32) {
        self.albedo.insert(*component, value);
    }
}
