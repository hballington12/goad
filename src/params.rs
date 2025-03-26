#[derive(Debug, PartialEq, Clone)]
pub struct Params {
    pub asymettry: Option<f32>,
    pub scat_cross: Option<f32>,
    pub ext_cross: Option<f32>,
    pub albedo: Option<f32>,
}

impl Params {
    pub fn new() -> Self {
        Self {
            asymettry: None,
            scat_cross: None,
            ext_cross: None,
            albedo: None,
        }
    }
}
