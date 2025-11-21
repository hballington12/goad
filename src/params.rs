use crate::result::GOComponent;
use rand_distr::num_traits::Pow;
use serde::Serialize;
use std::{
    collections::HashMap,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Debug, PartialEq, Clone, Serialize)]
pub struct Params {
    params: HashMap<(Param, GOComponent), f32>,
}

// all params must add linearly eg. asymmetry must be multiplied by scatt cross
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Serialize)]
pub enum Param {
    AsymmetryScatt, // asymmetry multiplied by scatt cross
    ScatCross,
    ExtCross,
    AlbedoExt, // albedo multiplied by ext cross
}

impl Params {
    pub fn set_param(&mut self, param: Param, component: GOComponent, value: f32) {
        self.params.insert((param, component), value);
    }
}

impl Div<f32> for Params {
    type Output = Params;
    fn div(self, rhs: f32) -> Self::Output {
        let mut params = self.clone();
        for ((param, component), value) in self.params.iter() {
            params.set_param(*param, *component, value / rhs);
        }
        params
    }
}

impl Pow<f32> for Params {
    type Output = Params;
    fn pow(self, rhs: f32) -> Self::Output {
        let mut params = self.clone();
        for (key, _) in self.params.iter() {
            params.params.entry(*key).or_default().pow(rhs);
        }
        params
    }
}

impl Mul<f32> for Params {
    type Output = Params;
    fn mul(self, rhs: f32) -> Self::Output {
        let mut params = self.clone();
        for (key, _) in self.params.iter() {
            *params.params.entry(*key).or_default() *= rhs;
        }
        params
    }
}

impl Mul<Params> for Params {
    type Output = Params;
    fn mul(self, other: Params) -> Self::Output {
        let mut params = self.clone();
        for (key, value) in other.params.iter() {
            *params.params.entry(*key).or_default() *= value;
        }
        params
    }
}

// care here because cannot just "add most things"
impl Add<Params> for Params {
    type Output = Params;
    fn add(self, other: Params) -> Self::Output {
        let mut params = self.clone();
        for (key, value) in other.params.iter() {
            *params.params.entry(*key).or_default() += value;
        }
        params
    }
}

// care here because cannot just "add most things"
impl Sub<Params> for Params {
    type Output = Params;
    fn sub(self, other: Params) -> Self::Output {
        let mut params = self.clone();
        for (key, value) in other.params.iter() {
            *params.params.entry(*key).or_default() -= value;
        }
        params
    }
}

impl Params {
    pub fn new() -> Self {
        Self {
            params: HashMap::new(),
        }
    }

    pub fn asymmetry(&self, component: &GOComponent) -> Option<f32> {
        if let (Some(asymmetry), Some(scatt_cross)) = (
            self.params.get(&(Param::AsymmetryScatt, *component)),
            self.params.get(&(Param::ScatCross, *component)),
        ) {
            Some(asymmetry / scatt_cross)
        } else {
            None
        }
    }

    pub fn albedo(&self, component: &GOComponent) -> Option<f32> {
        if let (Some(albedo_ext), Some(ext_cross)) = (
            self.params.get(&(Param::AlbedoExt, *component)),
            self.params.get(&(Param::ExtCross, *component)),
        ) {
            Some(albedo_ext / ext_cross)
        } else {
            None
        }
    }

    pub fn scatt_cross(&self, component: &GOComponent) -> Option<f32> {
        self.params.get(&(Param::ScatCross, *component)).copied()
    }

    pub fn ext_cross(&self, component: &GOComponent) -> Option<f32> {
        self.params.get(&(Param::ExtCross, *component)).copied()
    }
}
