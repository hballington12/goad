use crate::convergence::Convergeable;
use crate::result::GOComponent;
use rand_distr::num_traits::Pow;
use serde::ser::{SerializeMap, Serializer};
use serde::Serialize;
use std::{
    collections::HashMap,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Debug, PartialEq, Clone)]
pub struct Params {
    params: HashMap<(Param, GOComponent), f32>,
}

impl Serialize for Params {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.params.len()))?;
        for ((param, component), value) in &self.params {
            let key = format!("{:?}_{:?}", param, component);
            map.serialize_entry(&key, value)?;
        }
        map.end()
    }
}

// Params are stored as raw values. Weighted averaging is handled by Convergeable trait.
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Serialize)]
pub enum Param {
    Asymmetry, // raw asymmetry parameter g
    Albedo,    // raw single scattering albedo
    ScatCross,
    ExtCross,
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
        for (key, value) in self.params.iter() {
            params.set_param(key.0, key.1, value.powf(rhs));
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
        self.params.get(&(Param::Asymmetry, *component)).copied()
    }

    pub fn albedo(&self, component: &GOComponent) -> Option<f32> {
        self.params.get(&(Param::Albedo, *component)).copied()
    }

    pub fn scatt_cross(&self, component: &GOComponent) -> Option<f32> {
        self.params.get(&(Param::ScatCross, *component)).copied()
    }

    pub fn ext_cross(&self, component: &GOComponent) -> Option<f32> {
        self.params.get(&(Param::ExtCross, *component)).copied()
    }

    /// Returns a weighted version of Params for convergence tracking.
    /// - asymmetry becomes asymmetry * scat_cross
    /// - albedo becomes albedo * ext_cross
    /// - scat_cross and ext_cross stay the same
    pub fn to_weighted(&self) -> Self {
        let mut result = self.clone();
        for component in [GOComponent::Total, GOComponent::Beam, GOComponent::ExtDiff] {
            // Asymmetry weighted by ScatCross
            if let (Some(asym), Some(sc)) =
                (self.asymmetry(&component), self.scatt_cross(&component))
            {
                result.set_param(Param::Asymmetry, component, asym * sc);
            }
            // Albedo weighted by ExtCross
            if let (Some(alb), Some(ec)) = (self.albedo(&component), self.ext_cross(&component)) {
                result.set_param(Param::Albedo, component, alb * ec);
            }
        }
        result
    }

    /// Returns a Params struct containing the weights for each field.
    /// - asymmetry slot contains scat_cross
    /// - albedo slot contains ext_cross
    /// - scat_cross and ext_cross slots contain 1.0
    pub fn weights(&self) -> Self {
        let mut result = Params::new();
        for component in [GOComponent::Total, GOComponent::Beam, GOComponent::ExtDiff] {
            // Asymmetry weight is ScatCross
            if let Some(sc) = self.scatt_cross(&component) {
                result.set_param(Param::Asymmetry, component, sc);
            }
            // Albedo weight is ExtCross
            if let Some(ec) = self.ext_cross(&component) {
                result.set_param(Param::Albedo, component, ec);
            }
            // ScatCross and ExtCross weights are 1.0
            if self.scatt_cross(&component).is_some() {
                result.set_param(Param::ScatCross, component, 1.0);
            }
            if self.ext_cross(&component).is_some() {
                result.set_param(Param::ExtCross, component, 1.0);
            }
        }
        result
    }
}

impl Convergeable for Params {
    fn zero_like(&self) -> Self {
        Params::new()
    }

    fn weighted_add(&self, other: &Self, w1: f32, w2: f32) -> Self {
        let mut result = Params::new();
        let total_weight = w1 + w2;

        for component in [GOComponent::Total, GOComponent::Beam, GOComponent::ExtDiff] {
            // ScatCross: simple weighted average by count
            if let (Some(s1), Some(s2)) =
                (self.scatt_cross(&component), other.scatt_cross(&component))
            {
                result.set_param(
                    Param::ScatCross,
                    component,
                    (s1 * w1 + s2 * w2) / total_weight,
                );
            } else if let Some(s1) = self.scatt_cross(&component) {
                result.set_param(Param::ScatCross, component, s1);
            } else if let Some(s2) = other.scatt_cross(&component) {
                result.set_param(Param::ScatCross, component, s2);
            }

            // ExtCross: simple weighted average by count
            if let (Some(e1), Some(e2)) = (self.ext_cross(&component), other.ext_cross(&component))
            {
                result.set_param(
                    Param::ExtCross,
                    component,
                    (e1 * w1 + e2 * w2) / total_weight,
                );
            } else if let Some(e1) = self.ext_cross(&component) {
                result.set_param(Param::ExtCross, component, e1);
            } else if let Some(e2) = other.ext_cross(&component) {
                result.set_param(Param::ExtCross, component, e2);
            }

            // Asymmetry: weighted by ScatCross
            if let (Some(g1), Some(g2), Some(sc1), Some(sc2)) = (
                self.asymmetry(&component),
                other.asymmetry(&component),
                self.scatt_cross(&component),
                other.scatt_cross(&component),
            ) {
                let weight1 = sc1 * w1;
                let weight2 = sc2 * w2;
                let new_g = (g1 * weight1 + g2 * weight2) / (weight1 + weight2);
                result.set_param(Param::Asymmetry, component, new_g);
            } else if let Some(g1) = self.asymmetry(&component) {
                result.set_param(Param::Asymmetry, component, g1);
            } else if let Some(g2) = other.asymmetry(&component) {
                result.set_param(Param::Asymmetry, component, g2);
            }

            // Albedo: weighted by ExtCross
            if let (Some(a1), Some(a2), Some(ec1), Some(ec2)) = (
                self.albedo(&component),
                other.albedo(&component),
                self.ext_cross(&component),
                other.ext_cross(&component),
            ) {
                let weight1 = ec1 * w1;
                let weight2 = ec2 * w2;
                let new_a = (a1 * weight1 + a2 * weight2) / (weight1 + weight2);
                result.set_param(Param::Albedo, component, new_a);
            } else if let Some(a1) = self.albedo(&component) {
                result.set_param(Param::Albedo, component, a1);
            } else if let Some(a2) = other.albedo(&component) {
                result.set_param(Param::Albedo, component, a2);
            }
        }

        result
    }

    fn mul_elem(&self, other: &Self) -> Self {
        self.clone() * other.clone()
    }

    fn div_elem(&self, other: &Self) -> Self {
        let mut result = self.clone();
        for (key, value) in other.params.iter() {
            if let Some(self_val) = result.params.get_mut(key) {
                *self_val /= value;
            }
        }
        result
    }

    fn add_elem(&self, other: &Self) -> Self {
        self.clone() + other.clone()
    }

    fn sub_elem(&self, other: &Self) -> Self {
        self.clone() - other.clone()
    }

    fn scale(&self, scalar: f32) -> Self {
        self.clone() * scalar
    }

    fn sqrt_elem(&self) -> Self {
        self.clone().pow(0.5)
    }

    fn to_weighted(&self) -> Self {
        self.to_weighted()
    }

    fn weights(&self) -> Self {
        self.weights()
    }
}
