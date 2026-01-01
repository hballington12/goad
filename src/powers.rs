use crate::convergence::Convergeable;
use rand_distr::num_traits::Pow;
use serde::Serialize;
use std::{fmt, ops::*};

#[derive(Debug, Copy, Clone, PartialEq, Serialize)]
pub struct Powers {
    pub input: f32,       // near-field input power
    pub output: f32,      // near-field output power
    pub absorbed: f32,    // near-field absorbed power
    pub trnc_ref: f32,    // truncated power due to max reflections
    pub trnc_rec: f32,    // truncated power due to max recursions
    pub trnc_clip: f32,   // truncated power due to clipping
    pub trnc_energy: f32, // truncated power due to threshold beam power
    pub clip_err: f32,    // truncated power due to clipping error
    pub trnc_area: f32,   // truncated power due to area threshold
    pub trnc_cop: f32,    // truncated power due to cutoff power
    pub ext_diff: f32,    // external diffraction power
}

impl DivAssign<f32> for Powers {
    fn div_assign(&mut self, rhs: f32) {
        self.input /= rhs;
        self.output /= rhs;
        self.absorbed /= rhs;
        self.trnc_ref /= rhs;
        self.trnc_rec /= rhs;
        self.trnc_clip /= rhs;
        self.trnc_energy /= rhs;
        self.clip_err /= rhs;
        self.trnc_area /= rhs;
        self.trnc_cop /= rhs;
        self.ext_diff /= rhs;
    }
}

impl Div<f32> for Powers {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        Self {
            input: self.input / rhs,
            output: self.output / rhs,
            absorbed: self.absorbed / rhs,
            trnc_ref: self.trnc_ref / rhs,
            trnc_rec: self.trnc_rec / rhs,
            trnc_clip: self.trnc_clip / rhs,
            trnc_energy: self.trnc_energy / rhs,
            clip_err: self.clip_err / rhs,
            trnc_area: self.trnc_area / rhs,
            trnc_cop: self.trnc_cop / rhs,
            ext_diff: self.ext_diff / rhs,
        }
    }
}

impl Pow<f32> for Powers {
    type Output = Self;

    fn pow(self, rhs: f32) -> Self {
        Self {
            input: self.input.pow(rhs),
            output: self.output.pow(rhs),
            absorbed: self.absorbed.pow(rhs),
            trnc_ref: self.trnc_ref.pow(rhs),
            trnc_rec: self.trnc_rec.pow(rhs),
            trnc_clip: self.trnc_clip.pow(rhs),
            trnc_energy: self.trnc_energy.pow(rhs),
            clip_err: self.clip_err.pow(rhs),
            trnc_area: self.trnc_area.pow(rhs),
            trnc_cop: self.trnc_cop.pow(rhs),
            ext_diff: self.ext_diff.pow(rhs),
        }
    }
}

impl Mul<f32> for Powers {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        Self {
            input: self.input * rhs,
            output: self.output * rhs,
            absorbed: self.absorbed * rhs,
            trnc_ref: self.trnc_ref * rhs,
            trnc_rec: self.trnc_rec * rhs,
            trnc_clip: self.trnc_clip * rhs,
            trnc_energy: self.trnc_energy * rhs,
            clip_err: self.clip_err * rhs,
            trnc_area: self.trnc_area * rhs,
            trnc_cop: self.trnc_cop * rhs,
            ext_diff: self.ext_diff * rhs,
        }
    }
}

impl Mul for Powers {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            input: self.input * other.input,
            output: self.output * other.output,
            absorbed: self.absorbed * other.absorbed,
            trnc_ref: self.trnc_ref * other.trnc_ref,
            trnc_rec: self.trnc_rec * other.trnc_rec,
            trnc_clip: self.trnc_clip * other.trnc_clip,
            trnc_energy: self.trnc_energy * other.trnc_energy,
            clip_err: self.clip_err * other.clip_err,
            trnc_area: self.trnc_area * other.trnc_area,
            trnc_cop: self.trnc_cop * other.trnc_cop,
            ext_diff: self.ext_diff * other.ext_diff,
        }
    }
}

impl Add for Powers {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            input: self.input + other.input,
            output: self.output + other.output,
            absorbed: self.absorbed + other.absorbed,
            trnc_ref: self.trnc_ref + other.trnc_ref,
            trnc_rec: self.trnc_rec + other.trnc_rec,
            trnc_clip: self.trnc_clip + other.trnc_clip,
            trnc_energy: self.trnc_energy + other.trnc_energy,
            clip_err: self.clip_err + other.clip_err,
            trnc_area: self.trnc_area + other.trnc_area,
            trnc_cop: self.trnc_cop + other.trnc_cop,
            ext_diff: self.ext_diff + other.ext_diff,
        }
    }
}

impl Sub for Powers {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            input: self.input - other.input,
            output: self.output - other.output,
            absorbed: self.absorbed - other.absorbed,
            trnc_ref: self.trnc_ref - other.trnc_ref,
            trnc_rec: self.trnc_rec - other.trnc_rec,
            trnc_clip: self.trnc_clip - other.trnc_clip,
            trnc_energy: self.trnc_energy - other.trnc_energy,
            clip_err: self.clip_err - other.clip_err,
            trnc_area: self.trnc_area - other.trnc_area,
            trnc_cop: self.trnc_cop - other.trnc_cop,
            ext_diff: self.ext_diff - other.ext_diff,
        }
    }
}

impl AddAssign for Powers {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            input: self.input + other.input,
            output: self.output + other.output,
            absorbed: self.absorbed + other.absorbed,
            trnc_ref: self.trnc_ref + other.trnc_ref,
            trnc_rec: self.trnc_rec + other.trnc_rec,
            trnc_clip: self.trnc_clip + other.trnc_clip,
            trnc_energy: self.trnc_energy + other.trnc_energy,
            clip_err: self.clip_err + other.clip_err,
            trnc_area: self.trnc_area + other.trnc_area,
            trnc_cop: self.trnc_cop + other.trnc_cop,
            ext_diff: self.ext_diff + other.ext_diff,
        };
    }
}

impl Powers {
    pub fn new() -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            absorbed: 0.0,
            trnc_ref: 0.0,
            trnc_rec: 0.0,
            trnc_clip: 0.0,
            trnc_energy: 0.0,
            clip_err: 0.0,
            trnc_area: 0.0,
            trnc_cop: 0.0,
            ext_diff: 0.0,
        }
    }

    /// Returns a Powers struct with all fields set to 1.0 (for weights)
    pub fn ones() -> Self {
        Self {
            input: 1.0,
            output: 1.0,
            absorbed: 1.0,
            trnc_ref: 1.0,
            trnc_rec: 1.0,
            trnc_clip: 1.0,
            trnc_energy: 1.0,
            clip_err: 1.0,
            trnc_area: 1.0,
            trnc_cop: 1.0,
            ext_diff: 1.0,
        }
    }

    /// Returns the power unaccounted for.
    pub fn missing(&self) -> f32 {
        self.input
            - (self.output
                + self.absorbed
                + self.trnc_ref
                + self.trnc_rec
                // + self.trnc_clip
                + self.trnc_area
                + self.clip_err
                + self.trnc_cop
                + self.trnc_energy)
    }
}

impl fmt::Display for Powers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Powers:")?;
        writeln!(f, "  Input:            {:.6}", self.input)?;
        writeln!(f, "  Output:           {:.6}", self.output)?;
        writeln!(f, "  Absorbed:         {:.6}", self.absorbed)?;
        writeln!(f, "  Trunc. Refl:      {:.6}", self.trnc_ref)?;
        writeln!(f, "  Trunc. Rec:       {:.6}", self.trnc_rec)?;
        // writeln!(f, "  Trunc. Clip:   {:.6}", self.trnc_clip)?;
        writeln!(f, "  Clip Err:         {:.6}", self.clip_err)?;
        writeln!(f, "  Trunc. Energy:    {:.6}", self.trnc_energy)?;
        writeln!(f, "  Trunc. Area:      {:.6}", self.trnc_area)?;
        writeln!(f, "  Trunc. Cop:       {:.6}", self.trnc_cop)?;
        writeln!(f, "  Other:            {:.6}", self.missing())?;
        writeln!(f, "  External Diff:    {:.6}", self.ext_diff)
    }
}

impl Convergeable for Powers {
    fn zero_like(&self) -> Self {
        Powers::new()
    }

    fn weighted_add(&self, other: &Self, w1: f32, w2: f32) -> Self {
        // Simple weighted average by count for all power fields
        let total = w1 + w2;
        (*self * w1 + *other * w2) / total
    }

    fn mul_elem(&self, other: &Self) -> Self {
        *self * *other
    }

    fn div_elem(&self, other: &Self) -> Self {
        Self {
            input: self.input / other.input,
            output: self.output / other.output,
            absorbed: self.absorbed / other.absorbed,
            trnc_ref: self.trnc_ref / other.trnc_ref,
            trnc_rec: self.trnc_rec / other.trnc_rec,
            trnc_clip: self.trnc_clip / other.trnc_clip,
            trnc_energy: self.trnc_energy / other.trnc_energy,
            clip_err: self.clip_err / other.clip_err,
            trnc_area: self.trnc_area / other.trnc_area,
            trnc_cop: self.trnc_cop / other.trnc_cop,
            ext_diff: self.ext_diff / other.ext_diff,
        }
    }

    fn add_elem(&self, other: &Self) -> Self {
        *self + *other
    }

    fn sub_elem(&self, other: &Self) -> Self {
        *self - *other
    }

    fn scale(&self, scalar: f32) -> Self {
        *self * scalar
    }

    fn sqrt_elem(&self) -> Self {
        self.pow(0.5)
    }

    fn to_weighted(&self) -> Self {
        // Powers don't need special weighting
        self.clone()
    }

    fn weights(&self) -> Self {
        // All weights are 1.0 for powers
        Powers::ones()
    }
}
