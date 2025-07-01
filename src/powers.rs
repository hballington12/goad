//! Energy conservation tracking for electromagnetic simulations.
//!
//! This module provides comprehensive power tracking throughout electromagnetic
//! beam propagation to validate energy conservation and identify sources of
//! numerical error. It tracks all power contributions including absorbed power,
//! output power, and various truncation mechanisms that limit simulation accuracy.
//!
//! The power tracking system provides:
//! - Input and output power accounting
//! - Material absorption tracking
//! - Truncation power categorization by mechanism
//! - Conservation validation and error analysis
//! - Operator overloading for power combination
//! - Formatted output for debugging and validation
//!
//! # Power Budget Components
//!
//! - Input power: Total incident electromagnetic energy
//! - Output power: Far-field scattered power
//! - Absorbed power: Energy absorbed in materials
//! - Truncation categories: Power lost to various numerical limits
//! - Missing power: Unaccounted energy indicating numerical errors

use std::{fmt, ops::*};

/// Power conservation tracking for electromagnetic beam propagation.
/// 
/// **Context**: Energy conservation is a fundamental check for simulation accuracy.
/// In geometric optics with beam truncation, perfect conservation is impossible,
/// but tracking all power contributions enables validation and convergence analysis.
/// Each truncation mechanism contributes to the power budget.
/// 
/// **How it Works**: Tracks input power, output power reaching far field, absorbed
/// power in materials, and various truncation contributions. The missing() method
/// computes unaccounted power, which should approach zero for well-converged simulations.
#[derive(Debug, Copy, Clone, PartialEq)]
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
    /// Creates a new power tracking structure with zero initial values.
    /// 
    /// **Context**: Power tracking begins with zero values and accumulates
    /// contributions throughout the simulation.
    /// 
    /// **How it Works**: Initializes all power components to zero.
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

    /// Computes power missing from the conservation budget.
    /// 
    /// **Context**: Perfect power conservation would show zero missing power.
    /// Non-zero values indicate numerical errors, unconverged simulations,
    /// or untracked truncation mechanisms.
    /// 
    /// **How it Works**: Subtracts all tracked power components from input
    /// power to find the unaccounted remainder.
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
