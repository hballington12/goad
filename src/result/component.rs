//! GOComponent enum for distinguishing scattering contributions.

use serde::Serialize;

/// Geometric optics component type.
///
/// Distinguishes between total scattering, beam contribution only,
/// and external diffraction contribution only.
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Serialize)]
pub enum GOComponent {
    /// Total scattering (beam + external diffraction)
    Total,
    /// Beam contribution only
    Beam,
    /// External diffraction contribution only
    ExtDiff,
}

impl GOComponent {
    /// Returns the file extension suffix for the given GOComponent.
    pub fn file_extension(&self) -> &'static str {
        match self {
            GOComponent::Total => "",
            GOComponent::Beam => "beam",
            GOComponent::ExtDiff => "ext",
        }
    }
}
