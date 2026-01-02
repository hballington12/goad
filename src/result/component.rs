use serde::Serialize;

#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Serialize)]
pub enum GOComponent {
    Total,
    Beam,
    ExtDiff,
}

impl GOComponent {
    /// Returns the file extension for the given GOComponent.
    pub fn file_extension(&self) -> &'static str {
        match self {
            GOComponent::Total => "",
            GOComponent::Beam => "beam",
            GOComponent::ExtDiff => "ext",
        }
    }
}
