use super::component::GOComponent;
use super::mueller::MuellerMatrix;
use super::scatt_result::ScattResult1D;

/// Helper function to integrate over theta for a specific component with custom weighting
pub fn integrate_theta_weighted_component<F>(
    field_1d: &[ScattResult1D],
    component: GOComponent,
    weight_fn: F,
) -> f32
where
    F: Fn(f32, f32) -> f32, // (theta_radians, s11_value) -> weighted_value
{
    let sum: f32 = field_1d
        .iter()
        .map(|result| {
            let mueller = match component {
                GOComponent::Total => result.mueller_total,
                GOComponent::Beam => result.mueller_beam,
                GOComponent::ExtDiff => result.mueller_ext,
            };
            let s11 = mueller.s11();
            let theta_rad = result.bin.center.to_radians();
            let bin_width = result.bin.width().to_radians();
            weight_fn(theta_rad, s11) * bin_width
        })
        .sum();

    sum
}
