use ndarray::Array1;

/// Generate theta and phi combinations
pub fn generate_theta_phi_combinations() -> Vec<(f32, f32)> {
    let thetas = Array1::linspace(0.0, std::f32::consts::PI, 720).insert_axis(ndarray::Axis(1)); // Reshape to (50, 1)

    let phis = Array1::linspace(0.0, 2.0 * std::f32::consts::PI, 720).insert_axis(ndarray::Axis(0)); // Reshape to (1, 60)

    // Flatten the combinations of theta and phi into a 1D array of tuples
    thetas
        .iter()
        .flat_map(|&theta| phis.iter().map(move |&phi| (theta, phi)))
        .collect()
}
