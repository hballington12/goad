use super::Settings;

pub fn validate_config(config: &mut Settings) {
    assert!(
        config.beam_area_threshold_fac > 1e-5,
        "Beam area threshold factor must be greater than 1e-5"
    );
    assert!(config.wavelength > 0.0, "Wavelength must be greater than 0");
}
