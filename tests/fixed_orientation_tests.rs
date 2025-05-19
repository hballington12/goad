use goad::{
    bins::{self, BinningScheme},
    multiproblem::MultiProblem,
    orientation::Euler,
    problem::collect_mueller,
    settings,
};
use num_complex::Complex32;

use std::{fs::File, io::BufReader, path::Path};

// Tolerance for comparing Mueller matrix elements
const FRAC_TOL: f32 = 1e-4; // fractional error
const ABS_TOL: f32 = 1e2; // absolute error

pub fn compare_results(
    result: Vec<Vec<f32>>,
    reference: Vec<Vec<f32>>,
    frac_tolerance: f32,
    abs_tolerance: f32,
) -> Result<(), Box<dyn std::error::Error>> {
    for (r, ref_) in result.iter().zip(reference.iter()) {
        for (a, b) in r.iter().zip(ref_.iter()) {
            assert!(
                ((a - b) / a).abs() < frac_tolerance || (a - b).abs() < abs_tolerance,
                "value: {}, reference: {}, fractional error: {}, absolute error: {}",
                a,
                b,
                ((a - b) / a).abs(),
                (a - b).abs()
            );
        }
    }

    Ok(())
}

pub fn load_reference_mueller(filename: &str) -> Result<Vec<Vec<f32>>, Box<dyn std::error::Error>> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data")
        .join(filename);

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut data = Vec::new();

    // Read the file line by line
    for line in std::io::BufRead::lines(reader) {
        let line = line?;

        // Skip empty lines
        if line.trim().is_empty() {
            continue;
        }

        // Parse each value in the line, skipping the first 2 values
        let row: Vec<f32> = line
            .split_whitespace()
            .skip(2) // Skip theta, phi
            .filter_map(|s| s.parse::<f32>().ok())
            .collect();

        // Add the row to our data
        if !row.is_empty() {
            data.push(row);
        }
    }

    Ok(data)
}

#[test]
fn fixed_hex_30_30_30() {
    let mut settings = settings::load_default_config().unwrap();
    // Reduce binning for faster testing
    settings.binning = BinningScheme {
        scheme: bins::Scheme::Simple {
            num_theta: 19,
            num_phi: 19,
        },
    };
    settings.orientation = goad::orientation::Orientation {
        scheme: goad::orientation::Scheme::Discrete {
            eulers: vec![Euler::new(30.0, 30.0, 30.0)],
        },
        euler_convention: goad::orientation::EulerConvention::ZYZ,
    };

    let mut multiproblem = MultiProblem::new(settings);
    multiproblem.solve();

    let result = collect_mueller(&multiproblem.result.mueller);
    let reference = load_reference_mueller("fixed_hex_30_30_30_mueller_scatgrid").unwrap();
    compare_results(result, reference, FRAC_TOL, ABS_TOL).unwrap();
}

#[test]
fn fixed_hex_30_20_20() {
    let mut settings = settings::load_default_config().unwrap();
    // Reduce binning for faster testing
    settings.binning = BinningScheme {
        scheme: bins::Scheme::Simple {
            num_theta: 19,
            num_phi: 19,
        },
    };
    settings.orientation = goad::orientation::Orientation {
        scheme: goad::orientation::Scheme::Discrete {
            eulers: vec![Euler::new(30.0, 20.0, 20.0)],
        },
        euler_convention: goad::orientation::EulerConvention::ZYZ,
    };
    // Change the refractive index
    settings.particle_refr_index = vec![Complex32::new(1.3117, 0.1)];

    let mut multiproblem = MultiProblem::new(settings);
    multiproblem.solve();

    let result = collect_mueller(&multiproblem.result.mueller);
    let reference = load_reference_mueller("fixed_hex_30_20_20_mueller_scatgrid").unwrap();
    compare_results(result, reference, FRAC_TOL, ABS_TOL).unwrap();
}
