use std::{fs::File, io::BufReader, path::Path};

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
