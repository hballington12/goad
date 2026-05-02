use goad::settings::Settings;
use std::{fs::File, io::BufReader, path::Path};

pub const FRAC_TOL: f32 = 1e-4; // fractional error tolerance for Mueller comparisons
pub const ABS_TOL: f32 = 1e4; // absolute error tolerance for Mueller comparisons

/// One Mueller-matrix regression case. Each case knows its name (which is also
/// the reference filename in `tests/test_data/`) and how to build the
/// `Settings` to run. Both the `#[test]` and the `regen_mueller_refs` binary
/// pick cases up via `inventory::iter::<MuellerCase>`.
pub struct MuellerCase {
    pub name: &'static str,
    pub build_settings: fn() -> Settings,
}
inventory::collect!(MuellerCase);

/// Declare a Mueller regression case. Expands to:
///   - a `fn $name() -> Settings` builder
///   - an `inventory::submit!` registry entry
///   - (under `cfg(test)` only) a `#[test] fn test_$name()` that runs and compares
///
/// Usage:
///     mueller_case! {
///         name: my_case_name,
///         settings: |s| {
///             s.zones = ...;
///             s.orientation = ...;
///         }
///     }
#[macro_export]
macro_rules! mueller_case {
    (name: $name:ident, settings: |$s:ident| $body:block $(,)?) => {
        #[allow(dead_code)]
        fn $name() -> goad::settings::Settings {
            let mut $s = goad::settings::load_default_config().unwrap();
            $body
            $s
        }
        inventory::submit! {
            $crate::helpers::MuellerCase {
                name: stringify!($name),
                build_settings: $name,
            }
        }
        ::paste::paste! {
            #[cfg(test)]
            #[test]
            fn [<test_ $name>]() {
                use goad::result::MuellerMatrix;
                let mut mp = goad::multiproblem::MultiProblem::new(None, Some($name()))
                    .expect("Failed to create MultiProblem");
                mp.solve();
                let result: Vec<Vec<f32>> = mp
                    .result
                    .zones
                    .full_zone()
                    .expect("No full zone found")
                    .field_2d
                    .iter()
                    .map(|f| f.mueller_total.to_vec())
                    .collect();
                let reference =
                    $crate::helpers::load_reference_mueller(stringify!($name)).unwrap();
                $crate::helpers::compare_results(
                    result,
                    reference,
                    $crate::helpers::FRAC_TOL,
                    $crate::helpers::ABS_TOL,
                )
                .unwrap();
            }
        }
    };
}

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

    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let row: Vec<f32> = line
            .split_whitespace()
            .skip(2) // skip theta, phi
            .filter_map(|s| s.parse::<f32>().ok())
            .collect();
        if !row.is_empty() {
            data.push(row);
        }
    }

    Ok(data)
}
