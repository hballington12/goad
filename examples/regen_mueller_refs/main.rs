//! Regenerate Mueller-matrix regression reference files.
//!
//! Iterates the `MuellerCase` registry built up by `mueller_case! { ... }`
//! invocations in `tests/cases.rs` and writes a `mueller_scatgrid`-format file
//! per case into `tests/test_data/`. The output format is identical to what
//! goad's normal output writer produces, and is what `helpers::load_reference_mueller`
//! parses, so live-vs-reference comparison stays apples-to-apples.
//!
//! Usage:
//!     cargo run --release --example regen_mueller_refs                   # all cases
//!     cargo run --release --example regen_mueller_refs -- <name-substr>  # filter

use goad::{multiproblem::MultiProblem, result::MuellerMatrix};
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::PathBuf,
};

// Pull in helpers and case definitions from the test target so there is exactly
// one source of truth.
#[macro_use]
#[path = "../../tests/helpers.rs"]
mod helpers;
#[path = "../../tests/common/cases.rs"]
mod cases;

use helpers::MuellerCase;

fn main() {
    let filter = std::env::args().nth(1);
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data");
    fs::create_dir_all(&out_dir).expect("create test_data dir");

    let mut wrote_any = false;
    for case in inventory::iter::<MuellerCase> {
        if let Some(f) = &filter {
            if !case.name.contains(f) {
                continue;
            }
        }
        wrote_any = true;
        let settings = (case.build_settings)();
        let mut mp = MultiProblem::new(None, Some(settings))
            .expect("Failed to create MultiProblem");
        mp.solve();

        let zone = mp
            .result
            .zones
            .full_zone()
            .expect("regen case must produce a full zone");

        let path = out_dir.join(case.name);
        let mut w = BufWriter::new(File::create(&path).expect("create reference file"));
        for (bin, f2d) in zone.bins.iter().zip(&zone.field_2d) {
            write!(w, "{} {} ", bin.theta.center, bin.phi.center).unwrap();
            for v in f2d.mueller_total.to_vec() {
                write!(w, "{} ", v).unwrap();
            }
            writeln!(w).unwrap();
        }
        println!("wrote {}", path.display());
    }

    if !wrote_any {
        eprintln!(
            "no cases matched. registered cases:\n  {}",
            inventory::iter::<MuellerCase>()
                .into_iter()
                .map(|c| c.name)
                .collect::<Vec<_>>()
                .join("\n  ")
        );
        std::process::exit(1);
    }
}
