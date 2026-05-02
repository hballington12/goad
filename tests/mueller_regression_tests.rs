//! Mueller-matrix regression tests.
//!
//! Each case is declared once via `mueller_case! { ... }` in `tests/cases.rs`
//! and expands to a `#[test]` that compares the live result against a frozen
//! reference file in `tests/test_data/`. The same `cases.rs` is re-included by
//! `examples/regen_mueller_refs/main.rs`, which regenerates the reference
//! files when the physics intentionally changes.

#[macro_use]
pub mod helpers;

mod cases;
