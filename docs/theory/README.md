# Theory notes

- `stable_fraunhofer_edge_sum.typ` — derivation of the numerically stable edge sum for the Fraunhofer diffraction factor of a polygonal aperture. Covers:
  - equivalence with the current clamped implementation in `diff2.rs`
  - the two exact single-coordinate Green's-theorem forms
  - the small-argument series of the edge function
  - the analytic small-(kxx, kyy) limit (signed area plus first moments)
  - branch selection, error budget, and f32 validation results
- Compile with `typst compile stable_fraunhofer_edge_sum.typ`.
- Validation script: `fraunhofer_limit_check.py` in this directory (needs numpy; e.g. run with the venv in `examples/bullet_rosette_aggregate/.venv`).
