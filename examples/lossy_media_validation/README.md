# Lossy media validation (Phase B)

Validation harness for the lossy-media beam tracing rework. Theory:
`docs/theory/lossy_media_beam_tracing.typ`; plan:
`docs/theory/lossy_media_plan.md`.

- `lossy.py` — shared physics module:
  - reduced wave vector `w = N e + i A`, interface transfer (tangential continuity + principal square root), reflection
  - generalised Fresnel (Chang eq. 14) and the current `fresnel.rs` formulas for comparison
  - Macke `get_theta_t` port from `snell.rs`
  - stable polygon edge sum `A(p, q)` at complex argument, f32-emulated, in both the plain and the vertex-difference form
  - dense f64 quadrature reference with complex arguments
- `ramped_edge_sum_check.py` — task 1: complex-argument edge sum vs f64 quadrature over a sweep of ramp strengths up to the power-threshold bound; also the `D_eff = A(i mu, i nu)` identity.
- `lossy_wedge.py` — task 2: Chang section 4 lossy wedge, reproducing Figs. 3-7 from the transfer rule alone; checks the TIR angle, no-TIR transmission, and the standard-overestimates-loss ordering. Figures go to `figs/`.
- `fresnel_checks.py` — task 3: homogeneous-incidence reduction identity (general == current `fresnel.rs`), lossless TIR limit, Chang (16)/(18)/(20) expansion scalings, s-pol flux balance.
- `interface_chain_check.py` — task 4: dispersion and tangential-continuity invariants over random interface chains, slab reversibility, first-entry and exit-ramp structure, no-TIR property, and the Macke-vs-transfer direction comparison.
- `spectral_branch_check.py` — the branch adjudicator: hard-aperture ramped beam decomposed into real-tangential-wavenumber plane waves (branch-free elementary Fresnel), transmitted component-wise and recombined; compares the exact footprint reflection coefficient and power budget against the single-wave branch candidates (principal, other root, continuation rule, Chang eq. 12 rule). Result: the continuation rule (`Re + Im > 0`) is the physical branch, accurate to ~5e-4 away from the critical angle; breakdown window is +-1-2 deg of critical and ramp exponents near the power-threshold bound. Also verifies Dupertuis II `R + T = 1 + M` to 1e-16.
- The branch decision led to a production fix in `src/fresnel.rs` (`cos_theta_t` continuation-rule branch): before it, any `n' > 0` flipped the TIR phase discontinuously (hex at `n' = 1e-9` flipped the sign of S34/S11 over theta 40-63 deg); after it, differences vs lossless scale linearly with `n'` and lossless output is bit-identical.

Run with:

- `uv run python fresnel_checks.py`
- `uv run python interface_chain_check.py`
- `uv run python ramped_edge_sum_check.py`
- `uv run python lossy_wedge.py`
