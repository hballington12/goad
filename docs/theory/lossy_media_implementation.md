# Phase D implementation plan: lossy-media rework

Status: draft for review, 2026-08-17. Branch: `lossy-media` (off
`python`, carrying the TIR branch fix). Theory:
`lossy_media_beam_tracing.typ`; validation targets:
`examples/lossy_media_validation/`.

## Regression gates (run after every stage)

- G1, lossless bit-identity: a fixed-orientation hex run at `n' = 0`
  must be byte-identical to the pre-rework output recorded at stage 0.
- G2, forced homogeneity: with the beam wave state forced to
  `w = n_tilde k_hat` (a debug switch kept until stage 6), absorbing
  runs must reproduce the pre-rework absorbing outputs exactly. This
  separates plumbing bugs from intended physics changes.
- G3, harness parity: Rust unit tests pinned to golden values computed
  by the Python harness (transfer, Fresnel, complex-argument edge sum).
- `cargo test --lib --tests` throughout (untracked farfield example
  dirs do not compile against this base).

## Stage 0: baselines

- Record reference outputs (2D Mueller grids, results.json) for: hex
  `n' = 0`, `1e-9`, `0.01`, `0.1` at the C1 settings, into a
  reference-data directory used by G1/G2 integration tests.

## Stage 1: state plumbing, no behaviour change

- `Field` gains the wave state: apparent index `n_app: f32` and
  absorption vector `absorption: Vector3<f32>` (the pair `(N, A)`,
  equivalently `w = N e + iA` with `e = prop`). Constructors set the
  homogeneous values (`n_app = n.re`, `A = n.im * prop` when created
  inside a medium; `A = 0` in the lossless host). Debug assertion:
  dispersion residual `|w.w - n_tilde^2|` below tolerance.
- `Beam::refr_index` stays: the material index is still needed for
  Fresnel and the dispersion check.
- Polarisation basis (REVISED at implementation time): rather than
  changing the type of `e_perp` across ~140 call sites, the bivector
  is represented additively: `e_perp` remains the real part and an
  `e_perp_im: Vector3<f32>` (zero for homogeneous waves) joins it in
  stage 2, introduced together with the transfer that first produces
  nonzero values. Real-basis consumers are exact whenever
  `e_perp_im == 0`, which structurally provides the lossless fast
  path: the far-field hot path keeps its real rotation factors
  (`real_times_complex`, 40 vs 96 multiplies) with one dispatch per
  beam, not per bin. Benchmark at stage 6.
- Gate: G1 (nothing observable changed anywhere).

## Stage 2: interface transfer and Fresnel

- New `transfer` in `snell.rs` (replacing `get_theta_t` as the
  primary API): tangential continuity plus
  `gamma = sqrt(n2^2 - w_t.w_t)` with the continuation branch
  (`re + im > 0`). Returns transmitted `w2`, `gamma`, and the
  incident normal component `kappa_i`. Keep `get_theta_t` for the
  transition and delete it at stage 6.
- `fresnel.rs`: `refl`/`refr` variants taking `(kappa_i, kappa_t,
  n1, n2)` (Chang eq. 14 form). The homogeneous-incidence wrappers
  remain and are tested equal to the general form (the reduction
  identity).
- `beam.rs`:
  - propagation directions from `Re(w2)` (transmitted) and the
    mirrored `w` (reflected); reflected `A` mirrors too.
  - TIR restructured: no `(n2.re/n1.re).asin()` branch. A transmitted
    beam is spawned when its normal flux `Re(gamma)` exceeds zero
    (with the existing power threshold doing the real gating);
    lossless above-critical gives exactly `Re(gamma) = 0` and
    reproduces today's deterministic TIR behaviour, preserving
    `tir_count` semantics via the zero-flux test.
- Gates: G1, G2, G3 (transfer and Fresnel golden values, including
  the TIR continuity test at `n' -> 0`).

## Stage 3: propagation and the amplitude ramp

- `Field::propagate(displacement: Vector3, k)` with the vector
  displacement: phase `+= k N e . d` (f64 accumulate as today),
  amplitude `*= exp(-k A . d)`.
- Same change in the re-referencing sites: `earcut`,
  `remainders_to_beams`, and the sub-polygon centroid step in
  `intersections_to_beams`.
- `absorbed_power` keeps the current midpoint tally in this stage
  (upgraded in stage 5) so G2 stays exact.
- Gates: G1, G2.

## Stage 4: far field

- `diff2.rs::polygon_ft` over complex `(p, q)`: `cis` for complex
  argument (real exponential times the existing fast cis), complex
  `phi` with the same `DIFFERENCE_THRESHOLD` on `|omega|`, moments
  branch unchanged (analytic in `p, q`), `kappa = max(|p|, |q|)` in
  modulus. `EdgeData` unchanged.
- The mismatch `kxx, kyy` sourced from the exiting beam's `w` (its
  imaginary part is `k A` in the aperture frame). Beams with `A = 0`
  produce exactly the real code path.
- Karczewski matrix: evaluate at the complex direction `w / n_med`
  (analytic continuation), with the real-direction evaluation kept
  behind the same `A == 0` dispatch. Validate the continued form
  against far-field quadrature in the harness before enabling
  (extend `spectral_branch_check.py` with a vector-field check).
- FOV cone: unchanged (ramp-free width; conservative).
- Gates: G1, G3 (edge-sum golden values at complex argument).

## Stage 5: power bookkeeping

- `Beam::power()` and `csa()` with the `N` flux factor and the ramp
  integral `D_eff` via `polygon_ft(i mu, i nu)`; queue ordering and
  the power threshold use it.
- `powers` gains the `interface_mixed` counter (`M` from Dupertuis II
  eqs. 56/63-66) booked per interface, and `absorbed` moves to the
  `D_eff`-difference ledger. Reporting: `input = output + absorbed +
  interface_mixed + truncation` within tolerance.
- Gates: G1; G2 retired here (absorbed-power semantics intentionally
  change); energy-closure integration test added instead.

## Stage 6: cleanup and validation

- Remove `get_theta_t` and the forced-homogeneity switch; delete dead
  Fresnel paths.
- Port the remaining harness checks as Rust tests where cheap.
- Re-run C1 as a four-way comparison (ADDA, old, branch-fix-only,
  full rework) and update `goad-bench/parts/branch-fix-comparison`.
- Benchmark: far-field kernel cost for lossless (must be unchanged)
  and absorbing runs (complex-basis overhead reported).
- Convergence spot check: orientation-averaged absorbing hex, powers
  ledger closure.

## Risks and open points

- Complex `e_perp` touches many code paths; stage 1 is deliberately
  behaviour-neutral so G1 isolates mistakes early.
- f32 dispersion drift over long bounce chains: monitored by the
  debug assertion; renormalise `w` against the dispersion relation on
  transfer if drift is observed.
- The near-critical crossover window (theory note): no special
  handling; beams there inherit the continuation branch. Revisit only
  if validation shows artefacts.
- Distortion/roughness features change facet normals per intersection;
  `A` transfer uses the same perturbed normals, no extra work
  expected, but the distortion test suite must pass at G1.
