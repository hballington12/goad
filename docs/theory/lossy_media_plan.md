# Plan: beam tracing in lossy media

Status: draft for discussion, 2026-08-16.

Sources, all in agreement on the target physics:

- Chang, Walker, Hopcraft (2005), JQSRT 96, 327-341. Transcribed in
  `chang2005_absorbing_ray_tracing.md`. Cited below as "Chang (N)".
- The IITM PGOM Fortran reference code
  (`~/codes/IITM_PGOM_Computational_Program/PGOM_Aggregates`), which
  implements the Chang machinery in a polygon-beam tracer and cites
  Chang eq. (10) directly in `permittivity.f90`. Cited as "PGOM".
- The existing GOAD edge-sum note `stable_fraunhofer_edge_sum.typ`,
  whose outlook section already covers the linear-exponential aperture
  profile by analytic continuation.
- The GOAD reference paper (`~/Documents/work/post-doc/goad-bench/
  manuscript-submission/main.tex`), which fixes the notation the theory
  note should use and whose Appendix B provides ADDA absorbing
  benchmarks. Cited as "paper".

## 0. Notation

Adopt the paper's notation throughout the theory note; the beam tuple
`b = (F_b, k_hat_b, e_perp_b, S_b, phi_b, n_tilde_b)`, the propagation
factor `T(d)`, the near-field master equation
`S_R = F_R . R . T(d) . S_b`, the Karczewski entries `(a_b, c_b)`, and
the contour Fraunhofer factor `F_cal_b`. Two clashes with Chang must be
resolved up front:

- Chang's imaginary index `kappa` is the paper's `n'`. Use `n'`.
- Chang's apparent indices `(N, K)` collide with the paper's
  aperture-frame beam direction `K_hat = (K_x, K_y, K_z)`. Proposal:
  keep Chang's `N` for the apparent real index and rename Chang's `K`
  to `N'` (apparent imaginary index), mirroring the `n, n'` pairing.
  The pair `(N, N')` with inhomogeneity angle `alpha` then extends the
  beam tuple naturally.

## 1. Target physics

A beam inside a lossy medium is an inhomogeneous plane wave. Its state
beyond what GOAD currently tracks:

- Apparent indices `(N, K)` with `N^2 - K^2 = n^2 - kappa^2` and
  `N K cos(alpha) = n kappa` (Chang (3)), where `alpha` is the angle
  between the constant-phase normal `e` and the constant-amplitude
  normal `f`.
- The attenuation direction `f` is independent state. After entry from
  a lossless host, `f = -n_hat` of the entry facet, not the ray
  direction (PGOM `diff_refl.f90`).
- Equivalent formulation carried per beam: a real absorption vector `A`
  with amplitude proportional to `exp(-k0 A . r)` across the wavefront.
  `A` is what survives projection onto each new facet and is the
  natural GOAD representation (one covector per beam, no per-vertex
  state required until diffraction time).

Interface transfer rules:

- Generalised Snell: two laws, `N1 sin(theta_i) = N2 sin(theta_t)` and
  `K1 sin(psi_i) = K2 sin(psi_t)` (Chang (9)), with `N2` from the
  quartic Chang (10) via the invariants
  `Ns = N1 sin(theta_i)`, `Ks = K1 sin(psi_i)`, and the tangential
  azimuth `cos(phi)` between phase and attenuation directions. Closed
  form of the solution: PGOM `permittivity.f90:32-49`.
- Transmitted and reflected directions for both `e` and `f`:
  Chang (12).
- Fresnel: Chang (14), everything in terms of the complex normal
  wavevector component `k . n_hat = k0 (N cos(theta) + i K cos(psi))`.
  Reduces to the current `fresnel.rs` when the incident wave is
  homogeneous.
- No hard TIR in lossy media (Chang section 4). The solver form
  guarantees `sin(theta_t) <= 1`; the lossless evanescent case is a
  limit branch, not an angular test.

Path loss and power:

- Attenuation along the ray is `exp(-k0 K l cos(alpha))`, not
  `exp(-k0 kappa l)` (Chang (15)). The standard method overestimates
  absorption.
- Path length varies linearly across the beam cross-section, so the
  exit aperture carries an exact linear-exponential amplitude ramp.
- Flux factor is the apparent index `N`, not `n.re`
  (PGOM `refraction.f90:144-151`).
- Absorption-weighted effective cross-section
  `D_eff = integral over aperture of exp(-2 T(r)) dA`, computed exactly
  by a Stokes-theorem edge sum since `T` is linear
  (PGOM `refraction.f90:458-509`). Absorbed power per bounce is the
  difference of `D_eff` between successive interfaces.

Far field:

- The aperture integral with the ramp is the analytic continuation of
  the existing edge sum, `A_0(p - i mu, q - i nu)`, with `(mu, nu)` the
  projection of `A` onto the aperture plane
  (edge-sum note, outlook; PGOM `refraction.f90:511-720` implements
  exactly this with a complex-argument `phi` function).

## 2. Scope decision needed

Absorbing particle in a non-absorbing host is the primary target (ice
and dust at infrared wavelengths). The general two-lossy-media
formalism (Chang (10)) costs nothing extra at interior interfaces, so
we carry it. An absorbing *host* medium raises separate conceptual
issues in the far field (extinction of the scattered wave on the way
to infinity, cf. the Mie-in-absorbing-media literature) and is out of
scope for this round; the settings already accept a complex
`medium_refr_index`, so we must decide whether to reject `im != 0`
loudly or keep the current approximate behaviour with a warning.

## 3. Fallback requirement

The generalisation must collapse to the current physics symbolically
when `kappa = 0`, not just numerically:

- `K = 0`, `N = n`, `alpha` irrelevant, `A = 0`.
- Chang (10) reduces to real Snell; verified in the low-absorption
  expansions Chang (16)-(17).
- Chang (14) reduces to the current complex-cos-theta-t Fresnel,
  including the TIR phase behaviour (the current `fresnel.rs` is
  already the homogeneous-incident-wave special case).
- `D_eff` reduces to `area x cos(theta)`, the current `csa()`.
- The complex-argument edge sum at real argument is the current
  `polygon_ft` term by term.

Design intent: one general code path that collapses, rather than a
lossy/lossless branch at every site. A cheap `A == 0` fast path is
allowed only where profiling justifies it, and must be bit-compatible
or tolerance-tested against the general path.

The error budget for when the current physics suffices is quantified
in Chang section 3: direction and path loss deviate at second order in
`x = kappa / n`, Fresnel amplitudes at first order. Guidelines:
rigorous treatment needed for `x > 0.01` when polarisation matters;
the standard approach holds to about 1 percent for unpolarised
quantities up to `x ~ 0.07`.

## 4. Phases

### Phase A: theory note (typst)

Derive, in GOAD's beam-and-aperture language, self-contained:

1. Beam state: `(N, K, e, f)` or equivalently `(N, A)`; the relation
   between the two; what must be stored versus recomputed.
2. Interface transfer for `e`, `f`, `N`, `K` (Chang (9), (10), (12)),
   including the closed-form quartic solution and its branches
   (lossless, evanescent, normal incidence where `f_t = -n_hat`).
3. Generalised Fresnel (Chang (14)) and its reduction to `fresnel.rs`.
4. Per-beam attenuation: winding `A` through propagation
   (facet update `A <- A + (N_i - A.e) n_hat / (e.n_hat)`, PGOM
   `beam_clipping.f90:1873`), projection into the new facet plane,
   re-referencing on beam splitting.
5. Power: flux with `N`, `D_eff` edge sum, absorbed power as `D_eff`
   differences, and consistency with the optical theorem.
6. Far field: complex-argument edge sum. Redo the stability analysis
   of `stable_fraunhofer_edge_sum.typ` for complex `(p, q)`: the
   `DIFFERENCE_THRESHOLD` crossover, the small-argument moment limit,
   and the branch selection all need revisiting when the argument has
   an imaginary part (the difference form gains a real exponential
   ratio, which changes the cancellation analysis).
7. Symbolic fallback proofs (section 3 above), one per formula.

### Phase B: numerical validation harness (Python, before any Rust)

1. Extend `fraunhofer_limit_check.py`: f64 brute-force quadrature of
   the ramped aperture integral versus the complex-argument edge sum in
   f32, over a sweep of ramp strengths `|mu| r_max` from 0 to the
   strongest physically reachable value (set by `kappa` and particle
   size; establish that bound in Phase A).
2. Lossy wedge reproduction: implement Chang section 4 analytically
   (output angle, `(N, K)` curves, transmittances, attenuation
   coefficient `d`) and reproduce Figs. 3-7 qualitatively; this
   validates the interface-transfer chain end to end with no beam
   tracing involved.
3. Fresnel checks: Chang (18)-(20) low-absorption expansions as series
   checks against the general formula; Brewster-angle residual (20);
   lossless TIR phase agreement with current `fresnel.rs`.
4. Interface-chain invariant test: trace a closed polygon of
   interfaces and confirm `(Ns, Ks, cos phi)` bookkeeping conserves
   the tangential wavevector at every step.

### Phase C: cross-validation against PGOM and ADDA

STATUS (2026-08-17): C-ADDA (item 4) done and adopted as the canonical
benchmark; results in `goad-bench/parts/branch-fix-comparison/`
(three-way ADDA vs GOAD-old vs GOAD-fixed; the branch fix improves
sigma_abs to -0.04%, S34/S11 MAE from 0.141 to 0.097 at n' = 0.01; the
n' = 0.1 interior is optically dead so the fix is a no-op there and
the residual discrepancy is non-interior physics).

C-PGOM (items 1-3) attempted and PARKED, with these findings:

- The shipped example outputs are not reproducible: modern gfortran
  seeds RANDOM_NUMBER from entropy, so their single random orientation
  is unrecoverable, and PGOM aborts on an internal `eps > 1` error for
  some orientations (including our fixed test orientation for the
  non-absorbing config; the absorbing config completed).
- The controlled-comparison recipe (for a future attempt): patch
  main.f90 to hardcode BETA_ARR(1)/GAMMA_ARR(1) after the random
  draw, rebuild per index config with `make` (mpifort), run with
  `shape.obj` in cwd; output is 498 nonuniform theta rows,
  phi-integrated, cols theta + s11..s44.
- BLOCKER: the PGOM 8-column aggregate is built from interpenetrating
  columns (bbox overlaps of 100-300 units), generated internally by
  GENERATE_HEXAGG. PGOM handles overlap via projected-outline merging;
  GOAD's containment model requires non-overlapping closed volumes and
  silently books 100% of input power to `trnc_clip` on this geometry
  (found via a welded, planarized conversion of their obj export;
  polygon-soup welding and planarization scripts were needed and the
  recipe is recorded here). A true same-geometry comparison needs a
  mesh boolean union of the eight prisms, or PGOM run on a
  non-overlapping aggregate; both deferred.
- Side finding for GOAD robustness: a geometry whose closed shells
  overlap fails with all power in `trnc_clip` and no warning; the
  loader or powers report should flag this.

1. Run the PGOM 8-column aggregate example; it ships absorbing and
   non-absorbing reference outputs
   (`example_out_abs_8col.dat`, `example_out_nonabs_8col.dat`).
2. Run GOAD on the same geometry (`8col.obj`) at matched wavelength,
   index, and orientation set, before and after the lossy rework.
   Non-absorbing runs isolate baseline geometric differences;
   absorbing runs then isolate the physics delta.
3. Compare: phase function, extinction / scattering / absorption cross
   sections, single-scattering albedo, and the polarisation elements
   Chang flags as most sensitive (`S34`, `S43`).
4. Re-run the paper's Appendix B absorbing benchmarks against the
   existing ADDA references (hex column, `X = 63.4`, `n' = 0.01` and
   `0.1`, fixed orientation). The `n' = 0.1` case sits at
   `x = n'/n = 0.076`, right where Chang's error budget says the
   standard method starts to fail visibly, so a before/after
   comparison tests whether the rework recovers part of the reported
   backscatter disagreement (-45% and +45% in the paper's Table 1).

### Phase D: implementation plan (Rust, only after A-C hold)

Sketch, to be refined once the theory note is fixed:

- `Field`/`Beam`: add the absorption vector `A` (and apparent index
  `N` where `refr_index.re` is used as flux factor). Phase stays f64;
  `A` is f32 like the rest.
- `snell.rs`: replace `get_theta_t` with the `effectivem`-style solver
  taking `(Ns, Ks, cos phi, n2, kappa2)`; keep the Macke port around
  for regression comparison during transition.
- `fresnel.rs`: generalised coefficients from complex normal
  components.
- `beam.rs`: remove the `(n2.re/n1.re).asin()` TIR branches in favour
  of the general transfer (lossless TIR emerges as the limit); wind
  `A` in `propagate`, project on facet change, re-reference in
  `earcut` and `remainders_to_beams` alongside the existing phase
  re-referencing.
- `powers.rs` and `Beam::power()`: `N` as flux factor, `D_eff` in
  place of `csa` where absorption is on, absorbed power as `D_eff`
  differences replacing the current midpoint-intensity tally.
- `diff2.rs`: `polygon_ft` over complex `(p, q)`; `EdgeData` gains the
  aperture-plane ramp `(mu, nu)`; FOV logic reviewed (a strong ramp
  concentrates the effective aperture and widens the diffraction
  cone).
- Regression gate: full non-absorbing test suite must be unchanged
  within tolerance; new absorbing tests pinned to Phase B/C reference
  values.

## 5. Known open questions

- RESOLVED (Phase A): the beam carries the reduced wave vector
  `w = N e + i A`; per-vertex exponents are evaluated lazily. See
  `lossy_media_beam_tracing.typ`.
- RESOLVED (user decision, 2026-08-16): polarisation bases go complex
  (bivector), as PGOM carries; no real-basis approximation in the
  implementation target.
- Interaction with the external diffraction path
  (`BeamVariant::ExternalDiff`): those beams live in the host medium;
  with a lossless host they carry no ramp and are untouched.
- `absorbed_power` currently accumulates during tracing per
  intersection; the `D_eff`-difference ledger changes where that tally
  lives. Energy-conservation reporting (`powers`) must remain
  interpretable mid-transition.
- Distortion / roughness features: facet-level distortion changes the
  entry normal and therefore `A`; confirm no hidden assumption that
  attenuation is isotropic per shape.
- How far the f32 budget stretches: `T(r)` exponents are differences
  of large numbers for big particles; check whether the `A . r`
  evaluation needs the same f64 care as `phase`.
