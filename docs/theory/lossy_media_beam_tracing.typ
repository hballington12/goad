#set page(margin: 2.5cm)
#set text(size: 11pt)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Inhomogeneous-wave beam tracing for lossy media
  ]

  #v(0.5em)
  GOAD theory notes
]

= Motivation and scope

GOAD currently handles absorption with three approximations. A beam
transmitted into a lossy constituent takes its direction from a real
effective-index Snell's law (the Macke 1996 port in `snell.rs`), its
Fresnel coefficients from a complex $cos theta_t$ (`fresnel.rs`), and its
attenuation from a scalar factor $exp(-k n' d)$ applied along the
propagation direction over the centroid-to-centroid distance $d$
(`Field::propagate`). The beam cross-section keeps a uniform amplitude, and
total internal reflection is decided by the real-index critical angle.

The rigorous picture is that the field inside a lossy constituent is an
_inhomogeneous plane wave_: the planes of constant phase and the planes of
constant amplitude have different normals. The definitive ray-level
treatment is Chang, Walker and Hopcraft (2005), transcribed in
`chang2005_absorbing_ray_tracing.md` and cited below as "Chang (N)". The
IITM PGOM Fortran code implements the same physics in a polygon-beam
tracer and serves as a reference implementation. This note derives the
generalised beam tracing in GOAD's own notation (the beam tuple, the
propagation factor $T$, and the aperture edge sum of
`stable_fraunhofer_edge_sum.typ`), and proves that every generalised
formula collapses to the current implementation in the lossless limit.

The scope is an absorbing particle in a lossless host medium. All interior
interfaces may separate two lossy media; only the host is assumed
transparent, which is what makes the far-field limit unambiguous.

= Notation and the reduced wave vector

Material indices follow the reference paper: $tilde(n) = n + n' i$ per
constituent, wavenumber $k = 2 pi \/ lambda$ in the host. Chang's apparent
indices $(N, K)$ are renamed $(N, N')$ here, mirroring the $(n, n')$
pairing, because $K$ collides with the aperture-frame beam direction
$hat(bold(K))$ of the far-field formulae.

A time-harmonic plane wave in a lossy medium has the complex wave vector
$bold(k)_c = k (N hat(bold(e)) + i N' hat(bold(f)))$, where
$hat(bold(e))$ is the unit normal to planes of constant phase,
$hat(bold(f))$ the unit normal to planes of constant amplitude, and
$N, N' >= 0$ are the apparent refractive indices. Define the
_reduced wave vector_

$ bold(w) = bold(k)_c / k = N hat(bold(e)) + i bold(A), quad bold(A) = N' hat(bold(f)), $ <wdef>

so that the field of a wave with amplitude $E_0$ at a reference point
$bold(r)_c$ is

$ bold(E)(bold(r)) = bold(E)_0 exp[i k bold(w) dot (bold(r) - bold(r)_c)]
  = bold(E)_0
    underbrace(exp[i k N hat(bold(e)) dot (bold(r) - bold(r)_c)], "phase")
    thin
    underbrace(exp[-k bold(A) dot (bold(r) - bold(r)_c)], "amplitude"). $ <wavefield>

The real vector $bold(A)$ is the _absorption vector_: the gradient of the
log-amplitude in units of $k$. Everything in this note is a statement
about how $bold(w)$ is created, transported, and consumed. The angle
bookkeeping of Chang (the pairs $theta, psi$, the inhomogeneity angle
$alpha$, the azimuth $phi.alt$ between tangential directions) never needs
to be materialised: it is all contained in the components of $bold(w)$.

The beam tuple of the reference paper,
$b = (F_b, hat(bold(k))_b, hat(bold(e))_perp^b, bold(S)_b, phi.alt_b, tilde(n)_b)$,
is extended by the wave state:

$ b = (F_b, hat(bold(k))_b, hat(bold(e))_perp^b, bold(S)_b, phi.alt_b, tilde(n)_b, bold(w)_b), quad hat(bold(k))_b = "Re" thin bold(w)_b / |"Re" thin bold(w)_b|. $ <tuple>

Storing $bold(w)_b$ (one complex 3-vector, or equivalently the real pair
$N_b hat(bold(e))_b$ and $bold(A)_b$) is the entire cost of the
generalisation at the state level. The propagation direction is no longer
independent data; it is the direction of $"Re" thin bold(w)_b$.

= The inhomogeneous plane wave

Substituting the plane-wave ansatz into the Helmholtz equation for a
medium of index $tilde(n)$ gives the single complex dispersion relation

$ bold(w) dot bold(w) = tilde(n)^2. $ <dispersion>

Expanding with @wdef, the real and imaginary parts are

$ N^2 - |bold(A)|^2 = n^2 - n'^2, quad N thin hat(bold(e)) dot bold(A) = n n', $ <dispersion-parts>

which are Chang (3) with $bold(A) = N' hat(bold(f))$ and
$hat(bold(e)) dot hat(bold(f)) = cos alpha$. Two consequences used
throughout:

+ The attenuation per unit length _along the ray_ is
  $k thin bold(A) dot hat(bold(e)) = k n n' \/ N$, which is Chang (15).
  It is always at most $k n'$, since $N >= n$: the standard method
  overestimates path loss.

+ In a lossless medium ($n' = 0$) either $bold(A) = 0$ (homogeneous
  wave) or $bold(A) perp hat(bold(e))$ with $|bold(A)|^2 = N^2 - n^2$
  (a transversely decaying wave; the evanescent wave is the limiting
  case with $N$ below the tangential invariant). Both occur in GOAD:
  the incident illumination is the first, a beam exiting an absorbing
  particle is generically the second.

A homogeneous wave in a lossy medium is the special case
$hat(bold(f)) = hat(bold(e))$, for which @dispersion-parts forces
$N = n$, $N' = n'$ and @wavefield reduces to the current beam field with
$bold(w) = tilde(n) hat(bold(k))$. The current implementation is
therefore the general formalism with $bold(w)_b$ _forced_ to
$tilde(n)_b hat(bold(k))_b$ for every beam.

= Interface transfer

Let $hat(bold(n))$ be the interface normal oriented from medium 1 into
medium 2, and split any vector into tangential and normal parts,
$bold(w) = bold(w)_t + w_n hat(bold(n))$. Matching the total field phase
and amplitude over the interface plane (Chang (8), applied to the full
complex vector) gives the entire generalised Snell's law in one line:

$ bold(w)_(t,r) = bold(w)_(t,2) = bold(w)_(t,1) =: bold(w)_t. $ <tangential>

The normal component of the transmitted wave follows from @dispersion in
medium 2:

$ bold(w)_2 = bold(w)_t + gamma hat(bold(n)), quad gamma = sqrt(tilde(n)_2^2 - bold(w)_t dot bold(w)_t), $ <transfer>

with the branch of the square root fixed by the _continuation rule_ of
@sec-branch below: of the two roots $plus.minus gamma$, take the one
with $"Re" thin gamma + "Im" thin gamma > 0$. The reflected wave stays
in medium 1 with its normal component negated:

$ bold(w)_r = bold(w)_t - w_(n,1) hat(bold(n)). $ <reflected>

Note that @reflected mirrors the absorption vector as well as the
propagation direction; the reflected beam inherits the incident
inhomogeneity, which is Chang (12) for $hat(bold(f))_r$.

== Equivalence with Chang's closed form

Writing $bold(w)_t = N_s hat(bold(t))_e + i K_s hat(bold(t))_f$ with
$N_s = N_1 sin theta_i$, $K_s = N'_1 sin psi_i$ and
$hat(bold(t))_e dot hat(bold(t))_f = cos phi.alt$, the argument of the
root in @transfer is

$ gamma^2 = underbrace((n_2^2 - n_2'^2), M) - N_s^2 + K_s^2 + 2 i underbrace((n_2 n_2' - N_s K_s cos phi.alt), P). $ <gamma-sq>

Splitting $gamma = gamma_r + i gamma_i$ and using
$N_2^2 = |"Re" thin bold(w)_2|^2 = N_s^2 + gamma_r^2$ gives

$ N_2^2 = 1/2 [ N_s^2 + K_s^2 + M + sqrt((N_s^2 - K_s^2 - M)^2 + 4 P^2) ], quad N_2'^2 = N_2^2 - M, $ <quartic-solution>

which is the positive root of the quartic Chang (10) and the closed form
implemented in PGOM (`permittivity.f90`). The transmitted directions
$hat(bold(e))_2 = "Re" thin bold(w)_2 / N_2$ and
$hat(bold(f))_2 = "Im" thin bold(w)_2 / N_2'$ reproduce Chang (12). None
of these intermediates are needed by an implementation that carries
$bold(w)$; @transfer is one complex square root and one vector sum.

== Special cases

*First entry.* The incident illumination is homogeneous in a lossless
host: $bold(w)_1 = n_"med" hat(bold(k))_0$, so $bold(w)_t$ is real,
$K_s = 0$, and @quartic-solution reduces to the Liou form used at the
first interface in PGOM. The transmitted absorption vector is
$bold(A)_2 = gamma_i hat(bold(n))$: attenuation points _into the
surface_, not along the refracted ray. This is the single largest
qualitative difference from the current code, which sets
$bold(A) = n' hat(bold(e))$ implicitly.

*Interior interfaces.* The incident wave is already inhomogeneous;
@tangential and @transfer apply unchanged. This is where the current
complex-Snell treatment breaks down (see @sec-fresnel) and where Chang's
quartic is genuinely needed.

*Exit into the lossless host.* $tilde(n)_2 = n_"med"$ real, but
$bold(w)_t$ is complex because the interior wave carries attenuation.
By @dispersion-parts with $n' = 0$, the transmitted wave has
$bold(A)_2 perp hat(bold(e))_2$: a propagating wave with a purely
transverse amplitude ramp. This ramp is the exact boundary field that
the far-field stage must diffract (@sec-farfield); it is not an
artefact but the physical imprint of the oblique absorbing path on the
exiting wavefront.

== Branch selection <sec-branch>

The two roots $plus.minus gamma$ satisfy $gamma_r gamma_i = P$ with $P$
from @gamma-sq. When $P >= 0$ one root has both components
non-negative and there is nothing to choose. When $P < 0$, which is
_always_ the case at an exit face into a lossless host (there
$P = -N_s K_s cos phi.alt$), outgoing flux ($gamma_r > 0$) and
boundedness ($gamma_i > 0$) conflict, and no source in the literature
states a prescription: Dupertuis et al. leave the sign implicit,
Chang's "same sign as incident" rule (his eq. 12) produces a
$kappa_t$ that satisfies _neither_ root (it violates the medium-2
dispersion relation, and his wedge Figs. 5 and 7 inherit this), and
PGOM, like the current `fresnel.rs`, takes the principal root
throughout.

The correct prescription follows from what the inhomogeneous plane
wave _is_ for a finite beam: an analytic continuation of the beam's
spectrum off the real tangential axis. Continue
$gamma(k_t) = sqrt(tilde(n)_2^2 - k_t^2)$ from the beam's dominant
real tangential wavenumber to the complex $bold(w)_t$: below the
critical angle this lands on $gamma_r > 0$ (outgoing, with
$gamma_i < 0$ allowed: a transverse ramp tilt), above it on
$gamma_i > 0$ (decaying evanescent, with a small $gamma_r < 0$).
Uniformly: *take the root with $gamma_r + gamma_i > 0$*, equivalently
$arg gamma in (-pi \/ 4, 3 pi \/ 4]$. In code: take the principal
square root and negate it if the sum of its parts is negative.

This rule was validated against a branch-free spectral computation
(`spectral_branch_check.py`): a hard-aperture ramped beam decomposed
into real-$k$ plane waves, for which the elementary Fresnel problem
has no ambiguity, transmitted component-wise and recombined. Results
(ice-like $tilde(n) = 1.31 + 0.01 i$, 30-wavelength aperture): the
continuation root reproduces the exact footprint reflection
coefficient to $5 times 10^(-4)$ below critical and $4 times 10^(-4)$
in deep TIR, where the rejected root errs at order unity (it is the
complex conjugate phase) and Chang's rule errs at $5 times 10^(-2)$
to $2 times 10^(-1)$. The single-wave model itself breaks down only
(a) within roughly $plus.minus 1 degree$ of the critical angle, a
window set by the beam's spectral width, where every choice errs at
the $0.2$ level and the continuation root remains least bad, and
(b) at ramp exponents near the power-threshold bound, where accuracy
degrades to $tilde 0.2$. A third boundary of the same sliver: the
scalar rule is time-reversal consistent (a slab traversal
$1 arrow.r 2 arrow.r 1$ restores the wave exactly) if and only if the
incident normal component itself satisfies
$"Re" thin kappa + "Im" thin kappa > 0$; waves violating that have a
strongly anti-normal ramp at near-grazing incidence, within
$O(N'_s)$ of grazing, and lie inside the breakdown window already
counted.

Two consequences of the rule:

+ Above the critical angle the reflected-field coefficient satisfies
  $|r| > 1$ slightly (e.g. $1.026$ for ice at $n' = 0.01$,
  $theta_i = 60 degree$), which the spectral ground truth confirms is
  physical; the power bookkeeping is addressed in the energy section.
+ The principal root, in contrast, is _discontinuous_ with the
  lossless limit above critical: its $r$ is the complex conjugate of
  the correct one, so the TIR phase flips sign between $n' = 0$ and
  any $n' > 0$. This is a live artefact in the current code: a fixed
  orientation hex column at $n' = 10^(-9)$ (physically
  indistinguishable from lossless) flips the sign of $S_34 \/ S_11$
  across the TIR-dominated $40 degree$-$63 degree$ range with jumps
  up to $1.2$ and order-unity changes in $S_11$, identical at
  $n' = 10^(-6)$: the signature of a bistable branch, not of
  absorption. Every absorbing GOAD run to date, including the
  published ADDA benchmark comparisons, was computed on the wrong
  TIR branch.

*Total internal reflection.* From @quartic-solution,
$N_2 >= N_s$ always, with equality only when $P = 0$ and
$N_s^2 >= K_s^2 + M$. For a lossy medium 2 (or an inhomogeneous
incident wave with $P != 0$) the inequality is strict, $sin theta_(t,2) =
N_s \/ N_2 < 1$, and a transmitted wave always exists: _there is no
total internal reflection in lossy media_ (Chang section 4). The
classical TIR emerges as the lossless limit: $n'_2 = 0$, $K_s = 0$,
$N_s > n_2$ gives $gamma^2 < 0$, $gamma = i sqrt(N_s^2 - n_2^2)$, a pure
evanescent wave carrying no normal flux, and $|r| = 1$ from the Fresnel
formulae below. The beam-level consequence: the decision to spawn a
transmitted beam becomes an energy-threshold decision, not an angular
one. The current branches `theta_i > (n2.re/n1.re).asin()` in `beam.rs`
should be replaced by the general transfer, with the lossless
evanescent case suppressed by the existing power threshold (its
transmitted flux is exactly zero) or by an explicit
$"Re" thin gamma = 0$ test for determinism.

= Generalised Fresnel coefficients <sec-fresnel>

Define the complex normal components $kappa = bold(w) dot hat(bold(n))$
for the incident and transmitted waves,

$ kappa_i = bold(w)_1 dot hat(bold(n)) = N_1 cos theta_i + i N'_1 cos psi_i, quad kappa_t = gamma. $ <kappa-def>

Chang (14) gives the reflection and transmission coefficients entirely
in terms of these:

$ r_s = (kappa_i - kappa_t) / (kappa_i + kappa_t), quad
  t_s = (2 kappa_i) / (kappa_i + kappa_t), quad
  r_p = (tilde(n)_2^2 kappa_i - tilde(n)_1^2 kappa_t) / (tilde(n)_2^2 kappa_i + tilde(n)_1^2 kappa_t), quad
  t_p = (2 tilde(n)_1 tilde(n)_2 kappa_i) / (tilde(n)_2^2 kappa_i + tilde(n)_1^2 kappa_t). $ <fresnel-gen>

== Exact reduction to the current implementation

For a _homogeneous_ incident wave, $bold(w)_1 = tilde(n)_1 hat(bold(e))_1$,
so $kappa_i = tilde(n)_1 cos theta_i$ and $bold(w)_t dot bold(w)_t =
tilde(n)_1^2 sin^2 theta_i$; then

$ kappa_t = gamma = sqrt(tilde(n)_2^2 - tilde(n)_1^2 sin^2 theta_i) = tilde(n)_2 cos theta_t^c, $ <kappa-homog>

where $cos theta_t^c$ is precisely the complex cosine computed by
`fresnel.rs::cos_theta_t`. Substituting into @fresnel-gen and cancelling
$tilde(n)_1 tilde(n)_2$ from $r_p$ recovers the four formulae of the
reference paper's Appendix A term for term. Two corollaries:

+ The current Fresnel treatment is _exact_ at every below-critical
  interface whose incident wave is homogeneous, which includes the
  first entry into an absorbing particle, even between two lossy
  media. Above the critical angle the algebraic identity still holds
  but the principal-root $cos theta_t^c$ sits on the wrong branch for
  $n' > 0$ (@sec-branch), so exactness there requires the
  continuation-rule sign correction.

+ The error of the current treatment is confined to interior interfaces
  reached by beams that are already inhomogeneous. There, replacing
  $tilde(n)_2 cos theta_t^c$ by the true $gamma$ of @transfer, and
  $tilde(n)_1 cos theta_i$ by the true $kappa_i$, is the entire
  correction: the algebraic form of @fresnel-gen is unchanged.

== Polarisation basis

The exact field decomposition for inhomogeneous waves uses the complex
orthonormal pair of Chang (13),
$bold(h)_1 = hat(bold(n)) times bold(w) \/ sqrt((hat(bold(n)) times bold(w)) dot (hat(bold(n)) times bold(w)))$
and $bold(h)_2 = bold(w) times bold(h)_1 \/ tilde(n)$, which are
_bivectors_ (complex-valued vectors) whenever the tangential phase and
attenuation directions differ. The coefficients @fresnel-gen are the
scalars that act on the $(bold(h)_1, bold(h)_2)$ components.

GOAD's field structure assumes a real basis
$(hat(bold(e))_perp, hat(bold(e))_parallel)$. Since
$hat(bold(n)) times bold(w) = hat(bold(n)) times bold(w)_t$, the basis
$bold(h)_1$ is real exactly when the tangential parts of
$hat(bold(e))$ and $hat(bold(f))$ are parallel; in general the real-basis
description misses a component of relative size $O(|bold(A)| \/ N)$,
which is $O(n' \/ n)$. There are two consistent options:

+ _Exact_: promote $hat(bold(e))_perp$ (and the basis-rotation matrices)
  to complex vectors, as PGOM does (`VA`, `VB`, `VK` are complex). The
  rotation matrix construction of `Field::rotation_matrix` continues to
  apply with complex dot products.

+ _Real-basis approximation_: keep real bases built from
  $hat(bold(e))$, apply the scalar coefficients @fresnel-gen. The
  neglected terms are first order in $n' \/ n$, the same order as the
  Fresnel phase corrections themselves, so this option captures the
  dominant correction (the modified $kappa$ values) while leaving a
  residual of the same formal order with a smaller prefactor. Whether
  the prefactor is actually small must be established numerically in
  the validation phase, not assumed.

*Decision:* GOAD adopts option 1, the exact complex bivector bases,
following PGOM. The real-basis approximation is recorded above only as
an analysis tool (it quantifies which part of the correction the
$kappa$ substitution alone captures) and as a possible future fast
path, not as the implementation target.

= Beam propagation, splitting, and the amplitude ramp

== Propagation between interfaces

The current `Field::propagate(dist, k, n)` winds phase by
$k n d$ and multiplies amplitude by $exp(-k n' d)$, with $d$ the scalar
displacement along $hat(bold(k))_b$. The general rule, from @wavefield,
uses the vector displacement $bold(d) = bold(r)_c' - bold(r)_c$ between
the parent centroid and the child centroid:

$ phi.alt_b arrow.r phi.alt_b + k N thin hat(bold(e)) dot bold(d), quad
  bold(S)_b arrow.r bold(S)_b thin exp(-k bold(A) dot bold(d)). $ <propagate-gen>

The phase term is unchanged in structure (the scalar $d$ used today _is_
$hat(bold(e)) dot bold(d)$). The amplitude term is not: for an
inhomogeneous wave $bold(A) dot bold(d) != (bold(A) dot hat(bold(e)))
(hat(bold(e)) dot bold(d))$, so `propagate` must receive the full
displacement vector rather than a scalar distance. The same applies to
the phase re-referencing in `earcut` and `remainders_to_beams`.

Like the accumulated phase, the accumulated log-amplitude is a running
sum of increments; unlike the phase, its consumer is an exponential
rather than a trigonometric function, so it needs no modular reduction
and f32 relative error does not grow with the total. Carrying it as a
separate f64 accumulator alongside `phase` is nevertheless cheap and
removes any doubt; decide at implementation time.

== The ramp across a beam

Because $bold(S)_b$ is defined at the centroid and the field over $F_b$
is @wavefield, the amplitude across the beam polygon varies as
$exp(-k bold(A)_t dot (bold(r) - bold(r)_c))$, with $bold(A)_t$ the
projection of $bold(A)$ into the plane of $F_b$. No per-vertex state is
required: PGOM's per-vertex attenuation arrays (`TP`) are recoverable at
any time as $k bold(A) dot (bold(r)_j - bold(r)_c)$. GOAD should carry
only $bold(w)_b$ and evaluate vertex exponents lazily where they are
consumed (the effective cross-section and the far-field integral).

PGOM's facet update rule for its absorption vector,
$bold(A) arrow.l bold(A) + (N_i - bold(A) dot hat(bold(e))) hat(bold(n)) \/ (hat(bold(e)) dot hat(bold(n)))$,
is the combination "keep the tangential part, restore the dispersion
constraint $bold(A) dot hat(bold(e)) = n n' \/ N$ by adjusting the
normal part". This is algebraically identical to the transfer @transfer:
both determine $bold(A)_2$ from tangential continuity plus
@dispersion-parts. The two derivations agreeing is a useful
cross-validation of @transfer.

== Splitting

When a beam is clipped into sub-polygons or earcut into simple faces,
each child re-references $bold(S)$ and $phi.alt$ to its own centroid via
@propagate-gen with $bold(d)$ the centroid offset. The ramp is then
automatically consistent: parent and child describe the same field where
they overlap, because both are restrictions of the same inhomogeneous
plane wave. This is the same invariant the current code maintains for
phase, extended to amplitude.

= Power, effective cross-section, and the absorption ledger

== Energy flux

For the TE component the time-averaged Poynting vector of an
inhomogeneous wave is exactly

$ angle.l bold(S) angle.r = 1/2 sqrt(epsilon_0 / mu_0) thin N |E|^2 thin e^(-2 k bold(A) dot bold(r)) thin hat(bold(e)), $ <poynting-te>

that is, flux along $hat(bold(e))$ with the _apparent_ index $N$ as the
flux factor, not $n$. For the TM component there is an additional term
proportional to $"Im"[bold(E)^* (bold(E) dot bold(A))]$, of relative
size $O(|bold(A)| \/ N)$, which tilts the flux slightly out of
$hat(bold(e))$. The tracer adopts @poynting-te for both polarisations
and the TM correction joins the error budget (@sec-budget). The current
`Beam::power()` uses `refr_index.re`, i.e. $n$; the generalisation
replaces it by $N = |"Re" thin bold(w)|$.

== Effective cross-section

The power a beam carries through its polygon is the flux integrated over
the ramp:

$ P_b = 1/2 N_b thin ||bold(S)_b||^2 thin |hat(bold(e))_b dot hat(bold(n))_(F_b)| integral.double_(F_b) e^(-2 k bold(A)_t dot (bold(r) - bold(r)_c)) dif A, $ <deff>

which generalises `csa()` (the integral reduces to the polygon area when
$bold(A)_t = 0$). The integral is the Fourier transform of the polygon
indicator at a purely imaginary argument: in the aperture frame with
in-plane ramp components $(mu, nu) = 2 k (A_x, A_y)$,

$ integral.double_(F_b) e^(-(mu x + nu y)) dif x dif y = A(i mu, i nu), $ <deff-ft>

with $A(p, q)$ the edge-sum transform of
`stable_fraunhofer_edge_sum.typ` continued to complex arguments
(@sec-farfield). One evaluator therefore serves both the far field and
the power bookkeeping; PGOM's separate `PYEXP` contour form for
$D_"eff"$ is this same identity.

The beam-queue prioritisation and the power threshold should use
@deff, otherwise a strongly ramped beam is ranked by an amplitude its
bright edge does not represent.

== Interface energy bookkeeping for inhomogeneous incidence

A caveat discovered by the wedge validation (Phase B): when the incident
wave at an interface is inhomogeneous, $kappa_i$ is complex and
$1 - |r|^2$ is _not_ a flux ratio; $|r| > 1$ genuinely occurs in the
strongly absorbing regime (observed for $n' \/ n gt.tilde 0.3$ in the
former TIR zone of the Chang wedge), with energy balanced by
incident-reflected interference cross terms on the interface plane. The
per-beam fluxes downstream of the interface remain individually
well defined (the beams separate spatially), but the sum of reflected
and transmitted beam powers need not equal the incident beam power at
the interface itself. Consequences:

+ Per-interface power conservation checks must tolerate a discrepancy
  of order the cross terms, which are $O(("Im" kappa_i) / ("Re" kappa_i))$
  and negligible in GOAD's regime ($x lt.tilde 0.08$) below critical.
+ The `powers` ledger should book the discrepancy explicitly (its own
  counter, not `missing`) so it stays visible as absorption strengths
  grow. Dupertuis II supplies the exact per-interface balance
  $R + T = 1 + M$ (his eq. 60) with closed-form $M$; the counter is
  that $M$, which vanishes at first entry from a lossless host.
+ At above-critical bounces the single-wave power model errs more than
  the field model: the spectral ground truth for ice at $n' = 0.01$,
  $theta_i = 60 degree$ gives a true reflected power fraction $0.916$
  (with $0.027$ leaking through the sub-critical tail of the beam
  spectrum and the rest in interface interference), while $|r|^2$ of
  the correct branch is $1.053$. The coherent _field_ of the
  reflected beam is right to $4 times 10^(-4)$; the $tilde 5%$ power
  inflation per TIR bounce at this $n'$ is a model error to be booked
  against the $M$ counter, scaling linearly in $n'$.

== Absorption ledger

The power absorbed along one propagation leg is the difference between
@deff evaluated at the start and at the end of the leg (the amplitude
factors differ, and for oblique legs the ramp differs too). This
replaces the current midpoint tally
`absorbed_intensity * area * cos(theta_i) * n1.re` in
`intersections_to_beams`, which is the $bold(A) = n' hat(bold(e))$,
uniform-ramp special case. The ledger sums to the absorption
cross-section; extinction continues to come from the optical theorem on
the forward amplitude, which is unaffected because the host is lossless,
and energy-conservation checks compare the two as today.

= Far field with a ramped aperture <sec-farfield>

An exiting beam is a propagating wave in the host with a purely
transverse absorption vector (@dispersion-parts with $n' = 0$). Its
boundary field on the aperture polygon is the plane wave of the current
treatment multiplied by the real exponential ramp
$exp(-k bold(A)_t dot (bold(r) - bold(r)_c))$.

== Analytic continuation of the edge sum

In the aperture frame the Fraunhofer integral becomes

$ A_U (p, q) = integral.double_A e^(i (p x + q y)) e^(-(mu x + nu y)) dif x dif y = A(p + i mu, q + i nu), $ <continuation>

with $(mu, nu) = k (A_x, A_y)$ the in-plane ramp components. The stable
edge sum evaluates $A$ at complex arguments with no structural change:
the per-edge identity
$e^(i delta_j) phi(omega_j) = (e^(i delta_(j+1)) - e^(i delta_j)) \/ (i omega_j)$
is an identity of entire functions and holds for complex $delta, omega$;
the vertex factors become
$e^(i delta_j^r) e^(-delta_j^i)$, one `cis` and one real exponential per
vertex; and $phi$ keeps its Taylor branch with the threshold applied to
$|omega|$ (complex modulus). In the code, the continuation enters
through the wavevector mismatch: `kxx` and `kyy` gain imaginary parts
$k A_x$, $k A_y$, i.e. the mismatch is computed from $bold(w)$ of the
exiting beam rather than from $tilde(n) hat(bold(k))$.

== Stability at complex argument

The three error mechanisms of the real-argument analysis carry over with
one addition:

+ _Difference-form cancellation._ Unchanged in form: relative error
  $tilde epsilon_"f32" / |omega|$ near $omega = 0$; the existing
  threshold $|omega| < 0.3$ transfers with $|omega|$ the complex
  modulus, since the Taylor truncation bound $|omega|^6 / 5040$ is a
  bound in modulus.

+ _Small-argument limit._ The moment expansion
  $A = S + i (p M_x + q M_y) + O((kappa R)^2)$ is an expansion of an
  entire function and holds verbatim for complex $p, q$ with
  $kappa = max(|p|, |q|)$.

+ _Dynamic range (new)._ The vertex weights $e^(-delta_j^i)$ span a
  ratio $exp(2 k |bold(A)_t| r_max)$ across the polygon. The edge sum
  then computes a result that can be dominated by the bright side while
  f32 rounding is set by the largest terms; relative error grows
  linearly with the dynamic range ratio in the worst case. Physically a
  large ratio means the dim side of the beam is negligible, so the
  absolute error remains controlled, but this must be measured, not
  argued: the Phase B harness sweeps $k |bold(A)_t| r_max$ from 0 to
  the largest value reachable (bounded by the total attenuation a beam
  can suffer before the power threshold kills it, roughly
  $-ln("P"_"th")$ in exponent units).

== Field of view

The FOV cone half-angle currently scales with $lambda \/ r_max$. A
strong ramp concentrates the effective aperture into a region smaller
than $r_max$, widening the true diffraction lobe. A conservative
generalisation replaces $r_max$ by an effective radius
$r_"eff" = sqrt(A(i mu, i nu) \/ pi)$, monotonically shrinking with the
ramp. Whether the extra `polygon_ft` evaluation pays for itself is a
Phase D question; the conservative interim choice is to leave the FOV
cone at its ramp-free width, which can only over-include bins.

== Polarisation matrix

The Karczewski matrix is derived for a homogeneous plane wave through
the aperture. With a ramp, the aperture magnetic field acquires a term
$prop bold(A) times bold(E)$, out of phase with the electric surface
current, of relative size $|bold(A)| \/ N$. Since the Karczewski
entries are algebraic in the components of $hat(bold(K))$, the natural
exact treatment is the same analytic continuation as everywhere else:
evaluate them at the complex direction $bold(w) \/ n_"med"$ of the
exiting wave. The interim approximation keeps the real
$hat(bold(K)) = hat(bold(e))$ and accepts an $O(|bold(A)| \/ N)$ error
in the polarisation factor only (the scalar diffraction factor is exact
either way). Phase B validates the continued form against direct
quadrature of the vector integral before Phase D commits to either.

= Collapse to the current physics <sec-collapse>

Setting $n' = 0$ for every constituent (and the host):

+ _State._ @transfer preserves realness: real $bold(w)_t$, real
  $tilde(n)_2^2$, and below the critical angle real $gamma$, so
  $bold(A) = 0$ for every beam and $bold(w) = n hat(bold(k))$
  identically. Above the critical angle the continuation rule gives
  $gamma$ purely imaginary positive, the transmitted flux
  @poynting-te is zero, and no transmitted beam survives the power
  threshold: classical TIR, including the $|r| = 1$ phase-only
  reflection of the current code. Unlike the principal root, the
  continuation rule makes this limit continuous: $r(n' arrow.r 0^+)
  arrow.r r(0)$ on both sides of the critical angle (@sec-branch).
+ _Snell._ @transfer with real arguments is
  $sin theta_t = n_1 sin theta_i \/ n_2$; `get_theta_t` (Macke) is the
  homogeneous special case of @quartic-solution and agrees.
+ _Fresnel._ @kappa-homog shows @fresnel-gen equals the current
  formulae exactly whenever the incident wave is homogeneous, which
  with $bold(A) = 0$ is every interface.
+ _Propagation._ @propagate-gen with $bold(A) = 0$: amplitude factor 1,
  phase increment $k n hat(bold(e)) dot bold(d)$, the current winding.
+ _Power._ $N = n$, the ramp integral in @deff is the polygon area:
  `power()` and `csa()` unchanged.
+ _Far field._ $(mu, nu) = 0$: @continuation is the real edge sum,
  bit-identical code path if the complex evaluator short-circuits on a
  zero imaginary part, or tolerance-identical if it does not.

For absorbing media under the current approximation
($bold(w)$ forced to $tilde(n) hat(bold(k))$), the general formulae
reproduce the current behaviour exactly as well, which provides a second
regression axis: the new code with forced homogeneity must match the old
code on the absorbing benchmarks before the forcing is removed.

= Error budget and applicability <sec-budget>

From Chang section 3, with $x = n' \/ n$:

- Refraction angles and path loss: corrections are $O(x^2)$
  (Chang (16), (17)). Negligible below $x tilde 0.05$.
- Fresnel amplitudes: corrections are $O(x)$ (Chang (18)), and the
  induced ellipticity (the $S_34$-type elements) is directly
  proportional to $x$. This is the leading physical effect and the
  reason polarisation output is the most sensitive validation target.
- Reflectivities: $O(x^2)$ but with a large prefactor
  $4 n^3 \/ (n^2 - 1)^2$ (Chang (19)); at $n = 1.5$ a 1 percent
  deviation is reached already at $x tilde 0.034$.
- No-TIR transmission: qualitatively absent from the current code; the
  transmitted power in the former TIR zone is first order in $x$ and
  matters for beams trapped by repeated internal reflection.

Applicability guideline (Chang section 6): the current physics is
adequate for unpolarised quantities below $x tilde 0.07$ and for
polarised quantities below $x tilde 0.01$. The GOAD absorbing
benchmarks at $tilde(n) = 1.31 + 0.1 i$ sit at $x = 0.076$, outside
both.

Additional validity bounds from the spectral study (@sec-branch): the
single-inhomogeneous-wave beam model degrades within about
$plus.minus 1 degree$ of the critical angle (window set by the beam's
angular spectral width) and at ramp exponents approaching the
power-threshold bound; in both regimes the continuation root remains
the best available choice and the error is bounded near $0.2$ in the
reflection coefficient.

f32 concerns specific to this rework: the accumulated log-amplitude
(bounded by $-ln "P"_"th"$, no precision issue), the ramp dynamic range
in the edge sum (measured in Phase B), and $gamma$ near the branch
point (a $sqrt$ of a small complex number, benign in f32 but the
regression tests should pin the near-critical behaviour).

= Implementation implications

Deferred to `lossy_media_plan.md` Phase D; the state-level summary is:

- `Beam`/`Field` gain $bold(w)$ (equivalently $bold(A)$ plus $N$);
  `refr_index` remains the material index for Fresnel.
- `snell.rs` shrinks to @transfer; `fresnel.rs` keeps its algebraic
  form with $kappa$ inputs from $bold(w)$.
- `propagate`, `earcut`, `remainders_to_beams` take vector
  displacements.
- `power()`, `csa()`, the queue ordering, and the absorbed-power tally
  move to $N$ and the ramp integral @deff.
- `polygon_ft` continues to complex $(p, q)$; `EdgeData` unchanged;
  the mismatch `kxx`, `kyy` becomes complex, sourced from $bold(w)$.
- The polarisation basis goes complex (bivector), as in PGOM:
  `Field` carries a complex `e_perp`, and the basis-rotation
  construction extends with complex dot products.
