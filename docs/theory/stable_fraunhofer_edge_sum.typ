#set page(margin: 2.5cm)
#set text(size: 11pt)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    A numerically stable edge sum for the Fraunhofer \ diffraction factor of a polygonal aperture
  ]

  #v(0.5em)
  GOAD theory notes
]

= Motivation

GOAD maps each near-field beam to the far field by evaluating a Fraunhofer
diffraction integral over the beam's polygonal aperture. The implementation in
`diff2.rs` reduces this integral to a sum over the polygon edges. Each edge
term in the current form contains the factor $1 slash (p q)$, where $p$ and
$q$ are the transverse components of the scattering wavevector mismatch. The
factors are singular whenever the observation direction approaches the beam
propagation direction in either transverse coordinate. The code guards these
singularities by clamping $|p|$ and $|q|$ to a floor value
(`KXY_EPSILON` $= 10^(-3)$).

The clamp avoids division by zero but introduces two artefacts. First, the
result is frozen to a constant wrong value inside the clamp band, since the
formula no longer responds to the true $p$ or $q$. Second, the result jumps
discontinuously when $p$ or $q$ crosses the clamp boundary. Backscattering is
particularly sensitive: the retro-reflected beams that dominate the
$theta = 180 degree$ bin all have small $p$ and $q$, and the backscatter
Mueller element $S_11$ is a strongly cancelling coherent sum, which amplifies
per-beam errors of order $0.1%$ into errors of several percent. This is the
numerical reason behind the `ZONE_THETA_OFFSET` guard, which places the
forward and backward zone bins at $0.01 degree$ away from the exact poles.

This note derives an equivalent edge sum in which the singularity never
appears. The result is exact at $p = q = 0$, smooth through the poles, and
removes the need for `KXY_EPSILON`, the slope guards `m_adj` and `n_adj`, and
the `ZONE_THETA_OFFSET` workaround.

= Setting and notation

Work in the aperture frame, where the beam cross section is a planar polygon
in the $x y$ plane with vertices $v_j = (x_j, y_j)$ for $j = 0, dots, N - 1$,
taken cyclically so that $v_N = v_0$. Write the edge vectors as

$ Delta x_j = x_(j+1) - x_j, quad Delta y_j = y_(j+1) - y_j. $

The far-field contribution of the beam at observation direction $hat(k)$
factorises into an amplitude part and the scalar Fraunhofer factor. The
Fraunhofer factor is proportional to the two-dimensional Fourier transform of
the aperture indicator function,

$ A(p, q) = integral.double_A e^(i (p x + q y)) dif x dif y, $ <ftdef>

evaluated at the wavevector mismatch

$ p = k (hat(k)_"inc",x - hat(k)_x), quad q = k (hat(k)_"inc",y - hat(k)_y), $

where $k$ is the wavenumber and $hat(k)_"inc"$ is the beam propagation
direction in the aperture frame. In the code, $p$ and $q$ are `kxx` and
`kyy`. The full far-field contribution is

$ E(hat(k)) prop -i k / (2 pi) e^(i k hat(k) dot r_0) A(p, q), $

where $r_0$ is the aperture centre of mass in the lab frame. The prefactor
and phase correspond to `inv_denom` and `bvsk` in the code and are unchanged
by everything below. The task is to evaluate $A(p, q)$ stably for all $p, q$,
including the limit $p, q -> 0$ where $A -> "Area"$.

= The current edge sum and its singular structure

Let $m_j = Delta y_j slash Delta x_j$ be the edge slope and define

$ delta_j = p x_j + q y_j, quad delta_(1,j) = p + m_j q, quad omega_j = Delta x_j delta_(1,j) = p Delta x_j + q Delta y_j. $

The implementation accumulates, per edge, a term algebraically equal to

$ F_j = -i e^(i delta_j) (e^(i omega_j) - 1) (p - m_j q) / (2 p q delta_(1,j)), $ <currentform>

and the total is $sum_j F_j = -i A(p, q)$. The factor $1 slash (p q)$ in
every term is the problem: each term diverges as $p -> 0$ or $q -> 0$, even
though the sum remains finite because the divergences cancel between edges.
The cancellation is exact in real arithmetic but is destroyed by the clamps,
and is badly conditioned in `f32` even without them.

The origin of the double singularity is that @currentform is the average of
two distinct exact representations of $A$, as shown next. Averaging puts both
$1 slash p$ and $1 slash q$ into every term. Using either representation
alone, chosen appropriately, removes the smaller denominator entirely.

= Two exact single-coordinate forms

Green's theorem in the plane states, for differentiable $L$ and $M$,

$ integral.double_A ((diff M) / (diff x) - (diff L) / (diff y)) dif x dif y = integral.cont (L dif x + M dif y), $

with the contour traversed in the positive (counterclockwise) sense. Two
choices of potential reproduce the integrand $e^(i(p x + q y))$:

+ $L = 0$, $M = e^(i(p x + q y)) slash (i p)$, giving
  $ A(p, q) = 1 / (i p) integral.cont e^(i (p x + q y)) dif y. $ <formp>

+ $M = 0$, $L = -e^(i(p x + q y)) slash (i q)$, giving
  $ A(p, q) = -1 / (i q) integral.cont e^(i (p x + q y)) dif x. $ <formq>

Note that @formp contains no $1 slash q$ and @formq contains no
$1 slash p$.

On the straight edge from $v_j$ to $v_(j+1)$, parameterise
$r(t) = v_j + t (Delta x_j, Delta y_j)$ with $t in [0, 1]$. Then
$p x + q y = delta_j + t omega_j$ and $dif y = Delta y_j dif t$, so

$ integral_("edge" j) e^(i (p x + q y)) dif y = Delta y_j e^(i delta_j) integral_0^1 e^(i t omega_j) dif t = Delta y_j e^(i delta_j) phi(omega_j), $

with the entire function

$ phi(omega) = (e^(i omega) - 1) / (i omega), quad phi(0) = 1. $ <phidef>

The two stable edge sums are therefore

$ A(p, q) = 1 / (i p) sum_(j=0)^(N-1) Delta y_j e^(i delta_j) phi(omega_j), $ <stablep>

$ A(p, q) = -1 / (i q) sum_(j=0)^(N-1) Delta x_j e^(i delta_j) phi(omega_j). $ <stableq>

Both are exact for all $p, q$ where they are defined. Only $delta_j$ and
$omega_j$ are needed per edge. No slopes appear, so the guards on
near-vertical and near-horizontal edges (`m_adj`, `n_adj`, `DIFF_DMIN`) are
unnecessary in this formulation.

Averaging @stablep and @stableq and substituting
$Delta y_j = m_j Delta x_j$ recovers @currentform after a few lines of
algebra. This confirms the equivalence with the current implementation and
shows that the $1 slash (p q)$ structure is an artefact of the
symmetrisation, not intrinsic to the problem.

= Evaluation of $phi$ for small argument

$phi$ is entire, with Taylor series

$ phi(omega) = 1 + (i omega) / 2 - omega^2 / 6 - (i omega^3) / 24 + cal(O)(omega^4). $ <phiseries>

For $|omega| < omega_0 = 10^(-2)$ the truncation error of @phiseries is below
$omega_0^4 slash 120 approx 8 times 10^(-11)$, far below `f32` resolution, so
the series can be used there and the closed form
$phi(omega) = sin(omega) slash omega + i (1 - cos(omega)) slash omega$
elsewhere. This removes the remaining removable singularity at
$omega_j = 0$.

= The limit of small $p$ and $q$

When both $|p|$ and $|q|$ are small the prefactor $1 slash (i p)$ in
@stablep still amplifies rounding error, because the leading terms of the sum
cancel. The limit can instead be taken analytically. Expand
$e^(i delta_j) phi(omega_j)$ to first order in the small quantities
$delta_j$ and $omega_j$ and insert into @stablep:

$ A = 1 / (i p) sum_j Delta y_j (1 + i delta_j + (i omega_j) / 2 + cal(O)(2)) . $

Three elementary identities for a closed polygon evaluate the sums:

+ $sum_j Delta y_j = 0$, by telescoping.

+ $sum_j Delta y_j (y_j + Delta y_j / 2) = 1/2 sum_j (y_(j+1)^2 - y_j^2) = 0$, by telescoping.

+ $sum_j Delta y_j (x_j + Delta x_j / 2) = integral.cont x dif y = S$, the signed area, exact for straight edges.

Collecting terms, with
$delta_j + omega_j slash 2 = p (x_j + Delta x_j slash 2) + q (y_j + Delta y_j slash 2)$:

$ A(p, q) = S + i (p M_x + q M_y) + cal(O)((kappa R)^2) S, $ <arealimit>

where $kappa = max(|p|, |q|)$, $R$ is the maximum vertex distance from the
origin, and the first moments are

$ M_x = integral.double_A x dif A, quad M_y = integral.double_A y dif A, $

given for a polygon by the standard shoelace-type formulas

$ S = 1/2 sum_j (x_j y_(j+1) - x_(j+1) y_j), $

$ M_x = 1/6 sum_j (x_j + x_(j+1)) (x_j y_(j+1) - x_(j+1) y_j), quad M_y = 1/6 sum_j (y_j + y_(j+1)) (x_j y_(j+1) - x_(j+1) y_j). $

In GOAD the aperture is translated to its centre of mass before diffraction,
so $M_x approx M_y approx 0$ and @arealimit reduces to $A approx S$. The
moment terms are retained anyway since they are cheap and make @arealimit
correct to second order for any polygon origin. Note that $S$ is the signed
area. The sign encodes the winding direction of the aperture polygon and must
not be replaced by the absolute area, since @stablep and @stableq carry the
same sign convention.

= Branch selection and error budget

The evaluation strategy is:

+ If $kappa R < tau$: use the limit @arealimit. The neglected second-order
  term gives relative error $cal(O)(tau^2)$.

+ Otherwise, if $|p| >= |q|$: use @stablep, dividing only by the larger
  component. Rounding in the cancelling sum gives relative error of order
  $epsilon_"f32" dot cal(P) slash (|p| |S|)$ with $cal(P)$ the perimeter
  scale, which is largest at the branch boundary $|p| = tau slash R$.

+ Otherwise: use @stableq.

With $tau = 3 times 10^(-3)$ both error sources stay near or below
$10^(-4)$ relative, compared with $1.5 times 10^(-3)$ and discontinuous for
the clamped form in the small-$(p, q)$ regime. Any $tau$ in
$[10^(-3), 10^(-2)]$ is acceptable; the total error is insensitive to the
exact choice.

= Numerical validation

The formulas were validated against dense `float64` quadrature of @ftdef on a
centred hexagonal aperture of circumradius 3, with the edge sums evaluated in
`f32` arithmetic to mimic the Rust implementation. The script is
`fraunhofer_limit_check.py` in this directory.

#table(
  columns: (auto, auto, auto),
  align: (left, right, right),
  table.header([case], [clamped form error], [stable form error]),
  [generic $p, q$], [$0.001%$], [$0.0005%$],
  [$q$ inside clamp band], [$0.024%$], [$0.0003%$],
  [$q = 0$ exactly], [$0.036%$], [$0.0003%$],
  [both inside clamp band], [$0.152%$], [$0.00001%$],
  [backscatter pole ($p = 2 times 10^(-5)$, $q = -8 times 10^(-8)$)], [$0.152%$], [$approx 0%$],
  [$p = q = 0$ exactly], [$0.152%$], [$approx 0%$],
)

The clamped form additionally returns the identical value for every $(p, q)$
inside the clamp band and jumps at its boundary. In a full GOAD run on a
hexagonal column (`config/local.toml`, seed 6), this produced erratic jumps
of $-1.8%$ and $+3.0%$ in the backscatter $S_11$ within $0.0001 degree$ of
$theta = 180 degree$, against a smooth extrapolated limit. The stable form
removes the discontinuity, which makes evaluation at exactly
$theta = 0 degree$ and $theta = 180 degree$ legitimate and allows
`ZONE_THETA_OFFSET` to be set to zero.

= Implementation notes

- $p, q$ are the raw `kxx`, `kyy` from `calculate_kxx_kyy` with the clamps
  removed.
- The per-edge inputs are only $delta_j$ and $omega_j$; `EdgeData` needs the
  raw $Delta x_j$, $Delta y_j$ (without the `DIFF_DMIN` floor) and no longer
  needs `m_adj`, `n_adj`.
- The signed area and first moments can be precomputed once per beam in
  `EdgeData::from_vertices`.
- The final Fraunhofer factor keeps its existing prefactors:
  multiply $A(p, q)$ by $-i$, by `inv_denom` $= k slash (2 pi)$, and by the
  phase $e^(i dot "bvsk")$.
- The legacy geometric-optics path in `diff.rs` has its own copy of the
  clamped formulas and is not affected.
