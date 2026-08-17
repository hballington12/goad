"""Shared physics for the lossy-media Phase B validation harness.

Implements, in float64 (reference) and float32 (emulating the Rust
implementation), the machinery derived in
docs/theory/lossy_media_beam_tracing.typ:

  - the reduced wave vector w = N e + i A with dispersion w.w = n_tilde^2
  - interface transfer: tangential continuity + one principal square root
  - generalised Fresnel coefficients (Chang eq. 14) in kappa = w.n form
  - the current GOAD Fresnel (complex-Snell cos_theta_t), for comparison
  - the Macke get_theta_t port from snell.rs, for comparison
  - the stable polygon edge sum A(p, q) continued to complex arguments,
    in both the plain per-edge form and the vertex-difference form that
    the Rust code uses

Conventions: interface normal n_hat points from medium 1 INTO medium 2
(along transmission). kappa = w . n_hat. All indices are the full complex
material indices n_tilde = n + i n'.
"""

import numpy as np

# Thresholds mirroring diff2.rs / the theory notes.
DIFFERENCE_THRESHOLD = 0.3     # |omega| below which phi's Taylor branch is used
POLYGON_FT_SMALL_ARG = 3e-3    # kappa * r_max below which the moment limit is used


# =========================================================================
# Reduced wave vector: construction and interface transfer (float64)
# =========================================================================

def w_homogeneous(n_tilde, e_hat):
    """Reduced wave vector of a homogeneous wave: w = n_tilde * e_hat."""
    return n_tilde * np.asarray(e_hat, dtype=np.complex128)


def transfer(w1, n_hat, n2_tilde, branch="continuation"):
    """Transmitted reduced wave vector across an interface.

    w2 = w_t + gamma * n_hat with gamma^2 = n2_tilde^2 - w_t.w_t.
    Returns (w2, gamma). n_hat from medium 1 into medium 2.

    branch = "continuation" (default, the physical prescription): the
    root with Re(gamma) + Im(gamma) > 0, i.e. the analytic continuation
    from the beam's dominant real tangential wavenumber. Validated
    against the branch-free spectral ground truth in
    spectral_branch_check.py; see the theory note's branch section.
    branch = "principal": the raw principal square root (what
    fresnel.rs did before the branch fix, and what PGOM does); above
    the critical angle with an absorbing medium 1 this is the wrong
    (conjugate) branch.
    """
    n_hat = np.asarray(n_hat, dtype=np.float64)
    w1 = np.asarray(w1, dtype=np.complex128)
    w_t = w1 - (w1 @ n_hat) * n_hat
    gamma = np.sqrt(np.complex128(n2_tilde**2 - w_t @ w_t))
    if branch == "continuation" and gamma.real + gamma.imag < 0:
        gamma = -gamma
    return w_t + gamma * n_hat, gamma


def reflect(w1, n_hat):
    """Reflected reduced wave vector: negate the normal component."""
    n_hat = np.asarray(n_hat, dtype=np.float64)
    w1 = np.asarray(w1, dtype=np.complex128)
    return w1 - 2.0 * (w1 @ n_hat) * n_hat


def apparent_indices(w):
    """(N, N_prime, e_hat, A) decomposition of a reduced wave vector."""
    re = np.real(w)
    im = np.imag(w)
    N = np.linalg.norm(re)
    Np = np.linalg.norm(im)
    e_hat = re / N if N > 0 else re
    return N, Np, e_hat, im


# =========================================================================
# Fresnel coefficients
# =========================================================================

def fresnel_general(kappa_i, kappa_t, n1_tilde, n2_tilde):
    """Generalised Fresnel coefficients, Chang eq. (14).

    Returns (r_s, t_s, r_p, t_p). kappa values are the complex normal
    components w . n_hat of the incident and transmitted reduced wave
    vectors.
    """
    r_s = (kappa_i - kappa_t) / (kappa_i + kappa_t)
    t_s = 2.0 * kappa_i / (kappa_i + kappa_t)
    d_p = n2_tilde**2 * kappa_i + n1_tilde**2 * kappa_t
    r_p = (n2_tilde**2 * kappa_i - n1_tilde**2 * kappa_t) / d_p
    t_p = 2.0 * n1_tilde * n2_tilde * kappa_i / d_p
    return r_s, t_s, r_p, t_p


def fresnel_current(n1_tilde, n2_tilde, theta_i):
    """The fresnel.rs formulas: complex cos_theta_t from Snell, with the
    continuation-rule branch correction (as fixed in fresnel.rs).

    Returns (r_perp, r_par, t_perp, t_par) matching refl()/refr() in
    fresnel.rs: r_par = f11, r_perp = f22.
    """
    cti = np.cos(theta_i)
    ratio = n1_tilde / n2_tilde
    ctt = np.sqrt(np.complex128(1.0 - ratio**2 * np.sin(theta_i) ** 2))
    if ctt.real + ctt.imag < 0:
        ctt = -ctt
    r_par = (n2_tilde * cti - n1_tilde * ctt) / (n1_tilde * ctt + n2_tilde * cti)
    r_perp = (n1_tilde * cti - n2_tilde * ctt) / (n1_tilde * cti + n2_tilde * ctt)
    t_par = 2.0 * n1_tilde * cti / (n1_tilde * ctt + n2_tilde * cti)
    t_perp = 2.0 * n1_tilde * cti / (n1_tilde * cti + n2_tilde * ctt)
    return r_perp, r_par, t_perp, t_par


# =========================================================================
# Macke get_theta_t (snell.rs port), for the homogeneous-incidence check
# =========================================================================

def macke_theta_t(theta_i, m1, m2):
    """Direct port of snell.rs::get_theta_t (Macke 1996 rt_c.f90)."""
    if m1 == m2:
        return theta_i
    k1 = m1.imag / m1.real
    k2 = m2.imag / m2.real
    krel = (k2 - k1) / (1.0 + k1 * k2)
    nrel = m2.real / m1.real * (1.0 + k1 * k2) / (1.0 + k1 * k1)

    ref1 = nrel * nrel
    ref2 = krel * krel
    ref3 = (1.0 + ref2) ** 2
    ref6 = ref1 * ref3 / ((1.0 + krel * k2) ** 2)

    sintiq = np.sin(theta_i) ** 2
    ref4 = 1.0 - (1.0 - ref2) / ref1 / ref3 * sintiq
    ref5 = 2.0 * krel / ref1 / ref3 * sintiq

    q2 = np.sqrt(ref4 * ref4 + ref5 * ref5)
    g = np.arccos(ref4 / q2) / 2.0

    ref7 = (np.cos(g) - k2 * np.sin(g)) ** 2
    rnstar = np.sqrt(sintiq + ref6 * q2 * ref7)
    return np.arcsin(np.sin(theta_i) / rnstar)


# =========================================================================
# Polygon transform A(p, q) at complex argument
# =========================================================================

def hexagon(radius=3.0, rot=0.3):
    ang = rot + np.arange(6) * np.pi / 3
    return np.stack([radius * np.cos(ang), radius * np.sin(ang)], axis=1)


def polygon_moments(verts):
    """(signed_area, Mx, My, r_max) via the shoelace forms."""
    nv = len(verts)
    area = mx = my = 0.0
    r_max = 0.0
    for j in range(nv):
        x0, y0 = verts[j]
        x1, y1 = verts[(j + 1) % nv]
        cr = x0 * y1 - x1 * y0
        area += cr
        mx += (x0 + x1) * cr
        my += (y0 + y1) * cr
        r_max = max(r_max, np.hypot(x0, y0))
    return 0.5 * area, mx / 6.0, my / 6.0, r_max


def reference_quadrature(verts, p, q, n=300):
    """Dense float64 quadrature of A(p, q) with complex p, q.

    Fan-triangulates from the centroid and applies the midpoint rule on
    a barycentric grid, exactly as fraunhofer_limit_check.py does, but
    with complex arguments so the exponential ramp is included.
    """
    p = np.complex128(p)
    q = np.complex128(q)
    c = verts.mean(axis=0)
    nv = len(verts)
    total = 0.0 + 0.0j
    u, v = np.meshgrid(np.arange(n), np.arange(n))
    mask_low = (u + v) < n
    mask_up = (u + v) < (n - 1)
    uu1 = (u + 1 / 3)[mask_low] / n
    vv1 = (v + 1 / 3)[mask_low] / n
    uu2 = (u + 2 / 3)[mask_up] / n
    vv2 = (v + 2 / 3)[mask_up] / n
    uu = np.concatenate([uu1, uu2])
    vv = np.concatenate([vv1, vv2])
    for j in range(nv):
        a, b = verts[j], verts[(j + 1) % nv]
        ac, bc = a - c, b - c
        area = 0.5 * abs(ac[0] * bc[1] - ac[1] * bc[0])
        xs = c[0] + (a[0] - c[0]) * uu + (b[0] - c[0]) * vv
        ys = c[1] + (a[1] - c[1]) * uu + (b[1] - c[1]) * vv
        total += np.exp(1j * (p * xs + q * ys)).sum() * (area / (n * n))
    return total


def _phi_series_c64(omega):
    """phi(omega) = (e^{i omega} - 1)/(i omega) by Taylor series, complex64.

    Same retained orders as diff2.rs::phi: through omega^5, truncation
    |omega|^6 / 5040.
    """
    w = np.complex64(omega)
    iw = np.complex64(1j) * w
    return (
        np.complex64(1.0)
        + iw / np.float32(2.0)
        - w * w / np.float32(6.0)
        - iw * w * w / np.float32(24.0)
        + w * w * w * w / np.float32(120.0)
        + iw * w * w * w * w / np.float32(720.0)
    )


def _cis_c64(delta):
    """e^{i delta} for complex delta in float32 pieces: complex64 result."""
    dr = np.float32(np.real(delta))
    di = np.float32(np.imag(delta))
    mag = np.float32(np.exp(-di))
    return np.complex64(complex(mag * np.cos(dr), mag * np.sin(dr)))


def polygon_ft(verts, p, q, f32=True, difference_form=True,
               tau=POLYGON_FT_SMALL_ARG, thresh=DIFFERENCE_THRESHOLD):
    """Stable edge sum for A(p, q), continued to complex p, q.

    f32=True emulates the Rust float32 arithmetic (complex64 accumulation,
    float32 scalars). difference_form=True uses the vertex-difference
    identity that halves transcendental calls in diff2.rs, with the phi
    Taylor branch below |omega| = thresh; False uses the plain per-edge
    e^{i delta} phi(omega) form.
    """
    if f32:
        cplx = np.complex64
        real = np.float32
    else:
        cplx = np.complex128
        real = np.float64
    p = cplx(p)
    q = cplx(q)

    area, mx, my, r_max = polygon_moments(verts)
    kappa = max(abs(p), abs(q))
    if kappa * r_max < tau:
        return complex(cplx(area) + cplx(1j) * (p * cplx(mx) + q * cplx(my)))

    nv = len(verts)
    use_p = abs(p) >= abs(q)
    total = cplx(0)

    def cis(delta):
        if f32:
            return _cis_c64(delta)
        return np.exp(1j * np.complex128(delta))

    def phi(omega):
        if f32:
            if abs(omega) < thresh:
                return _phi_series_c64(omega)
            w = np.complex64(omega)
            return (cis(w) - np.complex64(1.0)) / (np.complex64(1j) * w)
        w = np.complex128(omega)
        if abs(w) < 1e-8:
            return np.complex128(1.0) + 1j * w / 2.0
        return (np.exp(1j * w) - 1.0) / (1j * w)

    xs = [real(v[0]) for v in verts]
    ys = [real(v[1]) for v in verts]
    deltas = [p * xs[j] + q * ys[j] for j in range(nv)]

    if difference_form:
        e = [cis(d) for d in deltas]
        for j in range(nv):
            jn = (j + 1) % nv
            dx = xs[jn] - xs[j]
            dy = ys[jn] - ys[j]
            omega = p * dx + q * dy
            if abs(omega) < thresh:
                factor = e[j] * phi(omega)
            else:
                factor = (e[jn] - e[j]) / (cplx(1j) * omega)
            weight = dy if use_p else dx
            total += factor * cplx(weight)
    else:
        for j in range(nv):
            jn = (j + 1) % nv
            dx = xs[jn] - xs[j]
            dy = ys[jn] - ys[j]
            omega = p * dx + q * dy
            weight = dy if use_p else dx
            total += cplx(weight) * cis(deltas[j]) * phi(omega)

    if use_p:
        return complex(total / (cplx(1j) * p))
    return complex(-total / (cplx(1j) * q))
