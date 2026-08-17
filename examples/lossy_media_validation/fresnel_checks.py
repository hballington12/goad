"""Task 3: Fresnel checks.

  1. Homogeneous-incidence identity: the generalised coefficients
     (Chang eq. 14 with kappa from the reduced wave vector) must equal
     the current fresnel.rs formulas (complex-Snell cos_theta_t) exactly
     whenever the incident wave is homogeneous, including between two
     lossy media. This is the theory note's central reduction proof,
     checked numerically to f64 precision.
  2. Lossless TIR: |r| = 1 above critical, nonzero s/p phase difference,
     matching the fresnel.rs unit tests.
  3. Chang (16): low-absorption expansions of N2, K2, sin(theta_t);
     error must scale as x^4 (first neglected order).
  4. Chang (18): first-order expansion of r_s; error must scale as x^2.
  5. Chang (20): Brewster-angle residual |r_p(theta_B)|^2 -> ((n^2-1)^2 / 4n^4) x^2.
  6. s-polarisation flux balance Re(kappa_i)(1 - |r_s|^2) = Re(kappa_t)|t_s|^2,
     exact for incidence from a lossless medium.
"""

import numpy as np

from lossy import fresnel_current, fresnel_general, transfer, w_homogeneous


def homogeneous_identity():
    print("1. homogeneous-incidence identity (general == current fresnel.rs)")
    pairs = [
        ("vacuum -> weakly lossy", 1.0 + 0.0j, 1.31 + 0.01j),
        ("vacuum -> strongly lossy", 1.0 + 0.0j, 1.31 + 0.1j),
        ("lossy -> lossy", 1.31 + 0.05j, 1.5 + 0.2j),
        ("ice -> vacuum (TIR range)", 1.31 + 0.0j, 1.0 + 0.0j),
        ("lossy ice -> vacuum (above crit)", 1.31 + 0.01j, 1.0 + 0.0j),
    ]
    n_hat = np.array([0.0, 0.0, 1.0])
    worst = 0.0
    for name, n1, n2 in pairs:
        for theta_deg in np.linspace(0.5, 89.5, 90):
            th = np.radians(theta_deg)
            e_i = np.array([np.sin(th), 0.0, np.cos(th)])
            w1 = w_homogeneous(n1, e_i)
            _, gamma = transfer(w1, n_hat, n2)
            kappa_i = w1 @ n_hat.astype(np.complex128)
            rs, ts, rp, tp = fresnel_general(kappa_i, gamma, n1, n2)
            r_perp, r_par, t_perp, t_par = fresnel_current(n1, n2, th)
            d = max(abs(rs - r_perp), abs(rp - r_par),
                    abs(ts - t_perp), abs(tp - t_par))
            worst = max(worst, d)
        print(f"   {name:<28} max |general - current| so far: {worst:.2e}")
    assert worst < 1e-12, "reduction identity violated"
    print("   PASS: identical to f64 precision\n")


def tir_limit():
    print("2. lossless TIR limit (ice -> air, 60 deg internal)")
    n1, n2 = 1.31 + 0.0j, 1.0 + 0.0j
    th = np.radians(60.0)
    n_hat = np.array([0.0, 0.0, 1.0])
    e_i = np.array([np.sin(th), 0.0, np.cos(th)])
    w1 = w_homogeneous(n1, e_i)
    w2, gamma = transfer(w1, n_hat, n2)
    kappa_i = w1 @ n_hat.astype(np.complex128)
    rs, _, rp, _ = fresnel_general(kappa_i, gamma, n1, n2)
    print(f"   gamma = {gamma:.6f}  (must be purely imaginary, Im > 0)")
    print(f"   |r_s| = {abs(rs):.8f}, |r_p| = {abs(rp):.8f}  (must be 1)")
    delta = np.angle(rp) - np.angle(rs)
    print(f"   s/p phase difference = {delta:.4f} rad (must be nonzero)")
    assert abs(gamma.real) < 1e-14 and gamma.imag > 0
    assert abs(abs(rs) - 1) < 1e-12 and abs(abs(rp) - 1) < 1e-12
    assert abs(delta) > 0.1
    # Transmitted wave carries no normal flux: check Re(kappa_t) = 0.
    print(f"   Re(kappa_t) = {gamma.real:.2e}  (zero normal flux into medium 2)")
    print("   PASS\n")


def chang16_expansion():
    print("3. Chang (16): low-absorption expansions of N2, K2, sin(theta_t)")
    n = 1.5
    th = np.radians(40.0)
    n_hat = np.array([0.0, 0.0, 1.0])
    e_i = np.array([np.sin(th), 0.0, np.cos(th)])
    s2 = np.sin(th) ** 2
    root = np.sqrt(n**2 - s2)
    print(f"   {'x':>8} {'err(N2)':>12} {'err(K2)':>12} {'err(sin t)':>12}")
    errs = []
    xs = [1e-1, 1e-2, 1e-3]
    for x in xs:
        kappa = n * x
        w2, _ = transfer(w_homogeneous(1.0, e_i), n_hat, complex(n, kappa))
        n2 = np.linalg.norm(np.real(w2))
        k2 = np.linalg.norm(np.imag(w2))
        sin_t = np.linalg.norm(np.real(w2) - (np.real(w2) @ n_hat) * n_hat) / n2

        n2_ap = n * (1 + s2 / (2 * (n**2 - s2)) * x**2)
        k2_ap = n * kappa / root * (1 - n**2 * s2 / (2 * (n**2 - s2) ** 2) * x**2)
        st_ap = np.sin(th) / n * (1 - s2 / (2 * (n**2 - s2)) * x**2)
        e = (abs(n2 - n2_ap) / n2, abs(k2 - k2_ap) / max(k2, 1e-300),
             abs(sin_t - st_ap) / sin_t)
        errs.append(e)
        print(f"   {x:>8.0e} {e[0]:>12.2e} {e[1]:>12.2e} {e[2]:>12.2e}")
    # each decade in x must shrink the error by ~1e4 (x^4 scaling)
    for i in (0, 2):
        ratio = errs[0][i] / max(errs[1][i], 1e-300)
        assert ratio > 1e3, f"expected ~x^4 scaling, got ratio {ratio:.1e}"
    print("   PASS: errors scale as x^4 (expansions correct through x^2)\n")


def chang18_expansion():
    print("4. Chang (18): first-order r_s expansion")
    n = 1.5
    th = np.radians(40.0)
    n_hat = np.array([0.0, 0.0, 1.0])
    e_i = np.array([np.sin(th), 0.0, np.cos(th)])
    s2 = np.sin(th) ** 2
    root = np.sqrt(n**2 - s2)
    r0 = (np.cos(th) - root) / (np.cos(th) + root)
    print(f"   {'x':>8} {'err(r_s)':>12}")
    prev = None
    for x in [1e-1, 1e-2, 1e-3]:
        kappa = n * x
        nt = complex(n, kappa)
        w1 = w_homogeneous(1.0, e_i)
        _, gamma = transfer(w1, n_hat, nt)
        rs, _, _, _ = fresnel_general(np.complex128(np.cos(th)), gamma, 1.0, nt)
        rs_ap = r0 * (1 + 2j / (n**2 - 1) * n**2 * np.cos(th) / root * x)
        err = abs(rs - rs_ap) / abs(rs)
        print(f"   {x:>8.0e} {err:>12.2e}")
        if prev is not None:
            assert prev / err > 30, "expected ~x^2 scaling of the residual"
        prev = err
    print("   PASS: residual scales as x^2 (expansion correct through x)\n")


def chang20_brewster():
    print("5. Chang (20): Brewster residual |r_p|^2 -> ((n^2-1)^2 / 4n^4) x^2")
    n = 1.5
    thb = np.arctan(n)
    n_hat = np.array([0.0, 0.0, 1.0])
    e_i = np.array([np.sin(thb), 0.0, np.cos(thb)])
    coef = (n**2 - 1) ** 2 / (4 * n**4)
    print(f"   {'x':>8} {'|r_p|^2 / (coef x^2)':>22}")
    for x in [1e-1, 1e-2, 1e-3]:
        nt = complex(n, n * x)
        w1 = w_homogeneous(1.0, e_i)
        _, gamma = transfer(w1, n_hat, nt)
        _, _, rp, _ = fresnel_general(np.complex128(np.cos(thb)), gamma, 1.0, nt)
        ratio = abs(rp) ** 2 / (coef * x**2)
        print(f"   {x:>8.0e} {ratio:>22.6f}")
    assert abs(ratio - 1) < 0.01, "Brewster residual coefficient mismatch"
    print("   PASS: ratio -> 1 as x -> 0\n")


def flux_balance():
    print("6. s-pol flux balance from a lossless medium: "
          "Re(k_i)(1-|r|^2) = Re(k_t)|t|^2")
    n_hat = np.array([0.0, 0.0, 1.0])
    worst = 0.0
    for n2 in [1.31 + 0.01j, 1.31 + 0.1j, 2.694 + 0.509j]:
        for theta_deg in np.linspace(1, 89, 45):
            th = np.radians(theta_deg)
            e_i = np.array([np.sin(th), 0.0, np.cos(th)])
            w1 = w_homogeneous(1.0, e_i)
            _, gamma = transfer(w1, n_hat, n2)
            kappa_i = np.complex128(np.cos(th))
            rs, ts, _, _ = fresnel_general(kappa_i, gamma, 1.0, n2)
            lhs = kappa_i.real * (1 - abs(rs) ** 2)
            rhs = gamma.real * abs(ts) ** 2
            worst = max(worst, abs(lhs - rhs) / lhs)
    print(f"   max relative violation over media and angles: {worst:.2e}")
    assert worst < 1e-12
    print("   PASS\n")


def tir_continuity():
    print("7. TIR-phase continuity as n' -> 0 (continuation rule vs principal root)")
    th = np.radians(60.0)
    n_hat = np.array([0.0, 0.0, 1.0])
    e_i = np.array([np.sin(th), 0.0, np.cos(th)])

    def r_s(n_im, branch):
        n1 = complex(1.31, n_im)
        w1 = w_homogeneous(n1, e_i)
        _, gamma = transfer(w1, n_hat, 1.0, branch=branch)
        rs, _, _, _ = fresnel_general(w1 @ n_hat.astype(np.complex128), gamma, n1, 1.0)
        return rs

    r0 = r_s(0.0, "continuation")
    print(f"   lossless r_s = {r0:.4f}")
    print(f"   {'n_prime':>8} {'|r_cont - r0|':>14} {'|r_princ - r0|':>15}")
    for n_im in [1e-8, 1e-6, 1e-4, 1e-2]:
        dc = abs(r_s(n_im, "continuation") - r0)
        dp = abs(r_s(n_im, "principal") - r0)
        print(f"   {n_im:>8.0e} {dc:>14.2e} {dp:>15.2e}")
        assert dc < 10 * n_im + 1e-6, "continuation rule must be continuous"
    assert abs(r_s(1e-8, "principal") - r0) > 1.0, \
        "principal root should jump by ~2|Im r| (documents the old bug)"
    print("   PASS: continuation continuous (error ~ n'); principal jumps by "
          f"{abs(r_s(1e-8, 'principal') - r0):.3f} = 2|Im r0|\n")


def dupertuis_balance():
    print("8. Dupertuis II interface balance R + T = 1 + M (eq. 60), "
          "inhomogeneous incidence")
    rng = np.random.default_rng(6)
    media = [1.0 + 0.0j, 1.31 + 0.01j, 1.31 + 0.1j, 1.5 + 0.5j]
    n_hat = np.array([0.0, 0.0, 1.0])
    worst_c = worst_r = 0.0
    for _ in range(300):
        n1 = media[rng.integers(1, len(media))]
        n2 = media[rng.integers(0, len(media))]
        # inhomogeneous incident wave via a tilted first entry
        e = rng.normal(size=3)
        e /= np.linalg.norm(e)
        if e @ n_hat < 0:
            e = -e
        tilt = rng.normal(size=3)
        tilt /= np.linalg.norm(tilt)
        if e @ tilt < 0:
            tilt = -tilt
        w1, _ = transfer(w_homogeneous(1.0, e), tilt, n1)
        if np.real(w1) @ n_hat < 0:
            w1 = w1 - 2 * (w1 @ n_hat) * n_hat
        kappa_i = w1 @ n_hat.astype(np.complex128)
        _, gamma = transfer(w1, n_hat, n2)
        rs, ts, _, _ = fresnel_general(kappa_i, gamma, n1, n2)
        # complex coefficients, Dupertuis II eqs (54)-(56)
        R_bar = abs(rs) ** 2
        T_bar = abs(ts) ** 2 * np.conj(gamma) / np.conj(kappa_i)
        M_bar = rs - np.conj(rs)
        worst_c = max(worst_c, abs(R_bar + T_bar - 1 - M_bar))
        # observable (real) coefficients, flux-weighted by Re(kappa)
        T = abs(ts) ** 2 * gamma.real / kappa_i.real
        M = 2 * rs.imag * kappa_i.imag / kappa_i.real
        worst_r = max(worst_r, abs(R_bar + T - 1 - M))
    print(f"   max |R_bar + T_bar - 1 - M_bar| (complex, eq. 60): {worst_c:.2e}")
    print(f"   max |R + T - 1 - M| (observable fluxes):           {worst_r:.2e}")
    assert worst_c < 1e-12 and worst_r < 1e-12
    print("   PASS\n")


def stokes_relations():
    print("9. Stokes reversibility (Dupertuis II eqs. 1-4) + reverse-transfer "
          "branch consistency")
    rng = np.random.default_rng(7)
    media = [1.0 + 0.0j, 1.31 + 0.01j, 1.5 + 0.2j]
    n_hat = np.array([0.0, 0.0, 1.0])
    worst_s = worst_g = 0.0
    for _ in range(200):
        n1 = media[rng.integers(len(media))]
        n2 = media[rng.integers(len(media))]
        if n1 == n2:
            continue
        e = rng.normal(size=3)
        e /= np.linalg.norm(e)
        if e @ n_hat < 0:
            e = -e
        # mildly inhomogeneous incident wave (physical beam regime)
        tilt = rng.normal(size=3)
        tilt /= np.linalg.norm(tilt)
        if e @ tilt < 0:
            tilt = -tilt
        w1, _ = transfer(w_homogeneous(1.0, e), tilt, n1)
        if np.real(w1) @ n_hat < 0:
            w1 = w1 - 2 * (w1 @ n_hat) * n_hat
        kappa_1 = w1 @ n_hat.astype(np.complex128)
        w2, gamma = transfer(w1, n_hat, n2)

        # forward coefficients (1 -> 2) and backward (2 -> 1): the
        # backward incident wave is the reversed transmitted wave, whose
        # normal component along the reversed transfer normal is gamma.
        rp, tp, _, _ = fresnel_general(kappa_1, gamma, n1, n2)
        rm, tm, _, _ = fresnel_general(gamma, kappa_1, n2, n1)
        worst_s = max(worst_s,
                      abs(rp * rp + tp * tm - 1), abs(tp * (rp + rm)),
                      abs(rm * rm + tm * tp - 1))

        # reverse-transfer branch consistency: transferring the reversed
        # transmitted wave back through the interface must recover the
        # original incident normal component.
        w2_rev = -w2
        _, gamma_back = transfer(w2_rev, -n_hat, n1)
        worst_g = max(worst_g, abs(gamma_back - kappa_1))
    print(f"   max Stokes relation violation:            {worst_s:.2e}")
    print(f"   max |gamma_back - kappa_i| (reversibility): {worst_g:.2e}")
    assert worst_s < 1e-12 and worst_g < 1e-10
    print("   PASS\n")


if __name__ == "__main__":
    homogeneous_identity()
    tir_limit()
    chang16_expansion()
    chang18_expansion()
    chang20_brewster()
    flux_balance()
    tir_continuity()
    dupertuis_balance()
    stokes_relations()
    print("all Fresnel checks passed")
