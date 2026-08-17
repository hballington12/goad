"""Task 4: interface-chain invariants of the reduced-wave-vector transfer.

  1. Random chains of lossy interfaces: dispersion w.w = n_tilde^2 and
     tangential continuity hold to f64 precision at every step.
  2. Slab reversibility: two parallel interfaces (medium 1 -> 2 -> 1)
     restore the original w exactly, for inhomogeneous incident waves.
  3. First entry from a lossless host: the absorption vector is normal
     to the entry facet (A parallel to -n_hat convention of PGOM).
  4. Exit into a lossless host: A perpendicular to the propagation
     direction (pure transverse ramp).
  5. No TIR in lossy media: sin(theta_t) < 1 strictly, across and beyond
     the lossless critical angle.
  6. Macke get_theta_t (snell.rs) vs the transfer rule, homogeneous
     incidence: quantifies how far the current direction computation is
     from the rigorous one.
"""

import numpy as np

from lossy import macke_theta_t, reflect, transfer, w_homogeneous

RNG = np.random.default_rng(6)

MEDIA = [1.0 + 0.0j, 1.31 + 0.01j, 1.31 + 0.1j, 1.5 + 0.5j, 2.694 + 0.509j]


def random_unit():
    v = RNG.normal(size=3)
    return v / np.linalg.norm(v)


def random_chain():
    print("1. random interface chains: dispersion + tangential continuity")
    worst_disp = worst_tang = 0.0
    for _ in range(200):
        n1 = MEDIA[RNG.integers(len(MEDIA))]
        w = w_homogeneous(n1, random_unit())
        for _ in range(8):
            n2 = MEDIA[RNG.integers(len(MEDIA))]
            n_hat = random_unit()
            # orient the normal along transmission
            if np.real(w) @ n_hat < 0:
                n_hat = -n_hat
            w2, _ = transfer(w, n_hat, n2)
            wr = reflect(w, n_hat)
            disp = abs(w2 @ w2 - n2**2)
            tang = np.linalg.norm(
                (w2 - (w2 @ n_hat) * n_hat) - (w - (w @ n_hat) * n_hat))
            disp_r = abs(wr @ wr - (w @ w))
            worst_disp = max(worst_disp, disp, disp_r)
            worst_tang = max(worst_tang, tang)
            w, n1 = w2, n2
    print(f"   max dispersion violation |w.w - n^2|: {worst_disp:.2e}")
    print(f"   max tangential discontinuity:         {worst_tang:.2e}")
    assert worst_disp < 1e-10 and worst_tang < 1e-12
    print("   PASS\n")


def slab_reversibility():
    """Reversibility holds wherever the scalar continuation rule is in
    its validity region. The rule picks the root with Re + Im > 0; time
    reversal returns the original wave iff the original normal
    component kappa_1 itself satisfies Re + Im > 0. Waves violating
    that (ramp strongly anti-normal at near-grazing incidence, i.e.
    within ~K_s of grazing) sit inside the single-wave model's known
    breakdown sliver, so reversibility is asserted only outside it and
    the excluded fraction is reported.
    """
    print("2. slab reversibility (1 -> 2 -> 1 through parallel faces)")
    n_hat = np.array([0.0, 0.0, 1.0])
    worst = 0.0
    excluded = 0
    total = 0
    worst_excluded = 0.0
    for _ in range(200):
        n1 = MEDIA[RNG.integers(len(MEDIA))]
        n2 = MEDIA[RNG.integers(len(MEDIA))]
        # build a genuinely inhomogeneous incident wave in medium 1 by
        # entering it from vacuum through a tilted facet
        e = random_unit()
        if e @ n_hat < 0:
            e = -e
        tilt = random_unit()
        if e @ tilt < 0:
            tilt = -tilt
        w0, _ = transfer(w_homogeneous(1.0, e), tilt, n1)
        if np.real(w0) @ n_hat < 0:
            w0 = reflect(w0, n_hat)  # ensure it heads into the slab
        total += 1
        kappa_1 = w0 @ n_hat.astype(np.complex128)
        w_in, _ = transfer(w0, n_hat, n2)
        w_back, _ = transfer(w_in, n_hat, n1)
        err = np.linalg.norm(w_back - w0)
        if kappa_1.real + kappa_1.imag > 1e-3:
            worst = max(worst, err)
        else:
            excluded += 1
            worst_excluded = max(worst_excluded, err)
    print(f"   max |w_out - w_in| in the validity region: {worst:.2e} "
          f"({total - excluded}/{total} draws)")
    print(f"   excluded (anti-normal-ramp near grazing, model breakdown "
          f"sliver): {excluded}/{total}, worst there {worst_excluded:.2e}")
    assert worst < 1e-12
    print("   PASS\n")


def first_entry():
    print("3. first entry from lossless host: A normal to the facet")
    n_hat = np.array([0.0, 0.0, 1.0])
    worst = 0.0
    for theta_deg in np.linspace(0, 89, 90):
        th = np.radians(theta_deg)
        e_i = np.array([np.sin(th), 0.0, np.cos(th)])
        w2, _ = transfer(w_homogeneous(1.0, e_i), n_hat, 1.31 + 0.1j)
        a = np.imag(w2)
        tang = np.linalg.norm(a - (a @ n_hat) * n_hat)
        worst = max(worst, tang)
    print(f"   max tangential |A| component: {worst:.2e} (must be 0)")
    assert worst < 1e-14
    print("   PASS\n")


def exit_transverse():
    print("4. exit into lossless host: A perpendicular to propagation")
    n_hat_in = np.array([0.0, 0.0, 1.0])
    worst = 0.0
    for theta_deg in np.linspace(1, 60, 60):
        th = np.radians(theta_deg)
        e_i = np.array([np.sin(th), 0.0, np.cos(th)])
        w_in, _ = transfer(w_homogeneous(1.0, e_i), n_hat_in, 1.31 + 0.1j)
        # exit through a tilted facet so the interior wave is oblique to it
        tilt = np.array([np.sin(0.4), 0.0, np.cos(0.4)])
        w_out, _ = transfer(w_in, tilt, 1.0)
        a = np.imag(w_out)
        e = np.real(w_out)
        if np.linalg.norm(a) > 0:
            worst = max(worst, abs(a @ e) / (np.linalg.norm(a) * np.linalg.norm(e)))
    print(f"   max |cos(angle(A, e))| on exit: {worst:.2e} (must be 0)")
    assert worst < 1e-12
    print("   PASS\n")


def no_tir():
    print("5. no TIR in lossy media (ice-like -> air, beyond lossless critical)")
    n_hat = np.array([0.0, 0.0, 1.0])
    crit = np.degrees(np.arcsin(1.0 / 1.31))
    print(f"   lossless critical angle: {crit:.2f} deg")
    for kappa in [0.0, 0.01, 0.1]:
        nt = complex(1.31, kappa)
        row = []
        for theta_deg in [40.0, 49.76, 60.0, 80.0]:
            th = np.radians(theta_deg)
            e_i = np.array([np.sin(th), 0.0, np.cos(th)])
            w1 = w_homogeneous(nt, e_i)  # homogeneous interior wave
            w2, _ = transfer(w1, n_hat, 1.0)
            re2 = np.real(w2)
            sin_t = np.linalg.norm(re2 - (re2 @ n_hat) * n_hat) / np.linalg.norm(re2)
            row.append(sin_t)
        marks = ", ".join(f"{s:.6f}" for s in row)
        print(f"   kappa = {kappa:>5}: sin(theta_t) at 40/49.76/60/80 deg = {marks}")
        if kappa > 0:
            assert all(s < 1.0 for s in row), "lossy medium must transmit"
    print("   PASS: sin(theta_t) = 1 only in the lossless case at/above critical\n")


def macke_comparison():
    print("6. Macke get_theta_t (snell.rs) vs transfer rule, homogeneous incidence")
    n_hat = np.array([0.0, 0.0, 1.0])
    pairs = [
        ("vacuum -> 1.31+0.01i", 1.0 + 0.0j, 1.31 + 0.01j),
        ("vacuum -> 1.31+0.1i", 1.0 + 0.0j, 1.31 + 0.1j),
        ("vacuum -> 1.5+0.5i", 1.0 + 0.0j, 1.5 + 0.5j),
        ("1.31+0.1i -> vacuum", 1.31 + 0.1j, 1.0 + 0.0j),
        ("1.31+0.01i -> 1.5+0.2i", 1.31 + 0.01j, 1.5 + 0.2j),
    ]
    for name, n1, n2 in pairs:
        worst = 0.0
        for theta_deg in np.linspace(0.5, 85, 170):
            th = np.radians(theta_deg)
            e_i = np.array([np.sin(th), 0.0, np.cos(th)])
            w2, _ = transfer(w_homogeneous(n1, e_i), n_hat, n2)
            re2 = np.real(w2)
            sin_t = np.linalg.norm(re2 - (re2 @ n_hat) * n_hat) / np.linalg.norm(re2)
            sin_t_macke = np.sin(macke_theta_t(th, n1, n2))
            worst = max(worst, abs(sin_t - sin_t_macke))
        print(f"   {name:<26} max |d sin(theta_t)| = {worst:.2e}")
    print("   (agreement here means snell.rs is exactly the homogeneous special case;")
    print("    any gap quantifies a pre-existing direction error in the current code)\n")


if __name__ == "__main__":
    random_chain()
    slab_reversibility()
    first_entry()
    exit_transverse()
    no_tir()
    macke_comparison()
    print("all interface-chain checks passed")
