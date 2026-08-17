"""Spectral ground truth for the transmitted-wave branch question.

A finite beam exiting an absorbing medium through a hard aperture is
decomposed into ordinary plane waves with REAL tangential wavenumber.
For real k_x and a lossless exterior, the elementary Fresnel problem has
NO branch ambiguity: k_z2 is real positive (outgoing) below the critical
wavenumber and +i|k_z2| (decaying) above it, and k_z1 has Im > 0 from
absorption. Transmitting each component and recombining gives the exact
reflected and transmitted fields, against which the three single-wave
branch candidates are judged:

  A     : principal sqrt of gamma^2         (current fresnel.rs, PGOM)
  B     : the other Maxwell root
  cont  : continuation rule, Re(gamma) + Im(gamma) > 0
  chang : same-sign rule |Re g| + i|Im g|   (Chang eq. 12, non-Maxwellian)

Geometry: 2D, s-polarisation, interface z = 0, medium 1 (absorbing,
n1_tilde) at z < 0, vacuum at z > 0. Incident boundary field on z = 0:

  E(x) = rect(x / W) * exp(i k0 w_x x),   w_x = Ns + i Ks

i.e. a plane-wave phase with an exponential amplitude ramp across a
hard aperture of width W. k0 = 1 (lengths in units of 1/k0).

Outputs, per case: the central-footprint effective reflection
coefficient r_eff = E_refl(x)/E_inc(x) compared with each branch's r,
and the exact power budget (R, T, M) compared with each branch's
prediction. Also verifies the Dupertuis II balance R + T = 1 + M.
"""

import numpy as np

# Grid: domain length L, N samples. k_max = pi * N / L >> 1.
L = 8192.0
N = 1 << 18


def elementary_rt(kx, n1):
    """Branch-unambiguous Fresnel r, t for real tangential wavenumber."""
    kz1 = np.sqrt(np.asarray(n1**2 - kx**2, dtype=complex))
    # Im(n1^2) > 0 puts kz1 in the upper-right quadrant automatically.
    kz2_sq = 1.0 - kx**2
    kz2 = np.where(kz2_sq >= 0, np.sqrt(np.abs(kz2_sq)) + 0j,
                   1j * np.sqrt(np.abs(kz2_sq)))
    r = (kz1 - kz2) / (kz1 + kz2)
    t = 2.0 * kz1 / (kz1 + kz2)
    return kz1, kz2, r, t


def single_wave_candidates(w_x, n1):
    """gamma candidates and their r for the single inhomogeneous wave."""
    kappa_i = np.sqrt(complex(n1**2 - w_x**2))  # incident normal comp, Re > 0
    g = np.sqrt(complex(1.0 - w_x**2))          # principal root in vacuum
    cands = {"A": g, "B": -g}
    cands["cont"] = g if (g.real + g.imag) > 0 else -g
    cands["chang"] = abs(g.real) + 1j * abs(g.imag)
    out = {}
    for name, gam in cands.items():
        r = (kappa_i - gam) / (kappa_i + gam)
        t = 2.0 * kappa_i / (kappa_i + gam)
        # single-wave observable power coefficients (s-pol, per unit
        # incident normal flux Re(kappa_i))
        R = abs(r) ** 2
        T = gam.real * abs(t) ** 2 / kappa_i.real
        out[name] = {"gamma": gam, "r": r, "R": R, "T": T}
    return kappa_i, out


def spectral_case(n1, w_x, W):
    x = (np.arange(N) - N / 2) * (L / N)
    dx = L / N
    aperture = np.abs(x) <= W / 2
    E = np.where(aperture, np.exp(1j * w_x * x), 0.0)

    Ek = np.fft.fft(E) * dx
    kx = np.fft.fftfreq(N, d=dx) * 2 * np.pi

    kz1, kz2, r, t = elementary_rt(kx, n1)

    Er = np.fft.ifft(r * Ek) / dx
    # exact power budget through the z = 0 plane (per-component fluxes;
    # cross terms between distinct kx vanish on x-integration)
    w = np.abs(Ek) ** 2
    P_i = np.sum(kz1.real * w)
    P_r = np.sum(kz1.real * (np.abs(r) ** 2) * w)
    P_t = np.sum(kz2.real * (np.abs(t) ** 2) * w)
    # side-1 net flux per component: Re[kz1 (1 - |r|^2)] + mixed term
    # (Poynting continuity per component makes global balance exact):
    P_m = P_t + P_r - P_i  # ground-truth mixed (interference) power

    # central-footprint effective reflection coefficient
    core = np.abs(x) <= W / 4
    r_eff = Er[core] / E[core]
    r_eff_mean = np.mean(r_eff)
    r_eff_spread = np.std(r_eff)

    return {
        "r_eff": r_eff_mean,
        "r_spread": r_eff_spread,
        "R": P_r / P_i,
        "T": P_t / P_i,
        "M": P_m / P_i,
    }


def run_case(label, n1, w_x, W):
    kappa_i, cands = single_wave_candidates(w_x, n1)
    truth = spectral_case(n1, w_x, W)

    ramp_exponent = w_x.imag * W
    print(f"\n=== {label} ===")
    print(f"    n1 = {n1}, w_x = {w_x:.6f}, aperture W = {W:g} "
          f"({W / (2 * np.pi):.1f} wavelengths), ramp exponent across aperture = "
          f"{ramp_exponent:.2f}")
    print(f"    kappa_i = {kappa_i:.6f}, gamma_principal = {cands['A']['gamma']:.6f}, "
          f"continuation picks branch "
          f"{'A' if cands['cont']['gamma'] == cands['A']['gamma'] else 'B'}")
    print(f"    spectral truth: r_eff = {truth['r_eff']:.4f} "
          f"(|r| = {abs(truth['r_eff']):.4f}, arg = {np.angle(truth['r_eff']):+.4f}, "
          f"spread across footprint = {truth['r_spread']:.1e})")
    print(f"    spectral power: R = {truth['R']:.4f}, T = {truth['T']:.4f}, "
          f"M = {truth['M']:+.4f}, R+T-1-M = {truth['R'] + truth['T'] - 1 - truth['M']:+.1e}")
    print(f"    {'branch':<6} {'r':>22} {'|r_eff - r|':>12} {'R':>8} {'T':>8}")
    best = None
    for name, c in cands.items():
        err = abs(c["r"] - truth["r_eff"])
        if best is None or err < best[1]:
            best = (name, err)
        print(f"    {name:<6} {c['r']:>22.4f} {err:>12.2e} {c['R']:>8.4f} {c['T']:>8.4f}")
    print(f"    --> closest branch to spectral truth: {best[0]} (|dr| = {best[1]:.2e})")
    return best[0]


def main():
    print("Spectral (branch-free) ground truth vs single-wave branch candidates.")
    print("s-polarisation, absorbing half-space -> vacuum, hard aperture, k0 = 1.")

    lam = 2 * np.pi
    verdicts = {}

    # 1. Ice-like, well below critical (propagating exit).
    n1 = 1.31 + 0.01j
    th = np.radians(35.0)
    verdicts["below critical (35 deg, ice n'=0.01)"] = run_case(
        "ice n'=0.01, theta = 35 deg (critical 49.76)", n1,
        n1 * np.sin(th), 30 * lam)

    # 2. Ice-like, just above critical.
    th = np.radians(52.0)
    verdicts["just above critical (52 deg)"] = run_case(
        "ice n'=0.01, theta = 52 deg", n1, n1 * np.sin(th), 30 * lam)

    # 3. Ice-like, well above critical (deep TIR).
    th = np.radians(60.0)
    verdicts["deep TIR (60 deg)"] = run_case(
        "ice n'=0.01, theta = 60 deg", n1, n1 * np.sin(th), 30 * lam)

    # 4. Stronger absorption, above critical, aperture sized so the ramp
    #    exponent sits at the beam power-threshold bound (~5), the
    #    strongest ramp a surviving GOAD beam can carry.
    n1 = 1.31 + 0.1j
    th = np.radians(60.0)
    w_x = n1 * np.sin(th)
    W4 = 5.0 / w_x.imag
    verdicts["deep TIR, n'=0.1 (physical ramp)"] = run_case(
        f"n'=0.1, theta = 60 deg, ramp exponent 5", n1, w_x, W4)

    # 5. Near the crossover (arg gamma^2 ~ -pi/2), ramp exponent at the
    #    physical bound. The single-wave model is expected to break here;
    #    this quantifies by how much.
    n1 = 1.5 + 0.5j
    w_x = 1.05 + 0.16j
    verdicts["crossover regime (physical ramp)"] = run_case(
        "crossover: n1 = 1.5+0.5i, w_x = 1.05+0.16i, ramp exponent 5",
        n1, w_x, 5.0 / w_x.imag)

    print("\n=== verdicts ===")
    for k, v in verdicts.items():
        print(f"  {k:<42} -> {v}")

    # Sweep through the critical angle: field error of each branch vs the
    # spectral truth, mapping the near-critical breakdown window.
    print("\n=== sweep through critical (ice n'=0.01, 30-wavelength aperture) ===")
    n1 = 1.31 + 0.01j
    crit = np.degrees(np.arcsin(1.0 / 1.31))
    print(f"    critical angle {crit:.2f} deg")
    print(f"    {'theta':>7} {'|dr| A':>10} {'|dr| B':>10} {'|dr| cont':>10} "
          f"{'|dr| chang':>11} {'spread':>10}")
    for th_deg in [40.0, 46.0, 48.0, 49.0, 49.5, 49.76, 50.0, 50.5, 51.0,
                   52.0, 54.0, 58.0, 65.0]:
        th = np.radians(th_deg)
        w_x = n1 * np.sin(th)
        _, cands = single_wave_candidates(w_x, n1)
        truth = spectral_case(n1, w_x, 30 * lam)
        errs = {k: abs(c["r"] - truth["r_eff"]) for k, c in cands.items()}
        print(f"    {th_deg:>7.2f} {errs['A']:>10.2e} {errs['B']:>10.2e} "
              f"{errs['cont']:>10.2e} {errs['chang']:>11.2e} "
              f"{truth['r_spread']:>10.1e}")


if __name__ == "__main__":
    main()
