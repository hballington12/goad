"""Task 2: analytic lossy wedge, reproducing Chang (2005) section 4.

Traces a single ray through a wedge of internal angle gamma = 30 deg and
index n_tilde = 1.5 + i kappa using only the reduced-wave-vector transfer
rule, and reproduces the content of Chang Figs. 3-7:

  fig3: output angle theta_o vs input angle theta_i
  fig4: apparent indices (N, K) of the output wave vs theta_i
  fig5: s and p transmittances (no path loss) vs theta_i
  fig6: dimensionless attenuation coefficient d, standard vs rigorous
  fig7: total transmittance including path loss, entry at L = 2 lambda/pi

This validates the interface-transfer chain end to end with no beam
tracing involved. Figures are written to ./figs/ and not displayed.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from lossy import transfer, fresnel_general

GAMMA = np.radians(30.0)
N_RE = 1.5
KAPPAS = [0.0, 0.1, 0.2, 0.5, 1.0, 1.5]
K0L = 4.0  # k0 * L for L = 2 lambda / pi

FIGS = Path(__file__).parent / "figs"

ENTRY_NORMAL = np.array([1.0, 0.0, 0.0])  # into the glass
EXIT_NORMAL = np.array([np.cos(GAMMA), np.sin(GAMMA), 0.0])  # out of the glass
EXIT_TANGENT = np.array([np.sin(GAMMA), -np.cos(GAMMA), 0.0])


def trace(theta_i, kappa):
    """One ray through the wedge. Returns a dict of everything plotted."""
    nt = complex(N_RE, kappa)
    e_i = np.array([np.cos(theta_i), np.sin(theta_i), 0.0])
    w_i = e_i.astype(np.complex128)

    # Entry face.
    w_in, g1 = transfer(w_i, ENTRY_NORMAL, nt)
    kap_i1 = w_i @ ENTRY_NORMAL.astype(np.complex128)
    rs1, _, rp1, _ = fresnel_general(kap_i1, g1, 1.0, nt)

    # Exit face: physical (continuation-rule) branch, plus the Chang
    # eq.-(12) same-sign rule purely for comparison with his Figs. 5/7.
    # The Chang-rule kappa_t satisfies neither Maxwell root of the
    # dispersion relation; see the theory note's branch section.
    w_out, g2 = transfer(w_in, EXIT_NORMAL, 1.0)
    kap_i2 = w_in @ EXIT_NORMAL.astype(np.complex128)
    rs2, _, rp2, _ = fresnel_general(kap_i2, g2, nt, 1.0)
    g2_chang = abs(g2.real) + 1j * abs(g2.imag)
    rs2_c, _, rp2_c, _ = fresnel_general(kap_i2, g2_chang, nt, 1.0)

    # Output angle: signed, from the exit-face normal.
    re_out = np.real(w_out)
    n_out = np.linalg.norm(re_out)
    sin_o = (re_out @ EXIT_TANGENT) / n_out
    cos_o = (re_out @ EXIT_NORMAL) / n_out
    theta_o = np.arctan2(sin_o, cos_o)

    # Apparent indices of the output wave (in air).
    n_app = n_out
    k_app = np.linalg.norm(np.imag(w_out))

    # Path loss: entry point at distance L from the apex along the entry
    # face; path length l to the exit face plane through the apex.
    e_in = np.real(w_in) / np.linalg.norm(np.real(w_in))
    a_in = np.imag(w_in)
    cos_exit = e_in @ EXIT_NORMAL
    l_over_L = np.sin(GAMMA) / cos_exit  # l / L
    d_rig = (a_in @ e_in) * l_over_L
    d_std = kappa * l_over_L

    return {
        "theta_o": theta_o,
        "n_app": n_app,
        "k_app": k_app,
        "ts": (1 - abs(rs1) ** 2) * (1 - abs(rs2) ** 2),
        "tp": (1 - abs(rp1) ** 2) * (1 - abs(rp2) ** 2),
        "ts_chang": (1 - abs(rs1) ** 2) * (1 - abs(rs2_c) ** 2),
        "tp_chang": (1 - abs(rp1) ** 2) * (1 - abs(rp2_c) ** 2),
        "d_rig": d_rig,
    }


def trace_standard(theta_i, kappa):
    """The 'standard method' ray: real-part index for direction and
    Fresnel, absorption only as a path-loss exponential.

    Direction from real Snell with n = 1.5; classical TIR at the exit
    face; Fresnel coefficients from the real index (independent of
    kappa); d_std = kappa * l / L along the standard path.
    """
    n = N_RE
    sin_t = np.sin(theta_i) / n
    theta_t = np.arcsin(sin_t)
    e_in = np.array([np.cos(theta_t), np.sin(theta_t), 0.0])

    cos_exit = e_in @ EXIT_NORMAL  # cos of internal incidence at exit face
    l_over_L = np.sin(GAMMA) / cos_exit
    d_std = kappa * l_over_L

    # Real-index Fresnel at both faces.
    def fresnel_real(n1, n2, cos_i):
        s2 = 1.0 - cos_i**2
        arg = 1.0 - (n1 / n2) ** 2 * s2
        if arg <= 0.0:
            return 1.0, 1.0  # TIR: |r_s| = |r_p| = 1
        cos_t = np.sqrt(arg)
        r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
        r_p = (n2 * cos_i - n1 * cos_t) / (n1 * cos_t + n2 * cos_i)
        return abs(r_s), abs(r_p)

    rs1, rp1 = fresnel_real(1.0, n, np.cos(theta_i))
    rs2, rp2 = fresnel_real(n, 1.0, cos_exit)

    return {
        "ts": (1 - rs1**2) * (1 - rs2**2),
        "tp": (1 - rp1**2) * (1 - rp2**2),
        "d_std": d_std,
    }


def main():
    FIGS.mkdir(exist_ok=True)
    thetas = np.radians(np.linspace(-89.0, 89.0, 500))
    deg = np.degrees(thetas)

    results = {k: [trace(t, k) for t in thetas] for k in KAPPAS}
    std_results = {k: [trace_standard(t, k) for t in thetas] for k in KAPPAS}

    def series(kappa, key):
        return np.array([r[key] for r in results[kappa]])

    def std_series(kappa, key):
        return np.array([r[key] for r in std_results[kappa]])

    # Numeric checks.
    theta_tir = np.degrees(np.arcsin(np.sqrt(N_RE**2 - 1) * np.sin(GAMMA) - np.cos(GAMMA)))
    print(f"lossless TIR input angle (Chang, unnumbered): {theta_tir:.2f} deg "
          f"(paper: -17.88)")

    i40 = np.argmin(abs(deg + 40.0))
    print("\nT = 1 - |r|^2 (Chang's metric, no path loss) at theta_i = -40 deg "
          "(lossless TIR zone):")
    print(f"  {'kappa':>6} {'physical branch':>16} {'Chang eq.12 rule':>17}   (Chang Fig 5 stays in [0,1])")
    for k in KAPPAS:
        ts = series(k, "ts")[i40]
        tc = series(k, "ts_chang")[i40]
        print(f"  {k:>6} {ts:>+16.4f} {tc:>+17.4f}")
    print("  The physical (continuation-rule) branch admits |r| > 1 here, so 1 - |r|^2")
    print("  goes negative: it is not a flux ratio for complex kappa_i (the energy books")
    print("  close via the Dupertuis II mixed term, R + T = 1 + M). Chang's same-sign")
    print("  rule stays in [0,1] and reproduces his Figs. 5/7, but its kappa_t violates")
    print("  the vacuum dispersion relation; overlay curves are for comparison only.")

    ok = True
    for k in KAPPAS[1:]:
        if not np.all(std_series(k, "d_std") >= series(k, "d_rig") - 1e-12):
            ok = False
    print(f"\nd_standard >= d_rigorous everywhere: {ok} "
          f"(standard method overestimates path loss)")
    d_max = std_series(1.5, "d_std").max()
    print(f"max d_standard for kappa = 1.5: {d_max:.3f} "
          f"(Chang Fig 6 dotted curve reaches ~2.4-2.5)")

    # fig3: output angle.
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for k in KAPPAS:
        ax.plot(deg, np.degrees(series(k, "theta_o")), label=f"$\\kappa$ = {k}")
    ax.set_xlabel(r"input angle $\theta_i$ [deg]")
    ax.set_ylabel(r"output angle $\theta_o$ [deg]")
    ax.axvline(theta_tir, color="gray", lw=0.6, ls=":")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "wedge_fig3_output_angle.png", dpi=150)
    plt.close(fig)

    # fig4: apparent indices of the output wave.
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for k in KAPPAS:
        axes[0].plot(deg, series(k, "n_app"), label=f"$\\kappa$ = {k}")
        axes[1].plot(deg, series(k, "k_app"), label=f"$\\kappa$ = {k}")
    axes[0].set_xlabel(r"$\theta_i$ [deg]")
    axes[0].set_ylabel(r"output apparent index $N$")
    axes[1].set_xlabel(r"$\theta_i$ [deg]")
    axes[1].set_ylabel(r"output apparent index $K$")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "wedge_fig4_apparent_indices.png", dpi=150)
    plt.close(fig)

    # fig5: transmittances without path loss. Solid: physical branch;
    # thin dashed: Chang same-sign rule (matches the paper's figure).
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for i, k in enumerate(KAPPAS):
        color = f"C{i}"
        axes[0].plot(deg, series(k, "ts"), color=color, label=f"$\\kappa$ = {k}")
        axes[0].plot(deg, series(k, "ts_chang"), color=color, ls="--", lw=0.8)
        axes[1].plot(deg, series(k, "tp"), color=color)
        axes[1].plot(deg, series(k, "tp_chang"), color=color, ls="--", lw=0.8)
    axes[0].set_xlabel(r"$\theta_i$ [deg]")
    axes[0].set_ylabel(r"$T_s$ (no path loss; dashed: Chang eq. 12 rule)")
    axes[1].set_xlabel(r"$\theta_i$ [deg]")
    axes[1].set_ylabel(r"$T_p$ (no path loss)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "wedge_fig5_transmittance.png", dpi=150)
    plt.close(fig)

    # fig6: attenuation coefficient, standard (dotted) vs rigorous (solid).
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for i, k in enumerate(KAPPAS[1:], start=1):
        color = f"C{i}"
        ax.plot(deg, series(k, "d_rig"), color=color, label=f"$\\kappa$ = {k}")
        ax.plot(deg, std_series(k, "d_std"), color=color, ls=":")
    ax.set_xlabel(r"input angle $\theta_i$ [deg]")
    ax.set_ylabel(r"attenuation coefficient $d$ (solid: rigorous, dotted: standard)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "wedge_fig6_attenuation.png", dpi=150)
    plt.close(fig)

    # fig7: total transmittance with path loss, k0 L = 4, log scale.
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for i, k in enumerate(KAPPAS):
        color = f"C{i}"
        ts_rig = series(k, "ts") * np.exp(-2 * K0L * series(k, "d_rig"))
        tp_rig = series(k, "tp") * np.exp(-2 * K0L * series(k, "d_rig"))
        ts_std = std_series(k, "ts") * np.exp(-2 * K0L * std_series(k, "d_std"))
        tp_std = std_series(k, "tp") * np.exp(-2 * K0L * std_series(k, "d_std"))
        ts_ch = series(k, "ts_chang") * np.exp(-2 * K0L * series(k, "d_rig"))
        tp_ch = series(k, "tp_chang") * np.exp(-2 * K0L * series(k, "d_rig"))
        axes[0].plot(deg, ts_rig, color=color, label=f"$\\kappa$ = {k}")
        axes[0].plot(deg, ts_std, color=color, ls=":")
        axes[0].plot(deg, ts_ch, color=color, ls="--", lw=0.8)
        axes[1].plot(deg, tp_rig, color=color)
        axes[1].plot(deg, tp_std, color=color, ls=":")
        axes[1].plot(deg, tp_ch, color=color, ls="--", lw=0.8)
    for ax in axes:
        ax.set_yscale("log")
        ax.set_ylim(1e-6, 1.5)
        ax.set_xlabel(r"$\theta_i$ [deg]")
    axes[0].set_ylabel(r"total $T_s$ (solid: physical, dotted: standard, dashed: Chang)")
    axes[1].set_ylabel(r"total $T_p$")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "wedge_fig7_total_transmittance.png", dpi=150)
    plt.close(fig)

    print(f"\nfigures written to {FIGS}/wedge_fig[3-7]_*.png")


if __name__ == "__main__":
    main()
