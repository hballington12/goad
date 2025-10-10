#!/usr/bin/env python3
"""
Example script demonstrating PHIPS detector DSCS convergence.

This script:
1. Uses Custom binning with PHIPS detector geometry (phips_bins_edges.toml)
2. Runs orientation averaging until PHIPS DSCS values converge at all 20 detectors
3. Tracks mean DSCS and SEM at each detector
4. Plots convergence behavior and final DSCS vs scattering angle
"""

import goad_py as goad
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
GEOMETRY_FILE = "../examples/data/hex.obj"
PHIPS_BINS_FILE = "../../phips_bins_edges.toml"
WAVELENGTH = 0.532  # microns
PARTICLE_REFR_INDEX = 1.31 + 0.0j  # Ice refractive index at 532nm
BATCH_SIZE = 24
MAX_ORIENTATIONS = 10_000
TOLERANCE = 0.25  # 25% relative tolerance


def main():
    print("=" * 60)
    print("PHIPS Detector DSCS Convergence Example")
    print("=" * 60)

    # Check that PHIPS bins file exists
    phips_bins_path = Path(PHIPS_BINS_FILE)
    if not phips_bins_path.exists():
        print(f"\nError: PHIPS bins file not found at {phips_bins_path}")
        print("Please run phips_detector_angles_edges.py first to generate the file.")
        return

    print(f"\nGeometry: {GEOMETRY_FILE}")
    print(f"PHIPS bins: {PHIPS_BINS_FILE}")
    print(f"Wavelength: {WAVELENGTH} μm")
    print(f"Batch size: {BATCH_SIZE}")
    print(f"Tolerance: {TOLERANCE * 100:.0f}% relative")

    # Create settings with Custom binning
    print("\nLoading PHIPS Custom binning scheme...")

    # Load the PHIPS bins edges file
    import toml

    with open(phips_bins_path, "r") as f:
        bins_data = toml.load(f)

    # Convert to format expected by BinningScheme.custom()
    custom_bins = bins_data["bins"]
    print(f"Loaded {len(custom_bins)} custom bins")

    # Create binning scheme
    binning = goad.BinningScheme.custom(custom_bins)

    # Create settings
    settings = goad.Settings(
        GEOMETRY_FILE,
        wavelength=WAVELENGTH,
        particle_refr_index_re=PARTICLE_REFR_INDEX.real,
        particle_refr_index_im=PARTICLE_REFR_INDEX.imag,
        binning=binning,
    )

    # Create PHIPS convergable
    convergable = goad.PHIPSConvergable(tolerance_type="relative", tolerance=TOLERANCE)

    # Create PHIPS convergence study
    phips_conv = goad.PHIPSConvergence(
        settings=settings,
        convergable=convergable,
        batch_size=BATCH_SIZE,
        max_orientations=MAX_ORIENTATIONS,
        min_batches=10,
    )

    # Run convergence
    print("\nStarting convergence study...")
    results = phips_conv.run()

    # Extract results
    phips_dscs = results.values["phips_dscs"]  # Array of 20 values
    phips_dscs_sem = results.sem_values["phips_dscs"]  # Array of 20 SEMs
    detector_angles = phips_conv.detector_centers  # 20 detector angles

    # Plot results
    print("\nGenerating plots...")

    # Figure 1: DSCS vs scattering angle with error bars
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Filter out NaN detectors
    valid_mask = ~np.isnan(phips_dscs)
    valid_angles = detector_angles[valid_mask]
    valid_dscs = phips_dscs[valid_mask]
    valid_sem = phips_dscs_sem[valid_mask]

    # Top panel: DSCS with error bars
    ax1.errorbar(
        valid_angles,
        valid_dscs,
        yerr=valid_sem,
        fmt="o-",
        capsize=5,
        capthick=2,
        markersize=6,
        label="Mean DSCS",
    )
    ax1.set_xlabel("Scattering angle θ (degrees)", fontsize=12)
    ax1.set_ylabel("DSCS (m²/sr)", fontsize=12)
    ax1.set_title(
        f"PHIPS Detector DSCS ({results.n_orientations} orientations)", fontsize=14
    )
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale("log")
    ax1.legend()

    # Bottom panel: Relative SEM
    with np.errstate(divide="ignore", invalid="ignore"):
        relative_sem = valid_sem / np.abs(valid_dscs) * 100

    ax2.plot(valid_angles, relative_sem, "o-", markersize=6, color="red")
    ax2.axhline(
        y=TOLERANCE * 100,
        color="green",
        linestyle="--",
        label=f"Target: {TOLERANCE * 100:.0f}%",
    )
    ax2.set_xlabel("Scattering angle θ (degrees)", fontsize=12)
    ax2.set_ylabel("Relative SEM (%)", fontsize=12)
    ax2.set_title("Convergence Quality", fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("phips_dscs_convergence.png", dpi=150, bbox_inches="tight")
    print(f"Saved: phips_dscs_convergence.png")

    # Figure 2: Convergence history (if available)
    if results.convergence_history:
        fig2, ax = plt.subplots(figsize=(10, 6))

        n_orients = [h[0] for h in results.convergence_history]
        worst_sem = [h[2] * 100 for h in results.convergence_history]  # Convert to %

        ax.plot(n_orients, worst_sem, "o-", markersize=4, linewidth=2)
        ax.axhline(
            y=TOLERANCE * 100,
            color="green",
            linestyle="--",
            label=f"Target: {TOLERANCE * 100:.0f}%",
        )
        ax.set_xlabel("Number of orientations", fontsize=12)
        ax.set_ylabel("Worst-case relative SEM (%)", fontsize=12)
        ax.set_title("Convergence Progress", fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.tight_layout()
        plt.savefig("phips_convergence_history.png", dpi=150, bbox_inches="tight")
        print(f"Saved: phips_convergence_history.png")

    # Summary statistics
    print("\n" + "=" * 60)
    print("Summary Statistics")
    print("=" * 60)
    print(f"Total orientations: {results.n_orientations}")
    print(f"Converged: {'Yes' if results.converged else 'No'}")
    print(f"Valid detectors: {np.sum(valid_mask)} / {len(phips_dscs)}")

    if np.any(valid_mask):
        print(
            f"\nDSCS range: [{np.min(valid_dscs):.4e}, {np.max(valid_dscs):.4e}] m²/sr"
        )
        print(f"Mean relative SEM: {np.mean(relative_sem):.2f}%")
        print(f"Max relative SEM: {np.max(relative_sem):.2f}%")

    # Save data to file
    output_file = "phips_dscs_results.txt"
    with open(output_file, "w") as f:
        f.write("# PHIPS Detector DSCS Results\n")
        f.write(f"# Geometry: {GEOMETRY_FILE}\n")
        f.write(f"# Wavelength: {WAVELENGTH} μm\n")
        f.write(f"# Orientations: {results.n_orientations}\n")
        f.write(f"# Converged: {results.converged}\n")
        f.write("#\n")
        f.write("# detector_idx  theta_deg  dscs_m2sr  sem_m2sr  rel_sem_%\n")

        for i in range(len(phips_dscs)):
            theta = detector_angles[i]
            if valid_mask[i]:
                dscs = phips_dscs[i]
                sem = phips_dscs_sem[i]
                rel_sem = sem / abs(dscs) * 100
                f.write(
                    f"{i:2d}  {theta:7.2f}  {dscs:.6e}  {sem:.6e}  {rel_sem:7.2f}\n"
                )
            else:
                f.write(f"{i:2d}  {theta:7.2f}  NaN  NaN  NaN\n")

    print(f"\nSaved data to: {output_file}")
    print("\nDone!")


if __name__ == "__main__":
    main()
