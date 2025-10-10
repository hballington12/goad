"""
Example: Converge only backscattering S11 and plot full angular distribution.

Demonstrates:
- Convergence on specific theta bin (backscattering only)
- Plotting full S11 distribution with SEM
- Highlighting the converged region
"""

import goad_py as goad
import matplotlib

matplotlib.use("Agg")  # Suppress plot display
import matplotlib.pyplot as plt
import numpy as np
import os

print("Backscattering Convergence Example")
print("=" * 50)

# Setup geometry path
current_dir = os.path.dirname(os.path.abspath(__file__))
geom_path = os.path.join(current_dir, "..", "..", "..", "examples", "data", "hex.obj")
geom_path = os.path.abspath(geom_path)

print(f"Geometry: {geom_path}")

# Create binning scheme
binning = goad.BinningScheme.simple(181, 181)
print("Binning: Simple 181x181 (181 theta bins)")

# Create settings
settings = goad.Settings(
    geom_path=geom_path,
    wavelength=0.532,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    binning=binning,
)

# Define convergence criteria - only backscattering bin (last bin, index 180)
convergables = [
    goad.Convergable(
        "S11", tolerance_type="relative", tolerance=0.10, theta_indices=[180]
    )
]

print(
    "\nConvergence criterion: S11 backscattering bin (θ=180°) with 10% relative tolerance"
)
print("Running convergence...")

# Run convergence
conv = goad.Convergence(
    settings=settings,
    convergables=convergables,
    batch_size=12,
    max_orientations=500,
    min_batches=10,
    mueller_1d=True,
    mueller_2d=False,
)

results = conv.run()

# Print summary
print("\n" + "=" * 50)
print("Convergence Results:")
print(f"Converged: {results.converged}")
print(f"Total orientations: {results.n_orientations}")

# Extract S11 data
s11_values = results.values["S11"]
s11_sem = results.sem_values["S11"]
theta_bins = np.linspace(0, 180, len(s11_values))

print(f"S11 shape: {s11_values.shape}")
print(f"Backscattering S11: {s11_values[180]:.4g} ± {s11_sem[180]:.4g}")

# Create plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot S11 with error bars
ax1.errorbar(
    theta_bins,
    s11_values,
    yerr=1.96 * s11_sem,
    fmt="o-",
    markersize=3,
    linewidth=1,
    capsize=2,
    label="S11 ± 95% CI",
    alpha=0.7,
)

# Highlight backscattering region
ax1.axvline(
    x=180, color="red", linestyle="--", linewidth=2, label="Backscattering (converged)"
)
ax1.plot(
    180, s11_values[180], "ro", markersize=10, label=f"θ=180°: {s11_values[180]:.4g}"
)

ax1.set_xlabel("Scattering Angle θ (degrees)")
ax1.set_ylabel("S11 (Phase Function)")
ax1.set_title(
    f"Mueller S11 Element (backscatter converged after {results.n_orientations} orientations)"
)
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_yscale("log")

# Plot relative SEM
relative_sem = np.where(s11_values != 0, s11_sem / np.abs(s11_values), 0) * 100
ax2.plot(
    theta_bins, relative_sem, "o-", markersize=3, linewidth=1, color="blue", alpha=0.7
)
ax2.axhline(y=10, color="k", linestyle="--", label="10% tolerance")
ax2.axvline(
    x=180, color="red", linestyle="--", linewidth=2, label="Backscattering (converged)"
)
ax2.plot(180, relative_sem[180], "ro", markersize=10)

ax2.set_xlabel("Scattering Angle θ (degrees)")
ax2.set_ylabel("Relative SEM (%)")
ax2.set_title("Convergence Quality Across Scattering Angles")
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_ylim(0, min(100, np.max(relative_sem[relative_sem < np.inf]) * 1.1))

plt.tight_layout()

# Save figure
output_file = "backscatter_convergence.png"
plt.savefig(output_file, dpi=150, bbox_inches="tight")
print(f"\nPlot saved to: {output_file}")

plt.close()
