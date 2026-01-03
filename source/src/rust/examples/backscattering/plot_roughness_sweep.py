"""Plot lidar ratio vs depolarization ratio for roughness sweep."""

import json
from pathlib import Path

import matplotlib.pyplot as plt

base_dir = Path(__file__).parent

# Load results
with open(base_dir / "roughness_sweep_results.json") as f:
    results = json.load(f)

# Extract data
sigmas = [r["sigma"] for r in results]
lidar_ratios = [r["lidar_ratio"] for r in results]
depol_ratios = [r["depol_ratio"] * 100 for r in results]  # convert to %

# Create plot
fig, ax = plt.subplots(figsize=(8, 6))

# Scatter plot with sigma as color
scatter = ax.scatter(
    depol_ratios,
    lidar_ratios,
    c=sigmas,
    cmap="viridis",
    s=100,
    edgecolors="black",
    linewidths=0.5,
)

# Add colorbar
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Roughness sigma", fontsize=12)

# Add labels for each point
for i, sigma in enumerate(sigmas):
    ax.annotate(
        f"{sigma}",
        (depol_ratios[i], lidar_ratios[i]),
        textcoords="offset points",
        xytext=(5, 5),
        fontsize=9,
    )

ax.set_xlabel("Depolarization Ratio (%)", fontsize=12)
ax.set_ylabel("Lidar Ratio (sr)", fontsize=12)
ax.set_title(
    "Backscattering Properties vs Surface Roughness\n(50 um hexagonal plate, 532 nm)",
    fontsize=14,
)

ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(base_dir / "roughness_sweep_plot.png", dpi=150)
print(f"Plot saved to: {base_dir / 'roughness_sweep_plot.png'}")

plt.show()
