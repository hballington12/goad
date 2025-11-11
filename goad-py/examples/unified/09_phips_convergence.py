#!/usr/bin/env python3
"""
Example 9: PHIPS Detector Convergence

Demonstrates converging on PHIPS detector bins (all except the first one).
PHIPS detectors measure scattering at 20 specific angles from 18° to 170°.
"""

from pathlib import Path

import goad_py as goad

print("=" * 80)
print("Example 9: PHIPS Detector Convergence")
print("=" * 80)

# Get paths relative to this script's location
script_dir = Path(__file__).parent
geom_file = script_dir / "../../../examples/data/hex.obj"
phips_bins_file = script_dir / "../../../phips_bins_edges.toml"

# Resolve to absolute paths
geom_file = str(geom_file.resolve())
phips_bins_file = str(phips_bins_file.resolve())

# Converge on all PHIPS detectors except the first one (indices 1-19)
# Detector 0 is at 18° (forward scatter) - we'll skip it
# Detectors 1-19 cover 26° to 170°
detector_indices = list(range(1, 20))  # [1, 2, 3, ..., 19]

print(
    f"\nConverging on {len(detector_indices)} PHIPS detectors (indices {detector_indices[0]}-{detector_indices[-1]})"
)
print(f"Detector angles: 26° to 170° (skipping 18° forward scatter)")
print(f"Target: 50% relative SEM (relaxed for demo)\n")

# Run PHIPS convergence
results = goad.run_convergence(
    geometry=geom_file,
    targets=[
        {
            "tolerance": 0.50,  # 50% relative SEM (relaxed for demo)
            "tolerance_type": "relative",
            "detector_indices": detector_indices,
        }
    ],
    mode=goad.PHIPSMode(bins_file=phips_bins_file),
    batch_size=24,  # Larger batches for faster convergence
    max_orientations=500,  # Reduced for demo
    min_batches=5,  # Reduced minimum batches
    mueller_1d=False,  # PHIPS mode doesn't support Mueller matrices
    log_file="phips_convergence.log",
)

# Print summary
print(results.summary())

# Save results
results.save("09_phips_convergence_results.json")

# Show detailed PHIPS output
print("\n" + "=" * 80)
print("PHIPS Detector DSCS Values (all 20 detectors)")
print("=" * 80)

import json

with open("09_phips_convergence_results.json", "r") as f:
    data = json.load(f)

phips_dscs = data.get("phips_dscs")
detector_angles = data.get("detector_angles")

if phips_dscs and detector_angles:
    print(
        f"\n{'Detector':<10} {'Angle':<10} {'DSCS Mean':<15} {'DSCS SEM':<15} {'Rel. SEM':<10}"
    )
    print("-" * 70)

    for i, (dscs_data, angle) in enumerate(zip(phips_dscs, detector_angles)):
        mean, sem = dscs_data[0], dscs_data[1]
        rel_sem = (sem / abs(mean) * 100) if mean != 0 else float("inf")

        # Mark which detectors were converged on
        marker = "✓" if i in detector_indices else "○"

        print(
            f"{marker} {i:<8} {angle:<10.1f} {mean:<15.4e} {sem:<15.4e} {rel_sem:<10.2f}%"
        )

print("\n" + "=" * 80)
print("Legend:")
print("  ✓ = Converged on this detector")
print("  ○ = Not a convergence target (but still computed)")
print("=" * 80)

print("\n✓ Example completed successfully!")
