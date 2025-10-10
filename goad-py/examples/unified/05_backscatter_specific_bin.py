#!/usr/bin/env python3
"""
Example 5: Convergence on specific theta bins.

This example demonstrates converging on a specific scattering angle
(backscattering at 180 degrees) rather than all angles.
"""

import goad_py as goad
from pathlib import Path

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 5: Backscattering Convergence (theta=180°)")
print("=" * 80)

# Run convergence on specific theta bin (backscattering)
results = goad.run_convergence(
    geometry=geom_file,
    targets=[
        {
            "variable": "S11",
            "tolerance": 0.50,  # 50% relative tolerance
            "tolerance_type": "relative",
            "theta_indices": [180],  # Only converge backscattering bin
        }
    ],
    batch_size=12,
    max_orientations=1000,
)

print(results.summary())
print(f"\nBackscattering S11 (theta=180°): {results.values['S11'][180]:.6e}")
print(f"Backscattering SEM: {results.sem_values['S11'][180]:.6e}")

results.save("05_backscatter_specific_bin_results.json")

print("\n✓ Backscatter convergence completed!")
