#!/usr/bin/env python3
"""
Example 6: Ensemble convergence (multiple geometries).

This example shows how to run convergence over an ensemble of geometries
by passing a directory instead of a single file. Each batch randomly
selects a geometry file.
"""

import goad_py as goad
from pathlib import Path

# Get path to geometry directory
current_dir = Path(__file__).parent
geom_dir = current_dir / ".." / ".." / ".." / "examples" / "data"

print("=" * 80)
print("Example 6: Ensemble Convergence - Multiple Geometries")
print("=" * 80)

# Run convergence on ensemble (directory instead of file)
results = goad.run_convergence(
    geometry=geom_dir,  # Directory automatically triggers ensemble mode
    targets=["asymmetry", "scatt"],
    tolerance=0.2,  # 20% relative tolerance
    tolerance_type="relative",
    batch_size=12,
    max_orientations=1000,
)

print(results.summary())
print(f"\nIs ensemble: {results.is_ensemble}")

results.save("06_ensemble_convergence_results.json")

print("\n✓ Ensemble convergence completed!")
