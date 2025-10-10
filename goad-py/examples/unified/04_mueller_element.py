#!/usr/bin/env python3
"""
Example 4: Mueller matrix element convergence.

This example shows how to converge on Mueller matrix elements (like S11)
across all scattering angles.
"""

import goad_py as goad
from pathlib import Path

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 4: Mueller S11 Element Convergence")
print("=" * 80)

# Run convergence on S11 Mueller element
results = goad.run_convergence(
    geometry=geom_file,
    targets="S11",
    tolerance=0.4,  # 40% relative tolerance
    tolerance_type="relative",
    batch_size=12,
    max_orientations=1000,
    min_batches=10,
)

print(results.summary())
print(f"\nS11 values shape: {results.values['S11'].shape}")
print(f"S11 SEM shape: {results.sem_values['S11'].shape}")

results.save("04_mueller_element_results.json")

print("\n✓ Mueller element convergence completed!")
