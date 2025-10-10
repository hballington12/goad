#!/usr/bin/env python3
"""
Example 3: Convergence on multiple integrated parameters.

This example demonstrates converging on multiple targets simultaneously,
ensuring all parameters meet their convergence criteria.
"""

import goad_py as goad
from pathlib import Path

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 3: Multiple Targets - Asymmetry and Scattering Cross-Section")
print("=" * 80)

# Run convergence on multiple targets
results = goad.run_convergence(
    geometry=geom_file,
    targets=["asymmetry", "scatt"],
    tolerance=0.2,  # 20% relative tolerance for both
    tolerance_type="relative",
    batch_size=12,
    max_orientations=1000,
)

print(results.summary())
results.save("03_multiple_targets_results.json")

print("\n✓ Multiple targets converged!")
print(f"  Asymmetry: {results.values['asymmetry']:.6f}")
print(f"  Scattering cross-section: {results.values['scatt']:.6e}")
