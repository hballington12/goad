#!/usr/bin/env python3
"""
Example 1: Simple convergence on a single integrated parameter.

This example demonstrates the simplest use of the unified convergence API:
converging on the asymmetry parameter for a single geometry.
"""

from pathlib import Path

import goad_py as goad

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 1: Simple Convergence - Asymmetry Parameter")
print("=" * 80)

# Run convergence with minimal configuration
results = goad.run_convergence(
    geometry=geom_file,
    targets="asymmetry",
    tolerance=0.1,  # 10% relative tolerance
    tolerance_type="relative",
    batch_size=12,
    max_orientations=1000,
    log_file="convergence.log",
)

# Print results
print(results.summary())

# Save to JSON
results.save("01_simple_asymmetry_results.json")

print("\n✓ Example completed successfully!")
