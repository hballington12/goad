#!/usr/bin/env python3
"""
Example 2: Convergence with custom interval binning.

This example shows how to use interval binning to focus angular resolution
where you need it most (e.g., fine resolution near forward scattering).
"""

from pathlib import Path

import goad_py as goad

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 2: Convergence with Interval Binning")
print("=" * 80)

# Create interval binning scheme with fine resolution near forward scattering
interval_binning = goad.BinningScheme.interval(
    thetas=[0.0, 30.0, 180.0],  # Split at 30 degrees
    theta_spacings=[1.0, 5.0],  # 1° spacing for 0-30°, 5° spacing for 30-180°
    phis=[0.0, 360.0],  # Full phi range
    phi_spacings=[5.0],  # 5° phi spacing
)

print("Created interval binning scheme:")
print("  Theta intervals: [0.0, 30.0, 180.0]")
print("  Theta spacings: [1.0°, 5.0°]")
print("  (Fine resolution in forward scattering, coarser elsewhere)")

# Run convergence with custom binning via StandardMode
results = goad.run_convergence(
    geometry=geom_file,
    targets="asymmetry",
    tolerance=0.1,
    tolerance_type="relative",
    mode=goad.StandardMode(binning=interval_binning),
    batch_size=12,
    max_orientations=1000,
)

print(results.summary())
results.save("02_interval_binning_results.json")

print("\n✓ Interval binning works with unified API!")
