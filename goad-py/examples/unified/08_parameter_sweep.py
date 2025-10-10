#!/usr/bin/env python3
"""
Example 8: Parameter sweep (multiple wavelengths).

This example demonstrates running convergence studies across multiple
parameter values (e.g., wavelength sweep) using run_convergence_sweep.
"""

import goad_py as goad
from pathlib import Path

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 8: Parameter Sweep - Multiple Wavelengths")
print("=" * 80)

# Define wavelengths to sweep
wavelengths = [0.532, 0.633]  # Green and red lasers

print(f"Running convergence for {len(wavelengths)} wavelengths:")
for wl in wavelengths:
    print(f"  - {wl:.3f} μm")

# Create configurations for each wavelength
configs = [
    goad.ConvergenceConfig(
        geometry=geom_file,
        mode=goad.StandardMode(n_theta=181, n_phi=181),
        convergence_targets=[
            {"variable": "asymmetry", "tolerance": 0.05, "tolerance_type": "relative"},
        ],
        wavelength=wl,
        batch_size=12,
        max_orientations=300,
        min_batches=5,
        mueller_1d=True,
    )
    for wl in wavelengths
]

# Run sweep
results_list = goad.run_convergence_sweep(configs)

# Print sweep results
print("\n" + "=" * 80)
print("Sweep Results:")
print("=" * 80)
for wl, result in zip(wavelengths, results_list):
    asymmetry = result.values["asymmetry"]
    asymmetry_sem = result.sem_values["asymmetry"]
    print(f"  λ={wl:.3f} μm: asymmetry = {asymmetry:.6f} ± {asymmetry_sem:.6f}")
    result.save(f"08_parameter_sweep_wl{wl:.3f}_results.json")

print("\n✓ Parameter sweep completed!")
