#!/usr/bin/env python3
"""
Example 7: Advanced configuration with ConvergenceConfig.

This example demonstrates using ConvergenceConfig directly for full control
over all parameters, including optical properties and convergence settings.
"""

import goad_py as goad
from pathlib import Path

# Get path to geometry file
current_dir = Path(__file__).parent
geom_file = current_dir / ".." / ".." / ".." / "examples" / "data" / "hex.obj"

print("=" * 80)
print("Example 7: Advanced Configuration with ConvergenceConfig")
print("=" * 80)

# Create detailed configuration
config = goad.ConvergenceConfig(
    geometry=geom_file,
    mode=goad.StandardMode(n_theta=181, n_phi=181),
    convergence_targets=[
        {"variable": "asymmetry", "tolerance": 0.005, "tolerance_type": "absolute"},
        {"variable": "scatt", "tolerance": 0.01, "tolerance_type": "relative"},
    ],
    wavelength=0.633,  # Red laser (633 nm)
    particle_refr_index_re=1.5,
    particle_refr_index_im=0.01,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    batch_size=24,
    max_orientations=1000,
    min_batches=10,
    mueller_1d=True,
    mueller_2d=False,
)

print("Configuration:")
print(f"  Wavelength: {config.wavelength} μm")
print(
    f"  Particle refractive index: {config.particle_refr_index_re} + {config.particle_refr_index_im}j"
)
print(f"  Batch size: {config.batch_size}")
print(f"  Max orientations: {config.max_orientations}")

# Run convergence with explicit configuration
conv = goad.UnifiedConvergence(config)
results = conv.run()

print(results.summary())
results.save("07_advanced_config_results.json")

print("\n✓ Advanced configuration example completed!")
