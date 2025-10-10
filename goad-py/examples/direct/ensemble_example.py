"""
Example script demonstrating EnsembleConvergence for multi-geometry averaging.

This script shows how to run convergence analysis over an ensemble of
particle geometries, randomly sampling from a directory each batch.
"""

import goad_py as goad
import os

print("Ensemble Convergence Example")
print("=" * 50)

# Get absolute path to geometry directory
current_dir = os.path.dirname(os.path.abspath(__file__))
geom_dir = os.path.join(current_dir, "..", "..", "..", "examples", "data")
geom_dir = os.path.abspath(geom_dir)

print(f"Using geometry directory: {geom_dir}")

# Create binning scheme
binning = goad.BinningScheme.simple(181, 181)
print("Created binning scheme: Simple 181x181")

# Create settings (geom_path will be overridden per batch by EnsembleConvergence)
# Use a valid geometry file for Settings initialization (will be overridden)
dummy_geom = os.path.join(geom_dir, "hex.obj")
settings = goad.Settings(
    geom_path=dummy_geom,  # Initial valid path, will be overridden per batch
    wavelength=0.532,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    binning=binning,
    beam_power_threshold=0.005,
)

# Define convergence criteria
convergables = [
    goad.Convergable(
        "asymmetry", tolerance_type="relative", tolerance=0.05
    ),  # 5% relative
    goad.Convergable("scatt", tolerance_type="relative", tolerance=0.05),  # 5% relative
]

# Create and run ensemble convergence
print("\nStarting ensemble convergence study...")
ensemble_conv = goad.EnsembleConvergence(
    settings=settings,
    convergables=convergables,
    geom_dir=geom_dir,
    batch_size=12,  # Smaller batches for faster testing
    max_orientations=1000,  # Reasonable limit for testing
    min_batches=5,
    mueller_1d=True,
    mueller_2d=False,
)

results = ensemble_conv.run()

# Print final results
print("\n" + "=" * 50)
print("Final Ensemble Results:")
print(f"Converged: {results.converged}")
print(f"Total orientations: {results.n_orientations}")
print(f"\nFinal values:")
for var, value in results.values.items():
    sem = results.sem_values[var]
    print(f"  {var}: {value:.6f} ± {sem:.6f}")

if results.warning:
    print(f"\nWarning: {results.warning}")
