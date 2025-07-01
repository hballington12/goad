"""
Sample script demonstrating the GOAD Python API with the new Results interface.

This script shows how to:
1. Create a problem with settings
2. Solve the problem
3. Access results through the unified Results object
"""

import goad_py as goad
import time

# Create settings for the simulation
print("Creating GOAD settings...")
settings = goad.Settings(
    geom_name="hex.obj",  # Geometry file to load
    euler=[0.0, 10.0, 0.0]  # Euler angles for orientation
)

# Create a problem instance
# Note: both settings and geom parameters are now optional
# If geom is None, it will be loaded from settings.geom_name
print("Creating GOAD problem...")
problem = goad.Problem(settings, geom=None)

# Solve the problem and measure computation time
print("Solving GOAD problem...")
start_time = time.time()
problem.py_solve()
end_time = time.time()
print(f"Solving took {end_time - start_time:.4f} seconds")

problem.py_print_stats()

# === NEW API: Access all results through the Results object ===
# Instead of calling individual getters on Problem, we now get the Results object
# and access all data through it. This provides better organization and discoverability.
print("\n" + "="*60)
print("ACCESSING RESULTS USING NEW API")
print("="*60)

results = problem.results

# Access the main Mueller matrix (2D scattering pattern)
mueller = results.mueller
print(f"\nMueller matrix shape: {len(mueller)}x{len(mueller[0]) if mueller else 0}")

# Access 1D integrated Mueller matrix and corresponding theta values
mueller_1d = results.mueller_1d
bins_1d = results.bins_1d
if mueller_1d and bins_1d:
    print(f"Mueller 1D s11 first value: {mueller_1d[0][0]}")
    print(f"First theta value: {bins_1d[0]} degrees")
else:
    print("1D Mueller matrix not computed")

# Access physical parameters
print(f"\nPhysical parameters:")
print(f"  Asymmetry parameter: {results.asymmetry}")
print(f"  Scattering cross section: {results.scat_cross}")
print(f"  Extinction cross section: {results.ext_cross}")
print(f"  Albedo: {results.albedo}")

# Access all parameters as a dictionary
params = results.params
print(f"\nAll parameters dict: {params}")

# Access bins information
bins = results.bins
print(f"\nBins information:")
print(f"  Number of bins: {len(bins)}")
if bins:
    print(f"  First bin (theta, phi): {bins[0]}")
    print(f"  Last bin (theta, phi): {bins[-1]}")

# Access power breakdown
powers = results.powers
print(f"\nPower breakdown:")
print(f"  Input: {powers['input']:.6f}")
print(f"  Output: {powers['output']:.6f}")
print(f"  Absorbed: {powers['absorbed']:.6f}")
print(f"  Missing: {powers['missing']:.6f}")

# Example of accessing other mueller matrices
mueller_beam = results.mueller_beam
mueller_ext = results.mueller_ext
print(f"\nOther Mueller matrices available:")
print(f"  Beam Mueller shape: {len(mueller_beam)}x{len(mueller_beam[0]) if mueller_beam else 0}")
print(f"  External diffraction Mueller shape: {len(mueller_ext)}x{len(mueller_ext[0]) if mueller_ext else 0}")
