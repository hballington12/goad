"""
Example script demonstrating MultiProblem Python API without config file dependencies.

This script shows how to:
1. Create binning schemes explicitly
2. Create orientation schemes explicitly
3. Use Settings with absolute geometry paths
4. Use MultiProblem for multi-orientation averaging with no external dependencies
"""

import goad_py as goad
import os

print("MultiProblem Example - Config-Free API")
print("=" * 50)

# Get absolute path to geometry file
current_dir = os.path.dirname(os.path.abspath(__file__))
geom_path = os.path.join(current_dir, "..", "examples", "data", "hex.obj")
geom_path = os.path.abspath(geom_path)

print(f"Using geometry file: {geom_path}")

# Example 1: Uniform orientation scheme with custom settings
print("\n1. Creating MultiProblem with uniform orientations...")

# Create custom binning scheme
binning = goad.BinningScheme.simple(181, 181)  # 181x181 angular grid
print("Created binning scheme: Simple 181x181")

# Create uniform orientation (100 random orientations)
uniform_orientation = goad.create_uniform_orientation(100)
print("Created uniform orientation: 100 random orientations")

# Create settings with all parameters explicit (no config file dependencies)
settings = goad.Settings(
    geom_path=geom_path,
    wavelength=0.532,  # 532nm green laser
    particle_refr_index_re=1.31,  # glass refractive index
    particle_refr_index_im=0.0,   # no absorption
    medium_refr_index_re=1.0,     # air/vacuum
    medium_refr_index_im=0.0,
    orientation=uniform_orientation,
    binning=binning,
    beam_power_threshold=0.005,
    cutoff=0.99,
    max_rec=10,
    max_tir=10
)

# Create and solve MultiProblem
print("Creating MultiProblem with uniform orientations...")
multi_problem = goad.MultiProblem(settings)
print(f"Number of orientations: {multi_problem.num_orientations}")

print("Solving MultiProblem (this averages over all orientations)...")
multi_problem.py_solve()

# Access averaged results
results = multi_problem.results
print(f"Averaged Mueller matrix shape: {len(results.mueller)}x{len(results.mueller[0])}")
print(f"Averaged asymmetry parameter: {results.asymmetry}")

# Example 2: Discrete orientation scheme with custom binning
print("\n2. Creating MultiProblem with discrete orientations...")

# Create custom binning scheme with higher resolution in forward direction
interval_binning = goad.BinningScheme.interval(
    thetas=[0.0, 30.0, 180.0],      # Split at 30 degrees
    theta_spacings=[1.0, 5.0],       # Fine resolution near forward, coarse elsewhere
    phis=[0.0, 360.0],              # Full phi range
    phi_spacings=[5.0]              # 5 degree phi spacing
)
print("Created interval binning scheme")

# Create specific Euler angles
euler1 = goad.Euler(0.0, 0.0, 0.0)      # No rotation
euler2 = goad.Euler(30.0, 30.0, 30.0)   # 30 degree rotations
euler3 = goad.Euler(45.0, 60.0, 90.0)   # Mixed rotations

# Create discrete orientation scheme
discrete_orientation = goad.create_discrete_orientation([euler1, euler2, euler3])
print("Created discrete orientation: 3 specific orientations")

# Create settings with discrete orientation and custom binning
settings2 = goad.Settings(
    geom_path=geom_path,
    wavelength=0.633,  # 633nm red laser
    particle_refr_index_re=1.5,   # different material
    particle_refr_index_im=0.01,  # slight absorption
    orientation=discrete_orientation,
    binning=interval_binning
)

# Create and solve MultiProblem
print("Creating MultiProblem with discrete orientations...")
multi_problem2 = goad.MultiProblem(settings2)
print(f"Number of orientations: {multi_problem2.num_orientations}")

print("Solving MultiProblem with specific orientations...")
multi_problem2.py_solve()

# Access averaged results
results2 = multi_problem2.results
print(f"Averaged Mueller matrix shape: {len(results2.mueller)}x{len(results2.mueller[0])}")
print(f"Averaged asymmetry parameter: {results2.asymmetry}")

# Example 3: Custom binning with specific angles
print("\n3. Creating MultiProblem with custom binning...")

# Create custom binning scheme with specific angle pairs
custom_bins = [
    (0.0, 0.0),     # Forward scattering
    (10.0, 0.0),    # Near-forward
    (30.0, 0.0),    # Small angle
    (90.0, 0.0),    # Side scattering
    (150.0, 0.0),   # Backscattering region
    (180.0, 0.0)    # Exact backscattering
]
custom_binning = goad.BinningScheme.custom(custom_bins)
print(f"Created custom binning with {len(custom_bins)} specific angles")

# Simple settings with default values
settings3 = goad.Settings(
    geom_path=geom_path,
    binning=custom_binning
    # All other parameters use defaults
)

multi_problem3 = goad.MultiProblem(settings3)
print(f"Number of orientations: {multi_problem3.num_orientations}")

print("Solving MultiProblem with custom binning...")
multi_problem3.py_solve()

results3 = multi_problem3.results
print(f"Custom binned results: {len(results3.mueller)} angle pairs")
