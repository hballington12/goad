"""
Example script demonstrating MultiProblem Python API with orientation schemes.

This script shows how to:
1. Create uniform and discrete orientation schemes
2. Use MultiProblem for multi-orientation averaging
3. Access averaged results through the Results object
"""

import goad_py as goad

print("MultiProblem Example - Orientation Averaging")
print("=" * 50)

# Example 1: Uniform orientation scheme
print("\n1. Creating MultiProblem with uniform orientations...")

# Create uniform orientation (100 random orientations)
uniform_orientation = goad.create_uniform_orientation(100)  # euler_convention is now optional
print(f"Created uniform orientation: {uniform_orientation}")

# Create settings with uniform orientation
settings = goad.Settings(
    geom_name="hex.obj",
    orientation=uniform_orientation
)

# Create and solve MultiProblem
print("Creating MultiProblem with uniform orientations...")
multi_problem = goad.MultiProblem(settings=settings)
print(f"Number of orientations: {multi_problem.num_orientations}")

print("Solving MultiProblem (this averages over all orientations)...")
multi_problem.py_solve()

# Access averaged results
results = multi_problem.results
print(f"Averaged Mueller matrix shape: {len(results.mueller)}x{len(results.mueller[0])}")
print(f"Averaged asymmetry parameter: {results.asymmetry}")

# Example 2: Discrete orientation scheme  
print("\n2. Creating MultiProblem with discrete orientations...")

# Create specific Euler angles
euler1 = goad.Euler(0.0, 0.0, 0.0)
euler2 = goad.Euler(30.0, 30.0, 30.0) 
euler3 = goad.Euler(45.0, 60.0, 90.0)

# Create discrete orientation scheme
discrete_orientation = goad.create_discrete_orientation([euler1, euler2, euler3])  # euler_convention is now optional
print(f"Created discrete orientation: {discrete_orientation}")

# Create settings with discrete orientation
settings2 = goad.Settings(
    geom_name="hex.obj", 
    orientation=discrete_orientation
)

# Create and solve MultiProblem
print("Creating MultiProblem with discrete orientations...")
multi_problem2 = goad.MultiProblem(settings=settings2)
print(f"Number of orientations: {multi_problem2.num_orientations}")

print("Solving MultiProblem with specific orientations...")
multi_problem2.py_solve()

# Access averaged results
results2 = multi_problem2.results
print(f"Averaged Mueller matrix shape: {len(results2.mueller)}x{len(results2.mueller[0])}")
print(f"Averaged asymmetry parameter: {results2.asymmetry}")

# Example 3: Using helper functions directly
print("\n3. Using helper functions for orientation schemes...")

# Create schemes directly (without full Orientation wrapper)
uniform_scheme = goad.uniform_orientation(50)
print(f"Uniform scheme: {uniform_scheme}")

discrete_scheme = goad.discrete_orientation([euler1, euler2])
print(f"Discrete scheme: {discrete_scheme}")

print("\n✓ MultiProblem examples completed successfully!")
print("\nKey benefits of MultiProblem:")
print("- Automatically averages results over multiple orientations")
print("- Supports both random uniform and specific discrete orientations") 
print("- Returns averaged Mueller matrices and physical parameters")
print("- Same Results API as single Problem for consistency")