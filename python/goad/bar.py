import numpy as np

import goad

euler = goad.Euler(0, 45, 0)
euler.alpha = 3
print(euler)

convention = goad.EulerConvention("xyz")
print(convention)

scheme = goad.Orientation.uniform(19)
binning = goad.BinningScheme.simple(num_theta=100, num_phi=100)
binning = goad.BinningScheme.interval(
    thetas=[0, 1, 180], theta_spacings=[0.1, 1], phis=[0, 360], phi_spacings=[2]
)
print(binning.thetas())
print(binning.phis())

orientation = goad.Orientation.discrete([euler, euler])
orientation = goad.Orientation.uniform(1)
settings = goad.Settings(
    geom_path="../../examples/data/hex.obj", orientation=orientation
)
print(settings)
mp = goad.MultiProblem(settings)
mp.solve()
print(mp.results.asymmetry)

mueller = mp.results.mueller
mueller_1d = mp.results.mueller_1d

if mueller_1d is not None:
    print(mueller_1d[:, 0])

# Test setter functionality

# Cache the results object to avoid getting fresh clones
results = mp.results

# Extract first row
first_row = mueller[0]
print(f"\nOriginal first row: {first_row}")

# Null test - set it to itself
results.mueller = mueller
mueller_copy = results.mueller
assert np.allclose(mueller, mueller_copy), "Null test failed!"
print("Null test passed ✓")

# Modify a value (change s11 of first bin)
mueller_modified = mueller.copy()
original_s11 = mueller_modified[0, 0]
mueller_modified[0, 0] = 999.0

# Set the modified array
results.mueller = mueller_modified

# Get it back and check
mueller_retrieved = results.mueller
print(f"Original s11: {original_s11}")
print(f"Modified s11: {mueller_retrieved[0, 0]}")
assert mueller_retrieved[0, 0] == 999.0, "Modification test failed!"
print("Modification test passed ✓")
