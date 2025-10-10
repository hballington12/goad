# Convergence Example Script
#
# This script demonstrates how to run a convergence analysis using GOAD.
# The convergence analysis tracks the variance of integrated scattering
# parameters and uses this to estimate the standard error of the mean (SEM).
# This provides a lower-bound estimate of the error in scattering parameters
# compute with GOAD. Typically, problems with large size parameters require
# more orientations and therefore take longer to converge. This code assists
# with ensuring your simulation runs for a sufficient number of orientations.
# GOAD is only an approximate scattering model, and therefore the true error
# is generally larger and can only be inferred by comparing with other methods.

import goad_py as goad
from goad_py import Convergence, Convergable

# Minimal setup
convergence = Convergence(
    # Default settings with path to geometry file
    settings = goad.Settings(geom_path = "../examples/data/hex.obj"),
    # Define convergence criteria for integrated scattering parameters
    convergables = [
        # Run until asymmetry has an absolute SEM of 0.01
        Convergable('asymmetry', 'absolute', 0.005),
        # Run until scattering cross section has a relative SEM of 0.005
        Convergable('scatt', 'relative', 0.01),
    ],
    # Set the number of samples per batch (Choose a number greater than number
    # of threads to allow parallelisation)
    batch_size = 100
)

# Run the convergence
results = convergence.run()
