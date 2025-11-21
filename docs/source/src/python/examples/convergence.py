# --8<-- [start:basic]
from goad import Settings
from goad.convergence.asymmetry import Asymmetry
from goad.convergence.base import Convergence
from goad.convergence.convergable import Tolerance

convergence = Convergence(
    Settings(geom_path="../../../examples/data/hex.obj", quiet=True),
    [Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.02)],
)
convergence.run()
# --8<-- [end:basic]

# --8<-- [start:multiple]
from goad import Settings
from goad.convergence.asymmetry import Asymmetry
from goad.convergence.base import Convergence
from goad.convergence.convergable import Tolerance
from goad.convergence.ext_cross import ExtCross

convergence = Convergence(
    Settings(
        geom_path="path/to/geometry.obj",
        quiet=True,
        particle_refr_index_im=0.001,
    ),
    [
        ExtCross(tolerance=Tolerance.RELATIVE, threshold=0.02),
        Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.02),
    ],
)
convergence.run()
# --8<-- [end:multiple]

# --8<-- [start:albedo]
from goad import Settings
from goad.convergence.albedo import Albedo
from goad.convergence.base import Convergence
from goad.convergence.convergable import Tolerance

convergence = Convergence(
    Settings(
        geom_path="path/to/geometry.obj",
        quiet=True,
        particle_refr_index_im=0.001,
    ),
    [Albedo(tolerance=Tolerance.RELATIVE, threshold=0.005)],
)
convergence.run()
# --8<-- [end:albedo]


# --8<-- [start:extcross]
from goad import Settings
from goad.convergence.base import Convergence
from goad.convergence.convergable import Tolerance
from goad.convergence.ext_cross import ExtCross

convergence = Convergence(
    Settings(
        geom_path="path/to/geometry.obj",
        quiet=True,
        particle_refr_index_im=0.001,
    ),
    [ExtCross(tolerance=Tolerance.RELATIVE, threshold=0.02)],
)
convergence.run()
# --8<-- [end:extcross]

# --8<-- [start:results]
import numpy as np
from goad import Settings
from goad.convergence.asymmetry import Asymmetry
from goad.convergence.base import Convergence
from goad.convergence.convergable import Tolerance

np.set_printoptions(threshold=4)
np.set_printoptions(precision=3)

convergence = Convergence(
    Settings(geom_path="./data/hex.obj", quiet=True, particle_refr_index_im=0.001),
    [Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.01)],
)
convergence.run()
results = convergence.results()

# Get the main integrated parameters
print(f"Asymmetry: {results.asymmetry:.4f}")
print(f"Scattering Cross Section: {results.scat_cross:.4f}")
print(f"Extinction Cross Section: {results.ext_cross:.4f}")
print(f"Absorption Cross Section: {results.ext_cross - results.scat_cross:.4f}")
print(f"Single Scattering Albdeo: {results.albedo:.4f}")

# Access the Mueller matrix
print(f"Theta bins: \n{results.bins_1d[:]}")
print(f"[Theta, Phi] bins: \n{results.bins[:]}")
print(f"Mueller matrix S11: {results.mueller_1d[:, 0]}")
print(f"Mueller matrix S12: {results.mueller_1d[:, 1]}")
# --8<-- [end:results]
