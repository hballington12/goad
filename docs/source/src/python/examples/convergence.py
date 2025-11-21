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
        geom_path="../../../examples/data/hex.obj",
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
