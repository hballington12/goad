from goad.convergence.asymmetry import Asymmetry
from goad.convergence.convergable import Tolerance
from goad.convergence.ext_cross import ExtCross

if __name__ == "__main__":
    from goad import Settings
    from goad.convergence.base import Convergence

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
