import numpy as np

from goad import Results
from goad.convergence.convergable import Convergable, ConvergenceTracker, Tolerance


class ExtCross(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__("Ext. Cross Section", tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results) -> None:
        if result.ext_cross is None or np.isnan(result.ext_cross):
            # print("INFO: Tried to update but scattering cross section is None or NaN")
            return

        self.i += 1
        self.tracker.update(value=result.ext_cross)

    @property
    def mean(self) -> float | None:
        return self.tracker.mean

    @property
    def sem(self) -> float | None:
        return self.tracker.sem


if __name__ == "__main__":
    from goad import Settings
    from goad.convergence.base import Convergence

    convergence = Convergence(
        Settings(
            geom_path="../../../examples/data/hex.obj",
            quiet=True,
            particle_refr_index_im=0.001,
        ),
        [ExtCross(tolerance=Tolerance.RELATIVE, threshold=0.02)],
    )
    convergence.run()
