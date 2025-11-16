import numpy as np

from goad import Results
from goad.convergence.base import Convergable, ConvergenceTracker, Tolerance


class ScattCross(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__(tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results, batch_size: int) -> None:
        if result.scat_cross is None or np.isnan(result.scat_cross):
            print("INFO: Tried to update but scattering cross section is None or NaN")
            return

        for _ in range(batch_size):  # dirty fix to sidestep varying batch sizes
            self.i += 1
            self.tracker.update(value=result.scat_cross)

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
        Settings(geom_path="../../../examples/data/hex.obj", quiet=True),
        [ScattCross(tolerance=Tolerance.RELATIVE, threshold=0.01)],
    )
    convergence.run()
