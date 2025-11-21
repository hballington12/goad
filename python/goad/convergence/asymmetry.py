import numpy as np

from goad import Results
from goad.convergence.convergable import Convergable, ConvergenceTracker, Tolerance


class Asymmetry(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__("Asymmetry", tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results) -> None:
        if result.asymmetry is None or np.isnan(result.asymmetry):
            return

        if result.scat_cross is None or np.isnan(result.scat_cross):
            return

        self.tracker.update(value=result.asymmetry, weight=result.scat_cross)

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
        [Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.02)],
    )
    convergence.run()
