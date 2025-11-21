import numpy as np

from goad import Results
from goad.convergence.convergable import Convergable, ConvergenceTracker, Tolerance


class ScattCross(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__("Scatt. Cross Section", tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results) -> None:
        if result.scat_cross is None or np.isnan(result.scat_cross):
            return

        self.tracker.update(value=result.scat_cross)

    @property
    def mean(self) -> float | None:
        return self.tracker.mean

    @property
    def sem(self) -> float | None:
        return self.tracker.sem
