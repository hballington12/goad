import numpy as np

from goad import Results
from goad.convergence.convergable import Convergable, ConvergenceTracker, Tolerance


class ExtCross(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__("Ext. Cross Section", tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results) -> None:
        if result.ext_cross is None or np.isnan(result.ext_cross):
            return

        self.i += 1
        self.tracker.update(value=result.ext_cross)

    @property
    def mean(self) -> float | None:
        return self.tracker.mean

    @property
    def sem(self) -> float | None:
        return self.tracker.sem
