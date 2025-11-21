import numpy as np

from goad import Results
from goad.convergence.convergable import Convergable, ConvergenceTracker, Tolerance


class Albedo(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__("Single Scatt. Albedo", tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results) -> None:
        if result.albedo is None or np.isnan(result.albedo):
            return

        if result.ext_cross is None or np.isnan(result.ext_cross):
            return

        self.tracker.update(value=result.albedo, weight=result.ext_cross)

    @property
    def mean(self) -> float | None:
        return self.tracker.mean

    @property
    def sem(self) -> float | None:
        return self.tracker.sem
