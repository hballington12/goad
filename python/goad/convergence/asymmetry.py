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
        # print(f"self.sem: {self.sem}")

    @property
    def mean(self) -> float | None:
        return self.tracker.mean

    @property
    def sem(self) -> float | None:
        return self.tracker.sem
