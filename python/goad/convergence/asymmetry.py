import numpy as np

from goad import Results
from goad.convergence.base import Convergable, ConvergenceTracker, Tolerance


class Asymmetry(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__(tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results, batch_size: int) -> None:
        if result.asymmetry is None or np.isnan(result.asymmetry):
            print("INFO: Tried to update but asymmetry parameter is None or NaN")
            return

        if result.scat_cross is None or np.isnan(result.scat_cross):
            print(
                "INFO: Tried to update but cannot weight asymmetry parameter because scat_cross is None or NaN"
            )
            return

        for _ in range(batch_size):  # dirty fix to sidestep varying batch sizes
            self.i += 1
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
        [Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.01)],
    )
    convergence.run()
