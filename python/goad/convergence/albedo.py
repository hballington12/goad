import numpy as np

from goad import Results
from goad.convergence.convergable import Convergable, ConvergenceTracker, Tolerance


class Albedo(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        super().__init__("Single Scatt. Albedo", tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results) -> None:
        if result.albedo is None or np.isnan(result.albedo):
            # print("INFO: Tried to update but albedo is None or NaN")
            return

        # batch size is always 1
        self.i += 1
        self.tracker.update(value=result.albedo)

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
        [Albedo(tolerance=Tolerance.RELATIVE, threshold=0.005)],
    )
    convergence.run()
