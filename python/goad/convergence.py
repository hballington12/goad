# convergence.py
from abc import ABC, abstractmethod
from enum import Enum

import numpy as np

from goad import MultiProblem, Results, Settings


class Tolerance(Enum):
    RELATIVE = "relative"
    ABSOLUTE = "absolute"


class Convergable(ABC):
    def __init__(self, tolerance: Tolerance, threshold: float):
        self.tolerance = tolerance
        self.threshold = threshold
        self.i = 0

    def is_converged(self) -> bool:
        mean = self.mean()
        sem = self.sem()
        if mean is None:
            return False
        if sem is None:
            return False
        print(
            f"i: {self.i}, mean: {mean:.4f}, sem: {sem:.4f}, threshold: {self.threshold:.4f}"
        )
        if self.tolerance == Tolerance.RELATIVE:
            return sem < self.threshold
        elif self.tolerance == Tolerance.ABSOLUTE:
            return sem < self.threshold
        else:
            raise ValueError(f"Unknown tolerance type: {self.tolerance}")

    @abstractmethod
    def update(self, result: Results):
        pass

    @abstractmethod
    def mean(self) -> float | None:
        pass

    @abstractmethod
    def sem(self) -> float | None:
        pass


class ConvergenceTracker:
    # https://en.wikipedia.org/wiki/Monte_Carlo_method#Determining_a_sufficiently_large_n
    # modified with weights
    def __init__(self):
        self.i = 0
        self.m = None  # weighted mean
        self.mm = None  # previous weighted mean
        self.d = None  # delta
        self.s = 0  # variance
        self.ss = 0  # variance of previous iteration
        self.w = 0  # mean weight
        self.ww = None  # previous mean weight
        self.dw = None  # delta weight

    def update(self, value: float, weight: float = 1.0):
        self.i += 1
        value *= weight
        if self.i == 1:
            self.m = value
            self.w = weight
        elif self.i > 1:
            self.mm = self.m
            self.ss = self.s
            self.d = value - self.mm
            self.m = self.mm + (self.d / self.i)
            self.s = self.ss + self.d**2 * ((self.i - 1) / self.i)
            self.ww = self.w
            self.dw = weight - self.ww
            self.w = self.ww + (self.dw / self.i)

    def sem(self):
        if self.i < 2:
            return None
        return (self.s / (self.i - 1) ** 2 / (self.w**2)) ** 0.5

    def mean(self):
        if self.m is None:
            return None
        return self.m / self.w


class Asymmetry(Convergable):
    def __init__(self, tolerance: Tolerance, threshold: float):
        super().__init__(tolerance, threshold)
        self.tracker = ConvergenceTracker()

    def update(self, result: Results):
        self.i += 1
        if result.asymmetry is None:
            print("INFO: Tried to update but asymmetry parameter is None")
            return

        if np.isnan(result.asymmetry):
            print("INFO: Tried to update but asymmetry parameter is NaN")
            return

        if result.scat_cross is None:
            print(
                "INFO: Tried to update but cannot weight asymmetry parameter because scat_cross is None"
            )
            return

        if np.isnan(result.scat_cross):
            print(
                "INFO: Tried to update but cannot weight asymmetry parameter because scat_cross is NaN"
            )
            return

        self.tracker.update(value=result.asymmetry, weight=result.scat_cross)

    def mean(self):
        return self.tracker.mean()

    def sem(self):
        return self.tracker.sem()


# A GOAD Convergence
class Convergence:
    def __init__(self, goad_settings: Settings, targets: list[Convergable]):
        self.goad_settings = goad_settings
        if len(targets) == 0:
            raise ValueError("You must provide at least one Convergable target")
        self.targets = targets
        # possibly add more config here later
        self.max_iterations = 1000

    def is_converged(self):
        return all(target.is_converged() for target in self.targets)

    def run(self):
        iterations = 0
        while not self.is_converged():
            iterations += 1
            result = self.iterate()
            self.update(result)

            if iterations > self.max_iterations:
                print("INFO: reached max iterations without converging.")
                return

    def update(self, result: Results):
        for target in self.targets:
            target.update(result)

    def iterate(self) -> Results:
        mp = MultiProblem(self.goad_settings)
        mp.solve()
        return mp.results


if __name__ == "__main__":
    convergence = Convergence(
        Settings(geom_path="../../examples/data/hex.obj", quiet=True),
        [Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.01)],
    )
    convergence.run()
