# convergence.py
import time
from abc import ABC, abstractmethod
from enum import Enum

import numpy as np

from goad import Geom, MultiProblem, Orientation, Results, Settings


class Tolerance(Enum):
    RELATIVE = "relative"
    ABSOLUTE = "absolute"


class Convergable(ABC):
    def __init__(self, tolerance: Tolerance, threshold: float) -> None:
        if threshold <= 0:
            raise ValueError("Threshold must be positive")
        self.tolerance = tolerance
        self.threshold = threshold
        self.i = 0

    def is_converged(self) -> bool:
        mean = self.mean
        sem = self.sem
        if mean is None:
            return False
        if sem is None:
            return False
        print(
            f"i: {self.i}, mean: {mean:.4f}, sem: {sem:.4f}, threshold: {self.threshold:.4f}"
        )
        if self.tolerance == Tolerance.RELATIVE:
            return sem / mean < self.threshold
        elif self.tolerance == Tolerance.ABSOLUTE:
            return sem < self.threshold
        else:
            raise ValueError(f"Unknown tolerance type: {self.tolerance}")

    @abstractmethod
    def update(self, result: Results, batch_size: int) -> None:
        pass

    @property
    @abstractmethod
    def mean(self) -> float | None:
        pass

    @property
    @abstractmethod
    def sem(self) -> float | None:
        pass


class ConvergenceTracker:
    # https://en.wikipedia.org/wiki/Monte_Carlo_method#Determining_a_sufficiently_large_n
    # modified with weights
    def __init__(self) -> None:
        self.i = 0
        self.m = None  # weighted mean
        self.mm = None  # previous weighted mean
        self.d = None  # delta
        self.s = 0  # variance
        self.ss = 0  # variance of previous iteration
        self.w = 0  # mean weight
        self.ww = None  # previous mean weight
        self.dw = None  # delta weight

    def update(self, value: float, weight: float = 1.0) -> None:
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

    @property
    def sem(self) -> float | None:
        if self.i < 2:
            return None
        return (self.s / (self.i - 1) ** 2 / (self.w**2)) ** 0.5

    @property
    def mean(self) -> float | None:
        if self.m is None:
            return None
        return self.m / self.w


# A GOAD Convergence
class Convergence:
    def __init__(self, goad_settings: Settings, targets: list[Convergable]) -> None:
        if len(targets) == 0:
            raise ValueError("You must provide at least one Convergable target")

        self.geometry = Geom.from_file(goad_settings.geom_path)
        self.targets = targets
        # possibly add more config here later
        self.max_orientations = 100000
        self.batch_size = 1  # start at 1, increase as needed
        self.sim_time = np.inf  # time per simulation
        self.refresh_rate = 30  # refresh rate in Hz
        self.min_orientations = 50  # stop early termination in lucky cases

        goad_settings.orientation = Orientation.uniform(self.batch_size)

        self.goad_settings = goad_settings

    def is_converged(self) -> bool:
        return all(target.is_converged() for target in self.targets)

    def run(self) -> None:
        iterations = 0
        while not self.is_converged() or iterations < self.min_orientations:
            iterations += self.batch_size
            result = self.iterate()
            self.update(result)

            if iterations > self.max_orientations:
                print("INFO: reached max iterations without converging.")
                break

    def update_batch_size(self) -> None:
        """Update the batch size to balance performance with repsonsiveness."""
        if np.isinf(self.sim_time):
            return
        self.batch_size = max(
            1, int(self.batch_size / self.sim_time / self.refresh_rate)
        )
        self.goad_settings.orientation = Orientation.uniform(self.batch_size)

    def update(self, result: Results) -> None:
        for target in self.targets:
            target.update(result, self.batch_size)
        self.update_batch_size()

    def iterate(self) -> Results:
        start = time.time()
        mp = MultiProblem(
            self.goad_settings, geom=self.geometry
        )  # using cached geometry
        mp.solve()
        self.sim_time = time.time() - start
        return mp.results
