# convergence.py
from abc import ABC, abstractmethod
from enum import Enum

from goad import Results


class Tolerance(Enum):
    RELATIVE = "relative"
    ABSOLUTE = "absolute"


class Convergable(ABC):
    def __init__(self, name: str, tolerance: Tolerance, threshold: float) -> None:
        if not name:
            raise ValueError("Name is required for Convergable")
        if threshold <= 0:
            raise ValueError("Threshold must be positive")
        self.name = name
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
