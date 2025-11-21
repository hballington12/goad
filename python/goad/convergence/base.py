# convergence.py
import time

import numpy as np
from rich.live import Live

from goad import Geom, MultiProblem, Orientation, Results, Settings
from goad.convergence.convergable import Convergable
from goad.convergence.display import ConvergenceDisplay


# A GOAD Convergence
class Convergence:
    def __init__(
        self,
        goad_settings: Settings,
        targets: list[Convergable],
        max_orientations: int = 100000,
        min_orientations: int = 50,
        refresh_rate: float = 30,
    ) -> None:
        if len(targets) == 0:
            raise ValueError("targets must be a non-empty list")

        self.geometry = Geom.from_file(goad_settings.geom_path)
        self.targets = targets
        # possibly add more config here later
        self.max_orientations = max_orientations
        self.sim_time = np.inf  # time per simulation
        self.refresh_rate = refresh_rate  # refresh rate in Hz
        self.min_orientations = (
            min_orientations  # stop early termination in lucky cases
        )
        self.i: int = 0
        self.result_m: None | Results = None
        self.result_mm: None | Results = None
        self.result_d: None | Results = None
        self.result_w: None | Results = None
        self.result_ww: None | Results = None
        self.result_dw: None | Results = None
        self.result_s: None | Results = None
        self.result_ss: None | Results = None

        goad_settings.orientation = Orientation.uniform(1)

        self.goad_settings = goad_settings

        self.display = ConvergenceDisplay(self.targets)

    def results(self) -> None | Results:
        return self.result_m

    def results_sem(self) -> None | Results:
        if self.result_m is None:
            return None
        return self.result_m**0.5 / (self.i - 1)

    def is_converged(self) -> bool:
        return all(target.is_converged() for target in self.targets)

    def run(self) -> None:
        self.i = 0

        # Build initial display
        initial_display = self.display.build(
            self.i, self.sim_time, self.min_orientations, self.max_orientations
        )

        # Start live display
        with Live(
            initial_display,
            console=self.display.console,
            refresh_per_second=self.refresh_rate,
            transient=False,
        ) as live:
            while True:
                if self.is_converged() and self.i > self.min_orientations:
                    print(f"Converged in {self.i} iterations")
                    break

                # Run iteration
                self.update(self.iterate())

                # Update display
                live.update(
                    self.display.build(
                        self.i,
                        self.sim_time,
                        self.min_orientations,
                        self.max_orientations,
                    )
                )

                if self.i > self.max_orientations:
                    print(f"Did not converge after {self.i} iterations")
                    break

            if self.result_m is not None:
                print(f"asymmetry is: {self.result_m.asymmetry}")

    def update(self, result: Results) -> None:
        for target in self.targets:
            target.update(result)
        self.inc_results(result)

    def inc_results(self, result: Results) -> None:
        self.i += 1
        if self.i == 1:
            self.result_m = result  # initliase mean
        elif self.i > 1:
            self.result_mm = self.result_m
            self.result_ss = self.result_s
            self.result_d = result - self.result_mm  # pyright: ignore[reportOperatorIssue]
            self.result_m = self.result_mm + (self.result_d / self.i)  # pyright: ignore[reportOptionalOperand]
            if self.i == 2:
                self.result_s = self.result_d**2 * ((self.i - 1) / self.i)
            else:
                self.result_s = self.result_ss + self.result_d**2 * (  # pyright: ignore[reportOptionalOperand]
                    (self.i - 1) / self.i
                )

    @staticmethod
    def inc_val(old: float, new: float, iterations: int) -> float:
        delta = new - old
        return old + delta / iterations

    def iterate(self) -> Results:
        start = time.time()
        mp = MultiProblem(
            self.goad_settings, geom=self.geometry
        )  # using cached geometry
        mp.solve()
        self.sim_time = time.time() - start
        return mp.results
