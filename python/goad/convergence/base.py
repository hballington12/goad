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
        self.batch_size = 1  # start at 1, increase as needed
        self.sim_time = np.inf  # time per simulation
        self.refresh_rate = refresh_rate  # refresh rate in Hz
        self.min_orientations = (
            min_orientations  # stop early termination in lucky cases
        )
        self.iterations: int = 0
        self.results: None | Results = None

        goad_settings.orientation = Orientation.uniform(self.batch_size)

        self.goad_settings = goad_settings

        self.display = ConvergenceDisplay(self.targets)

    def is_converged(self) -> bool:
        return all(target.is_converged() for target in self.targets)

    def run(self) -> None:
        self.iterations = 0

        # Build initial display
        initial_display = self.display.build(
            self.iterations,
            self.batch_size,
            self.sim_time,
            self.min_orientations,
            self.max_orientations,
        )

        # Start live display
        with Live(
            initial_display,
            console=self.display.console,
            refresh_per_second=self.refresh_rate,
            transient=False,
        ) as live:
            while not self.is_converged() or self.iterations < self.min_orientations:
                self.iterations += self.batch_size

                # Run iteration
                self.update(self.iterate())

                # Update display
                live.update(
                    self.display.build(
                        self.iterations,
                        self.batch_size,
                        self.sim_time,
                        self.min_orientations,
                        self.max_orientations,
                    )
                )

                if self.iterations > self.max_orientations:
                    break

    # def update_batch_size(self) -> None:
    #     """Update the batch size to balance performance with repsonsiveness."""
    #     if np.isinf(self.sim_time):
    #         return
    #     optimal_batch_size = max(
    #         1, int(self.batch_size / self.sim_time / self.refresh_rate)
    #     )
    #     #  weighted mean of current mean batch size and optimal based on last batch
    #     self.batch_size = int(
    #         (
    #             self.batch_size * (self.iterations - self.batch_size)
    #             + optimal_batch_size * self.batch_size
    #         )
    #         / self.iterations
    #     )
    #     self.goad_settings.orientation = Orientation.uniform(self.batch_size)

    def update(self, result: Results) -> None:
        for target in self.targets:
            target.update(result, self.batch_size)
        # self.update_batch_size()
        self.inc_results(result)

    def inc_results(self, result: Results) -> None:
        # logic
        print("inc-ing")

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
