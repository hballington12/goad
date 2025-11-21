# convergence.py
import time

import numpy as np
from rich.live import Live

from goad import Geom, MultiProblem, Orientation, Results, Settings
from goad.convergence.convergable import Convergable
from goad.convergence.display import ConvergenceDisplay


class Convergence:
    """
    GOAD convergence analyzer for orientation-averaged scattering.

    Manages iterative Monte Carlo simulations with incremental statistics tracking
    and live convergence monitoring. Uses Welford's online algorithm for numerical
    stability when computing running mean and variance.
    """

    def __init__(
        self,
        goad_settings: Settings,
        targets: list[Convergable],
        max_orientations: int = 100000,
        min_orientations: int = 50,
        refresh_rate: float = 30,
    ) -> None:
        """
        Create a new convergence analyzer.

        Args:
            `goad_settings`: GOAD settings object for the simulation
            `targets`: List of convergence targets to monitor
            `max_orientations`: Maximum number of orientations to simulate before stopping
            `min_orientations`: Minimum orientations before allowing early convergence
            `refresh_rate`: Display refresh rate in Hz

        Raises:
            ValueError: If targets list is empty
        """
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
        """
        Return the current mean results.

        Returns:
            Mean `Results` object, or `None` if no iterations completed
        """
        return self.result_m

    def results_sem(self) -> None | Results:
        """
        Return the standard error of the mean (SEM) for the results.

        Returns:
            SEM as a `Results` object, or `None` if insufficient iterations
        """
        if self.result_s is None:
            return None
        return self.result_s**0.5 / (self.i - 1)

    def _is_converged(self) -> bool:
        """
        Check if all convergence targets have been met.

        Returns:
            `True` if all targets are converged, `False` otherwise
        """
        return all(target.is_converged() for target in self.targets)

    def run(self) -> None:
        """
        Run the convergence loop until targets converge or max iterations reached.

        Executes Monte Carlo simulations iteratively, updating statistics and live
        display until convergence criteria are met or `max_orientations` is reached.
        """
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
                if self._is_converged() and self.i > self.min_orientations:
                    # print(f"Converged in {self.i} iterations")
                    break

                # Run iteration
                self._update(self._iterate())

                # Update display
                live.update(
                    self.display.build(
                        self.i,
                        self.sim_time,
                        self.min_orientations,
                        self.max_orientations,
                    )
                )

                # temporary print statement
                if self.results_sem() is not None:
                    print(f"self.sem (base.py): {self.results_sem().asymmetry}")

                if self.i > self.max_orientations:
                    print(f"Did not converge after {self.i} iterations")
                    break

    def _update(self, result: Results) -> None:
        """
        Update all convergence targets and running statistics with new results.

        Args:
            `result`: New `Results` object from the current iteration
        """
        for target in self.targets:
            target.update(result)
        self._inc_results(result)

    def _inc_results(self, result: Results) -> None:
        """
        Incrementally update mean and variance using Welford's online algorithm.

        Args:
            `result`: New `Results` object to incorporate into statistics
        """
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
                print(f"base.py: {self.result_s.asymmetry}")
                self.result_s = self.result_ss + self.result_d**2 * (  # pyright: ignore[reportOptionalOperand]
                    (self.i - 1) / self.i
                )

    def _iterate(self) -> Results:
        """
        Run a single GOAD simulation iteration.

        Creates and solves a `MultiProblem` with the cached geometry,
        timing the execution.

        Returns:
            `Results` object from the simulation
        """
        start = time.time()
        mp = MultiProblem(
            self.goad_settings, geom=self.geometry
        )  # using cached geometry
        mp.solve()
        self.sim_time = time.time() - start
        return mp.results


if __name__ == "__main__":
    from goad import Settings
    from goad.convergence.asymmetry import Asymmetry
    from goad.convergence.base import Convergence
    from goad.convergence.convergable import Tolerance

    convergence = Convergence(
        Settings(geom_path="../../../examples/data/hex.obj", quiet=True),
        [Asymmetry(tolerance=Tolerance.RELATIVE, threshold=0.04)],
    )
    convergence.run()
    results = convergence.results()
    error = convergence.results_sem()

    print(f"Asymmetry: {results.scat_cross:.4f} +/- {error.scat_cross:.4f}")
    print(f"Asymmetry: {results.asymmetry:.4f} +/- {error.asymmetry:.4f}")
