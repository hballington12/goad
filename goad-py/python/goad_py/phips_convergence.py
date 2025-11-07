"""
PHIPS-specific convergence extension for GOAD.

This module provides convergence tracking for PHIPS detector DSCS values,
which requires Custom binning with PHIPS detector geometry and post-processing
to compute mean DSCS at each of the 20 PHIPS detectors.
"""

import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rich.console import Console

from . import _goad_py as goad
from .convergence import ConvergenceResults
from .convergence_display import (
    ArrayConvergenceVariable,
    ConvergenceDisplay,
)


@dataclass
class PHIPSConvergable:
    """Convergence criteria for PHIPS detector DSCS values."""

    tolerance_type: str = "relative"  # 'relative' or 'absolute'
    tolerance: float = 0.25  # Default 25% relative tolerance
    detector_indices: Optional[List[int]] = (
        None  # Specific detectors to check (None = all)
    )

    def __post_init__(self):
        valid_types = {"relative", "absolute"}
        if self.tolerance_type not in valid_types:
            raise ValueError(
                f"Invalid tolerance_type '{self.tolerance_type}'. Must be one of {valid_types}"
            )

        if self.tolerance <= 0:
            raise ValueError(f"Tolerance must be positive, got {self.tolerance}")

        if self.detector_indices is not None:
            if not isinstance(self.detector_indices, list):
                raise ValueError("detector_indices must be a list of integers")
            if not all(0 <= idx < 20 for idx in self.detector_indices):
                raise ValueError("detector_indices must be in range [0, 19]")


class PHIPSConvergence:
    """
    Convergence study for PHIPS detector DSCS values.

    Requires Custom binning with PHIPS detector geometry (phips_bins_edges.toml).
    Computes mean DSCS at each of 20 PHIPS detectors and tracks convergence.
    """

    # PHIPS detector parameters (from phips_detector_angles.py)
    NUM_DETECTORS = 20
    THETA_START = 18.0  # degrees
    THETA_END = 170.0  # degrees
    DETECTOR_WIDTH = 7.0  # degrees (aperture)

    def __init__(
        self,
        settings: goad.Settings,
        convergable: PHIPSConvergable,
        batch_size: int = 24,
        max_orientations: int = 100_000,
        min_batches: int = 10,
        log_file: Optional[str] = None,
    ):
        """
        Initialize a PHIPS convergence study.

        Args:
            settings: GOAD settings with Custom binning scheme
            convergable: PHIPS convergence criteria
            batch_size: Number of orientations per iteration
            max_orientations: Maximum total orientations before stopping
            min_batches: Minimum number of batches before allowing convergence
            log_file: Optional path to log file for convergence progress
        """
        self.settings = settings
        # Enable quiet mode to suppress Rust progress bars
        self.settings.quiet = True
        self.convergable = convergable
        self.batch_size = batch_size
        self.max_orientations = max_orientations
        self.min_batches = min_batches

        # Validate inputs
        if batch_size <= 0:
            raise ValueError(f"batch_size must be positive, got {batch_size}")

        if max_orientations <= 0:
            raise ValueError(
                f"max_orientations must be positive, got {max_orientations}"
            )

        if min_batches <= 0:
            raise ValueError(f"min_batches must be positive, got {min_batches}")

        # Initialize tracking variables
        self.n_orientations = 0
        self.convergence_history = []

        # Batch-based statistics tracking
        self.batch_data = []  # List of batch statistics

        # PHIPS detector centers (20 detectors from 18° to 170°)
        self.detector_centers = np.linspace(
            self.THETA_START, self.THETA_END, self.NUM_DETECTORS
        )
        self.half_width = self.DETECTOR_WIDTH / 2.0

        # Accumulated PHIPS DSCS for final average
        self.phips_dscs_sum = None

        # Rich console
        self._console = Console()

        # Create display variable for PHIPS DSCS
        display_variable = ArrayConvergenceVariable(
            name="phips_dscs",
            tolerance=convergable.tolerance,
            tolerance_type=convergable.tolerance_type,
            indices=convergable.detector_indices,
        )

        # Initialize display system
        self._display = ConvergenceDisplay(
            variables=[display_variable],
            batch_size=self.batch_size,
            min_batches=self.min_batches,
            convergence_type=self._get_convergence_type(),
            console=self._console,
            log_file=log_file,
        )

    def _compute_phips_dscs_from_mueller2d(self, results: goad.Results) -> np.ndarray:
        """
        Compute mean DSCS at each of 20 PHIPS detectors from Custom binning results.

        Args:
            results: Results from MultiProblem with Custom binning

        Returns:
            Array of shape (20,) with mean DSCS per detector (NaN if no bins in detector)
        """
        # Get mueller_2d from Custom binning
        mueller_2d = np.array(results.mueller)  # Shape: (n_custom_bins, 16)
        bins_2d = results.bins  # List of (theta_center, phi_center) tuples

        # Extract theta angles from bin centers
        theta_angles = np.array([bin_tuple[0] for bin_tuple in bins_2d])

        # Extract S11 and convert to DSCS
        s11_values = mueller_2d[:, 0]
        k = 2 * np.pi / self.settings.wavelength
        dscs_conversion_factor = 1e-12 / k**2
        dscs_values = s11_values * dscs_conversion_factor

        # Compute mean DSCS for each detector
        detector_dscs = []
        for bin_center_theta in self.detector_centers:
            lower_bound = bin_center_theta - self.half_width
            upper_bound = bin_center_theta + self.half_width

            # Find custom bins within this detector's angular window
            indices = np.where(
                (theta_angles >= lower_bound) & (theta_angles < upper_bound)
            )[0]

            if len(indices) > 0:
                # Mean DSCS over bins in this detector window
                mean_dscs = np.mean(dscs_values[indices])
                detector_dscs.append(mean_dscs)
            else:
                # No bins in this detector window
                detector_dscs.append(np.nan)

        return np.array(detector_dscs)  # Shape: (20,)

    def _update_statistics(self, results: goad.Results, batch_size: int):
        """
        Update statistics with new batch results.

        Args:
            results: Results from a MultiProblem run
            batch_size: Number of orientations in this batch
        """
        # Compute PHIPS DSCS for this batch
        phips_dscs = self._compute_phips_dscs_from_mueller2d(results)

        # Store batch data
        batch_info = {
            "batch_size": batch_size,
            "phips_dscs": phips_dscs,  # Shape: (20,)
        }
        self.batch_data.append(batch_info)

        # Accumulate for final average
        if self.phips_dscs_sum is None:
            self.phips_dscs_sum = phips_dscs * batch_size
        else:
            self.phips_dscs_sum += phips_dscs * batch_size

        # Update total orientation count
        self.n_orientations += batch_size

    def _calculate_phips_mean_and_sem(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate mean and SEM arrays for PHIPS DSCS across detectors.

        Returns:
            Tuple of (mean_array, sem_array) where each is shape (20,)
        """
        if not self.batch_data:
            return np.full(self.NUM_DETECTORS, np.nan), np.full(
                self.NUM_DETECTORS, np.inf
            )

        # Extract batch values: shape (n_batches, 20)
        batch_arrays = np.array([batch["phips_dscs"] for batch in self.batch_data])
        batch_sizes = np.array([batch["batch_size"] for batch in self.batch_data])

        if len(batch_arrays) < 2:
            # Can't estimate variance with < 2 batches
            mean_array = batch_arrays[0]
            sem_array = np.full(self.NUM_DETECTORS, np.inf)
            return mean_array, sem_array

        # Calculate mean and SEM independently for each detector
        # Use nanmean to handle NaN values (detectors with no data)
        mean_array = np.average(
            batch_arrays, axis=0, weights=batch_sizes
        )  # Shape: (20,)

        # Variance between batches at each detector (ignoring NaNs)
        batch_means_variance = np.nanvar(batch_arrays, axis=0, ddof=1)  # Shape: (20,)

        # Scale up to estimate population variance
        avg_batch_size = np.mean(batch_sizes)
        estimated_population_variance = batch_means_variance * avg_batch_size

        # Calculate SEM for total sample
        total_n = np.sum(batch_sizes)
        sem_array = np.sqrt(
            estimated_population_variance / (total_n - 1)
        )  # Shape: (20,)

        return mean_array, sem_array

    def _check_convergence(self) -> bool:
        """
        Check if PHIPS DSCS values have converged.

        Returns:
            True if converged, False otherwise
        """
        if len(self.batch_data) < self.min_batches:
            return False

        mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()

        # Determine which detectors to check
        if self.convergable.detector_indices is not None:
            check_indices = self.convergable.detector_indices
        else:
            # Check all detectors that have data (not NaN)
            check_indices = np.where(~np.isnan(mean_dscs))[0]

        if len(check_indices) == 0:
            return False  # No valid detectors to check

        # Extract values for detectors to check
        mean_subset = mean_dscs[check_indices]
        sem_subset = sem_dscs[check_indices]

        # Check convergence based on tolerance type
        if self.convergable.tolerance_type == "relative":
            # Relative SEM
            with np.errstate(divide="ignore", invalid="ignore"):
                relative_sem = np.where(
                    mean_subset != 0, sem_subset / np.abs(mean_subset), np.inf
                )
            converged = np.all(relative_sem < self.convergable.tolerance)
        else:  # absolute
            converged = np.all(sem_subset < self.convergable.tolerance)

        return converged

    def _get_convergence_type(self) -> str:
        """Get the convergence type name for display."""
        class_name = self.__class__.__name__
        if class_name == "PHIPSEnsembleConvergence":
            return "PHIPS Ensemble"
        elif class_name == "PHIPSConvergence":
            return "PHIPS"
        else:
            return class_name

    def _get_detector_angles(self, variable: str) -> np.ndarray:
        """Get detector angles for PHIPS detectors."""
        return self.detector_centers

    def _get_phips_stats(self, variable: str) -> Tuple[float, float]:
        """Get mean and SEM for a single PHIPS detector (not used for array display)."""
        # This is not used since PHIPS uses array display, but required by interface
        return 0.0, 0.0

    def _update_convergence_history(self):
        """Update convergence history with current worst-case SEM."""
        mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()

        if len(mean_dscs) > 0:
            # Find worst-case detector
            if self.convergable.tolerance_type == "relative":
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative_sem = np.where(
                        mean_dscs != 0, sem_dscs / np.abs(mean_dscs), np.inf
                    )
                worst_sem = np.max(relative_sem)
            else:
                worst_sem = np.max(sem_dscs)

            self.convergence_history.append(
                (self.n_orientations, "phips_dscs", worst_sem)
            )

    def run(self) -> ConvergenceResults:
        """
        Run convergence study until criteria are met or max orientations reached.

        Returns:
            ConvergenceResults with PHIPS DSCS values
        """
        iteration = 0
        converged = False

        # Create Live context for smooth updating display
        with self._display.create_live_context() as live:
            # Show initial display before first batch
            initial_display = self._display.build_display(
                iteration=0,
                n_orientations=self.n_orientations,
                get_stats=self._get_phips_stats,
                get_array_stats=lambda var: self._calculate_phips_mean_and_sem(),
                get_bin_labels=self._get_detector_angles,
                power_ratio=None,
                geom_info=None,
            )
            live.update(initial_display)

            while not converged and self.n_orientations < self.max_orientations:
                iteration += 1

                # Create orientations for this batch
                orientations = goad.create_uniform_orientation(self.batch_size)
                self.settings.orientation = orientations

                # Run MultiProblem with error handling for bad geometries
                try:
                    mp = goad.MultiProblem(self.settings)
                    mp.py_solve()
                except Exception as e:
                    # Geometry loading failed (bad faces, degenerate geometry, etc.)
                    # For single-geometry convergence, we can't skip - must raise error
                    error_msg = (
                        f"Failed to initialize MultiProblem with geometry '{self.settings.geom_path}': {e}\n"
                        f"Please check geometry file for:\n"
                        f"  - Degenerate faces (area = 0)\n"
                        f"  - Non-planar geometry\n"
                        f"  - Faces that are too small\n"
                        f"  - Invalid mesh topology\n"
                        f"  - Geometry file corruption"
                    )
                    raise type(e)(error_msg) from e

                # Update statistics
                self._update_statistics(mp.results, self.batch_size)

                # Update convergence history
                self._update_convergence_history()

                # Check convergence
                converged = self._check_convergence()

                # Update live display
                display = self._display.build_display(
                    iteration=iteration,
                    n_orientations=self.n_orientations,
                    get_stats=self._get_phips_stats,
                    get_array_stats=lambda var: self._calculate_phips_mean_and_sem(),
                    get_bin_labels=self._get_detector_angles,
                    power_ratio=None,
                    geom_info=None,
                )
                live.update(display)

        # Compute final results
        mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()

        # Create results
        results = ConvergenceResults(
            converged=converged,
            n_orientations=self.n_orientations,
            values={"phips_dscs": mean_dscs},  # Array of 20 values
            sem_values={"phips_dscs": sem_dscs},  # Array of 20 SEMs
            mueller_1d=None,
            mueller_2d=None,
            convergence_history=self.convergence_history,
            warning=None
            if converged
            else f"Did not converge within {self.max_orientations} orientations",
        )

        # Print final summary
        if converged:
            print(f"\nConverged after {self.n_orientations} orientations.")
        else:
            print(
                f"\nWarning: Did not converge within {self.max_orientations} orientations"
            )

        return results


class PHIPSEnsembleConvergence(PHIPSConvergence):
    """
    Ensemble convergence study for PHIPS detector DSCS values.

    Combines PHIPS detector DSCS tracking with ensemble geometry averaging.
    Each batch randomly selects a geometry file and runs orientation averaging,
    allowing convergence of DSCS values averaged over both orientations and geometries.
    """

    def __init__(
        self,
        settings: goad.Settings,
        convergable: PHIPSConvergable,
        geom_dir: str,
        batch_size: int = 24,
        max_orientations: int = 100_000,
        min_batches: int = 10,
        log_file: Optional[str] = None,
    ):
        """
        Initialize a PHIPS ensemble convergence study.

        Args:
            settings: GOAD settings with Custom binning (geom_path will be overridden)
            convergable: PHIPS convergence criteria
            geom_dir: Directory containing .obj geometry files
            batch_size: Number of orientations per iteration
            max_orientations: Maximum total orientations before stopping
            min_batches: Minimum number of batches before allowing convergence
            log_file: Optional path to log file for convergence progress
        """
        # Discover all .obj files in directory
        geom_path = Path(geom_dir)
        if not geom_path.exists():
            raise ValueError(f"Geometry directory does not exist: {geom_dir}")

        if not geom_path.is_dir():
            raise ValueError(f"Path is not a directory: {geom_dir}")

        self.geom_files = sorted([f.name for f in geom_path.glob("*.obj")])

        if not self.geom_files:
            raise ValueError(f"No .obj files found in directory: {geom_dir}")

        self.geom_dir = str(geom_path.resolve())

        print(f"Found {len(self.geom_files)} geometry files in {self.geom_dir}")

        # Call parent constructor
        super().__init__(
            settings=settings,
            convergable=convergable,
            batch_size=batch_size,
            max_orientations=max_orientations,
            min_batches=min_batches,
            log_file=log_file,
        )

    def run(self) -> ConvergenceResults:
        """
        Run ensemble convergence study.

        Each batch iteration randomly selects a geometry file from the
        ensemble directory before running the orientation averaging.

        Returns:
            ConvergenceResults with ensemble-averaged PHIPS DSCS values
        """
        iteration = 0
        converged = False
        skipped_geometries = []  # Track skipped geometry files

        # Create Live context for smooth updating display
        with self._display.create_live_context() as live:
            # Show initial display before first batch
            initial_display = self._display.build_display(
                iteration=0,
                n_orientations=self.n_orientations,
                get_stats=self._get_phips_stats,
                get_array_stats=lambda var: self._calculate_phips_mean_and_sem(),
                get_bin_labels=self._get_detector_angles,
                power_ratio=None,
                geom_info=None,
            )
            live.update(initial_display)

            while not converged and self.n_orientations < self.max_orientations:
                iteration += 1

                # Randomly select a geometry file for this batch
                geom_file = random.choice(self.geom_files)
                geom_path = os.path.join(self.geom_dir, geom_file)

                # Create orientations for this batch
                orientations = goad.create_uniform_orientation(self.batch_size)

                # Update settings with selected geometry and orientations
                self.settings.geom_path = geom_path
                self.settings.orientation = orientations

                # Run MultiProblem
                try:
                    mp = goad.MultiProblem(self.settings)
                    mp.py_solve()
                except Exception as e:
                    # Geometry loading failed (bad faces, degenerate geometry, etc.)
                    print(f"\nWarning: Skipping geometry '{geom_file}': {e}")
                    skipped_geometries.append(geom_file)

                    # Check if all geometries have been skipped
                    if len(skipped_geometries) >= len(self.geom_files):
                        raise ValueError(
                            f"All {len(self.geom_files)} geometry files failed to load. "
                            "Please check geometry files for degenerate faces, non-planar geometry, "
                            "or faces that are too small."
                        )

                    # Skip this iteration without updating statistics
                    continue

                # Update statistics
                self._update_statistics(mp.results, self.batch_size)

                # Update convergence history
                self._update_convergence_history()

                # Check convergence
                converged = self._check_convergence()

                # Update live display with geometry info
                geom_info = f"Geom: {geom_file}"
                display = self._display.build_display(
                    iteration=iteration,
                    n_orientations=self.n_orientations,
                    get_stats=self._get_phips_stats,
                    get_array_stats=lambda var: self._calculate_phips_mean_and_sem(),
                    get_bin_labels=self._get_detector_angles,
                    power_ratio=None,
                    geom_info=geom_info,
                )
                live.update(display)

        # Compute final results
        mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()

        # Prepare warning message
        warning = None
        if not converged:
            warning = f"Did not converge within {self.max_orientations} orientations"

        # Add skipped geometries info to warning
        if skipped_geometries:
            skipped_msg = f"Skipped {len(skipped_geometries)} bad geometries"
            warning = f"{warning} | {skipped_msg}" if warning else skipped_msg

        # Create results
        results = ConvergenceResults(
            converged=converged,
            n_orientations=self.n_orientations,
            values={"phips_dscs": mean_dscs},
            sem_values={"phips_dscs": sem_dscs},
            mueller_1d=None,
            mueller_2d=None,
            convergence_history=self.convergence_history,
            warning=warning,
        )

        # Print final summary
        if converged:
            print(f"\nEnsemble converged after {self.n_orientations} orientations.")
        else:
            print(
                f"\nWarning: Did not converge within {self.max_orientations} orientations"
            )

        # Report skipped geometries if any
        if skipped_geometries:
            print(
                f"Note: Skipped {len(skipped_geometries)} geometry file(s) due to errors"
            )

        return results
