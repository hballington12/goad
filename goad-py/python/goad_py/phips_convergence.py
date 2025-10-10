"""
PHIPS-specific convergence extension for GOAD.

This module provides convergence tracking for PHIPS detector DSCS values,
which requires Custom binning with PHIPS detector geometry and post-processing
to compute mean DSCS at each of the 20 PHIPS detectors.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np
import os
import random
from pathlib import Path
from . import _goad_py as goad
from .convergence import ConvergenceResults


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
    ):
        """
        Initialize a PHIPS convergence study.

        Args:
            settings: GOAD settings with Custom binning scheme
            convergable: PHIPS convergence criteria
            batch_size: Number of orientations per iteration
            max_orientations: Maximum total orientations before stopping
            min_batches: Minimum number of batches before allowing convergence
        """
        self.settings = settings
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

    def _print_progress(self, converged: bool):
        """Print convergence progress."""
        mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()

        # Determine which detectors to display
        if self.convergable.detector_indices is not None:
            check_indices = self.convergable.detector_indices
        else:
            check_indices = np.where(~np.isnan(mean_dscs))[0]

        if len(check_indices) == 0:
            print(f"  PHIPS DSCS: No valid detectors")
            return

        # Find worst-case detector (highest relative SEM)
        mean_subset = mean_dscs[check_indices]
        sem_subset = sem_dscs[check_indices]

        with np.errstate(divide="ignore", invalid="ignore"):
            relative_sem = np.where(
                mean_subset != 0, sem_subset / np.abs(mean_subset), np.inf
            )

        worst_idx_in_subset = np.argmax(relative_sem)
        worst_detector_idx = check_indices[worst_idx_in_subset]
        worst_theta = self.detector_centers[worst_detector_idx]
        worst_mean = mean_subset[worst_idx_in_subset]
        worst_sem = sem_subset[worst_idx_in_subset]
        worst_rel_sem = relative_sem[worst_idx_in_subset]

        # Count converged detectors
        if self.convergable.tolerance_type == "relative":
            converged_mask = relative_sem < self.convergable.tolerance
        else:
            converged_mask = sem_subset < self.convergable.tolerance

        n_converged = np.sum(converged_mask)
        n_total = len(check_indices)

        status = "✓" if converged else "..."

        if self.convergable.tolerance_type == "relative":
            current_str = f"{worst_rel_sem * 100:.2f}%"
            target_str = f"{self.convergable.tolerance * 100:.2f}%"
        else:
            current_str = f"{worst_sem:.4g}"
            target_str = f"{self.convergable.tolerance:.4g}"

        print(
            f"  PHIPS DSCS: {n_converged}/{n_total} detectors converged | "
            f"Worst θ={worst_theta:.1f}°: {worst_mean:.4g} | "
            f"SEM: {current_str} (target: {target_str}) {status}"
        )

    def run(self) -> ConvergenceResults:
        """
        Run convergence study until criteria are met or max orientations reached.

        Returns:
            ConvergenceResults with PHIPS DSCS values
        """
        print(f"\nStarting PHIPS convergence study:")
        print(f"  Batch size: {self.batch_size}")
        print(f"  Max orientations: {self.max_orientations}")
        print(
            f"  Tolerance: {self.convergable.tolerance * 100:.1f}% ({self.convergable.tolerance_type})"
        )
        print(f"  Min batches: {self.min_batches}")

        converged = False

        while not converged and self.n_orientations < self.max_orientations:
            # Create orientations for this batch
            orientations = goad.create_uniform_orientation(self.batch_size)
            self.settings.orientation = orientations

            # Run MultiProblem
            mp = goad.MultiProblem(self.settings)
            mp.py_solve()

            # Update statistics
            self._update_statistics(mp.results, self.batch_size)

            # Check convergence
            converged = self._check_convergence()

            # Print progress
            min_required = self.min_batches * self.batch_size
            if self.n_orientations < min_required:
                print(
                    f"\nBatch {len(self.batch_data)} ({self.n_orientations}/{min_required} total orientations, min not reached):"
                )
            else:
                print(
                    f"\nBatch {len(self.batch_data)} ({self.n_orientations} total orientations, min {min_required} reached):"
                )
            self._print_progress(converged)

            # Store history
            mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()
            # Use worst-case SEM for history
            valid_mask = ~np.isnan(mean_dscs)
            if np.any(valid_mask):
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative_sem = sem_dscs[valid_mask] / np.abs(mean_dscs[valid_mask])
                worst_sem = np.max(relative_sem)
                self.convergence_history.append(
                    (self.n_orientations, "phips_dscs", worst_sem)
                )

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

        # Print summary
        print(f"\n{'=' * 60}")
        if converged:
            print(f"✓ Converged after {self.n_orientations} orientations")
        else:
            print(f"✗ Did not converge (reached {self.n_orientations} orientations)")
        print(f"{'=' * 60}")

        # Print detector summary
        print(f"\nPHIPS Detector DSCS Summary:")
        valid_mask = ~np.isnan(mean_dscs)
        for i in range(self.NUM_DETECTORS):
            theta = self.detector_centers[i]
            if valid_mask[i]:
                mean_val = mean_dscs[i]
                sem_val = sem_dscs[i]
                rel_sem = sem_val / abs(mean_val) * 100
                print(
                    f"  Detector {i:2d} (θ={theta:6.1f}°): {mean_val:.6e} ± {sem_val:.6e} ({rel_sem:.2f}%)"
                )
            else:
                print(f"  Detector {i:2d} (θ={theta:6.1f}°): No data")

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
        )

    def run(self) -> ConvergenceResults:
        """
        Run ensemble convergence study.

        Each batch iteration randomly selects a geometry file from the
        ensemble directory before running the orientation averaging.

        Returns:
            ConvergenceResults with ensemble-averaged PHIPS DSCS values
        """
        print(f"\nStarting PHIPS Ensemble convergence study:")
        print(f"  Geometry files: {len(self.geom_files)}")
        print(f"  Batch size: {self.batch_size}")
        print(f"  Max orientations: {self.max_orientations}")
        print(
            f"  Tolerance: {self.convergable.tolerance * 100:.1f}% ({self.convergable.tolerance_type})"
        )
        print(f"  Min batches: {self.min_batches}")

        converged = False

        while not converged and self.n_orientations < self.max_orientations:
            # Randomly select a geometry file for this batch
            geom_file = random.choice(self.geom_files)
            geom_path = os.path.join(self.geom_dir, geom_file)

            # Create orientations for this batch
            orientations = goad.create_uniform_orientation(self.batch_size)

            # Update settings with selected geometry and orientations
            self.settings.geom_path = geom_path
            self.settings.orientation = orientations

            # Run MultiProblem
            mp = goad.MultiProblem(self.settings)
            mp.py_solve()

            # Update statistics
            self._update_statistics(mp.results, self.batch_size)

            # Check convergence
            converged = self._check_convergence()

            # Print progress (with geometry info)
            min_required = self.min_batches * self.batch_size
            if self.n_orientations < min_required:
                print(
                    f"\nBatch {len(self.batch_data)} ({self.n_orientations}/{min_required} total orientations, min not reached) - Geometry: {geom_file}"
                )
            else:
                print(
                    f"\nBatch {len(self.batch_data)} ({self.n_orientations} total orientations, min {min_required} reached) - Geometry: {geom_file}"
                )
            self._print_progress(converged)

            # Store history
            mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()
            valid_mask = ~np.isnan(mean_dscs)
            if np.any(valid_mask):
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative_sem = sem_dscs[valid_mask] / np.abs(mean_dscs[valid_mask])
                worst_sem = np.max(relative_sem)
                self.convergence_history.append(
                    (self.n_orientations, "phips_dscs", worst_sem)
                )

        # Compute final results
        mean_dscs, sem_dscs = self._calculate_phips_mean_and_sem()

        # Create results
        results = ConvergenceResults(
            converged=converged,
            n_orientations=self.n_orientations,
            values={"phips_dscs": mean_dscs},
            sem_values={"phips_dscs": sem_dscs},
            mueller_1d=None,
            mueller_2d=None,
            convergence_history=self.convergence_history,
            warning=None
            if converged
            else f"Did not converge within {self.max_orientations} orientations",
        )

        # Print summary
        print(f"\n{'=' * 60}")
        if converged:
            print(f"✓ Ensemble Converged after {self.n_orientations} orientations")
        else:
            print(f"✗ Did not converge (reached {self.n_orientations} orientations)")
        print(f"  Geometries sampled: {len(self.geom_files)}")
        print(f"{'=' * 60}")

        # Print detector summary
        print(f"\nPHIPS Detector DSCS Summary (Ensemble-Averaged):")
        valid_mask = ~np.isnan(mean_dscs)
        for i in range(self.NUM_DETECTORS):
            theta = self.detector_centers[i]
            if valid_mask[i]:
                mean_val = mean_dscs[i]
                sem_val = sem_dscs[i]
                rel_sem = sem_val / abs(mean_val) * 100
                print(
                    f"  Detector {i:2d} (θ={theta:6.1f}°): {mean_val:.6e} ± {sem_val:.6e} ({rel_sem:.2f}%)"
                )
            else:
                print(f"  Detector {i:2d} (θ={theta:6.1f}°): No data")

        return results
