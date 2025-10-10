from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np
from . import _goad_py as goad
import os
import random
from pathlib import Path


@dataclass
class Convergable:
    """Represents a variable to monitor for convergence."""

    variable: str  # 'asymmetry', 'scatt', 'ext', 'albedo', or Mueller element like 'S11', 'S12', etc.
    tolerance_type: str = "relative"  # 'relative' or 'absolute'
    tolerance: float = 0.01
    theta_indices: Optional[List[int]] = (
        None  # For Mueller elements: specific theta bin indices to check (None = all bins)
    )

    def __post_init__(self):
        # Scalar integrated parameters
        valid_scalars = {"asymmetry", "scatt", "ext", "albedo"}
        # Mueller matrix elements (S11, S12, ..., S44)
        valid_mueller = {f"S{i}{j}" for i in range(1, 5) for j in range(1, 5)}
        valid_variables = valid_scalars | valid_mueller

        if self.variable not in valid_variables:
            raise ValueError(
                f"Invalid variable '{self.variable}'. Must be one of {valid_scalars} or Mueller element (S11-S44)"
            )

        valid_types = {"relative", "absolute"}
        if self.tolerance_type not in valid_types:
            raise ValueError(
                f"Invalid tolerance_type '{self.tolerance_type}'. Must be one of {valid_types}"
            )

        if self.tolerance <= 0:
            raise ValueError(f"Tolerance must be positive, got {self.tolerance}")

        # Validate theta_indices only for Mueller elements
        if self.theta_indices is not None:
            if not self.is_mueller():
                raise ValueError("theta_indices can only be used with Mueller elements")
            # Convert range to list if needed
            if isinstance(self.theta_indices, range):
                self.theta_indices = list(self.theta_indices)
            if not isinstance(self.theta_indices, list):
                raise ValueError("theta_indices must be a list or range of integers")

    def is_mueller(self) -> bool:
        """Check if this convergable is a Mueller matrix element."""
        return self.variable.startswith("S") and len(self.variable) == 3


@dataclass
class ConvergenceResults:
    """Results from a convergence study."""

    converged: bool
    n_orientations: int
    values: Dict[str, float]  # Final mean values for each tracked variable
    sem_values: Dict[str, float]  # Final SEM values for each tracked variable
    mueller_1d: Optional[np.ndarray] = None
    mueller_2d: Optional[np.ndarray] = None
    convergence_history: List[Tuple[int, str, float]] = (
        None  # (n_orientations, variable, sem)
    )
    warning: Optional[str] = None


class Convergence:
    """Runs multiple MultiProblems until convergence criteria are met."""

    def __init__(
        self,
        settings: goad.Settings,
        convergables: List[Convergable],
        batch_size: int = 24,
        max_orientations: int = 100_000,
        min_batches: int = 10,
        mueller_1d: bool = True,
        mueller_2d: bool = False,
    ):
        """
        Initialize a convergence study.

        Args:
            settings: GOAD settings for the simulation
            convergables: List of variables to monitor for convergence
            batch_size: Number of orientations per iteration
            max_orientations: Maximum total orientations before stopping
            min_batches: Minimum number of batches before allowing convergence
            mueller_1d: Whether to collect 1D Mueller matrices
            mueller_2d: Whether to collect 2D Mueller matrices
        """
        self.settings = settings
        self.convergables = convergables
        self.batch_size = batch_size
        self.max_orientations = max_orientations
        self.min_batches = min_batches
        self.mueller_1d = mueller_1d
        self.mueller_2d = mueller_2d

        # Validate inputs
        if not convergables:
            raise ValueError("Must specify at least one convergable")

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

        # Batch-based statistics tracking for rigorous SEM calculation
        self.batch_data = []  # List of batch statistics

        # Mueller matrix accumulation
        self.mueller_1d_sum = None
        self.mueller_2d_sum = None

    def _update_statistics(self, results: goad.Results, batch_size: int):
        """Update statistics with new batch results.

        Args:
            results: Results from a MultiProblem run (pre-averaged over batch_size orientations)
            batch_size: Number of orientations in this batch
        """
        # Check for None values indicating Custom binning
        if (
            results.asymmetry is None
            or results.scat_cross is None
            or results.ext_cross is None
            or results.albedo is None
        ):
            raise ValueError(
                "Received None values for integrated properties. "
                "This likely means Custom binning scheme is being used. "
                "Convergence requires Simple or Interval binning schemes."
            )

        # Store batch data for proper statistical analysis
        batch_info = {"batch_size": batch_size, "values": {}, "weights": {}}

        # Always store all 4 integrated parameters (for unified API output)
        batch_info["values"]["asymmetry"] = results.asymmetry
        batch_info["weights"]["asymmetry"] = results.scat_cross
        batch_info["values"]["scatt"] = results.scat_cross
        batch_info["weights"]["scatt"] = 1.0  # Equal weighting
        batch_info["values"]["ext"] = results.ext_cross
        batch_info["weights"]["ext"] = 1.0  # Equal weighting
        batch_info["values"]["albedo"] = results.albedo
        batch_info["weights"]["albedo"] = results.ext_cross + results.scat_cross

        # Always store ALL 16 Mueller elements (for unified API output with full SEM)
        if self.mueller_1d and results.mueller_1d is not None:
            mueller_1d_array = np.array(results.mueller_1d)  # Shape: (n_theta, 16)

            # Store all 16 Mueller elements (S11, S12, ..., S44)
            for row in range(1, 5):
                for col in range(1, 5):
                    element_name = f"S{row}{col}"
                    mueller_idx = (row - 1) * 4 + (col - 1)
                    mueller_element = mueller_1d_array[
                        :, mueller_idx
                    ]  # Shape: (n_theta,)

                    batch_info["values"][element_name] = mueller_element
                    batch_info["weights"][element_name] = 1.0  # Equal weighting

            # Store theta bins if not already stored (for display purposes)
            if "mueller_theta_bins" not in batch_info and results.bins_1d is not None:
                batch_info["mueller_theta_bins"] = np.array(results.bins_1d)

        self.batch_data.append(batch_info)

        # Update Mueller matrices if enabled
        if self.mueller_1d and results.mueller_1d is not None:
            mueller_1d_array = np.array(results.mueller_1d)
            if self.mueller_1d_sum is None:
                self.mueller_1d_sum = mueller_1d_array * batch_size
            else:
                self.mueller_1d_sum += mueller_1d_array * batch_size

        if self.mueller_2d and results.mueller is not None:
            mueller_2d_array = np.array(results.mueller)
            if self.mueller_2d_sum is None:
                self.mueller_2d_sum = mueller_2d_array * batch_size
            else:
                self.mueller_2d_sum += mueller_2d_array * batch_size

        # Update total orientation count
        self.n_orientations += batch_size

    def _calculate_mean_and_sem_array(
        self, variable: str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate mean and SEM arrays for Mueller matrix elements across theta bins.

        Args:
            variable: Mueller element name (e.g., 'S11')

        Returns:
            Tuple of (mean_array, sem_array) where each is shape (n_theta,)
        """
        if not self.batch_data:
            return np.array([]), np.array([])

        # Extract batch values (each is an array of theta values)
        batch_arrays = []
        batch_sizes = []

        for batch in self.batch_data:
            if variable in batch["values"]:
                batch_arrays.append(batch["values"][variable])  # Shape: (n_theta,)
                batch_sizes.append(batch["batch_size"])

        if not batch_arrays:
            return np.array([]), np.array([])

        # Stack batches: shape (n_batches, n_theta)
        batch_arrays = np.array(batch_arrays)
        batch_sizes = np.array(batch_sizes)
        n_theta = batch_arrays.shape[1]

        if len(batch_arrays) < 2:
            # Can't estimate variance with < 2 batches
            mean_array = batch_arrays[0]
            sem_array = np.full(n_theta, float("inf"))
            return mean_array, sem_array

        # Calculate mean and SEM independently for each theta bin
        # Mean: weighted by batch size
        mean_array = np.average(
            batch_arrays, axis=0, weights=batch_sizes
        )  # Shape: (n_theta,)

        # Variance between batches at each theta
        batch_means_variance = np.var(batch_arrays, axis=0, ddof=1)  # Shape: (n_theta,)

        # Scale up to estimate population variance
        avg_batch_size = np.mean(batch_sizes)
        estimated_population_variance = batch_means_variance * avg_batch_size

        # Calculate SEM for total sample
        total_n = np.sum(batch_sizes)
        sem_array = np.sqrt(
            estimated_population_variance / (total_n - 1)
        )  # Shape: (n_theta,)

        return mean_array, sem_array

    def _calculate_mean_and_sem(self, variable: str) -> Tuple[float, float]:
        """Calculate mean and standard error of the mean for a variable using batch data.

        Args:
            variable: Variable name

        Returns:
            Tuple of (mean, sem)
        """
        if not self.batch_data:
            return 0.0, float("inf")

        # Extract batch values and weights
        batch_values = []
        batch_weights = []
        batch_sizes = []

        for batch in self.batch_data:
            if variable in batch["values"]:
                batch_values.append(batch["values"][variable])
                batch_weights.append(batch["weights"][variable])
                batch_sizes.append(batch["batch_size"])

        if not batch_values:
            return 0.0, float("inf")

        batch_values = np.array(batch_values)
        batch_weights = np.array(batch_weights)
        batch_sizes = np.array(batch_sizes)

        # For weighted variables (asymmetry, albedo), use weighted statistics
        if variable in ["asymmetry", "albedo"]:
            # Calculate weighted mean across batches
            # Each batch contributes: weight * batch_size * value
            total_weighted_sum = np.sum(batch_weights * batch_sizes * batch_values)
            total_weight = np.sum(batch_weights * batch_sizes)
            weighted_mean = total_weighted_sum / total_weight

            # Calculate weighted variance between batches
            if len(batch_values) < 2:
                return weighted_mean, float(
                    "inf"
                )  # Cannot estimate variance with < 2 batches

            # For batch means, we need to account for the effective weight of each batch
            effective_weights = batch_weights * batch_sizes
            weighted_variance_batch_means = np.sum(
                effective_weights * (batch_values - weighted_mean) ** 2
            ) / np.sum(effective_weights)

            # Scale up to estimate population variance
            # Batch means have variance = population_variance / average_batch_size
            # So population_variance ≈ batch_means_variance * average_batch_size
            avg_batch_size = np.average(batch_sizes, weights=effective_weights)
            estimated_population_variance = (
                weighted_variance_batch_means * avg_batch_size
            )

            # Calculate SEM for the total sample (using n-1 for sample standard error)
            total_n = np.sum(batch_sizes)
            sem = np.sqrt(estimated_population_variance / (total_n - 1))

            return weighted_mean, sem

        else:
            # For unweighted variables (scatt, ext), use simple batch statistics
            # Calculate mean of batch means, weighted by batch size
            total_sum = np.sum(batch_sizes * batch_values)
            total_n = np.sum(batch_sizes)
            mean = total_sum / total_n

            # Calculate variance between batch means
            if len(batch_values) < 2:
                return mean, float("inf")

            batch_means_variance = np.var(batch_values, ddof=1)

            # Scale up to estimate population variance
            # Batch means have variance = population_variance / average_batch_size
            # So population_variance ≈ batch_means_variance * average_batch_size
            avg_batch_size = np.mean(batch_sizes)
            estimated_population_variance = batch_means_variance * avg_batch_size

            # Calculate SEM for the total sample (using n-1 for sample standard error)
            sem = np.sqrt(estimated_population_variance / (total_n - 1))

            return mean, sem

    def _check_convergence(self) -> Dict[str, bool]:
        """Check if all convergence criteria are met.

        Returns:
            Dict mapping variable names to convergence status
        """
        converged = {}

        for conv in self.convergables:
            if conv.is_mueller():
                # Mueller element - check theta bins (all or specific indices)
                mean_array, sem_array = self._calculate_mean_and_sem_array(
                    conv.variable
                )

                if len(mean_array) == 0:
                    converged[conv.variable] = False
                    continue

                # Select theta bins to check
                if conv.theta_indices is not None:
                    # Check only specified indices
                    indices = [i for i in conv.theta_indices if i < len(mean_array)]
                    if not indices:
                        converged[conv.variable] = False
                        continue
                    mean_subset = mean_array[indices]
                    sem_subset = sem_array[indices]
                else:
                    # Check all bins
                    mean_subset = mean_array
                    sem_subset = sem_array

                # Check convergence at selected theta bins
                if conv.tolerance_type == "relative":
                    # Relative tolerance: SEM / |mean| < tolerance
                    relative_sem = np.where(
                        mean_subset != 0,
                        sem_subset / np.abs(mean_subset),
                        sem_subset / conv.tolerance,
                    )
                    converged[conv.variable] = np.all(relative_sem < conv.tolerance)
                else:
                    # Absolute tolerance: SEM < tolerance
                    converged[conv.variable] = np.all(sem_subset < conv.tolerance)
            else:
                # Scalar variable
                mean, sem = self._calculate_mean_and_sem(conv.variable)

                # Calculate tolerance based on type
                if conv.tolerance_type == "relative":
                    # Relative tolerance: SEM / |mean| < tolerance
                    if mean != 0:
                        relative_sem = sem / abs(mean)
                        converged[conv.variable] = relative_sem < conv.tolerance
                    else:
                        # If mean is zero, use absolute comparison
                        converged[conv.variable] = sem < conv.tolerance
                else:
                    # Absolute tolerance: SEM < tolerance
                    converged[conv.variable] = sem < conv.tolerance

        return converged

    def _all_converged(self) -> bool:
        """Check if all variables have converged.

        Returns:
            True if all variables meet their convergence criteria and minimum batches completed
        """
        # Check minimum batches requirement first
        if len(self.batch_data) < self.min_batches:
            return False

        converged_status = self._check_convergence()
        return all(converged_status.values())

    def _print_progress(self, iteration: int):
        """Print convergence progress.

        Args:
            iteration: Current iteration number
        """
        # Calculate minimum required orientations
        min_required = self.min_batches * self.batch_size

        # Show progress with min orientations requirement
        if self.n_orientations < min_required:
            print(
                f"\nIteration {iteration} ({self.n_orientations}/{min_required} orientations, min not reached):"
            )
        else:
            print(
                f"\nIteration {iteration} ({self.n_orientations} orientations, min {min_required} reached):"
            )

        converged_status = self._check_convergence()

        for conv in self.convergables:
            if conv.is_mueller():
                # Mueller element - show worst-case theta bin
                mean_array, sem_array = self._calculate_mean_and_sem_array(
                    conv.variable
                )

                if len(mean_array) == 0:
                    print(f"  {conv.variable:<10}: No data yet")
                    continue

                # Get theta bins from results (assuming we have access to bins_1d)
                if hasattr(self, "settings") and hasattr(self.settings, "binning"):
                    # We'll get theta values from the first batch's mueller_1d if available
                    theta_bins = None
                    for batch in self.batch_data:
                        if "mueller_theta_bins" in batch:
                            theta_bins = batch["mueller_theta_bins"]
                            break
                    if theta_bins is None:
                        theta_bins = np.arange(len(mean_array))
                else:
                    theta_bins = np.arange(len(mean_array))

                # Calculate relative SEM for each theta
                if conv.tolerance_type == "relative":
                    relative_sem_array = np.where(
                        mean_array != 0, sem_array / np.abs(mean_array), float("inf")
                    )
                    worst_idx = np.argmax(relative_sem_array)
                    worst_sem = relative_sem_array[worst_idx]
                    target_str = f"{conv.tolerance * 100:.1f}%"
                    current_str = f"{worst_sem * 100:.2f}%"
                else:
                    worst_idx = np.argmax(sem_array)
                    worst_sem = sem_array[worst_idx]
                    target_str = f"{conv.tolerance}"
                    current_str = f"{worst_sem:.4g}"

                worst_theta = theta_bins[worst_idx]
                worst_mean = mean_array[worst_idx]

                # Count converged bins (either all or specified indices)
                if conv.theta_indices is not None:
                    # Only checking specific bins
                    indices = [i for i in conv.theta_indices if i < len(mean_array)]
                    if conv.tolerance_type == "relative":
                        converged_bins = np.sum(
                            relative_sem_array[indices] < conv.tolerance
                        )
                    else:
                        converged_bins = np.sum(sem_array[indices] < conv.tolerance)
                    total_bins = len(indices)
                    bin_desc = (
                        f"θ={[theta_bins[i] for i in indices]}"
                        if len(indices) <= 3
                        else f"{len(indices)} bins"
                    )
                else:
                    # Checking all bins
                    if conv.tolerance_type == "relative":
                        converged_bins = np.sum(relative_sem_array < conv.tolerance)
                    else:
                        converged_bins = np.sum(sem_array < conv.tolerance)
                    total_bins = len(mean_array)
                    bin_desc = f"{total_bins} bins"

                status = "✓" if converged_status[conv.variable] else "❌"

                # Print Mueller convergence info
                if conv.theta_indices is not None and len(conv.theta_indices) <= 3:
                    # For small number of specific bins, show them explicitly
                    print(
                        f"  {conv.variable:<10}: {converged_bins}/{total_bins} {bin_desc} | "
                        f"Worst θ={worst_theta:.1f}°: {worst_mean:.4g} | SEM: {current_str} (target: {target_str}) {status}"
                    )
                else:
                    # For many bins, use standard format
                    print(
                        f"  {conv.variable:<10}: {converged_bins}/{total_bins} bins converged | "
                        f"Worst θ={worst_theta:.1f}°: {worst_mean:.4g} | SEM: {current_str} (target: {target_str}) {status}"
                    )

                # Add worst SEM to convergence history
                self.convergence_history.append(
                    (self.n_orientations, conv.variable, worst_sem)
                )
            else:
                # Scalar variable
                mean, sem = self._calculate_mean_and_sem(conv.variable)

                # Calculate 95% CI
                ci_lower = mean - 1.96 * sem
                ci_upper = mean + 1.96 * sem

                # Format based on tolerance type
                if conv.tolerance_type == "relative":
                    if mean != 0:
                        relative_sem = sem / abs(mean)
                        target_str = f"{conv.tolerance * 100:.1f}%"
                        current_str = f"{relative_sem * 100:.2f}%"
                    else:
                        target_str = f"{conv.tolerance} (abs, mean=0)"
                        current_str = f"{sem:.4g}"
                else:
                    target_str = f"{conv.tolerance}"
                    current_str = f"{sem:.4g}"

                # Status indicator
                status = "✓" if converged_status[conv.variable] else "❌"

                # Print line with mean, SEM, CI, and convergence status
                print(
                    f"  {conv.variable:<10}: {mean:.6f} ± {sem:.6f} [{ci_lower:.6f}, {ci_upper:.6f}] | SEM: {current_str} (target: {target_str}) {status}"
                )

                # Add to convergence history
                self.convergence_history.append(
                    (self.n_orientations, conv.variable, sem)
                )

    def run(self) -> ConvergenceResults:
        """Run the convergence study.

        Returns:
            ConvergenceResults containing final values and convergence status
        """
        iteration = 0
        converged = False
        warning = None

        while not converged and self.n_orientations < self.max_orientations:
            iteration += 1

            # Determine batch size for this iteration
            remaining = self.max_orientations - self.n_orientations
            batch_size = min(self.batch_size, remaining)

            # Set batch size
            orientations = goad.create_uniform_orientation(batch_size)

            # Set the orientations for the settings
            self.settings.orientation = orientations

            mp = goad.MultiProblem(self.settings)
            mp.py_solve()

            # Update statistics
            self._update_statistics(mp.results, batch_size)

            # Print progress
            self._print_progress(iteration)

            # Check convergence
            converged = self._all_converged()

        # Prepare final results
        if converged:
            print(f"\nConverged after {self.n_orientations} orientations.")
        else:
            warning = f"Maximum orientations ({self.max_orientations}) reached without convergence"
            print(f"\nWarning: {warning}")

        # Calculate final values and SEMs
        final_values = {}
        final_sems = {}
        for conv in self.convergables:
            if conv.is_mueller():
                mean_array, sem_array = self._calculate_mean_and_sem_array(
                    conv.variable
                )
                final_values[conv.variable] = mean_array
                final_sems[conv.variable] = sem_array
            else:
                mean, sem = self._calculate_mean_and_sem(conv.variable)
                final_values[conv.variable] = mean
                final_sems[conv.variable] = sem

        # Prepare Mueller matrices with SEM
        mueller_1d = None
        mueller_1d_sem = None
        mueller_2d = None

        if self.mueller_1d and self.mueller_1d_sum is not None:
            mueller_1d = self.mueller_1d_sum / self.n_orientations

            # Compute SEM for all 16 Mueller elements
            # mueller_1d shape: (n_theta, 16)
            n_theta = mueller_1d.shape[0]
            mueller_1d_sem = np.zeros_like(mueller_1d)

            for row in range(1, 5):
                for col in range(1, 5):
                    element_name = f"S{row}{col}"
                    mueller_idx = (row - 1) * 4 + (col - 1)

                    # Calculate mean and SEM for this element across all theta bins
                    mean_array, sem_array = self._calculate_mean_and_sem_array(
                        element_name
                    )

                    if len(sem_array) > 0:
                        mueller_1d_sem[:, mueller_idx] = sem_array

            # Store mueller_1d_sem in final_values for unified API access
            final_values["mueller_1d_sem"] = mueller_1d_sem

        if self.mueller_2d and self.mueller_2d_sum is not None:
            mueller_2d = self.mueller_2d_sum / self.n_orientations

        return ConvergenceResults(
            converged=converged,
            n_orientations=self.n_orientations,
            values=final_values,
            sem_values=final_sems,
            mueller_1d=mueller_1d,
            mueller_2d=mueller_2d,
            convergence_history=self.convergence_history,
            warning=warning,
        )


class EnsembleConvergence(Convergence):
    """Runs convergence study over an ensemble of particle geometries.

    Each batch randomly samples from a directory of geometry files,
    allowing convergence analysis of orientation-averaged and
    geometry-averaged scattering properties.
    """

    def __init__(
        self,
        settings: goad.Settings,
        convergables: List[Convergable],
        geom_dir: str,
        batch_size: int = 24,
        max_orientations: int = 100_000,
        min_batches: int = 10,
        mueller_1d: bool = True,
        mueller_2d: bool = False,
    ):
        """
        Initialize an ensemble convergence study.

        Args:
            settings: GOAD settings for the simulation (geom_path will be overridden)
            convergables: List of variables to monitor for convergence
            geom_dir: Directory containing .obj geometry files
            batch_size: Number of orientations per iteration
            max_orientations: Maximum total orientations before stopping
            min_batches: Minimum number of batches before allowing convergence
            mueller_1d: Whether to collect 1D Mueller matrices
            mueller_2d: Whether to collect 2D Mueller matrices
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
            convergables=convergables,
            batch_size=batch_size,
            max_orientations=max_orientations,
            min_batches=min_batches,
            mueller_1d=mueller_1d,
            mueller_2d=mueller_2d,
        )

    def run(self) -> ConvergenceResults:
        """Run the ensemble convergence study.

        Each batch iteration randomly selects a geometry file from the
        ensemble directory before running the orientation averaging.

        Returns:
            ConvergenceResults containing final ensemble-averaged values
        """
        iteration = 0
        converged = False
        warning = None

        while not converged and self.n_orientations < self.max_orientations:
            iteration += 1

            # Randomly select a geometry file for this batch
            geom_file = random.choice(self.geom_files)
            geom_path = os.path.join(self.geom_dir, geom_file)

            # Determine batch size for this iteration
            remaining = self.max_orientations - self.n_orientations
            batch_size = min(self.batch_size, remaining)

            # Create orientations for this batch
            orientations = goad.create_uniform_orientation(batch_size)

            # Update settings with selected geometry and orientations
            self.settings.geom_path = geom_path
            self.settings.orientation = orientations

            # Run MultiProblem with selected geometry
            mp = goad.MultiProblem(self.settings)
            mp.py_solve()

            # Update statistics
            self._update_statistics(mp.results, batch_size)

            # Print progress (with geometry info)
            min_required = self.min_batches * self.batch_size
            if self.n_orientations < min_required:
                print(
                    f"\nIteration {iteration} ({self.n_orientations}/{min_required} orientations, min not reached) - Geometry: {geom_file}"
                )
            else:
                print(
                    f"\nIteration {iteration} ({self.n_orientations} orientations, min {min_required} reached) - Geometry: {geom_file}"
                )
            self._print_progress_without_header(iteration)

            # Check convergence
            converged = self._all_converged()

        # Prepare final results
        if converged:
            print(f"\nConverged after {self.n_orientations} orientations.")
        else:
            warning = f"Maximum orientations ({self.max_orientations}) reached without convergence"
            print(f"\nWarning: {warning}")

        # Calculate final values and SEMs
        final_values = {}
        final_sems = {}
        for conv in self.convergables:
            if conv.is_mueller():
                mean_array, sem_array = self._calculate_mean_and_sem_array(
                    conv.variable
                )
                final_values[conv.variable] = mean_array
                final_sems[conv.variable] = sem_array
            else:
                mean, sem = self._calculate_mean_and_sem(conv.variable)
                final_values[conv.variable] = mean
                final_sems[conv.variable] = sem

        # Prepare Mueller matrices with SEM
        mueller_1d = None
        mueller_1d_sem = None
        mueller_2d = None

        if self.mueller_1d and self.mueller_1d_sum is not None:
            mueller_1d = self.mueller_1d_sum / self.n_orientations

            # Compute SEM for all 16 Mueller elements
            # mueller_1d shape: (n_theta, 16)
            n_theta = mueller_1d.shape[0]
            mueller_1d_sem = np.zeros_like(mueller_1d)

            for row in range(1, 5):
                for col in range(1, 5):
                    element_name = f"S{row}{col}"
                    mueller_idx = (row - 1) * 4 + (col - 1)

                    # Calculate mean and SEM for this element across all theta bins
                    mean_array, sem_array = self._calculate_mean_and_sem_array(
                        element_name
                    )

                    if len(sem_array) > 0:
                        mueller_1d_sem[:, mueller_idx] = sem_array

            # Store mueller_1d_sem in final_values for unified API access
            final_values["mueller_1d_sem"] = mueller_1d_sem

        if self.mueller_2d and self.mueller_2d_sum is not None:
            mueller_2d = self.mueller_2d_sum / self.n_orientations

        return ConvergenceResults(
            converged=converged,
            n_orientations=self.n_orientations,
            values=final_values,
            sem_values=final_sems,
            mueller_1d=mueller_1d,
            mueller_2d=mueller_2d,
            convergence_history=self.convergence_history,
            warning=warning,
        )

    def _print_progress_without_header(self, iteration: int):
        """Print convergence progress without iteration header (already printed with geometry)."""
        converged_status = self._check_convergence()

        for conv in self.convergables:
            mean, sem = self._calculate_mean_and_sem(conv.variable)

            # Calculate 95% CI
            ci_lower = mean - 1.96 * sem
            ci_upper = mean + 1.96 * sem

            # Format based on tolerance type
            if conv.tolerance_type == "relative":
                if mean != 0:
                    relative_sem = sem / abs(mean)
                    target_str = f"{conv.tolerance * 100:.1f}%"
                    current_str = f"{relative_sem * 100:.2f}%"
                else:
                    target_str = f"{conv.tolerance} (abs, mean=0)"
                    current_str = f"{sem:.4g}"
            else:
                target_str = f"{conv.tolerance}"
                current_str = f"{sem:.4g}"

            # Status indicator
            status = "✓" if converged_status[conv.variable] else "❌"

            # Print line with mean, SEM, CI, and convergence status
            print(
                f"  {conv.variable:<10}: {mean:.6f} ± {sem:.6f} [{ci_lower:.6f}, {ci_upper:.6f}] | SEM: {current_str} (target: {target_str}) {status}"
            )

            # Add to convergence history
            self.convergence_history.append((self.n_orientations, conv.variable, sem))
