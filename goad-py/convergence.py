from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np
import goad_py as goad


@dataclass
class Convergable:
    """Represents a variable to monitor for convergence."""
    variable: str  # 'asymmetry', 'scatt', 'ext', or 'albedo'
    tolerance_type: str = 'relative'  # 'relative' or 'absolute'
    tolerance: float = 0.01

    def __post_init__(self):
        valid_variables = {'asymmetry', 'scatt', 'ext', 'albedo'}
        if self.variable not in valid_variables:
            raise ValueError(f"Invalid variable '{self.variable}'. Must be one of {valid_variables}")

        valid_types = {'relative', 'absolute'}
        if self.tolerance_type not in valid_types:
            raise ValueError(f"Invalid tolerance_type '{self.tolerance_type}'. Must be one of {valid_types}")

        if self.tolerance <= 0:
            raise ValueError(f"Tolerance must be positive, got {self.tolerance}")


@dataclass
class ConvergenceResults:
    """Results from a convergence study."""
    converged: bool
    n_orientations: int
    values: Dict[str, float]  # Final mean values for each tracked variable
    sem_values: Dict[str, float]  # Final SEM values for each tracked variable
    mueller_1d: Optional[np.ndarray] = None
    mueller_2d: Optional[np.ndarray] = None
    convergence_history: List[Tuple[int, str, float]] = None  # (n_orientations, variable, sem)
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
        mueller_2d: bool = False
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
            raise ValueError(f"max_orientations must be positive, got {max_orientations}")
        
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
        if (results.asymmetry is None or results.scat_cross is None or
            results.ext_cross is None or results.albedo is None):
            raise ValueError(
                "Received None values for integrated properties. "
                "This likely means Custom binning scheme is being used. "
                "Convergence requires Simple or Interval binning schemes."
            )

        # Store batch data for proper statistical analysis
        batch_info = {
            'batch_size': batch_size,
            'values': {},
            'weights': {}
        }

        # Store values and weights for tracked variables
        for conv in self.convergables:
            if conv.variable == 'asymmetry':
                batch_info['values']['asymmetry'] = results.asymmetry
                batch_info['weights']['asymmetry'] = results.scat_cross
            elif conv.variable == 'scatt':
                batch_info['values']['scatt'] = results.scat_cross
                batch_info['weights']['scatt'] = 1.0  # Equal weighting
            elif conv.variable == 'ext':
                batch_info['values']['ext'] = results.ext_cross
                batch_info['weights']['ext'] = 1.0  # Equal weighting
            elif conv.variable == 'albedo':
                batch_info['values']['albedo'] = results.albedo
                batch_info['weights']['albedo'] = results.ext_cross + results.scat_cross

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

    def _calculate_mean_and_sem(self, variable: str) -> Tuple[float, float]:
        """Calculate mean and standard error of the mean for a variable using batch data.

        Args:
            variable: Variable name

        Returns:
            Tuple of (mean, sem)
        """
        if not self.batch_data:
            return 0.0, float('inf')

        # Extract batch values and weights
        batch_values = []
        batch_weights = []
        batch_sizes = []

        for batch in self.batch_data:
            if variable in batch['values']:
                batch_values.append(batch['values'][variable])
                batch_weights.append(batch['weights'][variable])
                batch_sizes.append(batch['batch_size'])

        if not batch_values:
            return 0.0, float('inf')

        batch_values = np.array(batch_values)
        batch_weights = np.array(batch_weights)
        batch_sizes = np.array(batch_sizes)

        # For weighted variables (asymmetry, albedo), use weighted statistics
        if variable in ['asymmetry', 'albedo']:
            # Calculate weighted mean across batches
            # Each batch contributes: weight * batch_size * value
            total_weighted_sum = np.sum(batch_weights * batch_sizes * batch_values)
            total_weight = np.sum(batch_weights * batch_sizes)
            weighted_mean = total_weighted_sum / total_weight

            # Calculate weighted variance between batches
            if len(batch_values) < 2:
                return weighted_mean, float('inf')  # Cannot estimate variance with < 2 batches

            # For batch means, we need to account for the effective weight of each batch
            effective_weights = batch_weights * batch_sizes
            weighted_variance_batch_means = np.sum(effective_weights * (batch_values - weighted_mean)**2) / np.sum(effective_weights)

            # Scale up to estimate population variance
            # Batch means have variance = population_variance / average_batch_size
            # So population_variance ≈ batch_means_variance * average_batch_size
            avg_batch_size = np.average(batch_sizes, weights=effective_weights)
            estimated_population_variance = weighted_variance_batch_means * avg_batch_size

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
                return mean, float('inf')

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
            mean, sem = self._calculate_mean_and_sem(conv.variable)

            # Calculate tolerance based on type
            if conv.tolerance_type == 'relative':
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
        print(f"\nIteration {iteration} ({self.n_orientations} orientations):")

        converged_status = self._check_convergence()

        for conv in self.convergables:
            mean, sem = self._calculate_mean_and_sem(conv.variable)

            # Calculate 95% CI
            ci_lower = mean - 1.96 * sem
            ci_upper = mean + 1.96 * sem

            # Format based on tolerance type
            if conv.tolerance_type == 'relative':
                if mean != 0:
                    relative_sem = sem / abs(mean)
                    target_str = f"{conv.tolerance*100:.1f}%"
                    current_str = f"{relative_sem*100:.2f}%"
                else:
                    target_str = f"{conv.tolerance} (abs, mean=0)"
                    current_str = f"{sem:.4g}"
            else:
                target_str = f"{conv.tolerance}"
                current_str = f"{sem:.4g}"

            # Status indicator
            status = "✓" if converged_status[conv.variable] else "❌"

            # Print line with mean, SEM, CI, and convergence status
            print(f"  {conv.variable:<10}: {mean:.6f} ± {sem:.6f} [{ci_lower:.6f}, {ci_upper:.6f}] | SEM: {current_str} (target: {target_str}) {status}")

            # Add to convergence history
            self.convergence_history.append((self.n_orientations, conv.variable, sem))

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
            mean, sem = self._calculate_mean_and_sem(conv.variable)
            final_values[conv.variable] = mean
            final_sems[conv.variable] = sem

        # Prepare Mueller matrices
        mueller_1d = None
        mueller_2d = None
        if self.mueller_1d and self.mueller_1d_sum is not None:
            mueller_1d = self.mueller_1d_sum / self.n_orientations
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
            warning=warning
        )
