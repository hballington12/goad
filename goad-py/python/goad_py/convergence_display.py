"""
Rich-based convergence display for GOAD.

This module provides a unified, reusable display system for convergence tracking
that works with both standard convergence and PHIPS convergence modes.
"""

from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from rich.columns import Columns
from rich.console import Console, Group
from rich.live import Live
from rich.progress import BarColumn, Progress, TextColumn
from rich.spinner import Spinner
from rich.text import Text


class ConvergenceVariable:
    """
    Base class for convergence variables.

    Encapsulates the logic for calculating progress, formatting display values,
    and checking convergence for a single variable.
    """

    def __init__(
        self,
        name: str,
        tolerance: float,
        tolerance_type: str = "relative",
    ):
        """
        Initialize a convergence variable.

        Args:
            name: Variable name (e.g., "asymmetry", "S11", "phips_dscs")
            tolerance: Convergence tolerance threshold
            tolerance_type: "relative" or "absolute"
        """
        self.name = name
        self.tolerance = tolerance
        self.tolerance_type = tolerance_type

        if tolerance_type not in {"relative", "absolute"}:
            raise ValueError(f"Invalid tolerance_type '{tolerance_type}'")

        if tolerance <= 0:
            raise ValueError(f"Tolerance must be positive, got {tolerance}")

    def calculate_progress(self, mean: float, sem: float) -> float:
        """
        Calculate progress percentage (0-100) based on current SEM.

        Uses sqrt formula for smoother progression: sqrt(target/current) * 100

        Args:
            mean: Current mean value
            sem: Current SEM value

        Returns:
            Progress percentage (0-100)
        """
        if self.tolerance_type == "relative":
            if mean != 0:
                relative_sem = sem / abs(mean)
                if relative_sem > 0 and not np.isinf(relative_sem):
                    return min(100.0, np.sqrt(self.tolerance / relative_sem) * 100.0)
        else:
            # Absolute
            if sem > 0 and not np.isinf(sem):
                return min(100.0, np.sqrt(self.tolerance / sem) * 100.0)

        return 0.0

    def format_sem_info(self, mean: float, sem: float) -> str:
        """
        Format SEM info string for display.

        Args:
            mean: Current mean value
            sem: Current SEM value

        Returns:
            Formatted string like "[SEM: 19.7% / 20.0%]" or "[SEM: 0.0042 / 0.0050]"
        """
        if self.tolerance_type == "relative":
            if mean != 0:
                relative_sem = sem / abs(mean)
                current_str = f"{relative_sem * 100:.1f}%"
            else:
                current_str = f"{sem:.4g}"
            target_str = f"{self.tolerance * 100:.1f}%"
        else:
            current_str = f"{sem:.4g}"
            target_str = f"{self.tolerance:.4g}"

        return f"[SEM: {current_str} / {target_str}]"

    def is_converged(self, mean: float, sem: float) -> bool:
        """
        Check if this variable has converged.

        Args:
            mean: Current mean value
            sem: Current SEM value

        Returns:
            True if converged, False otherwise
        """
        if self.tolerance_type == "relative":
            if mean != 0:
                relative_sem = sem / abs(mean)
                return relative_sem < self.tolerance
            return False
        else:
            return sem < self.tolerance


class ArrayConvergenceVariable(ConvergenceVariable):
    """
    Convergence variable for array data (Mueller elements, PHIPS detectors).

    Tracks convergence across multiple bins/detectors, reporting worst-case progress.
    """

    def __init__(
        self,
        name: str,
        tolerance: float,
        tolerance_type: str = "relative",
        indices: Optional[List[int]] = None,
    ):
        """
        Initialize an array convergence variable.

        Args:
            name: Variable name
            tolerance: Convergence tolerance threshold
            tolerance_type: "relative" or "absolute"
            indices: Specific indices to check (None = all)
        """
        super().__init__(name, tolerance, tolerance_type)
        self.indices = indices

    def calculate_progress_array(
        self, mean_array: np.ndarray, sem_array: np.ndarray
    ) -> Tuple[float, int, float]:
        """
        Calculate progress for array data based on worst-case bin/detector.

        Args:
            mean_array: Array of mean values
            sem_array: Array of SEM values

        Returns:
            Tuple of (progress_percentage, worst_index, worst_sem)
        """
        if len(mean_array) == 0:
            return 0.0, -1, np.inf

        # Filter to specific indices if requested
        if self.indices is not None:
            mean_array = mean_array[self.indices]
            sem_array = sem_array[self.indices]

        # Find worst-case bin
        if self.tolerance_type == "relative":
            relative_sem_array = np.where(
                mean_array != 0, sem_array / np.abs(mean_array), float("inf")
            )
            worst_idx = np.argmax(relative_sem_array)
            worst_sem = relative_sem_array[worst_idx]
        else:
            worst_idx = np.argmax(sem_array)
            worst_sem = sem_array[worst_idx]

        # Calculate progress using worst SEM
        if worst_sem > 0 and not np.isinf(worst_sem):
            progress = min(100.0, np.sqrt(self.tolerance / worst_sem) * 100.0)
        else:
            progress = 0.0

        return progress, worst_idx, worst_sem

    def format_sem_info_array(self, worst_sem: float) -> str:
        """
        Format SEM info for array data.

        Args:
            worst_sem: Worst-case SEM (already relative if tolerance_type is relative)

        Returns:
            Formatted SEM string
        """
        if self.tolerance_type == "relative":
            current_str = f"{worst_sem * 100:.2f}%"
            target_str = f"{self.tolerance * 100:.1f}%"
        else:
            current_str = f"{worst_sem:.4g}"
            target_str = f"{self.tolerance:.4g}"

        return f"[SEM: {current_str} / {target_str}]"

    def count_converged(
        self, mean_array: np.ndarray, sem_array: np.ndarray
    ) -> Tuple[int, int]:
        """
        Count how many bins/detectors have converged.

        Args:
            mean_array: Array of mean values
            sem_array: Array of SEM values

        Returns:
            Tuple of (converged_count, total_count)
        """
        if len(mean_array) == 0:
            return 0, 0

        # Filter to specific indices if requested
        if self.indices is not None:
            mean_array = mean_array[self.indices]
            sem_array = sem_array[self.indices]

        # Check convergence for each bin
        if self.tolerance_type == "relative":
            relative_sem_array = np.where(
                mean_array != 0, sem_array / np.abs(mean_array), float("inf")
            )
            converged_mask = relative_sem_array < self.tolerance
        else:
            converged_mask = sem_array < self.tolerance

        converged_count = np.sum(converged_mask)
        total_count = len(mean_array)

        return int(converged_count), total_count

    def is_converged_array(self, mean_array: np.ndarray, sem_array: np.ndarray) -> bool:
        """
        Check if all bins/detectors have converged.

        Args:
            mean_array: Array of mean values
            sem_array: Array of SEM values

        Returns:
            True if all converged, False otherwise
        """
        converged_count, total_count = self.count_converged(mean_array, sem_array)
        return converged_count == total_count


class ConvergenceDisplay:
    """
    Rich-based display for convergence progress.

    Provides a unified, reusable display system that works with different
    convergence modes (standard, PHIPS, etc.).
    """

    def __init__(
        self,
        variables: List[ConvergenceVariable],
        batch_size: int,
        min_batches: int,
        convergence_type: str = "standard",
        console: Optional[Console] = None,
        log_file: Optional[str] = None,
    ):
        """
        Initialize convergence display.

        Args:
            variables: List of convergence variables to track
            batch_size: Number of orientations per batch
            min_batches: Minimum batches before convergence check
            convergence_type: Display string for convergence mode
            console: Optional Rich console (creates one if None)
            log_file: Optional path to log file for convergence progress
        """
        self.variables = variables
        self.batch_size = batch_size
        self.min_batches = min_batches
        self.convergence_type = convergence_type
        self._console = console or Console()

        # File logging
        self.log_file = log_file
        self._file_console = None
        self._file_handle = None
        if self.log_file:
            # Create a separate console for file output
            self._file_handle = open(self.log_file, "w")
            self._file_console = Console(file=self._file_handle, width=120)
            # Write header
            self._file_console.print(f"GOAD Convergence Log - {convergence_type}")
            self._file_console.print("=" * 120)

        # Progress tracking
        self._progress: Optional[Progress] = None
        self._progress_tasks: Dict[str, int] = {}

    def _initialize_progress(self):
        """Initialize Rich Progress instance and task IDs."""
        if self._progress is not None:
            return

        self._progress = Progress(
            TextColumn("[bold]{task.fields[variable]:<12}"),
            BarColumn(bar_width=25),
            TextColumn("[bold]{task.percentage:>3.0f}%"),
            TextColumn("[cyan]{task.fields[sem_info]}"),
            console=self._console,
            transient=False,
        )

        # Add tasks for each variable
        for var in self.variables:
            task_id = self._progress.add_task(
                "",
                total=100,
                variable=var.name,
                sem_info="[SEM: -- / --]",
            )
            self._progress_tasks[var.name] = task_id

    def build_display(
        self,
        iteration: int,
        n_orientations: int,
        get_stats: Callable[[str], Tuple[float, float]],
        get_array_stats: Optional[
            Callable[[str], Tuple[np.ndarray, np.ndarray]]
        ] = None,
        get_bin_labels: Optional[Callable[[str], np.ndarray]] = None,
        power_ratio: Optional[float] = None,
        geom_info: Optional[str] = None,
    ) -> Group:
        """
        Build rich display for current convergence state.

        Args:
            iteration: Current batch iteration number
            n_orientations: Total number of orientations processed
            get_stats: Callback to get (mean, sem) for a variable name
            get_array_stats: Optional callback to get (mean_array, sem_array) for array variables
            get_bin_labels: Optional callback to get bin labels (e.g., theta angles) for array variables
            power_ratio: Optional power ratio from solver
            geom_info: Optional geometry info (for ensemble mode)

        Returns:
            Rich Group containing the full display
        """
        # Initialize progress on first call
        self._initialize_progress()

        # Calculate minimum required orientations
        min_required = self.min_batches * self.batch_size

        # Build title with inline spinner
        spinner = Spinner("aesthetic", style="cyan")
        title_text = Text.assemble(
            ("GOAD: ", "bold cyan"),
            (f"[Convergence: {self.convergence_type}] ", "bold white"),
        )
        title = Columns([title_text, spinner], expand=False, padding=(0, 1))

        # Build batch header
        batch_str = f"[Batch: {iteration}/{self.min_batches}]"

        # Color orient based on whether min is reached
        orient_color = "green" if n_orientations >= min_required else "red"
        orient_text = Text(
            f"[Orient: {n_orientations}/{min_required} ({self.batch_size})]",
            style=orient_color,
        )

        # Power ratio (if available)
        if power_ratio is not None:
            power_color = self._get_power_ratio_color(power_ratio)
            power_text = Text(f"Power: {power_ratio:.3f}", style=power_color)
        else:
            power_text = Text("Power: N/A")

        # Build header
        header_parts = [batch_str, " ", orient_text, " ", power_text]

        if geom_info:
            header_parts.extend([" ", (f"[{geom_info}]", "dim")])

        header = Text.assemble(*header_parts)
        separator = Text("━" * 70)

        # Update progress for each variable
        progress_lines = []
        for var in self.variables:
            task_id = self._progress_tasks[var.name]

            if isinstance(var, ArrayConvergenceVariable):
                # Array variable (Mueller element, PHIPS detector)
                if get_array_stats is None:
                    raise ValueError(
                        f"get_array_stats callback required for array variable '{var.name}'"
                    )

                mean_array, sem_array = get_array_stats(var.name)

                if len(mean_array) == 0:
                    self._progress.update(
                        task_id, completed=0, sem_info="[SEM: -- / --]"
                    )
                    mean_str = "--"
                    progress_pct = 0
                    sem_info = "[SEM: -- / --]"
                else:
                    # Calculate progress and find worst bin
                    progress, worst_idx, worst_sem = var.calculate_progress_array(
                        mean_array, sem_array
                    )
                    sem_info = var.format_sem_info_array(worst_sem)

                    # Count converged bins
                    converged_count, total_count = var.count_converged(
                        mean_array, sem_array
                    )

                    # Get bin label (e.g., theta angle)
                    if get_bin_labels is not None:
                        bin_labels = get_bin_labels(var.name)
                        worst_label = bin_labels[worst_idx]
                        mean_str = f"({converged_count}/{total_count}) Worst: θ={worst_label:.0f}°"
                    else:
                        mean_str = (
                            f"({converged_count}/{total_count}) Worst: #{worst_idx}"
                        )

                    progress_pct = int(progress)

                    # Update progress
                    self._progress.update(
                        task_id, completed=progress, sem_info=sem_info
                    )
            else:
                # Scalar variable
                mean, sem = get_stats(var.name)

                progress = var.calculate_progress(mean, sem)
                sem_info = var.format_sem_info(mean, sem)
                progress_pct = int(progress)

                # Format mean value
                if not np.isinf(mean):
                    mean_str = f"{mean:.4f}"
                else:
                    mean_str = "--"

                # Update progress
                self._progress.update(task_id, completed=progress, sem_info=sem_info)

            # Render progress bar line
            task = self._progress.tasks[task_id]
            bar_width = 25
            filled = int((task.percentage / 100) * bar_width)
            bar = "█" * filled + "░" * (bar_width - filled)

            variable_text = f"{task.fields['variable']:<12}"
            mean_display = f"{mean_str:>35}"
            progress_text = f"{progress_pct}%"

            line = Text.assemble(
                f"{variable_text} ",
                (mean_display, "bold green"),
                f" [{bar}] {progress_text:>3} ",
                (sem_info, "cyan"),
            )
            progress_lines.append(line)

        # Build full display
        display = Group(
            separator,
            title,
            Text(""),  # Blank line
            header,
            separator,
            *progress_lines,
        )

        # Log to file if enabled
        if self._file_console is not None:
            self._file_console.print(display)

        return display

    def _get_power_ratio_color(self, power_ratio: float) -> str:
        """Get color for power ratio based on threshold."""
        if power_ratio >= 0.99:
            return "green"
        elif power_ratio >= 0.95:
            return "yellow"
        else:
            return "red"

    def create_live_context(self, refresh_per_second: float = 1.3) -> Live:
        """
        Create a Rich Live context for auto-updating display.

        Args:
            refresh_per_second: Refresh rate (default: 1.3 fps for smooth animation)

        Returns:
            Rich Live context manager
        """
        return Live(
            console=self._console,
            refresh_per_second=refresh_per_second,
            transient=False,
        )

    def close(self):
        """Close the log file if open."""
        if self._file_handle is not None:
            self._file_console.print("\n" + "=" * 120)
            self._file_console.print("End of convergence log")
            self._file_handle.close()
            self._file_handle = None
            self._file_console = None

    def __del__(self):
        """Cleanup when object is destroyed."""
        self.close()
