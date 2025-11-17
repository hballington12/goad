# python/goad/convergence/display.py
from typing import Dict, List

import numpy as np
from rich.columns import Columns
from rich.console import Console, Group
from rich.progress import BarColumn, Progress, TextColumn
from rich.spinner import Spinner
from rich.text import Text

from goad.convergence.base import Convergable
from goad.convergence.convergable import Tolerance


class ConvergenceDisplay:
    """Rich display for convergence tracking."""

    def __init__(self, targets: List[Convergable]):
        self.targets = targets
        self.console = Console()
        self._progress = None
        self._progress_tasks: Dict[str, int] = {}
        self._initialize_progress()

    def _initialize_progress(self):
        """Initialize Rich Progress bars for all targets."""
        self._progress = Progress(
            TextColumn("[bold]{task.fields[name]:<15}"),
            TextColumn("{task.fields[mean]:>10}"),
            BarColumn(bar_width=25),
            TextColumn("[bold]{task.percentage:>3.0f}%"),
            TextColumn("[cyan]{task.fields[sem_info]}"),
            console=self.console,
            transient=False,
        )

        for target in self.targets:
            task_id = self._progress.add_task(
                "",
                total=100,
                name=target.name,
                mean="--",
                sem_info="[SEM: -- / --]",
            )
            self._progress_tasks[target.name] = task_id

    def calculate_progress(self, target: Convergable) -> float:
        """Calculate progress percentage for a convergence target."""
        mean = target.mean
        sem = target.sem

        if mean is None or sem is None:
            return 0.0

        if target.tolerance == Tolerance.RELATIVE:
            if mean != 0:
                relative_sem = sem / abs(mean)
                if relative_sem > 0:
                    return min(100.0, np.sqrt(target.threshold / relative_sem) * 100.0)
        else:  # ABSOLUTE
            if sem > 0:
                return min(100.0, np.sqrt(target.threshold / sem) * 100.0)

        return 0.0

    def format_sem_info(self, target: Convergable) -> str:
        """Format SEM information for display."""
        mean = target.mean
        sem = target.sem

        if mean is None or sem is None:
            return "[SEM: -- / --]"

        if target.tolerance == Tolerance.RELATIVE:
            if mean != 0:
                relative_sem = sem / abs(mean)
                current_str = f"{relative_sem * 100:.2f}%"
            else:
                current_str = f"{sem:.4g}"
            target_str = f"{target.threshold * 100:.2f}%"
        else:
            current_str = f"{sem:.4g}"
            target_str = f"{target.threshold:.4g}"

        return f"[SEM: {current_str} / {target_str}]"

    def update_progress(self):
        """Update progress bars for all targets."""
        for target in self.targets:
            task_id = self._progress_tasks[target.name]

            progress = self.calculate_progress(target)
            sem_info = self.format_sem_info(target)
            mean = target.mean
            sem = target.sem
            if sem is None:
                sem = 0

            # Format mean value
            if mean is not None:
                mean_str = f"{mean:.4f} ± {sem:.4f}"
                if target.is_converged():
                    mean_str = f"[green]{mean_str}[/green]"
            else:
                mean_str = "--"

            self._progress.update(
                task_id, completed=progress, mean=mean_str, sem_info=sem_info
            )

    def build(
        self,
        iterations: int,
        batch_size: int,
        sim_time: float,
        min_orientations: int,
        max_orientations: int,
    ) -> Group:
        """Build the complete display."""
        # Update progress bars first
        self.update_progress()

        # Title with spinner
        spinner = Spinner("aesthetic", style="cyan")
        title_text = Text.assemble(
            ("GOAD: ", "bold cyan"),
            ("[Convergence] ", "bold white"),
        )
        title = Columns([title_text, spinner], expand=False, padding=(0, 1))

        # Status line
        orient_color = "green" if iterations >= min_orientations else "red"
        orient_text = Text(
            f"[Orientations: {iterations} (max {max_orientations})]", style=orient_color
        )

        # batch_text = Text(f"[Batch Size: {batch_size}]", style="yellow")

        if not np.isinf(sim_time):
            time_text = Text(f"[{sim_time:.3f} sec/orientation]", style="blue")
        else:
            time_text = Text("[Time/orientation: calculating...]", style="dim")

        min_status = (
            "✓"
            if iterations >= min_orientations
            else f"{iterations}/{min_orientations}"
        )
        min_text = Text(
            f"[Minimum orientations check: {min_status}]",
            style="green" if iterations >= min_orientations else "orange3",
        )

        header = Text.assemble(
            # orient_text, " ", batch_text, " ", time_text, " ", min_text
            orient_text,
            " ",
            time_text,
            " ",
            min_text,
        )

        separator = Text("━" * 80, style="dim")

        # Build display group
        return Group(
            separator,
            title,
            Text(""),
            header,
            separator,
            Text(""),
            self._progress.get_renderable(),
            Text(""),
            separator,
        )

    # def print_completion(self, iterations: int, converged: bool):
    #     """Print completion message."""
    #     if converged:
    #         self.console.print(
    #             f"[green]✓ Convergence achieved after {iterations} orientations.[/green]"
    #         )
    #         for target in self.targets:
    #             mean = target.mean
    #             sem = target.sem
    #             self.console.print(f"  • {target.name}: {mean:.6f} ± {sem:.6f}")
    #     else:
    #         self.console.print(
    #             "[yellow]INFO: Reached max iterations without full convergence.[/yellow]"
    #         )
