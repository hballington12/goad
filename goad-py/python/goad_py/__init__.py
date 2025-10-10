# Re-export everything from the compiled Rust module
from goad_py._goad_py import *

# Import Python modules
from .convergence import (
    Convergence,
    Convergable,
    ConvergenceResults,
    EnsembleConvergence,
)

from .phips_convergence import (
    PHIPSConvergence,
    PHIPSConvergable,
    PHIPSEnsembleConvergence,
)

# Import unified convergence API
from .unified_convergence import (
    # Main entry points
    run_convergence,
    run_convergence_sweep,
    # Classes
    UnifiedConvergence,
    UnifiedResults,
    ConvergenceConfig,
    # Modes
    ConvergenceMode,
    StandardMode,
    PHIPSMode,
)

__all__ = [
    # Legacy convergence API (still supported)
    "Convergence",
    "Convergable",
    "ConvergenceResults",
    "EnsembleConvergence",
    "PHIPSConvergence",
    "PHIPSConvergable",
    "PHIPSEnsembleConvergence",
    # Unified convergence API (recommended)
    "run_convergence",
    "run_convergence_sweep",
    "UnifiedConvergence",
    "UnifiedResults",
    "ConvergenceConfig",
    "ConvergenceMode",
    "StandardMode",
    "PHIPSMode",
]
