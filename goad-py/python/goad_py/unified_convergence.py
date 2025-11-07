"""
Unified convergence API for GOAD.

Provides a single entry point for all convergence types:
- Standard convergence (integrated parameters, Mueller elements)
- PHIPS detector convergence
- Single geometry or ensemble averaging

Features:
- Auto-detection of convergence mode
- Strict input validation
- Uniform output format (UnifiedResults)
- Support for parameter sweeps
- Full control over beam tracing, geometry transformations, and advanced optics
- Reproducible results via random seed control
"""

import json
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import toml

from . import _goad_py as goad
from .convergence import (
    Convergable,
    Convergence,
    ConvergenceResults,
    EnsembleConvergence,
)
from .phips_convergence import (
    PHIPSConvergable,
    PHIPSConvergence,
    PHIPSEnsembleConvergence,
)

# ============================================================================
# Convergence Modes
# ============================================================================


class ConvergenceMode(ABC):
    """Abstract base class for convergence modes."""

    @abstractmethod
    def validate_targets(self, targets: List[dict]) -> None:
        """Validate that targets are appropriate for this mode."""
        pass

    @abstractmethod
    def get_binning(self) -> goad.BinningScheme:
        """Get appropriate binning scheme for this mode."""
        pass

    @abstractmethod
    def get_mode_name(self) -> str:
        """Get mode name for logging/output."""
        pass


class StandardMode(ConvergenceMode):
    """Standard convergence mode for integrated parameters and Mueller elements."""

    # Valid targets for standard mode
    VALID_SCALAR_TARGETS = {"asymmetry", "scatt", "ext", "albedo"}
    VALID_MUELLER_TARGETS = {f"S{i}{j}" for i in range(1, 5) for j in range(1, 5)}

    def __init__(
        self,
        n_theta: int = 181,
        n_phi: int = 181,
        binning: Optional[goad.BinningScheme] = None,
    ):
        """
        Initialize standard mode.

        Args:
            n_theta: Number of theta bins for simple binning (default: 181 for 0-180 degrees)
            n_phi: Number of phi bins for simple binning (default: 181)
            binning: Optional custom BinningScheme (interval, custom, etc.)
                     If provided, n_theta and n_phi are ignored.
        """
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.custom_binning = binning

    def validate_targets(self, targets: List[dict]) -> None:
        """Validate that all targets are valid for standard mode."""
        valid_targets = self.VALID_SCALAR_TARGETS | self.VALID_MUELLER_TARGETS

        for target in targets:
            variable = target.get("variable")
            if variable is None:
                raise ValueError(f"Target missing 'variable' field: {target}")

            if variable not in valid_targets:
                raise ValueError(
                    f"Invalid target '{variable}' for StandardMode. "
                    f"Valid targets are: {sorted(valid_targets)}"
                )

            # Validate theta_indices only for Mueller elements
            if "theta_indices" in target:
                if variable not in self.VALID_MUELLER_TARGETS:
                    raise ValueError(
                        f"theta_indices can only be used with Mueller elements, "
                        f"but got variable='{variable}'"
                    )

    def get_binning(self) -> goad.BinningScheme:
        """Return binning scheme (custom if provided, otherwise simple)."""
        if self.custom_binning is not None:
            return self.custom_binning
        return goad.BinningScheme.simple(self.n_theta, self.n_phi)

    def get_mode_name(self) -> str:
        return "standard"


class PHIPSMode(ConvergenceMode):
    """PHIPS detector convergence mode."""

    def __init__(self, bins_file: str):
        """
        Initialize PHIPS mode.

        Args:
            bins_file: Path to PHIPS bins TOML file (e.g., phips_bins_edges.toml)

        Raises:
            FileNotFoundError: If bins_file doesn't exist
            ValueError: If bins_file is not a valid TOML file
        """
        self.bins_file = Path(bins_file)

        if not self.bins_file.exists():
            raise FileNotFoundError(
                f"PHIPS bins file not found: {self.bins_file}\n"
                f"Please generate it using phips_detector_angles_edges.py"
            )

        # Load and validate bins file
        try:
            with open(self.bins_file, "r") as f:
                bins_data = toml.load(f)

            if "bins" not in bins_data:
                raise ValueError(
                    f"PHIPS bins file missing 'bins' key: {self.bins_file}"
                )

            self.custom_bins = bins_data["bins"]

            if not isinstance(self.custom_bins, list) or len(self.custom_bins) == 0:
                raise ValueError(
                    f"PHIPS bins file has invalid 'bins' data: {self.bins_file}"
                )

        except Exception as e:
            raise ValueError(f"Failed to load PHIPS bins file {self.bins_file}: {e}")

    def validate_targets(self, targets: List[dict]) -> None:
        """
        Validate that targets are appropriate for PHIPS mode.

        PHIPS mode only supports convergence on DSCS values at detector bins.
        The 'variable' field is not used (always 'phips_dscs').
        """
        if len(targets) != 1:
            raise ValueError(
                f"PHIPSMode only supports a single convergence target, got {len(targets)}. "
                f"Use tolerance and tolerance_type to specify convergence criteria."
            )

        target = targets[0]

        # PHIPS doesn't use 'variable' field - it's always DSCS
        # But validate other fields
        tolerance_type = target.get("tolerance_type", "relative")
        if tolerance_type not in {"relative", "absolute"}:
            raise ValueError(
                f"Invalid tolerance_type '{tolerance_type}'. Must be 'relative' or 'absolute'."
            )

        tolerance = target.get("tolerance")
        if tolerance is None:
            raise ValueError("Target missing 'tolerance' field")

        if tolerance <= 0:
            raise ValueError(f"tolerance must be positive, got {tolerance}")

        # Validate detector_indices if present
        if "detector_indices" in target:
            detector_indices = target["detector_indices"]
            if not isinstance(detector_indices, list):
                raise ValueError("detector_indices must be a list of integers")
            if not all(isinstance(i, int) and 0 <= i < 20 for i in detector_indices):
                raise ValueError("detector_indices must be integers in range [0, 19]")

    def get_binning(self) -> goad.BinningScheme:
        """Return custom binning scheme from PHIPS bins file."""
        return goad.BinningScheme.custom(self.custom_bins)

    def get_mode_name(self) -> str:
        return "phips"


# ============================================================================
# Configuration
# ============================================================================


@dataclass
class BeamTracingConfig:
    """Beam tracing performance and accuracy parameters."""

    beam_power_threshold: float = 0.005
    beam_area_threshold_fac: float = 0.001
    cutoff: float = 0.999
    max_rec: int = 10
    max_tir: int = 10

    def __post_init__(self):
        """Validate beam tracing parameters."""
        if self.beam_power_threshold <= 0 or self.beam_power_threshold > 1:
            raise ValueError(
                f"beam_power_threshold must be in range (0, 1], got {self.beam_power_threshold}"
            )

        if self.beam_area_threshold_fac <= 0:
            raise ValueError(
                f"beam_area_threshold_fac must be positive, got {self.beam_area_threshold_fac}"
            )

        if self.cutoff < 0 or self.cutoff > 1:
            raise ValueError(f"cutoff must be between 0 and 1, got {self.cutoff}")

        if self.max_rec < 0:
            raise ValueError(f"max_rec must be non-negative, got {self.max_rec}")

        if self.max_tir < 0:
            raise ValueError(f"max_tir must be non-negative, got {self.max_tir}")


@dataclass
class GeometryTransformConfig:
    """Geometry transformation parameters.

    Attributes:
        scale: Problem scaling factor - scales the entire problem including geometry,
               wavelength, and beam area thresholds (default: 1.0)
        distortion: Geometry distortion factor (optional)
        geom_scale: Per-axis geometry scaling [x, y, z] - scales only the geometry
                    in each dimension independently (optional)
    """

    scale: float = 1.0
    distortion: Optional[float] = None
    geom_scale: Optional[List[float]] = None

    def __post_init__(self):
        """Validate geometry transformations."""
        if self.scale <= 0:
            raise ValueError(f"scale must be positive, got {self.scale}")

        if self.geom_scale is not None:
            if len(self.geom_scale) != 3:
                raise ValueError(
                    f"geom_scale must have exactly 3 values [x, y, z], "
                    f"got {len(self.geom_scale)}"
                )
            if any(s <= 0 for s in self.geom_scale):
                raise ValueError("All geom_scale values must be positive")


@dataclass
class AdvancedConfig:
    """Advanced optical calculation parameters."""

    mapping: Optional[Any] = (
        "ApertureDiffraction"  # String that will be converted to enum in __post_init__
    )
    coherence: bool = True
    fov_factor: Optional[float] = None

    def __post_init__(self):
        """Validate advanced optics settings and convert mapping string to enum."""
        # Convert string mapping to enum if needed
        if isinstance(self.mapping, str):
            if self.mapping == "ApertureDiffraction":
                self.mapping = goad.Mapping.ApertureDiffraction
            elif self.mapping == "GeometricOptics":
                self.mapping = goad.Mapping.GeometricOptics
            else:
                raise ValueError(
                    f"Invalid mapping '{self.mapping}'. Must be 'ApertureDiffraction' or 'GeometricOptics'"
                )

        if self.fov_factor is not None and self.fov_factor <= 0:
            raise ValueError(f"fov_factor must be positive, got {self.fov_factor}")


@dataclass
class ConvergenceConfig:
    """
    Unified configuration for all convergence types.

    This class provides a single configuration interface that works for:
    - Standard convergence (integrated parameters, Mueller elements)
    - PHIPS detector convergence
    - Single geometry or ensemble averaging

    Attributes:
        geometry: Path to .obj file or directory of .obj files (ensemble)
        mode: ConvergenceMode instance (StandardMode or PHIPSMode)
        convergence_targets: List of convergence target dicts

        wavelength: Wavelength in microns (default: 0.532)
        particle_refr_index_re: Real part of particle refractive index (default: 1.31)
        particle_refr_index_im: Imaginary part of particle refractive index (default: 0.0)
        medium_refr_index_re: Real part of medium refractive index (default: 1.0)
        medium_refr_index_im: Imaginary part of medium refractive index (default: 0.0)

        batch_size: Orientations per batch (default: 24)
        max_orientations: Maximum orientations (default: 100,000)
        min_batches: Minimum batches before convergence check (default: 10)

        beam_tracing: BeamTracingConfig instance for beam tracing parameters
        geometry_transform: GeometryTransformConfig instance for geometry transformations
        advanced_config: AdvancedConfig instance for advanced optical parameters
        seed: Random seed for reproducibility (optional)

        mueller_1d: Compute 1D Mueller matrix (default: True, standard mode only)
        output_dir: Output directory path (optional)
    """

    # Required fields
    geometry: Union[str, Path]
    mode: ConvergenceMode
    convergence_targets: List[dict]

    # Optical settings
    wavelength: float = 0.532
    particle_refr_index_re: float = 1.31
    particle_refr_index_im: float = 0.0
    medium_refr_index_re: float = 1.0
    medium_refr_index_im: float = 0.0

    # Convergence parameters
    batch_size: int = 24
    max_orientations: int = 100_000
    min_batches: int = 10

    # Beam tracing configuration
    beam_tracing: BeamTracingConfig = field(default_factory=BeamTracingConfig)

    # Geometry transformations
    geometry_transform: GeometryTransformConfig = field(
        default_factory=GeometryTransformConfig
    )

    # Advanced optics
    advanced_config: AdvancedConfig = field(default_factory=AdvancedConfig)

    # Random seed for reproducibility
    seed: Optional[int] = None

    # Mueller matrix output (only for StandardMode)
    mueller_1d: bool = True

    # Output options
    output_dir: Optional[str] = None
    log_file: Optional[str] = None

    def __post_init__(self):
        """Validate configuration after initialization."""
        # Convert geometry to Path
        self.geometry = Path(self.geometry)

        # Validate geometry exists
        if not self.geometry.exists():
            raise FileNotFoundError(f"Geometry path does not exist: {self.geometry}")

        # Validate mode is a ConvergenceMode instance
        if not isinstance(self.mode, ConvergenceMode):
            raise TypeError(
                f"mode must be a ConvergenceMode instance (StandardMode or PHIPSMode), "
                f"got {type(self.mode)}"
            )

        # Validate convergence targets
        if not isinstance(self.convergence_targets, list):
            raise TypeError("convergence_targets must be a list of dicts")

        if len(self.convergence_targets) == 0:
            raise ValueError("convergence_targets cannot be empty")

        # Let mode validate its targets
        self.mode.validate_targets(self.convergence_targets)

        # Validate numeric parameters
        if self.wavelength <= 0:
            raise ValueError(f"wavelength must be positive, got {self.wavelength}")

        if self.batch_size <= 0:
            raise ValueError(f"batch_size must be positive, got {self.batch_size}")

        if self.max_orientations <= 0:
            raise ValueError(
                f"max_orientations must be positive, got {self.max_orientations}"
            )

        if self.min_batches <= 0:
            raise ValueError(f"min_batches must be positive, got {self.min_batches}")

        # Validate Mueller options for PHIPS mode
        if isinstance(self.mode, PHIPSMode):
            if self.mueller_1d:
                raise ValueError(
                    "mueller_1d is not supported in PHIPSMode. "
                    "PHIPS custom binning does not compute integrated Mueller matrices."
                )

    def is_ensemble(self) -> bool:
        """Check if geometry is a directory (ensemble mode)."""
        return self.geometry.is_dir()

    def to_dict(self) -> dict:
        """Convert config to dictionary (for serialization)."""
        return {
            "geometry": str(self.geometry),
            "mode": self.mode.get_mode_name(),
            "convergence_targets": self.convergence_targets,
            "wavelength": self.wavelength,
            "particle_refr_index_re": self.particle_refr_index_re,
            "particle_refr_index_im": self.particle_refr_index_im,
            "medium_refr_index_re": self.medium_refr_index_re,
            "medium_refr_index_im": self.medium_refr_index_im,
            "batch_size": self.batch_size,
            "max_orientations": self.max_orientations,
            "min_batches": self.min_batches,
            "beam_tracing": {
                "beam_power_threshold": self.beam_tracing.beam_power_threshold,
                "beam_area_threshold_fac": self.beam_tracing.beam_area_threshold_fac,
                "cutoff": self.beam_tracing.cutoff,
                "max_rec": self.beam_tracing.max_rec,
                "max_tir": self.beam_tracing.max_tir,
            },
            "geometry_transform": {
                "scale": self.geometry_transform.scale,
                "distortion": self.geometry_transform.distortion,
                "geom_scale": self.geometry_transform.geom_scale,
            },
            "advanced_config": {
                "mapping": str(self.advanced_config.mapping),
                "coherence": self.advanced_config.coherence,
                "fov_factor": self.advanced_config.fov_factor,
            },
            "seed": self.seed,
            "mueller_1d": self.mueller_1d,
            "is_ensemble": self.is_ensemble(),
        }


# ============================================================================
# Results
# ============================================================================


@dataclass
class UnifiedResults:
    """
    Unified results container for all convergence types.

    Provides a consistent interface regardless of convergence mode,
    with JSON serialization support.
    """

    converged: bool
    n_orientations: int
    mode: str  # "standard" or "phips"
    is_ensemble: bool

    # Values and SEMs (format depends on mode)
    # For standard: Dict[str, float] or Dict[str, np.ndarray] for Mueller
    # For PHIPS: Dict['phips_dscs', np.ndarray] with shape (20,)
    values: Dict[str, Union[float, np.ndarray]]
    sem_values: Dict[str, Union[float, np.ndarray]]

    # Optional Mueller matrices (standard mode only)
    mueller_1d: Optional[np.ndarray] = None

    # Bins and detector info
    bins_1d: Optional[np.ndarray] = None  # Theta bins for standard mode
    detector_angles: Optional[np.ndarray] = None  # Detector angles for PHIPS mode

    # Metadata
    convergence_history: Optional[List[Tuple[int, str, float]]] = None
    warning: Optional[str] = None
    config: Optional[ConvergenceConfig] = None

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = []
        lines.append("=" * 60)
        lines.append("Convergence Results Summary")
        lines.append("=" * 60)
        lines.append(f"Mode: {self.mode}")
        lines.append(f"Ensemble: {self.is_ensemble}")
        lines.append(f"Converged: {self.converged}")
        lines.append(f"Total orientations: {self.n_orientations}")

        if self.warning:
            lines.append(f"\nWarning: {self.warning}")

        lines.append("\nConverged Values:")
        for var, value in self.values.items():
            # Skip mueller_1d_sem (it's metadata, not a result to display)
            if var == "mueller_1d_sem":
                continue

            sem = self.sem_values.get(var)

            if sem is None:
                # Skip if no SEM available
                continue

            if isinstance(value, np.ndarray):
                # Array value (Mueller element or PHIPS DSCS)
                lines.append(f"  {var}: array with shape {value.shape}")
                lines.append(f"    Mean: {np.nanmean(value):.6e}")
                lines.append(f"    Std:  {np.nanstd(value):.6e}")
                lines.append(
                    f"    SEM range: [{np.nanmin(sem):.6e}, {np.nanmax(sem):.6e}]"
                )
            else:
                # Scalar value
                lines.append(f"  {var}: {value:.6e} ± {sem:.6e}")

        lines.append("=" * 60)
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """
        Convert results to dictionary for serialization.

        Outputs all values as [mean, sem] tuples.
        Mode-specific outputs:
        - Standard mode: integrated parameters + mueller_1d + bins_1d
        - PHIPS mode: all 20 DSCS values + detector_angles (nulls for integrated params/mueller)
        """
        result = {
            "converged": bool(self.converged),  # Convert numpy bool to Python bool
            "n_orientations": int(
                self.n_orientations
            ),  # Convert numpy int to Python int
            "mode": self.mode,
            "is_ensemble": bool(self.is_ensemble),  # Convert numpy bool to Python bool
            "warning": self.warning,
        }

        # Mode-specific output
        if self.mode == "standard":
            # Standard mode: always output all integrated parameters + Mueller

            # Integrated parameters (always output, even if not converged on)
            for param in ["asymmetry", "scatt", "ext", "albedo"]:
                mean_val = self.values.get(param)
                sem_val = self.sem_values.get(param)

                if mean_val is not None and sem_val is not None:
                    # Convert to Python floats for JSON serialization
                    result[param] = [float(mean_val), float(sem_val)]
                else:
                    result[param] = None

            # Mueller 1D (if available)
            if self.mueller_1d is not None:
                # mueller_1d is shape (n_theta, 16)
                # Convert to list of [mean, sem] tuples per angle per element
                mueller_mean = self.mueller_1d
                mueller_sem = self.values.get("mueller_1d_sem")

                if mueller_sem is not None:
                    # Create [mean, sem] format for each position
                    mueller_output = []
                    for i in range(mueller_mean.shape[0]):
                        angle_data = []
                        for j in range(mueller_mean.shape[1]):
                            # Convert to Python floats for JSON serialization
                            angle_data.append(
                                [float(mueller_mean[i, j]), float(mueller_sem[i, j])]
                            )
                        mueller_output.append(angle_data)
                    result["mueller_1d"] = mueller_output
                else:
                    # No SEM available, just output mean
                    result["mueller_1d"] = mueller_mean.tolist()
            else:
                result["mueller_1d"] = None

            # Bins (if available)
            if self.bins_1d is not None:
                result["bins_1d"] = self.bins_1d.tolist()
            else:
                result["bins_1d"] = None

            # PHIPS fields are null in standard mode
            result["phips_dscs"] = None
            result["detector_angles"] = None

        elif self.mode == "phips":
            # PHIPS mode: all 20 DSCS values + detector angles

            # Integrated parameters are null in PHIPS mode
            result["asymmetry"] = None
            result["scatt"] = None
            result["ext"] = None
            result["albedo"] = None
            result["mueller_1d"] = None
            result["bins_1d"] = None

            # PHIPS DSCS values (always all 20, with [mean, sem] format)
            phips_mean = self.values.get("phips_dscs")
            phips_sem = self.sem_values.get("phips_dscs")

            if phips_mean is not None and phips_sem is not None:
                # Shape should be (20,)
                phips_output = []
                for i in range(len(phips_mean)):
                    # Convert numpy scalars to Python floats for JSON serialization
                    phips_output.append([float(phips_mean[i]), float(phips_sem[i])])
                result["phips_dscs"] = phips_output
            else:
                result["phips_dscs"] = None

            # Detector angles
            if self.detector_angles is not None:
                result["detector_angles"] = self.detector_angles.tolist()
            else:
                result["detector_angles"] = None

        # Add convergence history
        if self.convergence_history:
            result["convergence_history"] = [
                {"n_orientations": n, "variable": var, "sem": sem}
                for n, var, sem in self.convergence_history
            ]

        # Add config if present
        if self.config:
            result["config"] = self.config.to_dict()

        return result

    def save(self, path: str) -> None:
        """
        Save results to JSON file.

        Args:
            path: Output file path (will append .json if not present)
        """
        path = Path(path)
        if path.suffix != ".json":
            path = path.with_suffix(".json")

        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

        print(f"Results saved to: {path}")

    @classmethod
    def load(cls, path: str) -> "UnifiedResults":
        """
        Load results from JSON file.

        Args:
            path: Input file path

        Returns:
            UnifiedResults instance
        """
        with open(path, "r") as f:
            data = json.load(f)

        # New tuple format: each parameter is [mean, sem]
        # Extract values and sems from the new format
        values = {}
        sem_values = {}

        mode = data["mode"]

        if mode == "standard":
            # Standard mode: extract integrated parameters
            for param in ["asymmetry", "scatt", "ext", "albedo"]:
                param_data = data.get(param)
                if param_data is not None:
                    values[param] = param_data[0]  # mean
                    sem_values[param] = param_data[1]  # sem

            # Mueller 1D (if present)
            mueller_1d_data = data.get("mueller_1d")
            if mueller_1d_data is not None and len(mueller_1d_data) > 0:
                # mueller_1d is stored as list of [mean, sem] tuples per angle per element
                # Need to separate into mean array and sem array
                # But actually, we need to check if it's in the new format or old format
                if isinstance(mueller_1d_data[0][0], list):
                    # New format: [[mean, sem], [mean, sem], ...]
                    n_theta = len(mueller_1d_data)
                    n_elements = len(mueller_1d_data[0])
                    mueller_mean = np.zeros((n_theta, n_elements))
                    mueller_sem = np.zeros((n_theta, n_elements))
                    for i in range(n_theta):
                        for j in range(n_elements):
                            mueller_mean[i, j] = mueller_1d_data[i][j][0]
                            mueller_sem[i, j] = mueller_1d_data[i][j][1]
                    mueller_1d = mueller_mean
                    values["mueller_1d_sem"] = mueller_sem
                else:
                    # Old format: just the mean array
                    mueller_1d = np.array(mueller_1d_data)
            else:
                mueller_1d = None

            # Bins
            bins_1d = data.get("bins_1d")
            if bins_1d is not None:
                bins_1d = np.array(bins_1d)

            detector_angles = None

        elif mode == "phips":
            # PHIPS mode: extract DSCS values
            phips_data = data.get("phips_dscs")
            if phips_data is not None:
                # phips_dscs is list of [mean, sem] tuples (20 detectors)
                phips_mean = np.array([d[0] for d in phips_data])
                phips_sem = np.array([d[1] for d in phips_data])
                values["phips_dscs"] = phips_mean
                sem_values["phips_dscs"] = phips_sem

            # Detector angles
            detector_angles = data.get("detector_angles")
            if detector_angles is not None:
                detector_angles = np.array(detector_angles)

            mueller_1d = None
            bins_1d = None

        # Convert convergence history
        convergence_history = None
        if "convergence_history" in data and data["convergence_history"] is not None:
            convergence_history = [
                (h["n_orientations"], h["variable"], h["sem"])
                for h in data["convergence_history"]
            ]

        return cls(
            converged=data["converged"],
            n_orientations=data["n_orientations"],
            mode=data["mode"],
            is_ensemble=data["is_ensemble"],
            values=values,
            sem_values=sem_values,
            mueller_1d=mueller_1d,
            bins_1d=bins_1d,
            detector_angles=detector_angles,
            convergence_history=convergence_history,
            warning=data.get("warning"),
            config=None,  # Config not reconstructed from JSON
        )


# ============================================================================
# Unified Convergence
# ============================================================================


class UnifiedConvergence:
    """
    Unified convergence runner that handles all modes.

    Auto-selects appropriate convergence class based on configuration.
    """

    def __init__(self, config: ConvergenceConfig):
        """
        Initialize unified convergence.

        Args:
            config: ConvergenceConfig instance
        """
        self.config = config
        self._convergence = None
        self._setup()

    def _setup(self):
        """Create appropriate convergence instance based on config."""
        # Determine geometry path for Settings initialization
        # For ensemble mode, use the first .obj file found in the directory
        if self.config.is_ensemble():
            geom_path = Path(self.config.geometry)
            obj_files = sorted(geom_path.glob("*.obj"))
            if not obj_files:
                raise ValueError(
                    f"No .obj files found in directory: {self.config.geometry}"
                )
            geom_path_str = str(obj_files[0])
        else:
            geom_path_str = str(self.config.geometry)

        # Create GOAD settings with all parameters
        settings = goad.Settings(
            geom_path=geom_path_str,
            wavelength=self.config.wavelength,
            particle_refr_index_re=self.config.particle_refr_index_re,
            particle_refr_index_im=self.config.particle_refr_index_im,
            medium_refr_index_re=self.config.medium_refr_index_re,
            medium_refr_index_im=self.config.medium_refr_index_im,
            binning=self.config.mode.get_binning(),
            # Beam tracing parameters
            beam_power_threshold=self.config.beam_tracing.beam_power_threshold,
            beam_area_threshold_fac=self.config.beam_tracing.beam_area_threshold_fac,
            cutoff=self.config.beam_tracing.cutoff,
            max_rec=self.config.beam_tracing.max_rec,
            max_tir=self.config.beam_tracing.max_tir,
            # Geometry transformations
            scale=self.config.geometry_transform.scale,
            # Advanced configuration
            mapping=self.config.advanced_config.mapping,
            coherence=self.config.advanced_config.coherence,
        )

        # Set optional parameters if provided
        if self.config.seed is not None:
            settings.seed = self.config.seed

        if self.config.geometry_transform.distortion is not None:
            settings.distortion = self.config.geometry_transform.distortion

        if self.config.geometry_transform.geom_scale is not None:
            settings.geom_scale = self.config.geometry_transform.geom_scale

        if self.config.advanced_config.fov_factor is not None:
            settings.fov_factor = self.config.advanced_config.fov_factor

        # Create convergence instance based on mode
        if isinstance(self.config.mode, StandardMode):
            self._setup_standard(settings)
        elif isinstance(self.config.mode, PHIPSMode):
            self._setup_phips(settings)
        else:
            raise ValueError(f"Unknown mode type: {type(self.config.mode)}")

    def _setup_standard(self, settings):
        """Setup standard convergence."""
        # Convert target dicts to Convergable instances
        convergables = []
        for target_dict in self.config.convergence_targets:
            convergables.append(Convergable(**target_dict))

        # Select ensemble or single-geometry convergence
        if self.config.is_ensemble():
            self._convergence = EnsembleConvergence(
                settings=settings,
                convergables=convergables,
                geom_dir=str(self.config.geometry),
                batch_size=self.config.batch_size,
                max_orientations=self.config.max_orientations,
                min_batches=self.config.min_batches,
                mueller_1d=self.config.mueller_1d,
                mueller_2d=False,
                log_file=self.config.log_file,
            )
        else:
            self._convergence = Convergence(
                settings=settings,
                convergables=convergables,
                batch_size=self.config.batch_size,
                max_orientations=self.config.max_orientations,
                min_batches=self.config.min_batches,
                mueller_1d=self.config.mueller_1d,
                mueller_2d=False,
                log_file=self.config.log_file,
            )

    def _setup_phips(self, settings):
        """Setup PHIPS convergence."""
        # PHIPS only has one convergable target
        target_dict = self.config.convergence_targets[0]

        # Remove 'variable' if present (not used in PHIPS)
        target_dict = {k: v for k, v in target_dict.items() if k != "variable"}

        convergable = PHIPSConvergable(**target_dict)

        # Select ensemble or single-geometry convergence
        if self.config.is_ensemble():
            self._convergence = PHIPSEnsembleConvergence(
                settings=settings,
                convergable=convergable,
                geom_dir=str(self.config.geometry),
                batch_size=self.config.batch_size,
                max_orientations=self.config.max_orientations,
                min_batches=self.config.min_batches,
                log_file=self.config.log_file,
            )
        else:
            self._convergence = PHIPSConvergence(
                settings=settings,
                convergable=convergable,
                batch_size=self.config.batch_size,
                max_orientations=self.config.max_orientations,
                min_batches=self.config.min_batches,
                log_file=self.config.log_file,
            )

    def run(self) -> UnifiedResults:
        """
        Run convergence study.

        Returns:
            UnifiedResults instance with ALL parameters extracted (not just converged ones)
        """
        # Run convergence
        results = self._convergence.run()

        # Extract ALL integrated parameters (not just converged ones)
        all_values = {}
        all_sems = {}

        if isinstance(self.config.mode, StandardMode):
            # Standard mode: extract all 4 integrated parameters
            for param in ["asymmetry", "scatt", "ext", "albedo"]:
                if param in results.values:
                    # Already computed (was a convergence target)
                    all_values[param] = results.values[param]
                    all_sems[param] = results.sem_values[param]
                else:
                    # Not a convergence target - compute it now
                    mean, sem = self._convergence._calculate_mean_and_sem(param)
                    all_values[param] = mean
                    all_sems[param] = sem

            # Also include any Mueller elements that were converged on
            for key in results.values:
                if key.startswith("S"):  # Mueller element
                    all_values[key] = results.values[key]
                    all_sems[key] = results.sem_values[key]

            # Include mueller_1d_sem if available (for full SEM output)
            if "mueller_1d_sem" in results.values:
                all_values["mueller_1d_sem"] = results.values["mueller_1d_sem"]

        elif isinstance(self.config.mode, PHIPSMode):
            # PHIPS mode: extract all 20 DSCS values
            if "phips_dscs" in results.values:
                all_values["phips_dscs"] = results.values["phips_dscs"]
                all_sems["phips_dscs"] = results.sem_values["phips_dscs"]
            else:
                # This shouldn't happen in PHIPS mode, but handle it
                all_values["phips_dscs"] = None
                all_sems["phips_dscs"] = None

        # Get bins_1d or detector_angles based on mode
        bins_1d = None
        detector_angles = None

        if isinstance(self.config.mode, StandardMode):
            # Extract theta bins from Mueller 1D if available
            if results.mueller_1d is not None:
                # bins_1d should be theta values
                # For simple binning: np.linspace(0, 180, n_theta)
                # For custom binning: we need to extract from the binning scheme
                binning = self.config.mode.get_binning()
                if hasattr(binning, "get_theta_bins"):
                    bins_1d = binning.get_theta_bins()
                else:
                    # Fallback: assume uniform spacing
                    n_theta = results.mueller_1d.shape[0]
                    bins_1d = np.linspace(0, 180, n_theta)

        elif isinstance(self.config.mode, PHIPSMode):
            # PHIPS detector angles (18° to 170°, 20 detectors)
            detector_angles = np.array(
                [
                    18.0,
                    26.0,
                    34.0,
                    42.0,
                    50.0,
                    58.0,
                    66.0,
                    74.0,
                    82.0,
                    90.0,
                    98.0,
                    106.0,
                    114.0,
                    122.0,
                    130.0,
                    138.0,
                    146.0,
                    154.0,
                    162.0,
                    170.0,
                ]
            )

        # Convert to UnifiedResults
        unified = UnifiedResults(
            converged=results.converged,
            n_orientations=results.n_orientations,
            mode=self.config.mode.get_mode_name(),
            is_ensemble=self.config.is_ensemble(),
            values=all_values,
            sem_values=all_sems,
            mueller_1d=results.mueller_1d,
            bins_1d=bins_1d,
            detector_angles=detector_angles,
            convergence_history=results.convergence_history,
            warning=results.warning,
            config=self.config,
        )

        return unified


# ============================================================================
# Convenience Functions
# ============================================================================


def run_convergence(
    geometry: Union[str, Path],
    targets: Union[str, List[str], List[dict]],
    tolerance: float = 0.05,
    tolerance_type: str = "relative",
    mode: Union[str, ConvergenceMode] = "auto",
    **kwargs,
) -> UnifiedResults:
    """
    Run convergence study with unified interface.

    This is the primary entry point for most users.

    Args:
        geometry: Path to .obj file or directory of .obj files
        targets: What to converge on:
            - Single string: "asymmetry", "scatt", "S11", "phips_dscs"
            - List of strings: ["asymmetry", "scatt"]
            - List of dicts: [{"variable": "S11", "tolerance": 0.1, ...}]
        tolerance: Default tolerance (can be overridden per target)
        tolerance_type: "relative" or "absolute"
        mode: ConvergenceMode instance, or string "auto"/"standard"/"phips"
        **kwargs: Additional settings:
            # Optical settings
            - wavelength: Wavelength in microns (default: 0.532)
            - particle_refr_index_re: Real part of particle refractive index
            - particle_refr_index_im: Imaginary part of particle refractive index
            - medium_refr_index_re: Real part of medium refractive index
            - medium_refr_index_im: Imaginary part of medium refractive index

            # Convergence parameters
            - batch_size: Orientations per batch (default: 24)
            - max_orientations: Maximum orientations (default: 100,000)
            - min_batches: Minimum batches before convergence (default: 10)
            - mueller_1d: Compute 1D Mueller matrix (default: True, standard mode only)

            # Mode settings
            - phips_bins_file: Path to PHIPS bins TOML (required if mode="phips")
            - n_theta: Number of theta bins for standard mode (default: 181)
            - n_phi: Number of phi bins for standard mode (default: 181)

            # Beam tracing parameters
            - beam_power_threshold: Beam power threshold (default: 0.005)
            - beam_area_threshold_fac: Beam area threshold factor (default: 0.001)
            - cutoff: Ray power cutoff (default: 0.999)
            - max_rec: Max recursion depth (default: 10)
            - max_tir: Max TIR bounces (default: 10)

            # Geometry transformations
            - scale: Problem scaling factor - scales entire problem including geometry,
                     wavelength, and beam area thresholds (default: 1.0)
            - distortion: Geometry distortion factor (optional)
            - geom_scale: Per-axis geometry scaling [x, y, z] - scales only geometry
                          in each dimension independently (optional)

            # Advanced configuration
            - mapping: DSCS mapping scheme (default: goad.Mapping.ApertureDiffraction)
            - coherence: Enable coherent scattering (default: True)
            - fov_factor: Field of view factor (optional)

            # Reproducibility
            - seed: Random seed for orientations (optional)

            # Output options
            - log_file: Path to log file for convergence progress (optional)

    Returns:
        UnifiedResults object

    Examples:
        # Simple: converge asymmetry for single geometry
        results = run_convergence("hex.obj", "asymmetry", tolerance=0.01)

        # Multiple targets, ensemble
        results = run_convergence(
            "./test_obj",
            ["asymmetry", "scatt"],
            tolerance=0.05,
            batch_size=48
        )

        # PHIPS detectors
        results = run_convergence(
            "./test_obj",
            "phips_dscs",
            tolerance=0.25,
            mode="phips",
            phips_bins_file="phips_bins_edges.toml"
        )

        # Mueller elements with specific bins
        results = run_convergence(
            "hex.obj",
            [{"variable": "S11", "tolerance": 0.1, "theta_indices": [180]}],
            batch_size=12
        )

        # Advanced: custom beam tracing and geometry scaling
        results = run_convergence(
            "complex_particle.obj",
            "asymmetry",
            max_rec=200,
            max_tir=150,
            cutoff=0.0001,
            scale=2.0,
            seed=42
        )
    """
    # Extract beam tracing parameters from kwargs
    beam_tracing = BeamTracingConfig(
        beam_power_threshold=kwargs.pop("beam_power_threshold", 0.005),
        beam_area_threshold_fac=kwargs.pop("beam_area_threshold_fac", 0.001),
        cutoff=kwargs.pop("cutoff", 0.999),
        max_rec=kwargs.pop("max_rec", 10),
        max_tir=kwargs.pop("max_tir", 10),
    )

    # Extract geometry transform parameters
    geometry_transform = GeometryTransformConfig(
        scale=kwargs.pop("scale", 1.0),
        distortion=kwargs.pop("distortion", None),
        geom_scale=kwargs.pop("geom_scale", None),
    )

    # Extract advanced configuration parameters
    advanced_config = AdvancedConfig(
        mapping=kwargs.pop("mapping", "ApertureDiffraction"),
        coherence=kwargs.pop("coherence", True),
        fov_factor=kwargs.pop("fov_factor", None),
    )

    # Extract seed and log_file
    seed = kwargs.pop("seed", None)
    log_file = kwargs.pop("log_file", None)

    # Normalize targets to list of dicts
    target_dicts = _normalize_targets(targets, tolerance, tolerance_type)

    # Auto-detect or create mode
    if isinstance(mode, str):
        if mode == "auto":
            mode = _auto_detect_mode(target_dicts, kwargs)
        elif mode == "standard":
            n_theta = kwargs.pop("n_theta", 181)
            n_phi = kwargs.pop("n_phi", 181)
            mode = StandardMode(n_theta=n_theta, n_phi=n_phi)
        elif mode == "phips":
            phips_bins_file = kwargs.pop("phips_bins_file", None)
            if phips_bins_file is None:
                raise ValueError("phips_bins_file must be specified when mode='phips'")
            mode = PHIPSMode(bins_file=phips_bins_file)
        else:
            raise ValueError(
                f"Invalid mode string '{mode}'. Must be 'auto', 'standard', or 'phips'"
            )

    # Build config with new parameters
    config = ConvergenceConfig(
        geometry=geometry,
        mode=mode,
        convergence_targets=target_dicts,
        beam_tracing=beam_tracing,
        geometry_transform=geometry_transform,
        advanced_config=advanced_config,
        seed=seed,
        log_file=log_file,
        **kwargs,  # Remaining kwargs (wavelength, particle_refr_index_re, etc.)
    )

    # Run
    conv = UnifiedConvergence(config)
    return conv.run()


def run_convergence_sweep(
    configs: List[ConvergenceConfig], parallel: bool = False, n_jobs: int = -1
) -> List[UnifiedResults]:
    """
    Run multiple convergence studies (for parameter sweeps).

    Args:
        configs: List of ConvergenceConfig instances
        parallel: Run in parallel if True (not yet implemented)
        n_jobs: Number of parallel jobs (-1 = all cores)

    Returns:
        List of UnifiedResults

    Example:
        # Wavelength sweep
        wavelengths = [0.532, 0.633, 0.780]
        configs = [
            ConvergenceConfig(
                geometry="./test_obj",
                mode=PHIPSMode("phips_bins_edges.toml"),
                convergence_targets=[{"tolerance": 0.25, "tolerance_type": "relative"}],
                wavelength=wl
            )
            for wl in wavelengths
        ]
        results = run_convergence_sweep(configs)
    """
    if parallel:
        raise NotImplementedError("Parallel execution not yet implemented")

    results = []
    for i, config in enumerate(configs):
        print(f"\n{'=' * 60}")
        print(f"Running convergence {i + 1}/{len(configs)}")
        print(f"{'=' * 60}")

        conv = UnifiedConvergence(config)
        result = conv.run()
        results.append(result)

    return results


# ============================================================================
# Helper Functions
# ============================================================================


def _normalize_targets(
    targets: Union[str, List[str], List[dict]],
    default_tolerance: float,
    default_tolerance_type: str,
) -> List[dict]:
    """
    Normalize targets to list of dicts.

    Args:
        targets: Target specification (string, list of strings, or list of dicts)
        default_tolerance: Default tolerance value
        default_tolerance_type: Default tolerance type

    Returns:
        List of target dicts
    """
    # Convert single string to list
    if isinstance(targets, str):
        targets = [targets]

    # Convert to list of dicts
    target_dicts = []
    for target in targets:
        if isinstance(target, str):
            # String target - apply defaults
            target_dict = {
                "variable": target,
                "tolerance": default_tolerance,
                "tolerance_type": default_tolerance_type,
            }
            target_dicts.append(target_dict)
        elif isinstance(target, dict):
            # Dict target - fill in missing defaults
            target_dict = {
                "tolerance": default_tolerance,
                "tolerance_type": default_tolerance_type,
                **target,  # Override with provided values
            }
            target_dicts.append(target_dict)
        else:
            raise TypeError(
                f"Invalid target type: {type(target)}. Must be str or dict."
            )

    return target_dicts


def _auto_detect_mode(target_dicts: List[dict], kwargs: dict) -> ConvergenceMode:
    """
    Auto-detect convergence mode from targets.

    Args:
        target_dicts: List of target dictionaries
        kwargs: Additional keyword arguments

    Returns:
        ConvergenceMode instance
    """
    # Check if any target is "phips_dscs" or if phips_bins_file is specified
    has_phips_target = any(
        target.get("variable") == "phips_dscs" for target in target_dicts
    )
    has_phips_bins = "phips_bins_file" in kwargs

    if has_phips_target or has_phips_bins:
        # PHIPS mode
        phips_bins_file = kwargs.pop("phips_bins_file", None)
        if phips_bins_file is None:
            raise ValueError(
                "Auto-detected PHIPS mode but phips_bins_file not specified. "
                "Please provide phips_bins_file or set mode='phips' explicitly."
            )
        return PHIPSMode(bins_file=phips_bins_file)
    else:
        # Standard mode
        n_theta = kwargs.pop("n_theta", 181)
        n_phi = kwargs.pop("n_phi", 181)
        return StandardMode(n_theta=n_theta, n_phi=n_phi)
