"""
GOAD Python API Type Definitions

This file provides comprehensive type hints for the GOAD (Geometric Optics 
Approximation for Diffraction) Python bindings. GOAD simulates light scattering 
by arbitrary 3D geometries using geometric optics with diffraction corrections.

The main workflow is:
1. Create Settings with geometry path (all other parameters have sensible defaults)
2. Use Problem for single-orientation or MultiProblem for multi-orientation simulations  
3. Call py_solve() to run the computation
4. Access results through the .results property

Default behavior (minimal setup):
- Single random orientation (perfect for initial exploration)
- Interval angular binning with ~1000 bins (good coverage, reasonable performance)
- Wavelength: 532nm (green laser)
- Refractive index: 1.31 + 0.0i (typical glass particle in air)

Example (minimal setup):
    import goad_py as goad
    
    settings = goad.Settings("particle.obj")
    mp = goad.MultiProblem(settings)
    mp.py_solve()
    
    results = mp.results
    print(f"Scattering cross-section: {results.scat_cross}")
    print(f"Extinction cross-section: {results.ext_cross}")
    print(f"Asymmetry parameter: {results.asymmetry}")
"""

from typing import Optional, List, Dict, Any, Union, Tuple
import numpy as np

class Euler:
    """Euler angles for rotations."""
    alpha: float
    beta: float
    gamma: float
    
    def __init__(self, alpha: float, beta: float, gamma: float) -> None: ...
    def __repr__(self) -> str: ...

class EulerConvention:
    """Euler angle conventions."""
    XZX: 'EulerConvention'
    XYX: 'EulerConvention'
    YXY: 'EulerConvention'
    YZY: 'EulerConvention'
    ZYZ: 'EulerConvention'
    ZXZ: 'EulerConvention'
    XZY: 'EulerConvention'
    XYZ: 'EulerConvention'
    YXZ: 'EulerConvention'
    YZX: 'EulerConvention'
    ZYX: 'EulerConvention'
    ZXY: 'EulerConvention'

class Scheme:
    """Orientation scheme (uniform or discrete)."""
    pass

class Orientation:
    """Full orientation specification."""
    scheme: Scheme
    euler_convention: EulerConvention
    
    def __init__(self, scheme: Scheme, euler_convention: Optional[EulerConvention] = None) -> None: ...
    def __repr__(self) -> str: ...

class Geom:
    """Geometry representation."""
    pass

class Shape:
    """Shape within geometry."""
    pass

class Results:
    """Results from problem solving."""
    
    @property
    def bins(self) -> List[tuple[float, float]]: ...
    
    @property  
    def bins_1d(self) -> Optional[List[float]]: ...
    
    @property
    def mueller(self) -> List[List[float]]: ...
    
    @property
    def mueller_beam(self) -> List[List[float]]: ...
    
    @property
    def mueller_ext(self) -> List[List[float]]: ...
    
    @property
    def mueller_1d(self) -> List[List[float]]: ...
    
    @property
    def mueller_1d_beam(self) -> List[List[float]]: ...
    
    @property
    def mueller_1d_ext(self) -> List[List[float]]: ...
    
    @property
    def asymmetry(self) -> Optional[float]: ...
    
    @property
    def scat_cross(self) -> Optional[float]: ...
    
    @property
    def ext_cross(self) -> Optional[float]: ...
    
    @property
    def albedo(self) -> Optional[float]: ...
    
    @property
    def params(self) -> Dict[str, Optional[float]]: ...
    
    @property
    def powers(self) -> Dict[str, float]: ...

class BinningScheme:
    """Angular binning scheme for scattering calculations.
    
    Defines how to discretize the scattering sphere into angular bins
    for Mueller matrix and amplitude computations. Supports simple
    regular grids, custom intervals, and arbitrary bin arrangements.
    """
    
    def __init__(self, bins: List[tuple[float, float]]) -> None: ...
    
    @staticmethod
    def simple(num_theta: int, num_phi: int) -> 'BinningScheme':
        """Create a simple regular grid binning scheme.
        
        Args:
            num_theta: Number of bins in theta direction (0 to 180 degrees)
            num_phi: Number of bins in phi direction (0 to 360 degrees)
            
        Returns:
            BinningScheme with regular angular grid
            
        Raises:
            ValueError: If num_theta or num_phi is zero or negative
        """
        ...
    
    @staticmethod
    def interval(
        thetas: List[float],
        theta_spacings: List[float], 
        phis: List[float],
        phi_spacings: List[float]
    ) -> 'BinningScheme':
        """Create binning scheme with custom intervals.
        
        Args:
            thetas: Key theta angles in degrees
            theta_spacings: Spacing around each theta value
            phis: Key phi angles in degrees  
            phi_spacings: Spacing around each phi value
            
        Returns:
            BinningScheme with custom intervals
        """
        ...
    
    @staticmethod
    def custom(bins: List[tuple[float, float]]) -> 'BinningScheme':
        """Create binning scheme from explicit (theta, phi) pairs.
        
        Args:
            bins: List of (theta, phi) angle pairs in degrees
            
        Returns:
            BinningScheme with custom bin locations
        """
        ...

class Settings:
    """Simulation parameters and physical properties.
    
    Contains all parameters needed for a GOAD simulation including
    geometry path, optical properties, orientation scheme, and
    numerical settings. Most parameters have sensible defaults.
    """
    
    def __init__(
        self,
        geom_path: str,
        wavelength: float = 0.532,
        particle_refr_index_re: float = 1.31,
        particle_refr_index_im: float = 0.0,
        medium_refr_index_re: float = 1.0,
        medium_refr_index_im: float = 0.0,
        orientation: Optional[Orientation] = None,
        binning: Optional[BinningScheme] = None,
        beam_power_threshold: float = 0.005,
        beam_area_threshold_fac: float = 0.1,
        cutoff: float = 0.99,
        max_rec: int = 10,
        max_tir: int = 10,
        scale: float = 1.0,
        directory: str = "goad_run"
    ) -> None:
        """Initialize simulation settings.
        
        Args:
            geom_path: Path to geometry file (.obj format)
            wavelength: Incident wavelength in geometry units (default: 0.532)
            particle_refr_index_re: Real part of particle refractive index (default: 1.31)
            particle_refr_index_im: Imaginary part of particle refractive index (default: 0.0)
            medium_refr_index_re: Real part of medium refractive index (default: 1.0)
            medium_refr_index_im: Imaginary part of medium refractive index (default: 0.0)
            orientation: Orientation scheme (default: None → 1 random uniform orientation)
            binning: Angular binning scheme (default: None → interval binning: theta=[0°,30°,60°,120°,180°] with spacings [2°,5°,10°,5°], phi=[0°,90°,180°,270°,360°] with 15° spacings, creating ~1000 bins)
            beam_power_threshold: Ray termination threshold (default: 0.005)
            beam_area_threshold_fac: Area threshold factor (default: 0.1)
            cutoff: Power cutoff fraction 0-1 (default: 0.99)
            max_rec: Maximum recursion depth (default: 10)
            max_tir: Maximum total internal reflections (default: 10)
            scale: Geometry scaling factor (default: 1.0)
            directory: Output directory for results (default: "goad_run")
            
        Raises:
            ValueError: If wavelength <= 0 or cutoff not in [0,1]
            FileNotFoundError: If geometry file doesn't exist
        """
        ...
    
    @property
    def euler(self) -> List[float]: ...
    
    @euler.setter
    def euler(self, value: List[float]) -> None: ...
    
    @property
    def orientation(self) -> Orientation: ...
    
    @orientation.setter
    def orientation(self, value: Orientation) -> None: ...

class Problem:
    """Single orientation problem."""
    
    def __init__(self, settings: Optional[Settings] = None, geom: Optional[Geom] = None) -> None: ...
    
    def py_solve(self) -> None: ...
    
    def py_print_stats(self) -> None: ...
    
    @property
    def results(self) -> Results: ...

class MultiProblem:
    """Multi-orientation light scattering simulation for a single geometry.
    
    Computes orientation-averaged scattering properties by running multiple
    single-orientation simulations and averaging the results. Supports both
    random and systematic orientation sampling schemes. Results include
    Mueller matrices, cross-sections, and derived optical parameters.
    
    Example:
        orientations = goad.create_uniform_orientation(100)
        settings = goad.Settings("particle.obj", orientation=orientations)
        mp = goad.MultiProblem(settings)
        mp.py_solve()
        print(f"Scattering cross-section: {mp.results.scat_cross}")
    """
    
    def __init__(self, settings: Settings, geom: Optional[Geom] = None) -> None:
        """Initialize multi-orientation problem.
        
        Args:
            settings: Simulation parameters including orientation scheme
            geom: Geometry object (loaded from settings.geom_path if None)
            
        Raises:
            FileNotFoundError: If geometry file cannot be loaded
        """
        ...
    
    def py_solve(self) -> None:
        """Solve the multi-orientation scattering problem.
        
        Computes scattering properties averaged over all orientations using
        parallel processing. The Global Interpreter Lock (GIL) is released
        during computation to allow concurrent Python operations.
        
        Raises:
            RuntimeError: If computation fails
        """
        ...
    
    def py_writeup(self) -> None:
        """Write results to output files in the specified directory."""
        ...
    
    def py_reset(self) -> None:
        """Reset problem to initial state and regenerate orientations."""
        ...
    
    def py_regenerate_orientations(self) -> None:
        """Regenerate random orientations (useful for statistical sampling)."""
        ...
    
    @property
    def results(self) -> Results:
        """Access orientation-averaged simulation results.
        
        Returns the complete Results object containing Mueller matrices,
        amplitude matrices, power distributions, and derived parameters
        averaged over all orientations.
        """
        ...
    
    @property
    def num_orientations(self) -> int:
        """Number of orientations in the current simulation."""
        ...

# Helper functions
def uniform_orientation(num_orients: int) -> Scheme: ...

def discrete_orientation(eulers: List[Euler]) -> Scheme: ...

def create_uniform_orientation(num_orients: int, euler_convention: Optional[EulerConvention] = None) -> Orientation: ...

def create_discrete_orientation(eulers: List[Euler], euler_convention: Optional[EulerConvention] = None) -> Orientation: ...

def sum_as_string(a: int, b: int) -> str: ...

def goad_py_add() -> None: ...

# Convergence Analysis Classes

class Convergable:
    """Represents a variable to monitor for convergence.
    
    Defines convergence criteria for integrated scattering parameters
    including asymmetry parameter, scattering cross-section, extinction
    cross-section, and single-scattering albedo.
    """
    
    variable: str
    tolerance_type: str
    tolerance: float
    
    def __init__(
        self,
        variable: str,
        tolerance_type: str = 'relative',
        tolerance: float = 0.01
    ) -> None:
        """Initialize convergence criterion.
        
        Args:
            variable: Variable to monitor ('asymmetry', 'scatt', 'ext', 'albedo')
            tolerance_type: 'relative' or 'absolute' tolerance
            tolerance: Tolerance value (relative as fraction, absolute as value)
            
        Raises:
            ValueError: If variable name or tolerance_type is invalid
        """
        ...

class ConvergenceResults:
    """Results from a convergence study.
    
    Contains final convergence status, parameter values with uncertainties,
    and complete convergence history for analysis.
    """
    
    converged: bool
    n_orientations: int
    values: Dict[str, float]
    sem_values: Dict[str, float]
    mueller_1d: Optional[np.ndarray]
    mueller_2d: Optional[np.ndarray]
    convergence_history: List[Tuple[int, str, float]]
    warning: Optional[str]
    
    def __init__(
        self,
        converged: bool,
        n_orientations: int,
        values: Dict[str, float],
        sem_values: Dict[str, float],
        mueller_1d: Optional[np.ndarray] = None,
        mueller_2d: Optional[np.ndarray] = None,
        convergence_history: List[Tuple[int, str, float]] = None,
        warning: Optional[str] = None
    ) -> None: ...

class Convergence:
    """Runs multiple MultiProblems until convergence criteria are met.
    
    Implements statistical convergence analysis for scattering parameters
    using batch-based standard error estimation. Monitors multiple variables
    simultaneously and stops when all meet their convergence criteria.
    
    Example:
        convergence = Convergence(
            settings=goad.Settings("particle.obj"),
            convergables=[
                Convergable('asymmetry', 'absolute', 0.005),
                Convergable('scatt', 'relative', 0.01),
            ],
            batch_size=100
        )
        results = convergence.run()
    """
    
    def __init__(
        self,
        settings: Settings,
        convergables: List[Convergable],
        batch_size: int = 24,
        max_orientations: int = 100_000,
        min_batches: int = 10,
        mueller_1d: bool = True,
        mueller_2d: bool = False
    ) -> None:
        """Initialize convergence study.
        
        Args:
            settings: GOAD settings for the simulation
            convergables: List of variables to monitor for convergence
            batch_size: Number of orientations per iteration
            max_orientations: Maximum total orientations before stopping
            min_batches: Minimum number of batches before allowing convergence
            mueller_1d: Whether to collect 1D Mueller matrices
            mueller_2d: Whether to collect 2D Mueller matrices
            
        Raises:
            ValueError: If parameters are invalid or no convergables specified
        """
        ...
    
    def run(self) -> ConvergenceResults:
        """Run the convergence study.
        
        Executes batches of orientations until all convergence criteria
        are met or maximum orientations reached. Provides progress updates
        and rigorous statistical analysis.
        
        Returns:
            ConvergenceResults containing final values and convergence status
        """
        ...