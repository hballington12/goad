"""
Type stubs for goad_py module to help with editor syntax highlighting and autocompletion.
"""

from typing import Optional, List, Dict, Any, Union

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

class Settings:
    """Problem settings and configuration."""
    
    def __init__(
        self,
        wavelength: Optional[float] = None,
        beam_power_threshold: Optional[float] = None,
        beam_area_threshold_fac: Optional[float] = None,
        cutoff: Optional[float] = None,
        medium_refr_index_re: Optional[float] = None,
        medium_refr_index_im: Optional[float] = None,
        particle_refr_index_re: Optional[float] = None,
        particle_refr_index_im: Optional[float] = None,
        geom_name: Optional[str] = None,
        max_rec: Optional[int] = None,
        max_tir: Optional[int] = None,
        theta_res: Optional[int] = None,
        phi_res: Optional[int] = None,
        euler: Optional[List[float]] = None,
        orientation: Optional[Orientation] = None
    ) -> None: ...
    
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
    """Multi-orientation averaging problem."""
    
    def __init__(self, settings: Optional[Settings] = None, geom: Optional[Geom] = None) -> None: ...
    
    def py_solve(self) -> None: ...
    
    def py_writeup(self) -> None: ...
    
    def py_reset(self) -> None: ...
    
    def py_regenerate_orientations(self) -> None: ...
    
    @property
    def results(self) -> Results: ...
    
    @property
    def num_orientations(self) -> int: ...

# Helper functions
def uniform_orientation(num_orients: int) -> Scheme: ...

def discrete_orientation(eulers: List[Euler]) -> Scheme: ...

def create_uniform_orientation(num_orients: int, euler_convention: Optional[EulerConvention] = None) -> Orientation: ...

def create_discrete_orientation(eulers: List[Euler], euler_convention: Optional[EulerConvention] = None) -> Orientation: ...

def sum_as_string(a: int, b: int) -> str: ...

def goad_py_add() -> None: ...