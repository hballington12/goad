# GOAD-PY

Python bindings for GOAD (Geometric Optics Approximation with Diffraction) - a physical optics light scattering computation library.

## Installation

```bash
pip install goad-py
```

## Quick Start

```python
import goad_py

# Create a problem with minimal setup
settings = goad_py.Settings("path/to/geometry.obj")
mp = goad_py.MultiProblem(settings)
mp.py_solve()

# Access scattering data
results = mp.results
print(f"Scattering cross-section: {results.scat_cross}")
print(f"Extinction cross-section: {results.ext_cross}")
print(f"Asymmetry parameter: {results.asymmetry}")
```

### Convergence Analysis

For statistical error estimation, use the convergence analysis functionality:

```python
from goad_py import Convergence, Convergable

# Set up convergence analysis
convergence = Convergence(
    settings=goad_py.Settings(geom_path="path/to/geometry.obj"),
    convergables=[
        Convergable('asymmetry', 'absolute', 0.005),  # absolute SEM < 0.005
        Convergable('scatt', 'relative', 0.01),       # relative SEM < 1%
    ],
    batch_size=100
)

# Run until convergence
results = convergence.run()
print(f"Converged: {results.converged}")
print(f"Final values: {results.values}")
```

## Features

- Fast light scattering computations using physical optics
- Support for various 3D geometry formats
- Configurable wavelength, refractive index, and orientations
- Multi-orientation averaging capabilities
- Convergence analysis for statistical error estimation
- Efficient parallel computation with GIL release

## Documentation

- [Rust API Documentation](https://docs.rs/goad/0.1.0/goad/index.html)
- [GitHub Repository](https://github.com/hballington12/goad)

## License

GPL-3.0 License - see the LICENSE file in the main repository for details.