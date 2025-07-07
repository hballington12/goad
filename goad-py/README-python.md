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
problem = goad_py.Problem("path/to/geometry.obj")

# Solve and get results
results = problem.py_solve()

# Access scattering data
print(f"Number of angles: {results.num_angles}")
print(f"Scattering cross section: {results.sca_cross_section}")
```

## Features

- Fast light scattering computations using physical optics
- Support for various 3D geometry formats
- Configurable wavelength, refractive index, and orientations
- Multi-orientation averaging capabilities
- Efficient parallel computation with GIL release

## Documentation

- [Rust API Documentation](https://docs.rs/goad/0.1.0/goad/index.html)
- [GitHub Repository](https://github.com/hballington12/goad)

## License

GPL-3.0 License - see the LICENSE file in the main repository for details.