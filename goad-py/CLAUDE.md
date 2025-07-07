# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GOAD-PY is the Python binding for GOAD (Geometric Optics Approximation with Diffraction), a physical optics light scattering computation library. It uses PyO3 to create efficient Rust-Python bindings for simulating light scattering by arbitrary 3D geometries.

## Build and Development Commands

### Python Package Development

```bash
# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install development dependencies
pip install --upgrade pip maturin

# Build and install in development mode
maturin develop

# Build release wheel
maturin build --release

# Complete build and test workflow
./build_and_test.sh

# Test built wheels across Python versions
./test_wheels.sh
```

### Testing

```bash
# Run example scripts
python simple_example.py
python multiproblem_example.py

# Basic import test
python -c "import goad_py; print('✅ GOAD Python bindings imported successfully!')"
```

### Release Process

```bash
# Version management
./release.sh patch    # Bump patch version
./release.sh minor    # Bump minor version
./release.sh major    # Bump major version
./release.sh tag      # Create git tag for release

# Test publishing
./publish_test.sh     # Upload to TestPyPI

# Production release (automatic via GitHub Actions)
git push origin vX.Y.Z  # Push tag to trigger CI release
```

## Architecture and Code Structure

### Project Layout

```
goad-py/
├── src/lib.rs           # PyO3 bindings implementation
├── goad_py.pyi          # Python type stubs
├── Cargo.toml           # Rust dependencies and configuration
├── pyproject.toml       # Python package configuration
└── examples/            # Usage examples
```

### Key API Components

1. **Settings** - Configuration container for all simulation parameters
   - Geometry path (required)
   - Wavelength, refractive indices, orientation, binning scheme (optional with smart defaults)

2. **Problem** - Single orientation scattering calculation
   - Uses Settings for configuration
   - Results accessible after py_solve()

3. **MultiProblem** - Multi-orientation averaging
   - Parallel processing of multiple orientations
   - Statistical averaging of scattering properties

4. **BinningScheme** - Angular discretization control
   - Simple: Regular theta-phi grid
   - Interval: Variable resolution binning
   - Custom: Explicit angle pairs

### Python-Rust Bridge Architecture

- Direct PyO3 bindings without Python wrapper overhead
- GIL released during computations for parallel Python operations
- Efficient data transfer using numpy-compatible arrays
- Type safety through comprehensive .pyi stubs

### Default Behavior

When users create minimal Settings, these defaults apply:
- Wavelength: 532nm (green laser)
- Refractive indices: 1.31+0i (glass) in 1.0+0i (air)
- Single random orientation
- ~1000 bins with adaptive angular resolution

## Important Development Notes

1. **Path Handling**: Always convert relative paths to absolute internally
2. **Config-Free Design**: Unlike the Rust version, Python API doesn't use config files
3. **Type Stubs**: Keep goad_py.pyi synchronized with lib.rs changes
4. **Examples**: Update simple_example.py and multiproblem_example.py when API changes
5. **Multi-Platform**: Test on Linux, Windows, and macOS before releases
6. **Python Versions**: Support Python 3.8-3.12 using abi3 compatibility

## CI/CD Pipeline

GitHub Actions automatically:
1. Builds wheels for all platforms on push
2. Tests on Python 3.8-3.12
3. Publishes to TestPyPI from main branch
4. Publishes to PyPI on version tags (vX.Y.Z)

Uses trusted publishing (OIDC) - no manual API tokens needed.