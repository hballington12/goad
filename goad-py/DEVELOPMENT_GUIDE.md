# GOAD Python Development Guide

## Quick Setup & Testing

### 1. Build and Install
```bash
cd goad-py
./build_and_test.sh
```

### 2. Test the Bindings
```bash
source .venv/bin/activate
python test_multiproblem.py
```

### 3. Run Examples
```bash
python multiproblem_example.py
```

## Issues Fixed

### ✅ **Function Signature Bugs**
- **Problem**: `create_uniform_orientation()` missing required positional argument
- **Solution**: Added `#[pyo3(signature = (num_orients, euler_convention = None))]` to make `euler_convention` optional with default value

- **Problem**: `MultiProblem.__new__()` and `Problem.__new__()` missing required positional arguments
- **Solution**: Added `#[pyo3(signature = (settings = None, geom = None))]` to make both parameters optional with default values

### ✅ **Build Script Issues**  
- **Problem**: maturin.sh tries to install wrong Python version wheel
- **Solution**: Created `build_and_test.sh` that auto-detects Python version and handles the full build pipeline
- **Improvement**: Script creates venv if needed and installs correct wheel automatically

### ✅ **Editor Syntax Errors**
- **Problem**: Editor shows syntax errors for goad methods
- **Solution**: Created `goad_py.pyi` type stub file with complete type annotations
- **Benefit**: Provides autocompletion and eliminates false syntax errors

### ✅ **Slow Maturin Builds**
- **Improvement**: New script includes cargo build optimization and better feedback

## Development Workflow

### For Code Changes:
1. Make Rust changes
2. Run `./build_and_test.sh`
3. Test with `python test_multiproblem.py`

### For Editor Support:
- Ensure `goad_py.pyi` is in the same directory as your Python files
- Most editors will automatically pick up the type hints

### Common Commands:
```bash
# Full rebuild and test
./build_and_test.sh

# Quick test after install  
source .venv/bin/activate && python test_multiproblem.py

# Check Python can import module
python -c "import goad_py; print('OK')"
```

## API Usage Examples

### Basic MultiProblem:
```python
import goad_py as goad

# Uniform orientations
orientation = goad.create_uniform_orientation(100)
settings = goad.Settings(geom_name="hex.obj", orientation=orientation)
multi_problem = goad.MultiProblem(settings)
multi_problem.py_solve()
results = multi_problem.results

# Discrete orientations
eulers = [goad.Euler(0, 0, 0), goad.Euler(30, 30, 30)]
orientation = goad.create_discrete_orientation(eulers)
settings = goad.Settings(geom_name="hex.obj", orientation=orientation)
```

## Type Hints Available

The `goad_py.pyi` file provides complete type annotations for:
- All classes (Problem, MultiProblem, Results, Settings, etc.)
- All methods and properties
- Helper functions
- Optional parameters with defaults

This should eliminate most editor syntax errors and provide good autocompletion.