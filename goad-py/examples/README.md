# GOAD Python Examples

This directory contains examples demonstrating different ways to use GOAD's Python bindings.

## Directory Structure

```
examples/
├── unified/          # Unified convergence API examples (RECOMMENDED)
│   ├── 01_simple_asymmetry.py
│   ├── 02_interval_binning.py
│   ├── 03_multiple_targets.py
│   ├── 04_mueller_element.py
│   ├── 05_backscatter_specific_bin.py
│   ├── 06_ensemble_convergence.py
│   ├── 07_advanced_config.py
│   ├── 08_parameter_sweep.py
│   └── README.md
│
└── direct/           # Direct GOAD API examples (lower-level)
    ├── simple_example.py
    ├── multiproblem_example.py
    ├── convergence_example.py
    ├── ensemble_example.py
    ├── s11_convergence_example.py
    ├── backscatter_convergence_example.py
    ├── phips_convergence_example.py
    ├── phips_ensemble_convergence_example.py
    └── README.md
```

## Which Examples to Use?

### 🌟 New Users: Start with `unified/`

The unified API examples (`examples/unified/`) are **recommended for new code**:
- ✅ Simple, consistent interface
- ✅ Automatic validation with clear error messages
- ✅ Easy parameter sweeps
- ✅ Uniform output format (JSON serialization)

**Start here:** `unified/01_simple_asymmetry.py`

### 🔧 Advanced Users: `direct/` for fine-grained control

The direct API examples (`examples/direct/`) use lower-level classes:
- Access to intermediate results
- Maximum control over GOAD settings
- Custom convergence logic
- Integration with existing code

## Quick Start

### Unified API (Recommended)

```python
import goad_py as goad

# Simplest possible convergence
results = goad.run_convergence(
    geometry="hex.obj",
    targets="asymmetry",
    tolerance=0.01  # 1% relative tolerance
)

print(results.summary())
results.save("results.json")
```

### Direct API (Lower-level)

```python
import goad_py as goad

# Create settings
settings = goad.Settings("hex.obj")

# Create convergence criteria
convergables = [
    goad.Convergable('asymmetry', 'relative', 0.01)
]

# Run convergence
conv = goad.Convergence(settings, convergables, batch_size=24)
results = conv.run()
```

## Running Examples

From the repository root:

```bash
cd goad-py
source .venv/bin/activate

# Run unified examples
python examples/unified/01_simple_asymmetry.py

# Run direct examples
python examples/direct/simple_example.py
```

## Example Categories

### Unified API Examples

1. **Basic convergence** - Single parameter, multiple parameters, Mueller elements
2. **Custom binning** - Interval binning for focused resolution
3. **Ensemble averaging** - Multiple geometries
4. **Advanced configuration** - Full control with ConvergenceConfig
5. **Parameter sweeps** - Running multiple configurations

### Direct API Examples

1. **Simple usage** - Basic MultiProblem
2. **Custom binning** - Simple, interval, and custom schemes
3. **Convergence studies** - Using Convergence classes
4. **PHIPS detectors** - Custom detector geometry

## Documentation

- **UNIFIED_API.md** - Complete unified API reference
- **UNIFIED_API_SUMMARY.md** - Implementation details
- **unified/** and **direct/** - Individual README files in each directory

## Migration Guide

If you have existing code using the direct API, you can migrate to the unified API:

**Before (Direct API):**
```python
convergables = [Convergable('asymmetry', 'relative', 0.01)]
conv = Convergence(settings, convergables, batch_size=24)
results = conv.run()
```

**After (Unified API):**
```python
results = goad.run_convergence(
    geometry="hex.obj",
    targets="asymmetry",
    tolerance=0.01,
    batch_size=24
)
```

Both APIs remain supported - no breaking changes.

## Getting Help

- Check the README files in `unified/` and `direct/`
- Read `UNIFIED_API.md` for complete API documentation
- Run the examples to see working code
- Open an issue on GitHub if you find problems
