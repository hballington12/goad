# Unified Convergence API Examples

This directory contains examples demonstrating the unified convergence API for GOAD.

The unified API provides a single, consistent interface for all convergence types, with strict validation and easy-to-use defaults.

## Examples

### Basic Examples

**01_simple_asymmetry.py**
- Simplest use case: converge on asymmetry parameter
- Shows minimal configuration needed

**02_interval_binning.py**
- Using custom interval binning for focused angular resolution
- Demonstrates StandardMode with custom binning

**03_multiple_targets.py**
- Converging on multiple parameters simultaneously
- Shows how to use list of target strings

**04_mueller_element.py**
- Mueller matrix element (S11) convergence
- Demonstrates convergence across all scattering angles

**05_backscatter_specific_bin.py**
- Converging on specific theta bins only
- Shows using theta_indices parameter

**06_ensemble_convergence.py**
- Ensemble averaging over multiple geometries
- Automatically detected by passing a directory

### Advanced Examples

**07_advanced_config.py**
- Using ConvergenceConfig for full control
- Shows all available parameters
- Custom optical properties

**08_parameter_sweep.py**
- Running multiple convergences (wavelength sweep)
- Using run_convergence_sweep()

## Running Examples

From the repository root:

```bash
cd goad-py
source .venv/bin/activate
python examples/unified/01_simple_asymmetry.py
```

Or from the examples/unified directory:

```bash
cd goad-py/examples/unified
source ../../.venv/bin/activate
python 01_simple_asymmetry.py
```

## Quick Start

The simplest example:

```python
import goad_py as goad

results = goad.run_convergence(
    geometry="hex.obj",
    targets="asymmetry",
    tolerance=0.01
)

print(results.summary())
results.save("results.json")
```

## See Also

- **examples/direct/** - Examples using the direct GOAD API (lower-level)
- **UNIFIED_API.md** - Complete API documentation
- **unified_convergence_example.py** - All examples in one comprehensive script
