# Unified Convergence API

A single, consistent interface for all GOAD convergence studies.

## Overview

The unified convergence API provides a single entry point for running convergence studies with GOAD, supporting:

- **Standard mode**: Integrated scattering parameters (asymmetry, scattering/extinction cross-sections, albedo) and Mueller matrix elements (S11-S44)
- **PHIPS mode**: Custom detector DSCS values at 20 PHIPS detectors
- **Single geometry** or **ensemble averaging** (multiple geometries)
- **Parameter sweeps** for sensitivity studies
- **JSON serialization** for saving/loading results

## Quick Start

### Simple Example

```python
import goad_py as goad

# Run convergence on asymmetry parameter
results = goad.run_convergence(
    geometry="hex.obj",
    targets="asymmetry",
    tolerance=0.01,  # 1% relative tolerance
)

print(results.summary())
results.save("results.json")
```

### Multiple Targets

```python
results = goad.run_convergence(
    geometry="hex.obj",
    targets=["asymmetry", "scatt", "ext"],
    tolerance=0.05,  # 5% relative tolerance
    batch_size=24,
)
```

### Mueller Matrix Elements

```python
# Converge all theta bins
results = goad.run_convergence(
    geometry="hex.obj",
    targets="S11",
    tolerance=0.25,
)

# Or converge specific bins only (e.g., backscattering)
results = goad.run_convergence(
    geometry="hex.obj",
    targets=[{
        "variable": "S11",
        "tolerance": 0.10,
        "theta_indices": [180]  # Only backscattering
    }],
)
```

### Ensemble Convergence

```python
# Pass a directory instead of a file
results = goad.run_convergence(
    geometry="./geometries/",  # Directory with .obj files
    targets=["asymmetry", "scatt"],
    tolerance=0.05,
)
```

### PHIPS Detector Convergence

```python
results = goad.run_convergence(
    geometry="hex.obj",
    targets="phips_dscs",
    tolerance=0.25,
    mode="phips",
    phips_bins_file="phips_bins_edges.toml",
)
```

## Advanced Usage

### Using ConvergenceConfig

For full control over all parameters:

```python
from goad_py import ConvergenceConfig, StandardMode, UnifiedConvergence

config = ConvergenceConfig(
    geometry="hex.obj",
    mode=StandardMode(n_theta=181, n_phi=181),
    convergence_targets=[
        {"variable": "asymmetry", "tolerance": 0.005, "tolerance_type": "absolute"},
        {"variable": "scatt", "tolerance": 0.01, "tolerance_type": "relative"},
    ],
    wavelength=0.633,  # Red laser
    particle_refr_index_re=1.5,
    particle_refr_index_im=0.01,
    batch_size=24,
    max_orientations=10_000,
    min_batches=10,
    mueller_1d=True,
)

conv = UnifiedConvergence(config)
results = conv.run()
```

### Parameter Sweeps

```python
from goad_py import ConvergenceConfig, StandardMode, run_convergence_sweep

# Wavelength sweep
wavelengths = [0.532, 0.633, 0.780]
configs = [
    ConvergenceConfig(
        geometry="hex.obj",
        mode=StandardMode(),
        convergence_targets=[{"variable": "asymmetry", "tolerance": 0.05}],
        wavelength=wl,
    )
    for wl in wavelengths
]

results_list = run_convergence_sweep(configs)

for wl, result in zip(wavelengths, results_list):
    print(f"λ={wl}: asymmetry = {result.values['asymmetry']:.6f}")
```

## API Reference

### Functions

#### `run_convergence(geometry, targets, tolerance=0.05, tolerance_type="relative", mode="auto", **kwargs)`

Main entry point for convergence studies.

**Parameters:**
- `geometry` (str|Path): Path to .obj file or directory
- `targets` (str|list): What to converge on (e.g., "asymmetry", ["asymmetry", "scatt"], or list of dicts)
- `tolerance` (float): Default tolerance for convergence
- `tolerance_type` (str): "relative" or "absolute"
- `mode` (str|ConvergenceMode): "auto", "standard", "phips", or a ConvergenceMode instance
- `**kwargs`: Additional settings (wavelength, batch_size, etc.)

**Returns:** `UnifiedResults`

#### `run_convergence_sweep(configs, parallel=False)`

Run multiple convergence studies for parameter sweeps.

**Parameters:**
- `configs` (List[ConvergenceConfig]): List of configurations
- `parallel` (bool): Run in parallel (not yet implemented)

**Returns:** `List[UnifiedResults]`

### Classes

#### `StandardMode(n_theta=181, n_phi=181)`

Standard convergence mode for integrated parameters and Mueller elements.

**Valid targets:**
- Scalar: `asymmetry`, `scatt`, `ext`, `albedo`
- Mueller: `S11`, `S12`, ..., `S44`
- Mueller with specific bins: `{"variable": "S11", "theta_indices": [0, 180]}`

#### `PHIPSMode(bins_file)`

PHIPS detector convergence mode.

**Parameters:**
- `bins_file` (str): Path to PHIPS bins TOML file

**Valid targets:**
- Only supports a single target with tolerance specification
- Optionally specify `detector_indices` to converge specific detectors

#### `ConvergenceConfig`

Full configuration for convergence studies.

**Required fields:**
- `geometry` (str|Path): Geometry file or directory
- `mode` (ConvergenceMode): StandardMode or PHIPSMode
- `convergence_targets` (List[dict]): List of target specifications

**Optional fields:**
- `wavelength` (float): Default 0.532 μm
- `particle_refr_index_re/im` (float): Refractive index (default: 1.31+0j)
- `medium_refr_index_re/im` (float): Medium refractive index (default: 1.0+0j)
- `batch_size` (int): Orientations per batch (default: 24)
- `max_orientations` (int): Maximum orientations (default: 100,000)
- `min_batches` (int): Minimum batches before convergence (default: 10)
- `mueller_1d` (bool): Compute 1D Mueller matrix (default: True, standard mode only)
- `mueller_2d` (bool): Compute 2D Mueller matrix (default: False, standard mode only)

#### `UnifiedResults`

Results container with uniform interface.

**Attributes:**
- `converged` (bool): Whether convergence criteria were met
- `n_orientations` (int): Total orientations run
- `mode` (str): "standard" or "phips"
- `is_ensemble` (bool): Whether ensemble averaging was used
- `values` (dict): Final mean values
- `sem_values` (dict): Final SEM values
- `mueller_1d` (ndarray): 1D Mueller matrix (if computed)
- `mueller_2d` (ndarray): 2D Mueller matrix (if computed)
- `convergence_history` (list): Convergence progress
- `warning` (str): Warning message (if any)

**Methods:**
- `summary()`: Generate human-readable summary
- `save(path)`: Save to JSON file
- `load(path)`: Load from JSON file (class method)
- `to_dict()`: Convert to dictionary

## Validation

The API performs comprehensive input validation:

### Mode-Specific Validation

**StandardMode:**
- Rejects invalid variable names
- Rejects `theta_indices` on non-Mueller variables
- Validates theta_indices are within range

**PHIPSMode:**
- Requires valid PHIPS bins TOML file
- Only allows a single convergence target
- Rejects `mueller_1d` and `mueller_2d` options (not compatible with custom binning)

### Config Validation

- Geometry path must exist
- `mode` must be a ConvergenceMode instance
- `convergence_targets` cannot be empty
- Numeric parameters must be positive
- Mode-specific targets are validated

**The API throws errors immediately rather than failing silently.**

## Examples

See `unified_convergence_example.py` for comprehensive examples covering:

1. Simple convergence on single parameter
2. Multiple convergence targets
3. Mueller matrix element convergence
4. Backscatter convergence (specific theta bins)
5. Ensemble convergence
6. PHIPS detector convergence
7. PHIPS ensemble convergence
8. Advanced usage with ConvergenceConfig
9. Parameter sweeps
10. Save/load results

## Migration from Legacy API

The legacy convergence classes (`Convergence`, `EnsembleConvergence`, `PHIPSConvergence`, etc.) are still available, but the unified API is recommended for new code.

### Before (Legacy):
```python
from goad_py import Convergence, Convergable

convergables = [
    Convergable('asymmetry', 'relative', 0.01),
    Convergable('scatt', 'relative', 0.05),
]
conv = Convergence(settings, convergables, batch_size=24)
results = conv.run()
```

### After (Unified):
```python
from goad_py import run_convergence

results = run_convergence(
    geometry="hex.obj",
    targets=[
        {"variable": "asymmetry", "tolerance": 0.01},
        {"variable": "scatt", "tolerance": 0.05},
    ],
    batch_size=24,
)
```

## Future Work

- Parallel execution for parameter sweeps
- Additional output formats (HDF5, NetCDF)
- Plotting utilities integrated into UnifiedResults
- Progress callbacks for long-running convergences
