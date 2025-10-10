# Unified Convergence API - Implementation Summary

## What We Built

A unified, single-entry-point API for all GOAD convergence studies with:

✅ **Strict validation** - Errors are thrown immediately for invalid inputs
✅ **Mode-based architecture** - StandardMode and PHIPSMode with specific validation rules
✅ **Uniform output** - UnifiedResults with JSON serialization
✅ **Auto-detection** - Automatically selects appropriate mode based on targets
✅ **Single or ensemble** - Automatically handles file vs directory geometry input
✅ **Parameter sweeps** - Built-in support for running multiple configurations
✅ **Backward compatible** - Legacy API still available

## Files Created

1. **`python/goad_py/unified_convergence.py`** (~850 lines)
   - `ConvergenceMode` abstract base class
   - `StandardMode` - for integrated parameters and Mueller elements
   - `PHIPSMode` - for PHIPS detector DSCS
   - `ConvergenceConfig` - validated configuration dataclass
   - `UnifiedResults` - results container with JSON save/load
   - `UnifiedConvergence` - main convergence runner
   - `run_convergence()` - primary entry point
   - `run_convergence_sweep()` - parameter sweep support

2. **`python/goad_py/__init__.py`** - Updated exports

3. **`unified_convergence_example.py`** - Comprehensive examples

4. **`test_unified_validation.py`** - Validation tests (all passing ✓)

5. **`test_unified_quick.py`** - End-to-end tests (all passing ✓)

6. **`UNIFIED_API.md`** - Complete documentation

## Key Design Decisions

### 1. Mode Classes Enforce Validation

```python
class StandardMode:
    VALID_SCALAR_TARGETS = {"asymmetry", "scatt", "ext", "albedo"}
    VALID_MUELLER_TARGETS = {"S11", "S12", ..., "S44"}
    
class PHIPSMode:
    def __init__(self, bins_file: str):
        # Validates bins file exists and is valid TOML
        # Loads custom bins for later use
```

**Rationale:** Each mode knows what's valid for itself. StandardMode can't be used with PHIPS bins, and PHIPSMode can't compute integrated parameters.

### 2. ConvergenceConfig Validates in __post_init__

```python
@dataclass
class ConvergenceConfig:
    def __post_init__(self):
        # Validate geometry exists
        # Validate mode is correct type
        # Let mode validate its targets
        # Validate numeric parameters > 0
        # Check mode-specific constraints
```

**Rationale:** Fail fast with clear error messages rather than mysterious failures during convergence.

### 3. Single Entry Point with Smart Defaults

```python
run_convergence(
    geometry="hex.obj",
    targets="asymmetry",  # Can be string, list of strings, or list of dicts
    tolerance=0.05,       # Applied to all unless overridden
    mode="auto",          # Auto-detect from targets
)
```

**Rationale:** Easy for simple cases, but can scale to complex configurations.

### 4. Uniform Results Format

```python
@dataclass
class UnifiedResults:
    converged: bool
    n_orientations: int
    mode: str
    is_ensemble: bool
    values: Dict[str, Union[float, np.ndarray]]
    sem_values: Dict[str, Union[float, np.ndarray]]
    # ... plus save/load/summary methods
```

**Rationale:** Same interface regardless of mode, making parameter sweeps trivial.

## Usage Examples

### Simplest Case
```python
results = goad.run_convergence("hex.obj", "asymmetry", tolerance=0.01)
```

### Multiple Targets
```python
results = goad.run_convergence(
    "hex.obj", 
    ["asymmetry", "scatt"],
    tolerance=0.05
)
```

### Mueller Elements with Specific Bins
```python
results = goad.run_convergence(
    "hex.obj",
    [{"variable": "S11", "tolerance": 0.1, "theta_indices": [180]}]
)
```

### Ensemble
```python
results = goad.run_convergence("./geometries/", "asymmetry", tolerance=0.05)
```

### PHIPS
```python
results = goad.run_convergence(
    "./geometries/",
    "phips_dscs",
    tolerance=0.25,
    mode="phips",
    phips_bins_file="phips_bins_edges.toml"
)
```

### Parameter Sweep
```python
configs = [
    ConvergenceConfig(
        geometry="hex.obj",
        mode=StandardMode(),
        convergence_targets=[{"variable": "asymmetry", "tolerance": 0.05}],
        wavelength=wl
    )
    for wl in [0.532, 0.633, 0.780]
]
results_list = run_convergence_sweep(configs)
```

## Validation Features

### What Gets Validated

✓ Geometry path exists (file or directory)
✓ Mode is correct type (StandardMode or PHIPSMode)
✓ Targets are valid for the mode
✓ PHIPS bins file exists and is valid TOML
✓ Numeric parameters are positive
✓ Mueller options not used in PHIPS mode
✓ theta_indices only used with Mueller elements
✓ Empty convergence_targets rejected

### Error Examples

```python
# Invalid geometry
config = ConvergenceConfig(geometry="nonexistent.obj", ...)
# FileNotFoundError: Geometry path does not exist

# Invalid target for StandardMode
run_convergence("hex.obj", targets="invalid_target")
# ValueError: Invalid target 'invalid_target' for StandardMode

# Mueller options in PHIPS mode
ConvergenceConfig(mode=PHIPSMode(...), mueller_1d=True, ...)
# ValueError: mueller_1d not supported in PHIPSMode

# theta_indices on non-Mueller variable
run_convergence("hex.obj", [{
    "variable": "asymmetry", 
    "theta_indices": [180]
}])
# ValueError: theta_indices can only be used with Mueller elements
```

## Testing

All tests passing ✓

**Validation tests** (`test_unified_validation.py`):
- StandardMode creation and validation
- PHIPSMode validation
- ConvergenceConfig validation
- Invalid input rejection
- Config serialization

**End-to-end tests** (`test_unified_quick.py`):
- Simple convergence
- Multiple targets
- Mueller element convergence
- Save/load JSON
- Using ConvergenceConfig directly
- Summary and serialization

## Future Enhancements

1. **Parallel parameter sweeps** - Use multiprocessing for sweep execution
2. **Plotting integration** - `results.plot()` method for standard visualizations
3. **HDF5/NetCDF output** - For large datasets and Mueller matrices
4. **Progress callbacks** - For monitoring long convergences
5. **Convergence diagnostics** - Autocorrelation, effective sample size, etc.

## Migration Path

The legacy API (`Convergence`, `EnsembleConvergence`, `PHIPSConvergence`) remains available and functional. Users can migrate incrementally:

1. **Phase 1**: Use unified API for new scripts
2. **Phase 2**: Update documentation to recommend unified API
3. **Phase 3**: Add deprecation warnings to legacy API
4. **Phase 4**: Remove legacy API (major version bump)

No breaking changes in this release.
