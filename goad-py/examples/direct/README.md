# Direct GOAD API Examples

This directory contains examples using the direct GOAD API classes.

These examples use the lower-level convergence classes (`Convergence`, `EnsembleConvergence`, `PHIPSConvergence`) and provide more fine-grained control but require more setup.

## Examples

**simple_example.py**
- Basic MultiProblem usage
- Single orientation averaging

**multiproblem_example.py**
- Creating custom binning schemes (simple, interval, custom)
- Using discrete orientations
- Multiple binning examples

**convergence_example.py**
- Basic convergence study using Convergence class
- Shows Convergable configuration

**ensemble_example.py**
- Ensemble convergence using EnsembleConvergence class
- Multi-geometry averaging

**s11_convergence_example.py**
- S11 Mueller element convergence
- Includes plotting

**backscatter_convergence_example.py**
- Backscatter convergence on specific bins
- Demonstrates theta_indices usage

**phips_convergence_example.py**
- PHIPS detector DSCS convergence
- Custom binning for detector geometry

**phips_ensemble_convergence_example.py**
- PHIPS ensemble convergence
- Combines PHIPS detectors with ensemble averaging

## When to Use Direct API

Use the direct API when you need:
- Maximum control over GOAD settings
- Access to intermediate results
- Custom convergence logic
- Integration with existing code using these classes

## Recommended: Use Unified API Instead

For new code, we recommend using the **unified convergence API** (see `examples/unified/`):
- Simpler interface
- Automatic validation
- Consistent output format
- Easier parameter sweeps

## Running Examples

From the repository root:

```bash
cd goad-py
source .venv/bin/activate
python examples/direct/simple_example.py
```

## See Also

- **examples/unified/** - Recommended unified convergence API examples
- **UNIFIED_API.md** - Documentation for the unified API
