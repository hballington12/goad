# Convergence

## Basic Usage

When computing orientation averaged scattering, it's not usually known beforehand exactly how many orientations are required to converge on the desired result. GOAD's solution to this is called a `Convergence`, which uses Welford's algorithm to track the mean and variance of one or more prescribed convergence variables. The simulation runs until the convergence criteria are met, or some maximum number of orientations is reached. A simple example runs until the standard error in the mean asymmetry parameter has an error less than 2%:

> **Tip:** GOAD simulations can be compute-intensive. Always call `save()` after `solve()` to write results to disk for later analysis, avoiding the need to re-run expensive simulations.

{{code_block('examples/convergence', 'basic')}}

which produces the following output:

```console
⠋ GOAD: [Convergence]  [Elapsed: 2s]  [Status: RUNNING]  [2025-01-15 14:32:18]
  [Orientations: 158 (100|10000)] [0.010 sec/orientation]
  Asymmetry  7.7260e-01 ± 1.5400e-02 [ 1.99% /  2.00%] [████████████████████] 100%
Converged after 158 orientations
```

## Accessing Results

A GOAD `Convergence` class uses [Welford's algorithm](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm) to track the mean and variance of all scattering properties across each orientation. The example below shows how to access the mean results, and their corresponding errors:

{{code_block('examples/convergence', 'results')}}

which produces the following output:

```console
Asymmetry: 0.8338 +/- 0.0083
Scattering Cross Section: 172.7278
Extinction Cross Section: 207.7922
Absorption Cross Section: 35.0644
Single Scattering Albedo: 0.8313
Theta bins:
[[5.000e-02]
 [1.500e-01]
 [2.500e-01]
 ...
 [1.798e+02]
 [1.799e+02]
 [1.799e+02]]
[Theta, Phi] bins:
[[5.000e-02 3.750e+00]
 [5.000e-02 1.125e+01]
 [5.000e-02 1.875e+01]
 ...
 [1.799e+02 3.412e+02]
 [1.799e+02 3.488e+02]
 [1.799e+02 3.562e+02]]
Mueller matrix S11: [3.696e+07 3.667e+07 3.612e+07 ... 3.594e+03 3.607e+03 3.614e+03]
Mueller matrix S12: [ 17.73   -6.18  -12.768 ...   5.655   2.052   0.259]
```

It is important to note that the error here is only a best-case scenario estimate. It is the estimated error due to the Monte-Carlo orientation sampling. GOAD itself is an approximate method - the error in asymmetry parameter at size 60 is typically ~1% compared to more accurate methods like the discrete dipole approximation. For this reason, it doesn't make much sense to converge beyond a relative error of 0.1%. True error decreases with size, so you might want to converge to smaller thresholds then.

## Multiple Targets

It is possible to set multiple targets to converge on. The convergence will then run until all targets have converged. The following example runs until 2% error in the asymmetry parameter and 2% error in the extinction cross section for a particle with a modified imaginary part of the refractive index (see the [`Settings`](settings.md) class for full details on configuration options):

{{code_block('examples/convergence', 'multiple')}}

## Other Examples

### Single scattering albedo:

{{code_block('examples/convergence', 'albedo')}}

Albedo is of course just equal to 1 for non-absorbing particles, so it is not a useful parameter to converge on in those cases.

### Extinction Cross Section

{{code_block('examples/convergence', 'extcross')}}

### Backscatter Properties

For lidar applications, you may want to converge on backscatter properties. These are computed at exactly θ=180° using aperture diffraction. Note that backscatter properties tend to be more variable than integrated quantities, so a more lenient convergence threshold (e.g., 5%) is recommended:

```python
from goad import Convergence, Param, Settings

settings = Settings(geom_path="hexcol_2x30")
conv = Convergence(settings)
conv.add_target(Param.LidarRatio, 0.05)  # 5% relative error
conv.add_target(Param.DepolarizationRatio, 0.05)
conv.solve()

mean = conv.mean
print(f"Lidar Ratio: {mean.params.lidar_ratio}")
print(f"Depolarization Ratio: {mean.params.depolarization_ratio}")
print(f"Backscatter Cross Section: {mean.params.backscatter_cross}")
```

## Convergable Parameters

The following table lists the current convergable parameters and some recommendations for starting values:

| Parameter | Recommended Value | Description |
|-----------|----------|-------------|
| `Param.Asymmetry` | `0.01` | Asymmetry parameter, the integrated cosine-weighted scattering |
| `Param.ScatCross` | `0.01` | Scattering cross section, the integrated scattering |
| `Param.ExtCross` | `0.01` | Extinction cross section, the integrated scattering + absorption |
| `Param.Albedo` | `0.01` | Single scattering albedo, the ratio of scattering cross section to extinction cross section |
| `Param.BackscatterCross` | `0.05` | Backscatter cross section at θ=180°, computed via aperture diffraction |
| `Param.LidarRatio` | `0.05` | Lidar ratio, the ratio of extinction to backscatter cross section |
| `Param.DepolarizationRatio` | `0.05` | Linear depolarization ratio at backscatter, (S11-S22)/(S11+S22) |

## Saving Results

GOAD simulations can be expensive. To avoid re-running simulations when you want to process results differently, always save your results to disk:

```python
from goad import Convergence, Param, Settings

settings = Settings(geom_path="path/to/geometry.obj")
convergence = Convergence(settings)
convergence.add_target(Param.Asymmetry, 0.02)
convergence.solve()

# Save results to a directory
convergence.save("my_simulation_results")

# Or use default directory from settings
convergence.save()
```

This writes Mueller matrices, parameters, and other output files that can be loaded and analyzed later without re-running the simulation.

## Python API Reference

```python
from goad import Convergence, Param, Settings

# Create convergence solver
settings = Settings(geom_path="path/to/geometry.obj")
convergence = Convergence(settings)

# Add convergence targets (relative error thresholds)
convergence.add_target(Param.Asymmetry, 0.02)  # 2% relative SEM
convergence.add_target(Param.ScatCross, 0.01)  # 1% relative SEM

# Optional: set max orientations (default 100,000)
convergence.max_orientations = 5000

# Solve (supports Ctrl-C interruption)
convergence.solve()

# Save results to disk
convergence.save("output_directory")

# Access results
mean = convergence.mean  # Mean values
sem = convergence.sem    # Standard error of the mean
count = convergence.count  # Number of orientations computed
```
