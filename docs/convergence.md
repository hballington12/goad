# Convergence

## Basic Usage

When computing orientation averaged scattering, it's not usually known beforehand exactly how many orientations are required to converge on the desired result. GOAD's solution to this is called a `Convergence`, which uses Welford's algorithm to track the mean and variance of one or more prescribed convergence variables. The simulation runs until the convergence criteria are met, or some maximum number of orientations is reached. A simple example runs until the standard error in the mean asymmetry parameter has an error less than 2%:

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

## Convergable Parameters

The following table lists the current convergable parameters and some recommendations for starting values:

| Parameter | Recommended Value | Description |
|-----------|----------|-------------|
| `Param.Asymmetry` | `0.01` | Asymmetry parameter, the integrated cosine-weighted scattering |
| `Param.ScatCross` | `0.01` | Scattering cross section, the integrated scattering |
| `Param.ExtCross` | `0.01` | Extinction cross section, the integrated scattering + absorption |
| `Param.Albedo` | `0.01` | Single scattering albedo, the ratio of scattering cross section to extinction cross section |

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

# Access results
mean = convergence.mean  # Mean values
sem = convergence.sem    # Standard error of the mean
count = convergence.count  # Number of orientations computed
```
