# Convergence

## Basic Usage

When computing orientation averaged scattering, it's not usually known beforehand exactly how many orientations are required to converge on the desired result. GOAD's solution to this is called a `Convergence`, which uses Westford's algorithm to track the mean and variance of one or more prescribed convergence variables. The simulation runs until the convergence criteria are met, or some maximum number of orientations is reached. A simple example runs until the standard error in the mean asymmetry parameter has an error less than 2%:

{{code_block('examples/convergence', 'basic')}}

which produces the following output:

```bash
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GOAD: [Convergence]  ▰▱▱▱▱▱▱

[Orientations: 158 (max 10000)] [0.010 sec/orientation] [Minimum orientations check: ✓]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Asymmetry       0.7726 ± 0.0154 ━━━━━━━━━━━━━━━━━━━━━━━━━ 100% [SEM: 2.00% / 2.00%]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

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
| [`Asymmetry`](#basic-usage) | `Relative: 0.01` | Asymmetry parameter, the integrated cosine-weighted scattering |
| [`ScattCross`](#basic-usage) | `Relative: 0.01` | Scattering cross section, the integrated scattering |
| [`ExtCross`](#basic-usage) | `Relative: 0.01` | Extinction cross section, the integrated scattering + absorption |
| [`Albedo`](#basic-usage) | `Relative: 0.01` | Single scattering albedo, the ratio of scattering cross section to extinction cross section |
