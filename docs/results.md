# Results

After solving a GOAD problem, access the scattering results through the `results` property of the `MultiProblem` object. For information about checking the validity of the results, see [Checks](checks.md).

## Basic Usage

{{code_block('examples/results', 'basic')}}

## Mueller Matrices

The Mueller matrix describes the transformation of the Stokes vector during scattering. GOAD provides Mueller matrices in different forms and for different scattering components.

### 2D Mueller Matrix

Full angular distribution over theta and phi:

{{code_block('examples/results', 'mueller_2d')}}

The Mueller matrix is returned as a list of 16-element lists, where each element corresponds to `[s11, s12, s13, s14, s21, s22, s23, s24, s31, s32, s33, s34, s41, s42, s43, s44]` for each bin.

### 1D Mueller Matrix

Phi-integrated Mueller matrix (theta only):

{{code_block('examples/results', 'mueller_1d')}}

The 1D Mueller matrix integrates over all phi angles at each theta, providing an azimuthally-averaged scattering pattern.

### Mueller Matrix Components

GOAD separates scattering into beam and external diffraction components:

{{code_block('examples/results', 'mueller_components')}}

- **Beam component**: Direct scattering from ray tracing
- **External diffraction**: Diffraction around the particle exterior
- **Total**: Sum of beam and external diffraction components. If [coherence](settings.md#coherence) is disabled, the total is just the linear sum of beam and external diffraction components.

## Angular Bins

Access the angular coordinates for each bin:

{{code_block('examples/results', 'bins')}}

The bins correspond to the center values of each angular bin in the scattering calculation. It is also possible to directly access the theta and phi values from the binning scheme, which is useful if you need to access the bins without running the simulation.

## Integrated Parameters

GOAD computes several integrated optical parameters from the Mueller matrix if the [binning scheme](settings.md#angular-binning) is `simple` or `interval`.

### Scattering Cross Section

The total scattering cross section:

{{code_block('examples/results', 'scat_cross')}}

### Extinction Cross Section

The total extinction cross section (scattering + absorption):

{{code_block('examples/results', 'ext_cross')}}

### Asymmetry Parameter

The asymmetry parameter `g` (average cosine of scattering angle):

{{code_block('examples/results', 'asymmetry')}}

Values range from -1 (complete backscattering) to +1 (complete forward scattering).

### Single Scattering Albedo

The ratio of scattering to extinction:

{{code_block('examples/results', 'albedo')}}

Values range from 0 (pure absorption) to 1 (pure scattering).

## Power Budget

Track energy conservation throughout the simulation:

{{code_block('examples/results', 'powers')}}

The power dictionary contains:

- `input`: Incident beam power
- `output`: Total scattered power
- `absorbed`: Power absorbed by the particle
- `trnc_ref`: Power lost to reflection truncation
- `trnc_rec`: Power lost to recursion limit
- `trnc_clip`: Power lost to clipping
- `trnc_energy`: Power lost to energy cutoff
- `trnc_area`: Power lost to area threshold
- `trnc_cop`: Power lost to coplanar threshold
- `clip_err`: Error from clipping algorithm
- `ext_diff`: External diffraction power
- `missing`: Total unaccounted power

## Complete Example

{{code_block('examples/results', 'complete')}}

## Result Properties Reference

| Property | Type | Description |
|----------|------|-------------|
| [`bins`](#angular-bins) | `list[list[float]]` | 2D angular bins `[[theta, phi], ...]` |
| [`bins_1d`](#angular-bins) | `list[float] \| None` | 1D theta bins |
| [`mueller`](#2d-mueller-matrix) | `list[list[float]]` | 2D Mueller matrix (total) |
| [`mueller_beam`](#mueller-matrix-components) | `list[list[float]]` | 2D Mueller matrix (beam component) |
| [`mueller_ext`](#mueller-matrix-components) | `list[list[float]]` | 2D Mueller matrix (external diffraction) |
| [`mueller_1d`](#1d-mueller-matrix) | `list[list[float]]` | 1D Mueller matrix (total) |
| [`mueller_1d_beam`](#mueller-matrix-components) | `list[list[float]]` | 1D Mueller matrix (beam component) |
| [`mueller_1d_ext`](#mueller-matrix-components) | `list[list[float]]` | 1D Mueller matrix (external diffraction) |
| [`scat_cross`](#scattering-cross-section) | `float \| None` | Scattering cross section |
| [`ext_cross`](#extinction-cross-section) | `float \| None` | Extinction cross section |
| [`asymmetry`](#asymmetry-parameter) | `float \| None` | Asymmetry parameter g |
| [`albedo`](#single-scattering-albedo) | `float \| None` | Single scattering albedo |
| [`powers`](#power-budget) | `dict[str, float]` | Power budget dictionary |
