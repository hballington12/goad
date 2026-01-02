# Results

After solving a GOAD problem, access the scattering results through the `results` property of the `MultiProblem` object. For information about checking the validity of the results, see [Checks](checks.md).

## Basic Usage

{{code_block('examples/results', 'basic')}}

> GOAD organizes scattering results into [zones](#zones). Properties like `results.mueller`, `results.bins`, and `results.scat_cross` are convenience accessors that return data from the full zone. For zone-specific access (e.g., backscatter parameters), see the [Zones](#zones) section.

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

GOAD computes several integrated optical parameters from the Mueller matrix. These parameters are available directly on the results object for convenience, and are also accessible through zone-specific `params` dictionaries (see [Zones](#zones) for more details).

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

## Zones

GOAD organizes scattering results into zones. For configuration details, see the [Zones section](settings.md#zones) in Settings.

### Accessing Zones

Access individual zones from the results:

{{code_block('examples/results', 'zones')}}

### Zone-Specific Parameters

Each zone type provides different parameters:

{{code_block('examples/results', 'zone_params')}}

**Full zone** (requires full theta range 0-180°):

- `scatt_cross`: Scattering cross section
- `ext_cross`: Extinction cross section
- `asymmetry`: Asymmetry parameter g
- `albedo`: Single scattering albedo

**Backward zone** (θ = 180°):

- `backscatter_cross`: Backscattering cross section
- `lidar_ratio`: Extinction-to-backscatter ratio
- `depolarization_ratio`: Linear depolarization ratio

**Forward zone** (θ ≈ 0°):

- `ext_cross_optical_theorem`: Extinction cross section via optical theorem

### Zone Mueller Matrices

Each zone contains its own Mueller matrix data:

{{code_block('examples/results', 'zone_mueller')}}

## Power Budget

Track energy conservation throughout the simulation:

{{code_block('examples/results', 'powers')}}

The power dictionary contains:

- `input`: Incident beam power ie. the mean geometrical cross section
- `output`: Total near-field output power, the sum of scattering and absorption
- `absorbed`: Power absorbed by the particle
- `trnc_ref`: Power lost due to max total internal reflections reached
- `trnc_rec`: Power lost due to max recursions reached
- `trnc_clip`: Power lost during beam clipping
- `trnc_energy`: Power lost due to minimum energy cutoff
- `trnc_area`: Power lost due to minimum area cutoff
- `trnc_cop`: Power lost due to total cutoff power threshold
- `clip_err`: Error from clipping algorithm
- `ext_diff`: External diffraction power
- `missing`: Total unaccounted power

## Complete Example

{{code_block('examples/results', 'complete')}}

## Result Properties Reference

| Property | Type | Description |
|----------|------|-------------|
| [`bins`](#angular-bins) | `ndarray` | 2D angular bins `[[theta, phi], ...]` |
| [`bins_1d`](#angular-bins) | `ndarray` | 1D theta bins |
| [`mueller`](#2d-mueller-matrix) | `ndarray` | 2D Mueller matrix (total) |
| [`mueller_beam`](#mueller-matrix-components) | `ndarray` | 2D Mueller matrix (beam component) |
| [`mueller_ext`](#mueller-matrix-components) | `ndarray` | 2D Mueller matrix (external diffraction) |
| [`mueller_1d`](#1d-mueller-matrix) | `ndarray` | 1D Mueller matrix (total) |
| [`mueller_1d_beam`](#mueller-matrix-components) | `ndarray` | 1D Mueller matrix (beam component) |
| [`mueller_1d_ext`](#mueller-matrix-components) | `ndarray` | 1D Mueller matrix (external diffraction) |
| [`scat_cross`](#scattering-cross-section) | `float` | Scattering cross section |
| [`ext_cross`](#extinction-cross-section) | `float` | Extinction cross section |
| [`asymmetry`](#asymmetry-parameter) | `float` | Asymmetry parameter g |
| [`albedo`](#single-scattering-albedo) | `float` | Single scattering albedo |
| [`powers`](#power-budget) | `dict[str, float]` | Power budget dictionary |
| [`zones`](#zones) | `Zones` | Collection of all zones |
| [`full_zone`](#accessing-zones) | `Zone` | Full scattering zone (0-180°) |
| [`forward_zone`](#accessing-zones) | `Zone` | Forward scattering zone (θ ≈ 0°) |
| [`backward_zone`](#accessing-zones) | `Zone` | Backward scattering zone (θ = 180°) |

## Zone Properties Reference

| Property | Type | Description |
|----------|------|-------------|
| `name` | `str` | Zone name (e.g., "full", "forward", "backward") |
| `zone_type` | `ZoneType` | Zone type enum (`Full`, `Forward`, `Backward`, `Custom`) |
| `num_bins` | `int` | Number of angular bins in this zone |
| `bins` | `ndarray` | 2D angular bins for this zone |
| `bins_1d` | `ndarray` | 1D theta bins for this zone |
| `mueller` | `ndarray` | 2D Mueller matrix for this zone |
| `mueller_1d` | `ndarray` | 1D Mueller matrix for this zone |
| `params` | `dict[str, float]` | Zone-specific integrated parameters |
| `label` | `str \| None` | Optional user-defined label |
