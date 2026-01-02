# Settings

The `Settings` object configures a GOAD simulation. It controls physical parameters, numerical methods, and output options.

## Basic Usage

At a minimum, you must specify the path to a geometry file or directory containing geometry files:

{{code_block('examples/settings', 'basic')}}

The geometry defines the units of the problem. If your geometry file is in microns, then you should also specify the wavelength in microns. All faces in the geometry must be planar and have some non-zero area. GOAD will return with an error if there are faces with zero area (ie. extremely thin triangles), since it needs to compute normals of each face by a cross product of 2 non-colinear edge vectors. You can make geometries in the open-source [Blender](https://www.blender.org/) software, or use some example geometries straight from Python [here](https://github.com/hballington12/bpy-geometries).

If you specify a directory, GOAD will attempt to load all files with the `.obj` extension in the directory. It will then choose geometries at random for each orientation in the simulation. See [Orientation Distribution](#orientation-distribution) for more details.

## Physical Parameters

### Wavelength

The wavelength of incident light in micrometers:

{{code_block('examples/settings', 'wavelength')}}

**Default:** `0.532` (532 nm, green laser)

### Refractive Indices

Specify the complex refractive index for both the particle and surrounding medium:

{{code_block('examples/settings', 'refractive')}}

**Defaults:**

- `particle_refr_index_re`: `1.31` (typical glass)
- `particle_refr_index_im`: `0.0`
- `medium_refr_index_re`: `1.0` (vacuum/air)
- `medium_refr_index_im`: `0.0`

### Particle Scaling

Scales the entire problem, including geometry and wavelength. Does not change the physics, only used for improving clipping algorithm accuracy. The default value is usually sufficient.

**Parameter:** `scale`  
**Default:** `1.0`

## Orientation Distribution

Define how particle orientations are sampled:

{{code_block('examples/settings', 'orientation')}}

**Default:** `Orientation.uniform(num_orients=1)` with `EulerConvention('ZYZ')`

For discrete orientations:

{{code_block('examples/settings', 'orientation_discrete')}}

The output will be an average over the orientations. Results from individual orientations are not stored. See [Results](results.md) for more details.

### Seed

Seeds the random number generator for reproducibility.

**Parameter:** `seed`  
**Default:** `None`

## Zones

Zones define a set of query points to evaluate the far-field scattering at. By default, GOAD creates three zones:

1. **Full zone** - The full scattering sphere from θ=0° to θ=180°. Used for computing integrated parameters like asymmetry, scattering cross-section, and albedo.
2. **Forward zone** - A single point at θ=0.01° (slightly off-axis for numerical stability). Used for computing extinction cross-section via the optical theorem.
3. **Backward zone** - A single point at θ=180°. Used for computing backscatter cross-section, lidar ratio, and depolarization ratio.

You can add additional zones as needed for your application. For backscattering-only applications, you can exclude the full zone and compute only at forward and backward scattering for faster computations.

>Note: GOAD automatically determines the zone type based on the theta range of your binning scheme. If the theta range covers 0° to 180°, it is classified as a Full zone and integrated parameters (asymmetry, scattering cross-section, etc.) will be computed. If the range is partial, it is classified as a Custom zone and these parameters will not be computed.

Different parameters are computed depending on the zone type - see [Integrated Parameters](results.md#integrated-parameters) for more details.

Each zone is specified by a binning scheme and an optional label. See [Binning Schemes](#binning-schemes) for more info about binning schemes.

{{code_block('examples/settings', 'zones')}}

**Default:** A single full zone with interval binning (high resolution at forward and backward angles).

## Binning Schemes

Control the angular resolution of scattering calculation in the far-field. As particle size increases, the width of peaks in the scattering decreases. GOAD currently requires the user to choose a sufficiently fine binning scheme to resolve the peaks. For `phi` angles, a relatively course binning scheme can be used, but for `theta` angles, care should be taken to ensure sufficient resolution, otherwise the integrated parameters lose accuracy. The compute time approximately scales with the number of bins, which for `simple` and `interval` binning schemes is just the product of the number of bins in each dimension.

### Simple Binning

Uniform spacing in theta and phi:

{{code_block('examples/settings', 'binning')}}

**Default:** `BinningScheme.interval(thetas=[0, 5, 175, 179, 180], theta_spacings=[0.1, 2.0, 0.5, 0.1], phis=[0, 360], phi_spacings=[7.5])`

### Interval Binning

Variable resolution for different angular regions:

{{code_block('examples/settings', 'binning_interval')}}

This example uses 1° resolution for forward scattering (0-90°) and 2° for backward scattering (90-180°). `interval` binning schemes are useful if the user is only interested in a specific angular scattering range, eg 6°-25°.

### Custom Binning

Specify arbitrary bin edges:

{{code_block('examples/settings', 'binning_custom')}}

GOAD will compute the bin centres automatically. GOAD will not compute the 1D mueller matrix or integrated parameters if using a `custom` binning scheme.

## Mapping Method

Choose how near-field results map to the far-field:

{{code_block('examples/settings', 'mapping')}}

**Options:**

- `'ad'`: Aperture Diffraction (default, more accurate). Suitable for fixed-orientation computations.
- `'go'`: Geometric Optics (faster, suitable for very large particles). Generally not suitable for fixed-orientation computations.

**Default:** `Mapping('ad')`

## Beam Tracing Parameters

Beams traced in the near-field if they pass the following checks, in order:

- Beam power is below the threshold
- Beam area is below the threshold
- Beam number of recursion is above the threshold

If the beam is from a total internal reflection event, the beam is traced even if the recursion is above the threshold, as long as the total internal reflection count is below the threshold.

### Beam Thresholds

Control when beams are truncated during ray tracing:

{{code_block('examples/settings', 'thresholds')}}

**Defaults:**

- `beam_power_threshold`: `0.005` (discard beams below 0.5% of incident power).
- `beam_area_threshold_fac`: `0.1` (factor × λ² determines the physical area threshold, below which beams are discarded. It scales with λ² following the applicability of geometric optics).
- `cutoff`: `0.99` (trace 99% of energy in the near field, then map. You generally want to use a value of at least 0.95 here, unless you have a good reason to do otherwise).

Lower thresholds and higher cutoff increases accuracy but slows computation.

### Recursion Limits

Limit internal beams bounces:

{{code_block('examples/settings', 'recursion')}}

**Defaults:**

- `max_rec`: `10` (maximum internal reflections)
- `max_tir`: `10` (maximum total internal reflections)

Increase these for complex internal ray paths, but expect slower performance. High numbers of total internal reflections recommended for backscattering computations, ie. `20`.

## Output Options

### Directory

Specify output directory for simulation data:

**Parameter:** `directory`  
**Default:** `"goad-run"` 

### Coherence

Enable coherent beam addition (phase tracking):

**Parameter:** `coherence`  
**Default:** `True`

Enables coherent beam addition with phase tracking for interference effects. If coherence is enabled, GOAD traces the phase of each beam, **and** combines amplitude matrices in the far-field with interference. If coherence is disabled, GOAD traces the phase of each beam, but combines the far-field contributions of beams in the far-field only by a linear summation of Mueller matrices. Coherence should be enabled for backscattering computations.

### Verbosity

Suppress console output:

**Parameter:** `quiet`  
**Default:** `False`

Set to `True` to silence progress messages.

## Complete Example

{{code_block('examples/settings', 'advanced')}}

## Parameter Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| [`geom_path`](#basic-usage) | `str` | **Required** | Path to geometry file |
| [`wavelength`](#wavelength) | `float` | `0.532` | Wavelength in geometry units |
| [`particle_refr_index_re`](#refractive-indices) | `float` | `1.31` | Real part of particle refractive index |
| [`particle_refr_index_im`](#refractive-indices) | `float` | `0.0` | Imaginary part of particle refractive index |
| [`medium_refr_index_re`](#refractive-indices) | `float` | `1.0` | Real part of medium refractive index |
| [`medium_refr_index_im`](#refractive-indices) | `float` | `0.0` | Imaginary part of medium refractive index |
| [`orientation`](#orientation-distribution) | `Orientation` | `Orientation.uniform(1)` | Orientation distribution |
| [`zones`](#zones) | `list[ZoneConfig]` | Single full zone with interval binning | Zone configurations for far-field evaluation |
| [`mapping`](#mapping-method) | `Mapping` | `Mapping('ad')` | Near-to-far field mapping method |
| [`beam_power_threshold`](#beam-thresholds) | `float` | `0.005` | Beam power truncation threshold |
| [`beam_area_threshold_fac`](#beam-thresholds) | `float` | `0.1` | Beam area truncation factor |
| [`cutoff`](#beam-thresholds) | `float` | `0.99` | Energy tracking cutoff |
| [`max_rec`](#recursion-limits) | `int` | `10` | Maximum internal reflections |
| [`max_tir`](#recursion-limits) | `int` | `10` | Maximum total internal reflections |
| [`scale`](#particle-scaling) | `float` | `1.0` | Geometry scaling factor |
| [`seed`](#seed) | `int` | `None` | Seed for random number generator |
| [`directory`](#directory) | `str` | `"goad_run"` | Output directory path |
| [`coherence`](#coherence) | `bool` | `True` | Enable coherent beam addition |
| [`quiet`](#verbosity) | `bool` | `False` | Suppress console output |
