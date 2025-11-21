# Checks

There are a few ways to check that the results from a GOAD solve are valid. Or rather, there are a few ways to measure the quality of results in a best case scenario, since GOAD itself is an approximate method.

## Prerequisite

GOAD is only valid in the geometric optics regime, which is when the wavelength is much smaller than the size of the geometry. The size parameter `X` is defined as `dπ/λ`, where `d` is (loosely) defined as the maximum dimension of the geometry and `λ` is the wavelength. GOAD is only valid for `X > 20`, with an accuracy that increases  with `X`.

## Energy Conservation

Your first check should be to ensure a reasonable input to output energy conservation. You can do this by inspecting the [powers](results.md#power-budget) in a goad [Results](results.md) object. Anything under 95% conservation usually indicates some issue with the solve. Possible causes include:

- The [beam thresholds](settings.md#beam-thresholds) are poorly set, or the cutoff power is too low. Generally you can inspect each component of the powers object to determine how much of the input power is lost to each truncation criteria. 
- The [recursion limits](settings.md#recursion-limits) are too low. Generally you will need at least 5 maximum recursions and total internal reflections. The recommended values are 10 each, but you may need to increase to 20 or more in some cases.
- The geometry contains extremely small faces compared to the overall geometry size. This can sometimes lead to issues with the beam clipping algorithm. A reasonable rule-of-thumb is to ensure that the smallest face has an area no less than 1% of the largest face area.

## Binning Resolution

The second check is usually to ensure that the binning resolution (especially in theta) is fine enough to resolve the forward scattering peak. In the geometric optics regime, in which GOAD is valid, the width of the forward scattering peak in radians is approximately `λ/d`. Therefore if you are computing the field starting from theta = 0°, there should be at least 3 bins until `λ/d`, otherwise the asymmetry parameter and scattering cross sections will be inaccurately computed. As an example:

{{code_block('examples/checks', 'binning_resolution')}}

## Number of Orientations

When running an orientation averaging problem, there is no easy way to know beforehand how many orientations are needed to converge. As a general rule, the number of required orientations:

- Increases with the particle size
- Decreases with particle complexity
- Decreases with number of particles in the ensemble (if using an ensemble)
- Decreases with wavelength

In the geometric optics regime where GOAD operates, you typically need somewhere in the range of 1000 to 10000 orientations to converge on the integrated scattering parameters to within 1% error, which is the point at which the real error from GOAD vs. a numerically exact solution like the DDA becomes the leading error term. The recommended way to run an orientation averaging problem where you don't need to specify the number of orientations beforehand is via a GOAD [Convergence](convergence.md).
