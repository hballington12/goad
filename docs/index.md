# Welcome to GOAD

Geometric Optics with Aperture Diffraction (GOAD) is a code for simulating light scattering from large particles. It approximates the near-field scattering for an incident plane wave for large particles, and then uses aperture diffraction theory to map the near-field to the far-field. It computes the Mueller matrix and integrated optical scattering parameters. The core is written in Rust, with bindings to Python.

## When is GOAD applicable?

You can usually use GOAD when the following conditions are met:

- The overall particle size `d` is much larger than the wavelength `λ`.
- The field of interest is in the far-field zone, ie. at a distance `r` where `r >> λ` and `r >> d`.

## Example

{{code_block('examples/multiproblem', 'multiproblem')}}

## Recommended Usage

For orientation averaging problems, it is recommended to run GOAD with the [`Convergence`](convergence.md) class.

## Contents

### User Guide
- [Settings](settings.md)
- [Results](results.md)
- [Checks](checks.md)
- [Convergence](convergence.md)
