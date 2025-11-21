# Welcome to GOAD

Geometric Optics with Aperture Diffraction (GOAD) is a code for simulating light scattering from large particles. It approximates the near-field scattering for an incident plane wave for large particles, and then uses aperture diffraction theory to map the near-field to the far-field. It computes the Mueller matrix and integrated optical scattering parameters. The core is written in Rust, with bindings to Python.

## Quickstart

To get started with an orientation averaging problem, see the [`Convergence`](convergence.md) class, which is the recommended way to run GOAD. Then head over to [`Settings`](settings.md) for more config options, and then to [`Results`](results.md) for more information on the output.

## When is GOAD applicable?

You can usually use GOAD when the following conditions are met:

- The overall particle size `d` is much larger than the wavelength `λ`.
- The field of interest is in the far-field zone, ie. at a distance `r` where `r >> λ` and `r >> d`.

## Example

{{code_block('examples/multiproblem', 'multiproblem')}}

## Contents

### User Guide
- [Settings](settings.md)
- [Results](results.md)
- [Checks](checks.md)
- [Convergence](convergence.md)
