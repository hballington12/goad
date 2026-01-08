<div align="center">

<!-- badges: start -->
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![minimum rustc 1.85](https://img.shields.io/badge/rustc-1.85+-red.svg)
[![Rust](https://github.com/hballington12/goad/actions/workflows/rust.yml/badge.svg)](https://github.com/hballington12/goad/actions/workflows/rust.yml)
<!-- badges: end -->

</div>

# GOAD - Geometric Optics with Aperture Diffraction

GOAD is a Rust-based physical-optics hybrid light scattering model based on geometric optics with aperture diffraction. It computes the 2D Mueller matrix by using geometric optics and a polygon clipping algorithm to compute the electric field on the particle surface. The surface field is then mapped to the far-field on the basis of the electromagnetic equivalence theorem, which takes the form of a vector surface integral diffraction equation. Green's theorem is used to reduce the surface integral to a line integral around the contours of outgoing beam cross sections, which leads to fast computations compared to some other methods. Compared to the [PBT](https://github.com/hballington12/pbt) method, GOAD uses a beam clipping algorithm instead of ray backtracing on a meshed geometry, which makes the computation more accurate and faster if the particle has smooth planar surfaces.

<div align="center">

> **📖 Reference Paper**
> If you use this code in your work, please cite:
> [A Light Scattering Model for Large Particles with Surface Roughness](https://doi.org/10.1016/j.jqsrt.2024.109054)
> *H. Ballington, E. Hesse*
> [JQSRT, 2024](https://www.journals.elsevier.com/journal-of-quantitative-spectroscopy-and-radiative-transfer)

</div>

---

For documentation, installation instructions, and usage examples, see the [GOAD Wiki](https://hballington12.github.io/goad/).
