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

## 📚 Contents

- [GOAD - Geometric Optics with Aperture Diffraction](#goad---geometric-optics-with-aperture-diffraction)
  - [📚 Contents](#-contents)
  - [🚀 Quickstart](#-quickstart)
  - [✨ Features](#-features)
  - [🛠️ Installation](#️-installation)
  - [⚙️ Usage](#️-usage)
    - [Configuration](#configuration)
    - [Command-Line Arguments](#command-line-arguments)
  - [▶️ Running the Simulation](#️-running-the-simulation)
  - [🛠️ Testing](#️-testing)
  - [🤝 Contributing](#-contributing)
  - [📜 License](#-license)

---

## 🚀 Quickstart

1. Run the `setup.sh` script from the project root to compile the code and initialise settings:

    ```sh
    ./setup.sh
    ```

2. Execute the binary located at `./target/release/goad`:

    ```sh
    ./target/release/goad [OPTIONS]
    ```

For more information, see the [quickstart guide](https://docs.rs/goad/0.1.0/goad/_quickstart/index.html) in the docs.

---

## ✨ Features

- **Full Mueller Matrix Output**: Rigorous vector diffraction theory for computation of all Mueller matrix elements.
- **Extensive Geometry Possibilities**: GOAD is built with the flexibility to extend beyond simple convex polyhedral geometries, such as concavities, inclusions, layered media, negative refractive indices, and surrounding mediums other than air.
- **Fixed and Multiple Orientation Scattering**: Rapid computation of 2D scattering patterns in fixed orientation at arbitrary scattering angles, as well as fast orientation-averaged scattering computations for radiative transfer and remote sensing applications.

---

## 🛠️ Installation

Before building the project, ensure you have Rust's package manager, Cargo, installed.
You can install Rust and Cargo by following the instructions on the [official Rust website](https://doc.rust-lang.org/cargo/getting-started/installation.html).

On Linux and macOS:

```sh
curl https://sh.rustup.rs -sSf | sh
```

Clone the repository and build the project:

```sh
git clone git@github.com:hballington12/goad.git
cd goad
cargo build --release
```

After building, the binary will be in the `target/release` directory.

---

## ⚙️ Usage

### Configuration

The application uses a default configuration file (`config/default.toml`).
**To customise:**

- Copy it to `config/local.toml` and edit as needed.
- Options in config files are overridden by command line arguments, which are in turn overridden by environment variables.

**Example ways to set the wavelength:**

1. Edit `config/local.toml`:

    ```toml
    wavelength = 0.532
    ```

2. Use a command line argument:

    ```sh
    goad -- -wavelength 0.532
    ```

3. Use an environment variable:

    ```sh
    export GOAD_wavelength=0.532
    goad
    ```

### Command-Line Arguments

```sh
GOAD - Geometric Optics with Aperture Diffraction
Harry Ballington

Usage: goad [OPTIONS]

Options:
  -w, --w <W>
          Wavelength in units of the geometry. Should be larger than the smallest feature in the geometry
      --bp <BP>
          Minimum beam power threshold for propagation. Beams with less power than this will be truncated
      --baf <BAF>
          Minimum area factor threshold for beam propagation. The actual area threshold is wavelength² × factor. Prevents geometric optics from modeling sub-wavelength beams
      --cop <COP>
          Total power cutoff fraction (0.0-1.0). Simulation stops when this fraction of input power is accounted for. Set to 1.0 to disable and trace all beams to completion
      --rec <REC>
          Maximum recursion depth for beam tracing. Typical values: 8-15. Higher values rarely improve results when reasonable beam power thresholds are set
      --tir <TIR>
          Maximum allowed total internal reflections. Prevents infinite TIR loops by truncating beams after this many TIR events
  -g, --geo <GEO>
          Path to geometry file (.obj format). Contains all input shapes for the simulation
      --ri0 <RI0>
          Surrounding medium refractive index. Format: "re+im" (e.g., "1.3117+0.0001i")
  -r, --ri <RI>...
          Particle refractive indices, space-separated. Each shape in the geometry is assigned a refractive index. If fewer values than shapes are provided, the first value is reused
      --distortion <DISTORTION>
          Distortion factor for the geometry. Applies distortion sampled from a Gaussian distribution. Default: sigma = 0.0 (no distortion). Sigma is the standard deviation of the facet theta tilt (in radians)
      --geom-scale <GEOM_SCALE>...
          Geometry scale factors for each axis (x, y, z). Format: "x y z" (e.g., "1.0 1.0 1.0"). Default: "1.0 1.0 1.0" (no scaling)
      --uniform <UNIFORM>
          Use uniform random orientation scheme. The value specifies the number of random orientations
      --discrete <DISCRETE>...
          Use discrete orientation scheme with specified Euler angles (degrees). Format: alpha1,beta1,gamma1 alpha2,beta2,gamma2 ...
      --euler <EULER>
          Specify Euler angle convention for orientation. Valid values: XYZ, XZY, YXZ, YZX, ZXY, ZYX, etc. Default: ZYZ
      --simple <SIMPLE> <SIMPLE>
          Use simple equal-spacing binning scheme. Format: <num_theta_bins> <num_phi_bins>
      --interval
          Enable interval binning scheme with variable spacing. Allows fine binning in regions of interest like forward/backward scattering
      --theta <THETA> <THETA> <THETA>...
          Theta angle bins for interval binning (degrees). Format: start step1 mid1 step2 mid2 ... stepN end Example: 0 1 10 2 180 = 0° to 10° in 1° steps, then 10° to 180° in 2° steps
      --phi <PHI> <PHI> <PHI>...
          Phi angle bins for interval binning (degrees). Format: start step1 mid1 step2 mid2 ... stepN end Example: 0 2 180 = 0° to 180° in 2° steps
      --custom <CUSTOM>
          Path to custom binning scheme file. Contains a list of (theta, phi) bin pairs in TOML format. Overrides other binning parameters
  -s, --seed <SEED>
          Random seed for reproducibility. Omit for a randomized seed
      --dir <DIR>
          Output directory for simulation results. If not specified, a directory in the format 'run00001' will be created automatically
  -h, --help
          Print help
  -V, --version
          Print version

EXAMPLES:
    # Run with a specific wavelength and geometry file
    goad -w 0.5 --geo geometry.obj

    # Run with a specific refractive index and random orientations
    goad --ri 1.31+0.01i --uniform 100

    # Run over discrete orientations with an interval binning scheme
    goad --discrete="-30.0,20.0,1.0 -40.0,13.0,12.1" --interval \
         --theta 0 1 10 2 180 --phi 0 2 180

    # Run inside a medium other than air
    goad --ri0 1.5+0.0i

    # Run with multiple shapes with different refractive indices
    goad --ri 1.31+0.0i 1.5+0.1i --geo geometries.obj

    # Save output to a specific directory
    goad --dir /path/to/output
```

---

## ▶️ Running the Simulation

```sh
cargo run --release -- [OPTIONS]
```

or

```sh
./target/release/goad -w 0.532 --geo ./examples/data/cube.obj
```

---

## 🛠️ Testing

To run the tests:

```sh
cargo test
```

---

## 🤝 Contributing

Contributions are welcome!
Please open an issue or submit a pull request on [GitHub](https://github.com/hballington12/goad).

---

## 📜 License

This project is licensed under the GNU General Public License.
See the [LICENSE](LICENSE) file for details.
