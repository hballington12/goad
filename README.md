# GOAD - Geometric Optics with Aperture Diffraction

GOAD is a Rust-based application for simulating geometric optics with aperture diffraction. It computes the 2D mueller matrix by using geometric optics and a polygon clipping algorithm to approximate the electric field on the particle surface. The surface field is then mapped to the far-field on the basis of the electromagnetic equivalence theorem, which takes the form of a vector surface integral diffraction equation. Green's theorem is used to reduce the surface integral to a line integral around the contours of outgoing beam cross sections, which is important for computational efficiency.

## Features

- **Full Mueller Matrix Output**: Rigorous vector diffraction theory for computation of all mueller matrix elements.
- **Extensive Geometry Possibilities**: GOAD is built with the flexibility to extend beyond simple convex polyhedral goemetries, such as concavities, inclusions, layered media, negative refractive indices, and surrounding mediums other than air.
- **Fixed and Multiple Orientation Scattering**: Rapid computation of 2D scattering patterns in fixed orientation at arbitrary scattering angles, as well as fast orientation-averaged scattering computations for radiative transfer and remote sensing application.

## Installation

To install and run the project, ensure you have Rust installed. Clone the repository and build the project using Cargo:

```sh
git clone <repository-url>
cd pbt
cargo build --release
```

## Usage

### Configuration

The application uses a configuration file (`config/default`) to set various simulation parameters. You can override these parameters using a local configuration file (`config/local`) or through environment variables.

### Command-Line Arguments

You can also override configuration values using command-line arguments. Here are some of the available options:

```sh
USAGE:
    pbt [OPTIONS]

OPTIONS:
    -w, --wavelength <WAVELENGTH>        Wavelength in units of the geometry.
        --bp <BEAM_POWER_THRESHOLD>      Minimum absolute beam power threshold for new beams to propagate.
        --baf <BEAM_AREA_THRESHOLD_FAC>  Minimum area factor for new beams to propagate.
        --cop <TOTAL_POWER_CUTOFF>       Cutoff power for beam propagation.
        --rec <MAX_REC>                  Maximum number of recursions before a beam is truncated.
        --tir <MAX_TIR>                  Maximum number of total internal reflections before a beam is truncated.
    -g, --geo <GEOMETRY_FILE>            File path to the input geometry.
        --ri0 <MEDIUM_REFRACTIVE_INDEX>  Refractive index of the surrounding medium.
    -r, --ri <PARTICLE_REFRACTIVE_INDEX> Refractive index of the particles.
        --orient <ORIENTATION_SCHEME>    Orientation scheme for the simulation.
    -s, --seed <SEED>                    Random seed for the simulation.
```

### Running the Simulation

To run the simulation, use the following command:

```sh
cargo run --release -- [OPTIONS]
```

### Example

```sh
cargo run --release -- -w 0.5 --bp 1e-3 --geo ./examples/data/geometry.obj
```

## Testing

To run the tests, use the following command:

```sh
cargo test
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
