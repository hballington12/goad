# GOAD - Geometric Optics with Aperture Diffraction

GOAD is a Rust-based physical-optics hybrid light scattering model based on geometric optics with aperture diffraction. It computes the 2D mueller matrix by using geometric optics and a polygon clipping algorithm to approximate the electric field on the particle surface. The surface field is then mapped to the far-field on the basis of the electromagnetic equivalence theorem, which takes the form of a vector surface integral diffraction equation. Green's theorem is used to reduce the surface integral to a line integral around the contours of outgoing beam cross sections, which is important for computational efficiency.

## Quickstart

- Simply run the `setup.sh` script from the project root to compile the code and initialise some settings.
- Then, you can execute the binary which should be located at `./target/release/goad`

## Features

- **Full Mueller Matrix Output**: Rigorous vector diffraction theory for computation of all mueller matrix elements.
- **Extensive Geometry Possibilities**: GOAD is built with the flexibility to extend beyond simple convex polyhedral goemetries, such as concavities, inclusions, layered media, negative refractive indices, and surrounding mediums other than air.
- **Fixed and Multiple Orientation Scattering**: Rapid computation of 2D scattering patterns in fixed orientation at arbitrary scattering angles, as well as fast orientation-averaged scattering computations for radiative transfer and remote sensing application.

## Installation

Before building the project ensure you have Rust's package manager, Cargo installed. You can install Rust and Cargo by following the instructions on the [official Rust website](https://doc.rust-lang.org/cargo/getting-started/installation.html).

On Linux and macOS, you can usually install Cargo using:

```sh
curl https://sh.rustup.rs -sSf | sh
```

Once Cargo is installed, you clone you can clone the repository and build the project using Cargo:

```sh
git clone git@github.com:hballington12/goad.git
cd goad
cargo build --release
```

After building the project, the binary will be located in the `target/release` directory. You can run it directly from there or use Cargo to manage the execution.

```sh
./target/release/goad [OPTIONS]
```

## Usage

### Configuration

The application uses a default configuration file (`config/default.toml`) to set various simulation parameters. Users should make a local copy of this file called `config/local.toml` and customise it as needed. Options set in these config files are overridden by any command line arguments, which can be viewed with:

```sh
goad -- --help
```

Command line arguments are overridden by environment variables, which are prefixed by `GOAD_`. For example, to set the wavelength of incident light, the user can either:

1. Modify the wavelength entry in `config/local.toml`

    ```sh
    # Wavelength of the light source in dimensions of the geometry file
    wavelength = 0.532

    # ... other configuration options ...
    ```

2. Set the wavelength using a command line argument:

    ```sh
    goad -- -wavelength 0.532
    ```

3. Set the wavelength using an environment variable:

    ```sh
    export GOAD_wavelength=0.532
    goad
    ```

### Command-Line Arguments

Configuration values in the config files can be overridden using command-line arguments. Custom far-field binning can be set in the config file. Here are some of the available command line arguments:

```sh
GOAD - Geometric Optics with Aperture Diffraction

Usage: goad [OPTIONS] [COMMAND]

Commands:
  uniform   Solve the problem by averaging over a uniform distribution of angles. Example: `uniform 100`
  discrete  Solve the problem by averaging over a discrete set of angles (in degrees). Example: `discrete 0,0,0 20,30,40`
  help      Print this message or the help of the given subcommand(s)

Options:
  -w, --w <W>        Wavelength in units of the geometry
      --bp <BP>      Minimum absolute beam power threshold for new beams to propagate
      --baf <BAF>    Minimum area factor for new beams to propagate. The actual area threshold is calculated as `wavelength^2 * factor`
      --cop <COP>    Cutoff power. The total acceptable output power per orientation before beam propagation is terminated. Once this threshold is reached, the near-field simulation will stop
      --rec <REC>    The maximum number of recursions before a beam is truncated
      --tir <TIR>    The maximum number of total internal reflections before a beam is truncated
  -g, --geo <GEO>    File path to the input geometry. All input shapes should be defined in this file. Currently, only the Wavefront .obj format is supported
      --ri0 <RI0>    The refractive index of the surrounding medium
  -r, --ri <RI>...   The refractive index of the particle/s, separated by spaces. If multiple values are provided, each shape in the geometry will be assigned a refractive index. If fewer values are provided than the number of shapes, the first value will be used for the remaining shapes
  -s, --seed <SEED>  Random seed for the simulation
  -h, --help         Print help
  -V, --version      Print version
```

### Running the Simulation

To run the simulation, use the following command:

```sh
cargo run --release -- [OPTIONS]
```

### Example

```sh
cargo run --release -- -w 0.532 --geo ./examples/data/cube.obj
```

Or, equivalently:

```sh
./target/release/goad -w 0.532 --geo ./examples/data/cube.obj
```

## Testing

To check that everything is running as it should be on your system, you can run the tests with the following command:

```sh
cargo test
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This project is licensed under the GNU General Public License. See the [LICENSE](LICENSE) file for details.
