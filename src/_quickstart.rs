//! # Quickstart Guide
//!
//! ## Compiling the code:
//! All Rust crates are built using `cargo`, the Rust package manager. You can
//! compile this crate with:
//! ```bash
//! cargo build --release
//! ```
//! The executable binary will then be located at `target/release/goad`. You can
//! then execute this from the command line, or a shell script, such as:
//! ```bash
//! ./target/release/goad
//! ```
//!
//! ## Setting up the configuration:
//! By default, if you try to run the GOAD binary, the computation will run with
//! some default settings defined in `config/default.toml`. There are three typical
//! ways that you might want to customise the configuration to your application:
//!
//! 1. **Configuration file**: The file `config/default.toml` contains the default
//! configuration options for the code. Make a copy of this and rename it to
//! `local.toml`. You can then place a `local.toml` file in one of two locations:
//!     1. In the `config/` directory.
//!     2. In a directory where you execute the `goad` binary from.
//!
//!     GOAD will looks for a configuration file in the following order:
//!     1. A `local.toml` file in the current working directory.
//!     2. A `local.toml` file in the `config/` directory of this crate.
//!     3. The `default.toml` file in the `config/` directory of this crate.
//!
//! Once GOAD finds a `config.toml` file, other `config.toml` files are ignored.
//!
//! 2. **Command Line Arguments**: You can specify command line arguments to
//! override any values defined in a configuration file. For instance, you might
//! run goad with:
//! ```bash
//! /path/to/goad/crate/target/release/goad --uniform 10000 --geo my_particle.obj
//! ```
//! which would run a simulation with geometry given in `my_particle.obj`, averaged
//! over 10000 uniformly distributed orientations. For a complete list of command
//! line arguments, see the help information:
//! ```bash
//! goad --help
//! ```
//!
//! 3. **Environment Variables**: You can set environment variables to override
//! options specified in a configuration file **and** those specified on the
//! command line. This is done by prepending `GOAD_` to an environment variable
//! with the same name is the command line argument you wish to fix. For example:
//! ```bash
//! export GOAD_wavelength=0.532
//! ```
//! would set the wavelength of any simulations you run to 0.532.
//!
//! ## Running the code
//! Running the code is as simple as executing the binary:
//! ```bash
//! /path/to/goad/crate/target/release/goad
//! ```
//! The output, by default, is placed in a directory called `goad_run` (you can
//! set the directory name in the config file, or with the `--dir` option). GOAD
//! will overwrite any existing directory called `goad_run`. The output typically
//! contains at least some of the following files:
//! - `mueller_scatgrid` - unnormalised Mueller matrix. First 2 columns are polar
//! scattering angle theta and azimuthal scattering angle phi. The remaining
//! columns are the 16 elements of the Mueller matrix.
//! - `mueller_scatgrid_1d` - unnormalised Mueller matrix, integrated over the
//! azimuthal scattering angle. First column is polar scattering angle and the
//! remainging columns ar ethe 16 elements of the Mueller matrix.
//! - `results.dat` - contains integrated scattering parameters, energy conservation,
//! and other useful information.
