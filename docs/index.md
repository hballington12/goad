# Welcome to GOAD

Geometric Optics with Aperture Diffraction (GOAD) is a code for simulating light scattering from large particles. It approximates the near-field scattering for an incident plane wave for large particles, and then uses aperture diffraction theory to map the near-field to the far-field. It computes the Mueller matrix and integrated optical scattering parameters. The core is written in Rust, with bindings to Python.

## Example

{{code_block('home/example','example',['scan_csv','filter','group_by','collect'])}}

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
