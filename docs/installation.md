# Installation

GOAD is hosted on the Python Package Index (PyPI), at [https://pypi.org/project/goad-py/](https://pypi.org/project/goad-py/).

## Setup

=== "venv"

    ```bash
    python3.11 -m venv .venv
    source .venv/bin/activate
    pip install goad-py
    ```

=== "conda"

    ```bash
    conda create -n goad-backscatter python=3.11 -y
    conda activate goad-backscatter
    pip install goad-py
    ```

=== "uv"

    ```bash
    uv venv --python 3.11
    source .venv/bin/activate
    uv add goad-py
    ```
