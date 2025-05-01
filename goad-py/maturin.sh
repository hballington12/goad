#!/bin/bash
# BLENDER_PYTHON=/snap/blender/5673/4.3/python/bin/python3.11
# $BLENDER_PYTHON -m pip install --upgrade maturin
pip install --upgrade maturin

cd ..
cargo build --release
cd -

maturin develop --release
maturin build --release

source .venv/bin/activate
pip install ../target/wheels/goad_py*.whl --force-reinstall
