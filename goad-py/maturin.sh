BLENDER_PYTHON=/snap/blender/5673/4.3/python/bin/python3.11
# $BLENDER_PYTHON -m pip install --upgrade maturin

maturin build --release --interpreter $BLENDER_PYTHON
