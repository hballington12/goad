"""Generate a roughened hexagonal plate for backscattering simulation."""

from pathlib import Path

import bpy  # must import first to make bmesh available
from bpy_geometries import HexagonalColumn, Roughened

output_dir = Path(__file__).parent

# 50 micron plate with aspect ratio 0.1
# radius = 25 microns, length = 5 microns
plate = Roughened(
    HexagonalColumn(length=5.0, radius=25.0, output_dir=output_dir),
    max_edge_length=15.0,
    displacement_sigma=0.1,
    merge_distance=1.0,
)

filepath = plate.generate()
print(f"Generated geometry: {filepath}")
