"""Generate a single hexagonal bullet arm for a three-bullet rosette aggregate."""
import bpy  # must import first to make bmesh available

from pathlib import Path
from bpy_geometries import HexagonalBullet

output_dir = Path(__file__).parent

# 10 micron long arm, radius 1.5 microns, indented tip inset by 2 microns
bullet = HexagonalBullet(
    length=10.0,
    radius=1.5,
    indentation_factor=0.3,
    inset=2.0,
    output_dir=output_dir,
)

filepath = bullet.generate()
print(f"Generated geometry: {filepath}")
