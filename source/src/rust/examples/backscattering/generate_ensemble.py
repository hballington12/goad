"""Generate 5 roughened hexagonal plates for ensemble averaging demonstration."""

from pathlib import Path

import bpy  # must import first to make bmesh available
from bpy_geometries import HexagonalColumn, Roughened

output_dir = Path(__file__).parent / "ensemble"
output_dir.mkdir(exist_ok=True)

# Generate 5 particles with same statistical parameters but different realizations
for i in range(5):
    plate = Roughened(
        HexagonalColumn(length=5.0, radius=25.0, output_dir=output_dir),
        max_edge_length=15.0,
        displacement_sigma=0.1,
        merge_distance=1.0,
    )
    filepath = plate.generate()
    print(f"Generated particle {i + 1}: {filepath}")

print(f"\nAll particles saved to: {output_dir}")
