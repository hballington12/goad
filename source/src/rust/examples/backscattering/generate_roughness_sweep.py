"""Generate particles for roughness sweep: 5 particles per roughness level."""

from pathlib import Path

import bpy  # must import first to make bmesh available
from bpy_geometries import HexagonalColumn, Roughened

base_dir = Path(__file__).parent

roughness_levels = [0.25 * i for i in range(1, 21)]  # 0.25, 0.5, 0.75, ..., 5.0

for sigma in roughness_levels:
    # Create directory for this roughness level
    sigma_str = str(sigma).replace(".", "p")
    output_dir = base_dir / f"roughness_{sigma_str}"
    output_dir.mkdir(exist_ok=True)

    print(f"\n{'=' * 60}")
    print(f"Generating 5 particles with sigma={sigma}")
    print(f"{'=' * 60}")

    for i in range(5):
        plate = Roughened(
            HexagonalColumn(length=5.0, radius=25.0, output_dir=output_dir),
            max_edge_length=15.0,
            displacement_sigma=sigma,
            merge_distance=1.0,
        )
        filepath = plate.generate()
        print(f"  Particle {i + 1}: {Path(filepath).name}")

print(f"\nGenerated {len(roughness_levels) * 5} particles total")
