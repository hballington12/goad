# Backscattering Case Study

This example demonstrates computing backscatter properties for ice crystals using GOAD, which is particularly relevant for lidar remote sensing applications.

## Setup

Create a working directory and set up a Python 3.11 environment:

```bash
mkdir -p examples/backscattering
cd examples/backscattering
```

=== "venv"

    ```bash
    python3.11 -m venv .venv
    source .venv/bin/activate
    pip install bpy==4.5.3
    pip install git+https://github.com/hballington12/bpy-geometries.git
    pip install goad-py
    ```

=== "conda"

    ```bash
    conda create -n goad-backscatter python=3.11 -y
    conda activate goad-backscatter
    pip install bpy==4.5.3
    pip install git+https://github.com/hballington12/bpy-geometries.git
    pip install goad-py
    ```

=== "uv"

    ```bash
    uv venv --python 3.11
    source .venv/bin/activate
    uv pip install bpy==4.5.3
    uv pip install git+https://github.com/hballington12/bpy-geometries.git
    uv pip install goad-py
    ```

## Geometry Generation

The [bpy-geometries](https://github.com/hballington12/bpy-geometries) module can be used to create basic ice crystal shapes with modifiers. Available base shapes include `HexagonalColumn`, `Bullet`, `Rosette`, and `Aggregate`. Modifiers like `Roughened` can be applied to add surface roughness.

We'll create a roughened hexagonal plate - a 50 micron plate with aspect ratio 0.1:

- Radius: 25 microns
- Length: 5 microns (50 × 0.1)
- Max edge length: 15 microns (keep this large - at least 20% of particle size and several times the wavelength; too small means slower computation and less accurate results)
- Displacement sigma: 0.1 (subtle roughness)

Create `generate_plate.py`:

```python
"""Generate a roughened hexagonal plate for backscattering simulation."""

import bpy  # must import first to make bmesh available

from pathlib import Path
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
```

Run it:

```bash
python generate_plate.py
```

This produces an OBJ file with the roughened plate geometry.

![Roughened hexagonal plate](images/plate-render.png)

> `roughened_edge15.0_sigma0p1_merge1p0_hexagonal_column_l5.0_r25.0_6904c8.obj`

## Running the Simulation

Create `run_simulation.py`:

```python
"""Run GOAD backscattering simulation on the roughened hexagonal plate."""

from pathlib import Path
from goad import Convergence, Param, Settings

output_dir = Path(__file__).parent

geometry_file = (
    output_dir
    / "roughened_edge15.0_sigma0p1_merge1p0_hexagonal_column_l5.0_r25.0_6904c8.obj"
)

settings = Settings(
    geom_path=str(geometry_file),
    wavelength=0.532,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    zones=[],
    max_tir=20,
    directory=str(output_dir),
    beam_area_threshold_fac=0.01,
    beam_power_threshold=0.01,
    cutoff=0.999,
)

convergence = Convergence(settings)
convergence.add_target(Param.LidarRatio, 0.1)
convergence.add_target(Param.DepolarizationRatio, 0.1)
convergence.solve()
convergence.save()
```

Key settings:

- **wavelength**: 0.532 µm (green lidar)
- **particle_refr_index_re**: 1.31 (ice at 532nm)
- **zones**: empty list disables the main scattering zone, leaving just forward and backward for greater speed
- **max_tir**: 20 (high value for backscattering accuracy)
- **convergence targets**: 10% error on lidar ratio and depolarization ratio

Run it:

```bash
python run_simulation.py
```

Output:

```
⠴ GOAD: [Convergence]  [Elapsed: 42s]  [Status: RUNNING]  [2026-01-03 10:17:38]
  [Orientations: 1731 (10|100000)] [0.025 sec/orientation]
  LidarRatio   3.2150e2 ± 2.8452e1   [ 8.85% / 10.00%] [████████████████████] 100%
Simulation complete!
```

## Results

The simulation saves results to `results.json`. Key backscattering parameters:

| Parameter | Value |
|-----------|-------|
| Lidar Ratio | 18.0 sr |
| Depolarization Ratio | 20.3% |
| Backscatter Cross Section | 120.7 µm² |
| Extinction Cross Section | 2177.2 µm² |

## Ensemble Averaging

Backscattering properties are sensitive to particle shape details. Different realizations of the same statistical roughness parameters can give different results. To get robust estimates, we ensemble average over multiple particle geometries.

Generate 5 particles with identical statistical parameters:

```python
"""Generate 5 roughened hexagonal plates for ensemble averaging."""

from pathlib import Path

import bpy
from bpy_geometries import HexagonalColumn, Roughened

output_dir = Path(__file__).parent / "ensemble"
output_dir.mkdir(exist_ok=True)

for i in range(5):
    plate = Roughened(
        HexagonalColumn(length=5.0, radius=25.0, output_dir=output_dir),
        max_edge_length=15.0,
        displacement_sigma=0.1,
        merge_distance=1.0,
    )
    filepath = plate.generate()
    print(f"Generated particle {i+1}: {filepath}")
```

![Ensemble of 5 roughened plates](images/plate-ensemble-render.png)

We can run the particles individually like before, or we can pass the directory and let GOAD run an average over the ensemble:

```python
settings = Settings(
    geom_path=str(ensemble_dir),  # pass directory instead of single file
    # ... other settings
)
```

> Note: The ensemble average is not a multiple scattering simulation. GOAD selects a particle from the ensemble at random, then selects a random orientation, then runs the simulation. This repeats as the simulation converges. The result is an orientation-averaged, ensemble-averaged single scattering computation.

**Individual particle results:**

| Particle | Lidar Ratio (sr) | Depolarization (%) | Backscatter Cross (µm²) |
|----------|------------------|-------------------|------------------------|
| 1 | 21.2 | 28.4 | 99.4 |
| 2 | 21.3 | 25.1 | 99.6 |
| 3 | 23.6 | 23.3 | 92.3 |
| 4 | 22.9 | 27.5 | 94.3 |
| 5 | 19.4 | 22.5 | 110.8 |

**Ensemble results (3 runs):**

| Run | Lidar Ratio (sr) | Depolarization (%) | Backscatter Cross (µm²) |
|-----|------------------|-------------------|------------------------|
| 1 | 19.4 | 22.9 | 109.7 |
| 2 | 18.8 | 23.4 | 114.0 |
| 3 | 22.2 | 25.9 | 95.7 |

## Roughness Sweep

To explore how surface roughness affects backscattering properties, we can run a parameter sweep over different roughness levels. Here we vary the displacement sigma from 0.25 to 5.0 in steps of 0.25, generating 5 particles for each level.

```python
roughness_levels = [0.25 * i for i in range(1, 21)]  # 0.25, 0.5, ..., 5.0

for sigma in roughness_levels:
    output_dir = base_dir / f"roughness_{sigma}"
    output_dir.mkdir(exist_ok=True)
    
    for i in range(5):
        plate = Roughened(
            HexagonalColumn(length=5.0, radius=25.0, output_dir=output_dir),
            max_edge_length=15.0,
            displacement_sigma=sigma,
            merge_distance=1.0,
        )
        plate.generate()
```

Running ensemble simulations for each roughness level (with 2.5% convergence tolerance):

![Roughness sweep results](images/roughness_sweep_plot.png)
