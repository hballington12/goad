"""Run GOAD backscattering simulation with ensemble averaging over all particles."""

import json
from pathlib import Path

from goad import Convergence, Param, Settings

base_dir = Path(__file__).parent
ensemble_dir = base_dir / "ensemble"
results_dir = base_dir / "ensemble_results"
results_dir.mkdir(exist_ok=True)

# Pass the directory containing all OBJ files for ensemble averaging
settings = Settings(
    geom_path=str(ensemble_dir),
    wavelength=0.532,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    zones=[],
    max_tir=20,
    directory=str(results_dir),
    beam_area_threshold_fac=0.01,
    beam_power_threshold=0.01,
    cutoff=0.999,
)

convergence = Convergence(settings)
convergence.add_target(Param.LidarRatio, 0.1)
convergence.add_target(Param.DepolarizationRatio, 0.1)
convergence.add_target(Param.BackscatterCross, 0.1)
convergence.solve()
convergence.save()

# Read and display results
with open(results_dir / "results.json") as f:
    result = json.load(f)

backward = next(z for z in result["zones"] if z["zone_type"] == "Backward")
params = backward["params"]

print(f"\nEnsemble Results:")
print(f"  Lidar Ratio: {params['LidarRatio_Total']:.2f} sr")
print(f"  Depolarization: {params['DepolarizationRatio_Total'] * 100:.1f}%")
print(f"  Backscatter Cross Section: {params['BackscatterCross_Total']:.2f} um^2")
