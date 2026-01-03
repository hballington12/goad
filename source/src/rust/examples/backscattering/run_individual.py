"""Run GOAD backscattering simulation on each particle individually."""

import json
from pathlib import Path

from goad import Convergence, Param, Settings

base_dir = Path(__file__).parent
ensemble_dir = base_dir / "ensemble"
results_dir = base_dir / "individual_results"
results_dir.mkdir(exist_ok=True)

# Get all OBJ files in ensemble directory
obj_files = sorted(ensemble_dir.glob("*.obj"))

results = []

for i, geometry_file in enumerate(obj_files):
    print(f"\n{'=' * 60}")
    print(f"Running particle {i + 1}/{len(obj_files)}: {geometry_file.name}")
    print(f"{'=' * 60}")

    # Create output directory for this particle
    particle_dir = results_dir / f"particle_{i + 1}"
    particle_dir.mkdir(exist_ok=True)

    settings = Settings(
        geom_path=str(geometry_file),
        wavelength=0.532,
        particle_refr_index_re=1.31,
        particle_refr_index_im=0.0,
        zones=[],
        max_tir=20,
        directory=str(particle_dir),
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

    # Read results
    with open(particle_dir / "results.json") as f:
        result = json.load(f)

    # Extract key parameters
    backward = next(z for z in result["zones"] if z["zone_type"] == "Backward")
    params = backward["params"]

    results.append(
        {
            "particle": i + 1,
            "file": geometry_file.name,
            "lidar_ratio": params["LidarRatio_Total"],
            "depol_ratio": params["DepolarizationRatio_Total"],
            "backscatter_cross": params["BackscatterCross_Total"],
        }
    )

    print(f"  Lidar Ratio: {params['LidarRatio_Total']:.2f} sr")
    print(f"  Depolarization: {params['DepolarizationRatio_Total'] * 100:.1f}%")

# Summary
print(f"\n{'=' * 60}")
print("SUMMARY - Individual Particle Results")
print(f"{'=' * 60}")
print(f"{'Particle':<10} {'Lidar Ratio (sr)':<18} {'Depol (%)':<12}")
print("-" * 40)
for r in results:
    print(
        f"{r['particle']:<10} {r['lidar_ratio']:<18.2f} {r['depol_ratio'] * 100:<12.1f}"
    )

# Save summary
with open(results_dir / "summary.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to: {results_dir}")
