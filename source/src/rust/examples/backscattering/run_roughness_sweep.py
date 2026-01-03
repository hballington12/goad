"""Run GOAD simulations for roughness sweep."""

import json
from pathlib import Path

from goad import Convergence, Param, Settings

base_dir = Path(__file__).parent

roughness_levels = [0.25 * i for i in range(1, 21)]  # 0.25, 0.5, 0.75, ..., 5.0
results = []

for sigma in roughness_levels:
    sigma_str = str(sigma).replace(".", "p")
    ensemble_dir = base_dir / f"roughness_{sigma_str}"
    results_dir = base_dir / f"roughness_{sigma_str}_results"
    results_dir.mkdir(exist_ok=True)

    print(f"\n{'=' * 60}")
    print(f"Running ensemble for sigma={sigma}")
    print(f"{'=' * 60}")

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
    convergence.add_target(Param.LidarRatio, 0.025)
    convergence.add_target(Param.DepolarizationRatio, 0.025)
    convergence.add_target(Param.BackscatterCross, 0.025)
    convergence.solve()
    convergence.save()

    # Read results
    with open(results_dir / "results.json") as f:
        result = json.load(f)

    backward = next(z for z in result["zones"] if z["zone_type"] == "Backward")
    params = backward["params"]

    results.append(
        {
            "sigma": sigma,
            "lidar_ratio": params["LidarRatio_Total"],
            "depol_ratio": params["DepolarizationRatio_Total"],
            "backscatter_cross": params["BackscatterCross_Total"],
        }
    )

    print(f"  Lidar Ratio: {params['LidarRatio_Total']:.2f} sr")
    print(f"  Depolarization: {params['DepolarizationRatio_Total'] * 100:.1f}%")
    print(f"  Backscatter Cross: {params['BackscatterCross_Total']:.2f} um^2")

# Save all results
with open(base_dir / "roughness_sweep_results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{'=' * 60}")
print("SUMMARY")
print(f"{'=' * 60}")
print(
    f"{'Sigma':<8} {'Lidar Ratio (sr)':<18} {'Depol (%)':<12} {'BS Cross (um^2)':<16}"
)
print("-" * 54)
for r in results:
    print(
        f"{r['sigma']:<8} {r['lidar_ratio']:<18.2f} {r['depol_ratio'] * 100:<12.1f} {r['backscatter_cross']:<16.2f}"
    )

print(f"\nResults saved to: {base_dir / 'roughness_sweep_results.json'}")
