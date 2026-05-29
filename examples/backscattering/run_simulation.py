"""Run GOAD backscattering simulation on the roughened hexagonal plate."""

from pathlib import Path

from goad import Convergence, Geom, Param, Settings

output_dir = Path(__file__).parent

# Path to the generated geometry
geometry_file = (
    output_dir
    / "roughened_edge15.0_sigma0p1_merge1p0_hexagonal_column_l5.0_r25.0_6904c8.obj"
)

# Load geometry with the ice refractive index
geoms = Geom.from_file(str(geometry_file), [1.31 + 0j])

# Set up simulation settings
# - wavelength 0.532 microns (green laser)
# - empty zones for faster computation (just forward and backward scattering)
# - high max_tir for backscattering
settings = Settings(
    wavelength=0.532,
    zones=[],
    max_tir=20,
    directory=str(output_dir),
    beam_area_threshold_fac=0.01,
    beam_power_threshold=0.01,
    cutoff=0.999,
)

# Set up convergence solver
convergence = Convergence(settings, geoms)

# Converge on lidar ratio and depolarisation ratio within 10% error
convergence.add_target(Param.LidarRatio, 0.1)
convergence.add_target(Param.DepolarizationRatio, 0.1)

# Run the simulation
convergence.solve()

# Save results
convergence.save()

print("Simulation complete!")
print(f"Results saved to: {output_dir}")
