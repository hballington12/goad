"""Run a fixed-orientation GOAD simulation on the large bullet rosette aggregate."""
from pathlib import Path
from goad import Geom, MultiProblem, Settings

output_dir = Path(__file__).parent
geometry_file = output_dir / "aggregate-2p9mm_tri.obj"

geoms = Geom.from_file(str(geometry_file), [1.39 + 0j])

settings = Settings(
    wavelength=0.2,
    seed=42,
    directory=str(output_dir / "output"),
)

mp = MultiProblem(settings, geoms)
mp.solve()
mp.save()

print(mp.results.asymmetry, mp.results.scat_cross, mp.results.ext_cross, mp.results.albedo)
