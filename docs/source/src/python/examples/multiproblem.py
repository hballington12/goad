# --8<-- [start:multiproblem]
from goad import Geom, MultiProblem, Settings

# Load geometry with its per-shape refractive index, then run with default settings
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
mp = MultiProblem(Settings(), geoms)
mp.solve()
mp.save("my_results")  # Save results to disk for later analysis
# --8<-- [end:multiproblem]
