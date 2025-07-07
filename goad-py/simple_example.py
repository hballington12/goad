import goad_py as goad
import os

# Path to example geometry
current_dir = os.path.dirname(os.path.abspath(__file__))
geom_path = os.path.join(current_dir, "..", "examples", "data", "hex.obj")
geom_path = os.path.abspath(geom_path)

# Basic setup
settings = goad.Settings(geom_path)
mp = goad.MultiProblem(settings)
mp.py_solve()

# Results
results = mp.results
print(f"Scattering cross-section: {results.scat_cross:.6f}")
print(f"Extinction cross-section: {results.ext_cross:.6f}")
print(f"Absorption cross-section: {results.ext_cross - results.scat_cross:.6f}")
print(f"Single scattering albedo: {results.albedo:.6f}")
print(f"Asymmetry parameter: {results.asymmetry:.6f}")