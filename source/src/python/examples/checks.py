# --8<-- [start:binning_resolution]
from goad import BinningScheme, MultiProblem, Settings

geom_path = "hex.obj"
wavelength = 0.532
# hex.obj is a column with length 10 and radius 5
# wavelength is 0.532, so we expect the forward scattering peak to
# have width in degrees ≈ 0.532/10 * 180/π = 3°
# A reasonable binning scheme would have at least 3 bins to cover this range
# ie. 180 theta bins (from 0° to 180°)
binning = BinningScheme.simple(num_theta=180, num_phi=60)
settings = Settings(
    geom_path=geom_path, wavelength=wavelength, binning=binning, seed=42
)
mp = MultiProblem(settings)
mp.solve()
print(f"Asymmetry: {mp.results.asymmetry}")  # prints: "0.8129"

# Example of what not to do: Only 1 bin to cover the forward scattering peak:
binning = BinningScheme.simple(num_theta=60, num_phi=60)
settings = Settings(
    geom_path=geom_path, wavelength=wavelength, binning=binning, seed=42
)
mp = MultiProblem(settings)
mp.solve()
print(f"Asymmetry: {mp.results.asymmetry}")  # prints: "0.8531" (5% error)

# Check that 1 degree binning is sufficient by increasing the resolution without
# significant change in asymmetry:
binning = BinningScheme.simple(num_theta=360, num_phi=60)
settings = Settings(
    geom_path=geom_path, wavelength=wavelength, binning=binning, seed=42
)
mp = MultiProblem(settings)
mp.solve()
print(f"Asymmetry: {mp.results.asymmetry}")  # prints: "0.8108"

# --8<-- [end:binning_resolution]
