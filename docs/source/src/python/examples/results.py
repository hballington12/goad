# --8<-- [start:basic]
from goad import MultiProblem, Settings

# Solve the problem and access results
mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()
results = mp.results
# --8<-- [end:basic]

# --8<-- [start:mueller_2d]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Get the 2D Mueller matrix
mueller = mp.results.mueller
print(f"Number of bins: {len(mueller)}")
print(f"Mueller matrix elements per bin: {len(mueller[0])}")  # 16 elements
# --8<-- [end:mueller_2d]

# --8<-- [start:mueller_1d]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Get the 1D phi-integrated Mueller matrix
mueller_1d = mp.results.mueller_1d
print(f"Number of theta bins: {len(mueller_1d)}")
# --8<-- [end:mueller_1d]

# --8<-- [start:mueller_components]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Access different scattering components
mueller_total = mp.results.mueller  # Total scattering
mueller_beam = mp.results.mueller_beam  # Beam component
mueller_ext = mp.results.mueller_ext  # External diffraction

# Same for 1D Mueller matrices
mueller_1d_total = mp.results.mueller_1d
mueller_1d_beam = mp.results.mueller_1d_beam
mueller_1d_ext = mp.results.mueller_1d_ext
# --8<-- [end:mueller_components]

# --8<-- [start:bins]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Get 2D bins (theta, phi pairs)
bins_2d = mp.results.bins
for theta, phi in bins_2d[:5]:  # First 5 bins
    print(f"Theta: {theta}°, Phi: {phi}°")

# Get 1D bins (theta only)
bins_1d = mp.results.bins_1d
if bins_1d:
    print(f"Theta values: {bins_1d}")
# --8<-- [end:bins]

# --8<-- [start:scat_cross]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

scat_cross = mp.results.scat_cross
print(f"Scattering cross section: {scat_cross}")
# --8<-- [end:scat_cross]

# --8<-- [start:ext_cross]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

ext_cross = mp.results.ext_cross
print(f"Extinction cross section: {ext_cross}")
# --8<-- [end:ext_cross]

# --8<-- [start:asymmetry]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

g = mp.results.asymmetry
print(f"Asymmetry parameter: {g}")
# --8<-- [end:asymmetry]

# --8<-- [start:albedo]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

albedo = mp.results.albedo
print(f"Single scattering albedo: {albedo}")
# --8<-- [end:albedo]

# --8<-- [start:powers]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

powers = mp.results.powers
print(f"Input power: {powers['input']}")
print(f"Output power: {powers['output']}")
print(f"Absorbed power: {powers['absorbed']}")
print(f"Missing power: {powers['missing']}")

# Check energy conservation
total_accounted = (
    powers["output"]
    + powers["absorbed"]
    + powers["trnc_ref"]
    + powers["trnc_rec"]
    + powers["trnc_clip"]
    + powers["trnc_energy"]
    + powers["trnc_area"]
)
print(f"Energy conservation error: {powers['input'] - total_accounted}")
# --8<-- [end:powers]

# --8<-- [start:zones]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Access all zones
zones = mp.results.zones
print(f"Available zones: {zones}")

# Access specific zones by type
full_zone = mp.results.full_zone
forward_zone = mp.results.forward_zone
backward_zone = mp.results.backward_zone

print(f"Full zone: {full_zone.name}, bins: {full_zone.num_bins}")
print(f"Forward zone: {forward_zone.name}, bins: {forward_zone.num_bins}")
print(f"Backward zone: {backward_zone.name}, bins: {backward_zone.num_bins}")
# --8<-- [end:zones]

# --8<-- [start:zone_params]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Full zone parameters (requires full theta coverage)
full_zone = mp.results.full_zone
full_params = full_zone.params
print(f"Scattering cross section: {full_params['scatt_cross']}")
print(f"Extinction cross section: {full_params['ext_cross']}")
print(f"Asymmetry parameter: {full_params['asymmetry']}")
print(f"Single scattering albedo: {full_params['albedo']}")

# Backward zone parameters (lidar-relevant)
backward_zone = mp.results.backward_zone
back_params = backward_zone.params
print(f"Backscatter cross section: {back_params['backscatter_cross']}")
print(f"Lidar ratio: {back_params['lidar_ratio']}")
print(f"Depolarization ratio: {back_params['depolarization_ratio']}")

# Forward zone parameters (optical theorem)
forward_zone = mp.results.forward_zone
fwd_params = forward_zone.params
print(f"Extinction (optical theorem): {fwd_params['ext_cross_optical_theorem']}")
# --8<-- [end:zone_params]

# --8<-- [start:zone_mueller]
from goad import MultiProblem, Settings

mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()

# Access Mueller matrix for a specific zone
zone = mp.results.full_zone
mueller = zone.mueller  # 2D Mueller matrix for this zone
mueller_1d = zone.mueller_1d  # 1D phi-integrated Mueller matrix
bins = zone.bins  # Angular bins for this zone
bins_1d = zone.bins_1d  # 1D theta bins

print(f"Zone '{zone.name}' has {zone.num_bins} bins")
print(f"Mueller matrix shape: {mueller.shape}")
# --8<-- [end:zone_mueller]

# --8<-- [start:complete]
import numpy as np
from goad import BinningScheme, MultiProblem, Orientation, Settings, ZoneConfig

# Configure and solve
settings = Settings(
    geom_path="path/to/geometry.obj",
    wavelength=0.532,
    orientation=Orientation.uniform(num_orients=100),
    zones=[ZoneConfig(BinningScheme.simple(num_theta=180, num_phi=48))],
)
mp = MultiProblem(settings)
mp.solve()

# Access all results
results = mp.results

# Extract scattering phase function (S11 element)
mueller_1d = np.array(results.mueller_1d)
s11 = mueller_1d[:, 0]  # First column is S11
theta = np.array(results.bins_1d)

# Print integrated parameters
print(f"Scattering cross section: {results.scat_cross:.6f}")
print(f"Extinction cross section: {results.ext_cross:.6f}")
print(f"Asymmetry parameter: {results.asymmetry:.4f}")
print(f"Single scattering albedo: {results.albedo:.4f}")

# Access zone-specific parameters
backward = results.backward_zone
print(f"Lidar ratio: {backward.params['lidar_ratio']:.4f}")
print(f"Backscatter cross section: {backward.params['backscatter_cross']:.6f}")

# Check power budget
powers = results.powers
efficiency = powers["output"] / powers["input"]
print(f"Scattering efficiency: {efficiency:.4f}")
print(f"Missing power fraction: {powers['missing'] / powers['input']:.2e}")
# --8<-- [end:complete]
