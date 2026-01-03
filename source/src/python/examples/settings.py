# --8<-- [start:basic]
from goad import MultiProblem, Settings

# Basic settings with minimal configuration
settings = Settings(geom_path="path/to/geometry.obj")
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:basic]

# --8<-- [start:wavelength]
from goad import MultiProblem, Settings

# Configure wavelength (in micrometers)
settings = Settings(
    geom_path="path/to/geometry.obj",
    wavelength=0.532,  # 532 nm
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:wavelength]

# --8<-- [start:refractive]
from goad import MultiProblem, Settings

# Configure refractive indices for particle and medium
settings = Settings(
    geom_path="path/to/geometry.obj",
    particle_refr_index_re=1.5,  # Real part of particle refractive index
    particle_refr_index_im=0.01,  # Imaginary part (absorption)
    medium_refr_index_re=1.33,  # Water as medium
    medium_refr_index_im=0.0,
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:refractive]

# --8<-- [start:orientation]
from goad import EulerConvention, MultiProblem, Orientation, Settings

# Configure particle orientation distribution
settings = Settings(
    geom_path="path/to/geometry.obj",
    orientation=Orientation.uniform(
        num_orients=100, euler_convention=EulerConvention("ZYZ")
    ),
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:orientation]

# --8<-- [start:orientation_discrete]
from goad import Euler, EulerConvention, MultiProblem, Orientation, Settings

# Configure discrete orientations
orients = Orientation.discrete(
    eulers=[Euler(0, 0, 0), Euler(45, 90, 0)], euler_convention=EulerConvention("ZYZ")
)
settings = Settings(geom_path="path/to/geometry.obj", orientation=orients)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:orientation_discrete]

# --8<-- [start:zones]
from goad import BinningScheme, MultiProblem, Settings, ZoneConfig

# Default: single full zone with interval binning (high-res forward/back)
settings = Settings(geom_path="path/to/geometry.obj")

# Custom full zone with simple binning
settings = Settings(
    geom_path="path/to/geometry.obj",
    zones=[ZoneConfig(BinningScheme.simple(180, 48))],
)

# Labeled zone
settings = Settings(
    geom_path="path/to/geometry.obj",
    zones=[ZoneConfig(BinningScheme.simple(90, 24), label="coarse")],
)

# Backscatter-only (no full zone, just forward + backward)
settings = Settings(geom_path="path/to/geometry.obj", zones=[])

mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:zones]

# --8<-- [start:binning]
from goad import BinningScheme, MultiProblem, Settings, ZoneConfig

# Configure angular binning for scattering output
settings = Settings(
    geom_path="path/to/geometry.obj",
    zones=[ZoneConfig(BinningScheme.simple(num_theta=180, num_phi=48))],
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:binning]

# --8<-- [start:binning_interval]
from goad import BinningScheme, MultiProblem, Settings, ZoneConfig

# Use variable angular resolution
settings = Settings(
    geom_path="path/to/geometry.obj",
    zones=[
        ZoneConfig(
            BinningScheme.interval(
                thetas=[0, 90, 180],
                theta_spacings=[1, 2],  # 1° steps up to 90°, then 2° steps
                phis=[0, 360],
                phi_spacings=[2],
            )
        )
    ],
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:binning_interval]

# --8<-- [start:binning_custom]
from goad import BinningScheme, MultiProblem, Settings, ZoneConfig

# Specify arbitrary bin edges
binning = BinningScheme.custom(
    bins=[
        [[0, 10], [0, 360]],  # Forward scattering cone
        [[10, 170], [0, 360]],  # Side scattering
        [[170, 180], [0, 360]],  # Backscattering cone
    ]
)
settings = Settings(geom_path="path/to/geometry.obj", zones=[ZoneConfig(binning)])
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:binning_custom]

# --8<-- [start:mapping]
from goad import Mapping, MultiProblem, Settings

# Configure near-to-far field mapping method
settings = Settings(
    geom_path="path/to/geometry.obj",
    mapping=Mapping("ad"),  # 'ad' for Aperture Diffraction, 'go' for Geometric Optics
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:mapping]

# --8<-- [start:thresholds]
from goad import MultiProblem, Settings

# Configure beam tracing thresholds
settings = Settings(
    geom_path="path/to/geometry.obj",
    beam_power_threshold=1e-6,  # Stop tracking beams below this power
    beam_area_threshold_fac=1e-3,  # Stop tracking beams smaller than this fraction
    cutoff=1e-10,  # Global energy cutoff
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:thresholds]

# --8<-- [start:recursion]
from goad import MultiProblem, Settings

# Configure ray tracing limits
settings = Settings(
    geom_path="path/to/geometry.obj",
    max_rec=10,  # Maximum internal reflections
    max_tir=5,  # Maximum total internal reflections
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:recursion]

# --8<-- [start:advanced]
from goad import BinningScheme, Mapping, MultiProblem, Orientation, Settings, ZoneConfig

# Complete configuration example
settings = Settings(
    geom_path="path/to/geometry.obj",
    wavelength=0.532,
    particle_refr_index_re=1.5,
    particle_refr_index_im=0.01,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    orientation=Orientation.uniform(num_orients=100),
    zones=[ZoneConfig(BinningScheme.simple(num_theta=180, num_phi=48))],
    mapping=Mapping("ad"),
    beam_power_threshold=1e-6,
    beam_area_threshold_fac=1e-3,
    cutoff=0.999,
    max_rec=10,
    max_tir=5,
    coherence=False,
    quiet=False,
    directory="output/",
)
mp = MultiProblem(settings)
mp.solve()
# --8<-- [end:advanced]
