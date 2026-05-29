# --8<-- [start:basic]
from goad import Geom, MultiProblem, Settings

# Basic settings with minimal configuration
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings()
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:basic]

# --8<-- [start:containment_tree]
from goad import Geom  # noqa: E402

geom = Geom.from_file("path/to/multi_shape.obj", [1.31 + 0j])[0]
print(geom.containment_tree())
# --8<-- [end:containment_tree]

# --8<-- [start:shape_refr_index]
from goad import Geom  # noqa: E402

geom = Geom.from_file("path/to/multi_shape.obj", [1.31 + 0j])[0]

# Override the refractive index of individual shapes
geom.set_refr_index(0, 1.5 + 0.0j)  # Set shape 0 to 1.5
geom.set_refr_index(1, 1.33 + 0.01j)  # Set shape 1 to 1.33 + 0.01j

print(geom.containment_tree())
# --8<-- [end:shape_refr_index]

# --8<-- [start:wavelength]
from goad import Geom, MultiProblem, Settings  # noqa: E402

# Configure wavelength (in micrometers)
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
    wavelength=0.532,  # 532 nm
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:wavelength]

# --8<-- [start:refractive]
from goad import Geom, MultiProblem, Settings  # noqa: E402

# Configure refractive indices for particle and medium
geoms = Geom.from_file("path/to/geometry.obj", [1.5 + 0.01j])
settings = Settings(
    medium_refr_index=1.33 + 0.0j,
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:refractive]

# --8<-- [start:orientation]
from goad import EulerConvention, Geom, MultiProblem, Orientation, Settings  # noqa: E402

# Configure particle orientation distribution
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
    orientation=Orientation.uniform(
        num_orients=100, euler_convention=EulerConvention("ZYZ")
    ),
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:orientation]

# --8<-- [start:orientation_discrete]
from goad import (  # noqa: E402
    Euler,
    EulerConvention,
    Geom,
    MultiProblem,
    Orientation,
    Settings,
)

# Configure discrete orientations
orients = Orientation.discrete(
    eulers=[Euler(0, 0, 0), Euler(45, 90, 0)], euler_convention=EulerConvention("ZYZ")
)
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(orientation=orients)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:orientation_discrete]

# --8<-- [start:zones]
from goad import BinningScheme, Geom, MultiProblem, Settings, ZoneConfig  # noqa: E402

geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])

# Default: single full zone with interval binning (high-res forward/back)
settings = Settings()

# Custom full zone with simple binning
settings = Settings(
    zones=[ZoneConfig(BinningScheme.simple(180, 48))],
)

# Labeled zone
settings = Settings(
    zones=[ZoneConfig(BinningScheme.simple(90, 24), label="coarse")],
)

# Backscatter-only (no full zone, just forward + backward)
settings = Settings(zones=[])

mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:zones]

# --8<-- [start:binning]
from goad import BinningScheme, Geom, MultiProblem, Settings, ZoneConfig  # noqa: E402

# Configure angular binning for scattering output
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
    zones=[ZoneConfig(BinningScheme.simple(num_theta=180, num_phi=48))],
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:binning]

# --8<-- [start:binning_interval]
from goad import BinningScheme, Geom, MultiProblem, Settings, ZoneConfig  # noqa: E402

# Use variable angular resolution
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
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
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:binning_interval]

# --8<-- [start:binning_custom]
from goad import BinningScheme, Geom, MultiProblem, Settings, ZoneConfig  # noqa: E402

# Specify arbitrary bin edges
binning = BinningScheme.custom(
    bins=[
        [[0, 10], [0, 360]],  # Forward scattering cone
        [[10, 170], [0, 360]],  # Side scattering
        [[170, 180], [0, 360]],  # Backscattering cone
    ]
)
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(zones=[ZoneConfig(binning)])
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:binning_custom]

# --8<-- [start:mapping]
from goad import Geom, Mapping, MultiProblem, Settings  # noqa: E402

# Configure near-to-far field mapping method
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
    mapping=Mapping("ad"),  # 'ad' for Aperture Diffraction, 'go' for Geometric Optics
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:mapping]

# --8<-- [start:thresholds]
from goad import Geom, MultiProblem, Settings  # noqa: E402

# Configure beam tracing thresholds
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
    beam_power_threshold=1e-6,  # Stop tracking beams below this power
    beam_area_threshold_fac=1e-3,  # Stop tracking beams smaller than this fraction
    cutoff=1e-10,  # Global energy cutoff
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:thresholds]

# --8<-- [start:recursion]
from goad import Geom, MultiProblem, Settings  # noqa: E402

# Configure ray tracing limits
geoms = Geom.from_file("path/to/geometry.obj", [1.31 + 0j])
settings = Settings(
    max_rec=10,  # Maximum internal reflections
    max_tir=5,  # Maximum total internal reflections
)
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:recursion]

# --8<-- [start:advanced]
from goad import (  # noqa: E402
    BinningScheme,
    Geom,
    Mapping,
    MultiProblem,
    Orientation,
    Settings,
    ZoneConfig,
)

# Complete configuration example
geoms = Geom.from_file("path/to/geometry.obj", [1.5 + 0.01j])
settings = Settings(
    medium_refr_index=1.0 + 0.0j,
    wavelength=0.532,
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
mp = MultiProblem(settings, geoms)
mp.solve()
# --8<-- [end:advanced]
