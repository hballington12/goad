# --8<-- [start:imports]
from pathlib import Path

from goad import (
    BinningScheme,
    Euler,
    EulerConvention,
    Mapping,
    MultiProblem,
    Orientation,
    Settings,
    ZoneConfig,
)

# --8<-- [end:imports]


# --8<-- [start:orientations]
EULER_CONVENTION = EulerConvention("ZYZ")

# Six orientations of a hexagonal column, sweeping the gamma rotation.
EULERS = [
    (0.0, 30.0, 0.0),
    (0.0, 30.0, 10.0),
    (0.0, 30.0, 20.0),
    (0.0, 30.0, 30.0),
    (0.0, 30.0, 40.0),
    (0.0, 30.0, 50.0),
]
# --8<-- [end:orientations]


# --8<-- [start:zones]
# A 2D zone over the full sphere, fine enough to render as an image.
zones = [
    ZoneConfig(
        BinningScheme.interval(
            thetas=[0, 90, 180],
            theta_spacings=[0.5, 1],
            phis=[0, 360],
            phi_spacings=[2],
        )
    )
]
# --8<-- [end:zones]


# --8<-- [start:settings]
GEOM_PATH = "hex.obj"
OUTPUT_ROOT = Path("runs")

settings = Settings(
    geom_path=GEOM_PATH,
    wavelength=0.532,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    zones=zones,
    mapping=Mapping("ad"),
    beam_power_threshold=0.001,
    beam_area_threshold_fac=0.001,
    cutoff=0.99999,
    max_rec=10,
    max_tir=20,
    seed=None,
    directory=str(OUTPUT_ROOT),
    coherence=True,
    quiet=False,
)
# --8<-- [end:settings]


# --8<-- [start:loop]
OUTPUT_ROOT.mkdir(exist_ok=True)

for i, (alpha, beta, gamma) in enumerate(EULERS):
    settings.orientation = Orientation.discrete(
        eulers=[Euler(float(alpha), float(beta), float(gamma))],
        euler_convention=EULER_CONVENTION,
    )

    out_dir = OUTPUT_ROOT / f"orient_{i:04d}"
    out_dir.mkdir(exist_ok=True)

    mp = MultiProblem(settings)
    mp.solve()
    mp.save(str(out_dir))
# --8<-- [end:loop]
