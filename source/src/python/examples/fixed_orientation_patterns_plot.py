"""Render 2D rectangular S11 scattering patterns from the per-orientation
runs produced by `fixed_orientation_patterns.py`. Hidden script — invoked
during docs build, output PNGs ship with the site, raw run data does not.

Adapted from /Users/ixguard/Documents/work/post-doc/frozen-droplets/plotting/
droxtal-bench/plot_2d_rectangular.py.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Layout of where the loop script writes its output, and where rendered
# images should go for the docs. All paths are anchored off this script's
# location so the plotter can be invoked from any cwd.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[4]  # docs/source/src/python/examples -> repo root
RUN_ROOT = SCRIPT_DIR / "runs"
RENDER_DIR = SCRIPT_DIR / "renders"  # disposable runs/ won't clobber these
ZONE_DIR = "full"  # display_name auto-inferred for a 0–180 theta interval zone
IMG_OUT = REPO_ROOT / "docs/examples/images/fixed-orientation-patterns"

EULERS = [(0.0, 30.0, gamma) for gamma in (0, 10, 20, 30, 40, 50)]

S11_COL = 2  # mueller_scatgrid columns: theta, phi, S11, S12, ..., S44


def load_2d_mueller(filepath: Path):
    """Load mueller_scatgrid (theta, phi, S11..S44 row-major) and return
    (theta_unique, phi_in_[-180,180], S11_grid)."""
    data = np.loadtxt(filepath)
    theta = data[:, 0]
    phi = data[:, 1]
    s11 = data[:, S11_COL]

    theta_unique = np.unique(theta)
    phi_unique = np.unique(phi)
    s11_grid = s11.reshape((len(theta_unique), len(phi_unique)))

    # phi: [0, 360] -> [-180, 180] for symmetric display
    phi_t = np.where(phi_unique > 180, phi_unique - 360, phi_unique)
    sort_idx = np.argsort(phi_t)
    return theta_unique, phi_t[sort_idx], s11_grid[:, sort_idx]


THETA_TICKS = [0, 90, 180]
PHI_TICKS = [-180, -90, 0, 90, 180]
# Upper-right corner inset for the orientation render. Axes-fraction units;
# size = 25% of the axes; centre is fixed (same in every panel).
INSET_BOUNDS = [0.74, 0.74, 0.25, 0.25]  # [x0, y0, w, h] in transAxes


def render_panel(ax, theta, phi, s11, vmin, vmax, render_path=None):
    th_grid, ph_grid = np.meshgrid(theta, phi, indexing="ij")
    mesh = ax.pcolormesh(
        th_grid,
        ph_grid,
        np.maximum(s11, 1e-12),
        cmap="inferno",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="nearest",
        rasterized=True,
    )
    ax.set_xlim(0, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(THETA_TICKS)
    ax.set_yticks(PHI_TICKS)
    ax.set_box_aspect(1)  # square subplot box regardless of data ranges

    if render_path is not None and Path(render_path).exists():
        inset = ax.inset_axes(INSET_BOUNDS, transform=ax.transAxes)
        img = plt.imread(render_path)
        inset.imshow(img)
        inset.set_xticks([])
        inset.set_yticks([])
        inset.patch.set_alpha(0.0)  # transparent inset background
        for s in inset.spines.values():
            s.set_visible(False)
    return mesh


def main():
    IMG_OUT.mkdir(parents=True, exist_ok=True)

    # Load all six orientations.
    loaded = []
    for i, euler in enumerate(EULERS):
        path = RUN_ROOT / f"orient_{i:04d}" / ZONE_DIR / "mueller_scatgrid"
        if not path.exists():
            raise FileNotFoundError(
                f"missing {path} -- run fixed_orientation_patterns.py first"
            )
        theta, phi, s11 = load_2d_mueller(path)
        render_path = RENDER_DIR / f"orient_{i:04d}.png"
        loaded.append((euler, theta, phi, s11, render_path))

    # Common log-scale color limits so frames are directly comparable.
    s_max = max(d[3].max() for d in loaded)
    vmin = s_max * 1e-7
    vmax = s_max

    # 2 x 3 composite grid for the docs page. Each panel is a square box.
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True)
    mesh = None
    for ax, (_, theta, phi, s11, render_path) in zip(axes.flat, loaded):
        mesh = render_panel(ax, theta, phi, s11, vmin, vmax, render_path=render_path)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\theta$ (deg)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\phi$ (deg)")

    # Manually positioned colorbar spanning the full vertical extent of the
    # rightmost column. The two panels in that column are forced to a square
    # box via set_box_aspect, so we read their final positions after a draw
    # and place a colorbar axes at that exact y-extent — avoids the
    # `fig.colorbar(..., ax=axes)` shrinkage trap.
    fig.canvas.draw()
    top_pos = axes[0, -1].get_position()
    bot_pos = axes[-1, -1].get_position()
    cbar_x = top_pos.x1 + 0.02
    cbar_w = 0.02
    cbar_y = bot_pos.y0
    cbar_h = top_pos.y1 - bot_pos.y0
    cax = fig.add_axes([cbar_x, cbar_y, cbar_w, cbar_h])
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(r"$S_{11}$")

    out = IMG_OUT / "patterns_grid.png"
    fig.savefig(out, dpi=120, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out}")

    # Also render each orientation individually (handy for the GIF later).
    for i, (_, theta, phi, s11, render_path) in enumerate(loaded):
        fig, ax = plt.subplots(figsize=(5, 5))
        mesh_i = render_panel(ax, theta, phi, s11, vmin, vmax, render_path=render_path)
        ax.set_xlabel(r"$\theta$ (deg)")
        ax.set_ylabel(r"$\phi$ (deg)")
        fig.canvas.draw()
        pos = ax.get_position()
        cax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.025, pos.y1 - pos.y0])
        cbar = fig.colorbar(mesh_i, cax=cax)
        cbar.set_label(r"$S_{11}$")
        frame = IMG_OUT / f"frame_{i:04d}.png"
        fig.savefig(frame, dpi=120, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print(f"wrote {frame}")


if __name__ == "__main__":
    main()
