"""Assemble the per-orientation frame_*.png renders into a looping animated
GIF with proper transparency, using ffmpeg's two-pass palettegen/paletteuse
filters. Hidden script — invoked during docs build, output GIF ships with
the site."""

import shutil
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[4]
IMG_DIR = REPO_ROOT / "docs/examples/images/fixed-orientation-patterns"

FRAME_DURATION_S = 0.4  # per-frame duration in seconds
FRAMERATE = 1.0 / FRAME_DURATION_S


def find_ffmpeg() -> str:
    """Prefer imageio-ffmpeg's bundled static binary so we don't depend on
    a system install. Fall back to whatever's on PATH."""
    try:
        import imageio_ffmpeg  # type: ignore

        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        binary = shutil.which("ffmpeg")
        if binary is None:
            raise SystemExit(
                "ffmpeg not found. install imageio-ffmpeg or put ffmpeg on PATH"
            )
        return binary


def main():
    ffmpeg = find_ffmpeg()

    frames = sorted(IMG_DIR.glob("frame_*.png"))
    if not frames:
        raise SystemExit(
            f"no frame_*.png found in {IMG_DIR} -- run the plotter first"
        )

    out = IMG_DIR / "animation.gif"
    palette = IMG_DIR / "palette.png"

    # Pass 1: build a global palette across all frames, reserving one slot
    # for transparency so RGBA inputs map cleanly to a 1-bit-alpha GIF.
    subprocess.run(
        [
            "ffmpeg", "-y",
            "-framerate", str(FRAMERATE),
            "-i", str(IMG_DIR / "frame_%04d.png"),
            "-vf", "palettegen=reserve_transparent=1",
            str(palette),
        ],
        check=True,
    )

    # Pass 2: encode the GIF using that palette. `alpha_threshold=128`
    # threshold-binarises the input alpha; `dither=bayer` keeps file size
    # down without obvious banding on the inferno colormap.
    subprocess.run(
        [
            "ffmpeg", "-y",
            "-framerate", str(FRAMERATE),
            "-i", str(IMG_DIR / "frame_%04d.png"),
            "-i", str(palette),
            "-lavfi", "paletteuse=alpha_threshold=128:dither=bayer:bayer_scale=5",
            "-loop", "0",
            str(out),
        ],
        check=True,
    )

    palette.unlink(missing_ok=True)
    print(f"wrote {out} ({len(frames)} frames, {FRAME_DURATION_S * 1000:.0f} ms each)")


if __name__ == "__main__":
    main()
