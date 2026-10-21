#!/usr/bin/env python3
"""Assemble the Heatwave animation GIF from Octave-exported PNG frames.

GNU Octave writes one PNG per computed time level into
``docs/results/frames/frame_001.png`` ... using the gnuplot graphics
toolkit (headless PNG export).  Octave-side animated GIF export is not
reliable on every build, so this script does the assembly with Pillow.
Every GIF frame is therefore an actual computed solver state, unmodified
apart from GIF palette quantization.

Usage (from the repository root, with Pillow installed):

    python tools/make_gif.py
    python tools/make_gif.py --frames docs/results/frames \
        --out docs/results/animation.gif --duration 80
"""

import argparse
from pathlib import Path

from PIL import Image


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Assemble animation.gif from frame_*.png files."
    )
    parser.add_argument(
        "--frames",
        default="docs/results/frames",
        help="directory containing frame_*.png (default: docs/results/frames)",
    )
    parser.add_argument(
        "--out",
        default="docs/results/animation.gif",
        help="output GIF path (default: docs/results/animation.gif)",
    )
    parser.add_argument(
        "--duration",
        type=int,
        default=80,
        help="frame duration in milliseconds (default: 80)",
    )
    args = parser.parse_args()

    frame_dir = Path(args.frames)
    files = sorted(frame_dir.glob("frame_*.png"))
    if len(files) < 2:
        raise SystemExit(
            f"need at least 2 frames in {frame_dir}, found {len(files)}"
        )

    frames = [
        f.convert("P", palette=Image.ADAPTIVE, colors=256) for f in
        (Image.open(p) for p in files)
    ]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_path,
        save_all=True,
        append_images=frames[1:],
        duration=args.duration,
        loop=0,
        optimize=True,
        disposal=2,
    )
    print(f"wrote {out_path} with {len(frames)} frames "
          f"({args.duration} ms/frame, {args.duration * len(frames) / 1000:.1f} s loop)")


if __name__ == "__main__":
    main()
