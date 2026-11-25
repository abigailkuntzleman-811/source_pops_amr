#!/usr/bin/env python3
"""Stack the three UpSet plots into one figure, laid out like Figure S2:
(A) Neanderthal, (B) Denisovan, (C) all archaic sequence, one above another,
with bold panel letters at the top left.

Run 01_upset_shared_unique.py for all three types first.

Usage:
    python figure_code/upset/02_combine_panels.py

Output: $FIGURE_DIR/upset_combined.png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import FIGURE_DIR  # noqa: E402

PANELS = [("neanderthal", "A"), ("denisovan", "B"), ("all", "C")]   # same order as Figure S2
DPI = 300   # the panels are saved at 300 dpi, so this keeps them at full resolution


def main():
    images = []
    for archaic, _ in PANELS:
        path = FIGURE_DIR / f"upset_{archaic}.png"
        if not path.exists():
            sys.exit(f"Missing {path}\nRun: python figure_code/upset/01_upset_shared_unique.py {archaic}")
        images.append(mpimg.imread(path))

    width = max(img.shape[1] for img in images)
    heights = [img.shape[0] for img in images]
    gap = int(0.06 * width)   # white space between panels, room for the letters
    fig = plt.figure(figsize=(width / DPI, (sum(heights) + gap * len(images)) / DPI), dpi=DPI)
    grid = fig.add_gridspec(len(images), 1, height_ratios=[h + gap for h in heights], hspace=0)
    for row, (img, (_, letter)) in enumerate(zip(images, PANELS)):
        ax = fig.add_subplot(grid[row])
        ax.imshow(img, interpolation="none")
        ax.set_xlim(0, width)
        ax.set_ylim(img.shape[0], -gap)    # gap above each panel for its letter
        ax.axis("off")
        ax.text(0, -gap * 0.25, letter, fontsize=24, fontweight="bold", va="bottom", ha="left")

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    out = FIGURE_DIR / "upset_combined.png"
    fig.savefig(out, dpi=DPI, facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
