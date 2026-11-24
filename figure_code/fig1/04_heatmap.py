#!/usr/bin/env python3
"""Figure 1: heatmap of pairwise shared archaic ancestry between populations.

Reads the pairwise sharing files written by 03_pairwise_sharing.sh, builds a
symmetric population-by-population matrix of mean shared bp per haplotype pair,
and draws the lower-triangle heatmap used in the paper. Diagonal entries are
within-population sharing.

Usage:
    python figure_code/fig1/04_heatmap.py figure        # Figure 1: A Neanderthal, B Denisovan
    python figure_code/fig1/04_heatmap.py neanderthal   # single panels
    python figure_code/fig1/04_heatmap.py denisovan
    python figure_code/fig1/04_heatmap.py all

Output:
    $FIGURE_DIR/shared_heatmap.png                  (figure; the name used in the paper)
    $FIGURE_DIR/<archaic>_shared_heatmap.png        (single panels)
    $FIGURE_DIR/<archaic>_shared_heatmap_values.csv
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import FIGURE_DIR, PAIRWISE_DIR  # noqa: E402

# Display order for the heatmap axes: East Asian, European, South Asian,
# then admixed American.
SUB_POPS = [
    "CHB", "JPT", "CHS", "CDX", "KHV",
    "CEU", "TSI", "GBR", "FIN", "IBS",
    "GIH", "PJL", "BEB", "STU", "ITU",
    "MXL", "PEL", "PUR", "CLM",
]


# Colours sampled from the colour bar of Figure 1 as published (light yellow to
# dark green), which is more muted than matplotlib's YlGn.
CMAP = LinearSegmentedColormap.from_list("fig1", [
    "#fdfecb", "#f9fcc6", "#e0eeb1", "#bedb9d", "#96c287",
    "#71a76b", "#528050", "#3b6540", "#1c3422",
])


# get_averages is copied unchanged from Chapter_1.ipynb (cell 8), so the
# numbers are identical to the original analysis.
#returns # of pairs, average of all pairs, # of pairs that have no sharing, average of pairs conditioned on sharing
def get_averages(filename):
    with open(filename, 'r') as file:
        integers = []
        no_share = 0
        for line in file:
            line = line.strip().split()
            id1 = line[0]
            id2 = line[1]
            if len(line) == 3:
                num_intersect = int(line[2])
            if 'sorted' in id1 and 'sorted' in id2:
                integers.append(num_intersect)
                if num_intersect == 0:
                    no_share += 1
        nonzero = [x for x in integers if x > 0]
        avg_cond_on_share = sum(nonzero) / len(nonzero) if nonzero else 0

        # #pairs, average of all pairs, # of pairs that have no sharing, average of pairs conditioned on sharing
        return len(integers), sum(integers)/len(integers), no_share, avg_cond_on_share


def build_matrix(archaic: str) -> np.ndarray:
    """Symmetric matrix of mean pairwise sharing, in SUB_POPS order."""
    directory = PAIRWISE_DIR / archaic
    if not directory.is_dir():
        sys.exit(
            f"Missing input directory: {directory}\n"
            f"Run: bash fig1/03_pairwise_sharing.sh {archaic}"
        )

    n = len(SUB_POPS)
    matrix = np.zeros((n, n))
    missing = []

    # Off-diagonal: between-population comparisons.
    for i in range(n):
        for j in range(i + 1, n):
            # The filename uses whichever population came first in the pairwise
            # run, so try both orderings.
            candidates = [
                directory / f"{SUB_POPS[i]}_{{{SUB_POPS[j]}}}.txt",
                directory / f"{SUB_POPS[j]}_{{{SUB_POPS[i]}}}.txt",
            ]
            path = next((p for p in candidates if p.is_file()), None)
            if path is None:
                missing.append(f"{SUB_POPS[i]}-{SUB_POPS[j]}")
                continue

            _, mean_all, _, _ = get_averages(path)
            matrix[i, j] = matrix[j, i] = mean_all

    # Diagonal: within-population comparisons.
    for i, pop in enumerate(SUB_POPS):
        path = directory / f"all_{pop}.txt"
        if path.is_file():
            _, mean_all, _, _ = get_averages(path)
            matrix[i, i] = mean_all
        else:
            missing.append(f"{pop}-{pop}")

    if missing:
        print(f"warning: {len(missing)} comparisons missing, left at zero: "
              f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}",
              file=sys.stderr)

    return matrix


def sci_label(x, _pos):
    """3500000 -> 3.5 x 10^6 (same style as Figures 3 and S3)."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = round(x / 10 ** e, 2)
    if m == 10:
        m, e = 1, e + 1
    return rf"${m:g} \times 10^{{{e}}}$"


def draw_panel(ax, matrix: np.ndarray, fig, letter: str = None, title: str = None) -> None:
    """Lower-triangle heatmap (diagonal = within-population sharing) with its colour bar."""
    lower_mask = np.tri(matrix.shape[0], matrix.shape[1], k=-1)
    masked = np.ma.array(matrix, mask=lower_mask.T)
    im = ax.imshow(masked, cmap=CMAP, vmin=np.nanmin(matrix), vmax=np.nanmax(matrix))
    ax.set_xticks(range(len(SUB_POPS)), labels=SUB_POPS,
                  rotation=45, ha="right", rotation_mode="anchor", fontsize=12)
    ax.set_yticks(range(len(SUB_POPS)), labels=SUB_POPS, fontsize=12)
    # colour bar exactly as tall as the heatmap
    cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.5)
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(sci_label))
    cbar.ax.tick_params(labelsize=11)
    if letter:
        ax.text(-0.22, 1.0, letter, transform=ax.transAxes, fontsize=26, va="top")
    if title:
        ax.set_title(title, fontsize=16)


def plot_comparisons(matrix: np.ndarray, archaic: str, out_path: Path) -> None:
    """One panel, titled with the archaic type."""
    fig, ax = plt.subplots()
    draw_panel(ax, matrix, fig, title=archaic.capitalize())
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_figure(neanderthal: np.ndarray, denisovan: np.ndarray, out_path: Path) -> None:
    """Figure 1 as published: (A) Neanderthal, (B) Denisovan, side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    draw_panel(axes[0], neanderthal, fig, letter="A")
    draw_panel(axes[1], denisovan, fig, letter="B")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_values(matrix: np.ndarray, archaic: str) -> None:
    # the underlying values, so the figure can be checked without rerunning
    csv_path = FIGURE_DIR / f"{archaic}_shared_heatmap_values.csv"
    np.savetxt(csv_path, matrix, delimiter=",", header=",".join(SUB_POPS), comments="")
    print(f"wrote {csv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archaic", choices=["figure", "all", "neanderthal", "denisovan"])
    args = parser.parse_args()

    if args.archaic == "figure":
        mats = {a: build_matrix(a) for a in ("neanderthal", "denisovan")}
        out_path = FIGURE_DIR / "shared_heatmap.png"
        plot_figure(mats["neanderthal"], mats["denisovan"], out_path)
        print(f"wrote {out_path}")
        for a, m in mats.items():
            save_values(m, a)
        return

    matrix = build_matrix(args.archaic)
    out_path = FIGURE_DIR / f"{args.archaic}_shared_heatmap.png"
    plot_comparisons(matrix, args.archaic, out_path)
    print(f"wrote {out_path}")
    save_values(matrix, args.archaic)


if __name__ == "__main__":
    main()
