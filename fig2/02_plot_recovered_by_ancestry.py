#!/usr/bin/env python3
"""Recovered archaic sequence split by continental ancestry background.

For each admixed American population, plots how much archaic sequence is
recovered on European-ancestry haplotypes (circles) versus Indigenous
American-ancestry haplotypes (triangles), as more individuals are sampled.
Reproduces Chapter_1.ipynb cell 25.

Prerequisite: fig2/01_recovered_by_ancestry.py for each population and
ancestry, which requires the local ancestry inference files (LAI_DIR).

Usage:
    python fig2/02_plot_recovered_by_ancestry.py
    python fig2/02_plot_recovered_by_ancestry.py --ancestries EUR AMR AFR

Output:
    results/figures/real_recovered_from_<ancestries>.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import AMR_POPS, FIGURE_DIR, RECOVERED_DIR  # noqa: E402

COLORS = {
    "CLM": "#A3D9FF",
    "MXL": "#DA3E52",
    "PEL": "#96E6B3",
    "PUR": "#7E6B8F",
}

# Marker per continental ancestry, in the order they appear in the legend.
MARKERS = {"EUR": "o", "AMR": "^", "AFR": "s"}


def read_curve(path: Path):
    """Read the two-line 'x: ...' / 'y: ...' format written by step 02."""
    with open(path) as f:
        x = np.array([float(v) for v in f.readline()[3:].strip().split(", ")])
        y = np.array([float(v) for v in f.readline()[3:].strip().split(", ")])
    return x, y


def sci_formatter(val, _pos):
    scaled = val / 1e8
    if float(scaled).is_integer():
        return rf"{int(scaled)} × $10^{{8}}$"
    return rf"{scaled:.1f} × $10^{{8}}$"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ancestries", nargs="+", default=["EUR", "AMR"],
                        help="continental ancestries to plot (default: EUR AMR)")
    args = parser.parse_args()

    unknown = [a for a in args.ancestries if a not in MARKERS]
    if unknown:
        sys.exit(f"No marker defined for: {', '.join(unknown)}")

    plt.figure(figsize=(10, 6))

    # Group legend entries by ancestry, as in the paper: all EUR, then all AMR.
    handles_by_ancestry = {a: [] for a in args.ancestries}
    missing = []

    for ancestry in args.ancestries:
        for pop in AMR_POPS:
            path = RECOVERED_DIR / f"{pop}_only_{ancestry}_archaic_recovered.txt"
            if not path.exists():
                missing.append(f"{pop}/{ancestry}")
                continue

            x, y = read_curve(path)
            handle = plt.scatter(x, y, color=COLORS[pop],
                                 marker=MARKERS[ancestry], s=40, alpha=0.9,
                                 label=f"{pop} ({ancestry})")
            handles_by_ancestry[ancestry].append(handle)

    if missing:
        print(f"warning: missing curves for {', '.join(missing)};\n"
              f"         run fig2/01_recovered_by_ancestry.py (needs LAI_DIR)",
              file=sys.stderr)

    handles = [h for a in args.ancestries for h in handles_by_ancestry[a]]
    if not handles:
        sys.exit(f"No input found in {RECOVERED_DIR}")

    plt.xlabel("Number of individuals sampled", fontsize=12)
    plt.ylabel("Total archaic sequence recovered (bp)", fontsize=12)
    plt.gca().yaxis.set_major_formatter(FuncFormatter(sci_formatter))
    plt.ylim(bottom=0)
    plt.legend(handles=handles, title="Population / Source", frameon=False)
    plt.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / f"real_recovered_from_{'_'.join(args.ancestries).lower()}.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
