#!/usr/bin/env python3
"""Figure S2: archaic sequence recovered vs number of individuals sampled.

Three panels from the curves written by 01_recovered_from_hmmix.py:
(A) Neanderthal, (B) Denisovan, (C) all archaic. Same style as the notebook
(Chapter_1.ipynb, plot_recovered), stacked into one figure.

Usage:
    python figure_code/figS2_hmmix_recovered/02_plot_recovered_from_hmmix.py

Output: $FIGURE_DIR/combined_hmmix_recovered.png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import FIGURE_DIR, RECOVERED_DIR  # noqa: E402

POPS = ["CLM", "MXL", "PEL", "PUR"]
COLORS = {"CLM": "#A3D9FF", "MXL": "#DA3E52", "PEL": "#96E6B3", "PUR": "#7E6B8F"}
PANELS = [  # (archaic, panel letter, y label, power of ten for the tick labels)
    ("neanderthal", "A", "Neanderthal sequence recovered (bp)", 8),
    ("denisovan", "B", "Denisovan sequence recovered (bp)", 7),
    ("all", "C", "Total archaic sequence recovered (bp)", 8),
]


def read_curve(path):
    with open(path) as f:
        x = np.array([float(v) for v in f.readline()[3:].split(", ")])
        y = np.array([float(v) for v in f.readline()[3:].split(", ")])
    return x, y


def main():
    fig, axes = plt.subplots(3, 1, figsize=(10, 18))
    for ax, (archaic, letter, ylabel, exp) in zip(axes, PANELS):
        for pop in POPS:
            path = RECOVERED_DIR / f"{pop}_total_{archaic}_recovered_from_hmmix_pyranges.txt"
            if not path.exists():
                sys.exit(f"Missing {path}\nRun: python figure_code/figS2_hmmix_recovered/01_recovered_from_hmmix.py {archaic}")
            x, y = read_curve(path)
            ax.plot(x, y, color=COLORS[pop], label=pop)
            ax.scatter(x, y, color=COLORS[pop], s=30, alpha=0.9)
        scale = 10 ** exp
        ax.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _p, s=scale, e=exp: rf"{v / s:g} × $10^{{{e}}}$"))
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Number of individuals sampled", fontsize=20)
        ax.set_ylabel(ylabel, fontsize=20)
        ax.tick_params(axis="both", labelsize=14)
        ax.legend(title="Population", frameon=False, fontsize=18, title_fontsize=18,
                  loc="upper left" if archaic == "denisovan" else "best")
        ax.text(-0.12, 1.02, letter, transform=ax.transAxes, fontsize=24, fontweight="bold")
    fig.tight_layout()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "combined_hmmix_recovered.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
