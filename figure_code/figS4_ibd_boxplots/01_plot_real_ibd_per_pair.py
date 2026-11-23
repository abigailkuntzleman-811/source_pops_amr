#!/usr/bin/env python3
"""Figure S4: distribution of IBD per haplotype pair in the real data.

Reads the per-pair table written by figure_code/fig3/01_within_population_ibd.py
(pileup-masked IBD, KING unrelated set, every within-population haplotype pair
including those with no IBD), so S4 uses exactly the same data as Figure 3.
Same calculation as combo_dict in Chapter_1.ipynb.

As before, pairs with no IBD are left out of the boxes (a log axis cannot show
them) and counted in the table.

Usage:
    python figure_code/figS4_ibd_boxplots/01_plot_real_ibd_per_pair.py

Output:
    $FIGURE_DIR/IBD_boxplots_excludes0s.png
    $TABLE_DIR/IBD_boxplots_excludes0s.csv
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import FIGURE_DIR, IBD_DIR, TABLE_DIR  # noqa: E402

# same order and colours as before (and as Figure S5)
LABELS = ["IBS", "CLM", "MXL", "PEL", "PUR"]
COLORS = ["#1f4e79", "#A3D9FF", "#DA3E52", "#96E6B3", "#7E6B8F"]


def main():
    path = IBD_DIR / "within_population_pairs.csv.gz"
    if not path.exists():
        sys.exit(f"Missing {path}\nRun: python figure_code/fig3/01_within_population_ibd.py")
    pairs = pd.read_csv(path, usecols=["pop", "cM"])
    pairs = pairs[pairs["pop"].isin(LABELS)]
    values = {lab: pairs.loc[pairs["pop"] == lab, "cM"].to_numpy() for lab in LABELS}

    rows = []
    for lab in LABELS:
        v = values[lab]
        rows.append({
            "Population": lab,
            "Total pairs": len(v),
            "Pairs with 0 IBD": int((v == 0).sum()),
            "Mean IBD across all pairs": v.mean(),
            "Mean IBD excluding 0-IBD pairs": v[v > 0].mean() if (v > 0).any() else np.nan,
            "Median IBD excluding 0-IBD pairs": np.median(v[v > 0]) if (v > 0).any() else np.nan,
        })
    table = pd.DataFrame(rows)
    print(table.round(3).to_string(index=False))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE_DIR / "IBD_boxplots_excludes0s.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    box = ax.boxplot([values[lab][values[lab] > 0] for lab in LABELS], patch_artist=True)
    ax.set_xticks(range(1, len(LABELS) + 1))
    ax.set_xticklabels(LABELS)
    for patch, color in zip(box["boxes"], COLORS):
        patch.set_facecolor(color)

    # log axis (values span 2 to a few hundred cM) with plain-number labels
    ax.set_yscale("log")
    lo, hi = ax.get_ylim()
    nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
    ax.yaxis.set_major_locator(FixedLocator([t for t in nice if lo <= t <= hi]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_ylabel("Total IBD (cM) per Pair")

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "IBD_boxplots_excludes0s.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}\nwrote {TABLE_DIR / 'IBD_boxplots_excludes0s.csv'}")


if __name__ == "__main__":
    main()
