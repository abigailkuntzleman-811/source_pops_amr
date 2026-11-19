#!/usr/bin/env python3
"""Figure 6 diagnostic: the original IBD calculation next to the fixed one.

Left panel:  original get_ibd_dict (tskit min_span = 1 Mb applied to raw,
             unjoined pieces), converted to cM.
Right panel: fixed (sim_ibd.ibd_per_individual: pieces joined, 2 cM cutoff).
Shared y axis, so the change in scale is visible.

Usage:
    python figure_code/fig6/compare_original_vs_fixed.py --reps 5
    python figure_code/fig6/compare_original_vs_fixed.py --pop CLM

Output:
    $FIGURE_DIR/<POP>_IBD_original_vs_fixed.png
    $TABLE_DIR/<POP>_IBD_original_vs_fixed.csv
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURE_DIR, SIM_DIR, TABLE_DIR  # noqa: E402
from sim_ibd import (MODELS, admix_props, ibd_per_individual,  # noqa: E402
                     ibd_per_individual_original)

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "plot_donor_ibd", Path(__file__).resolve().parent / "01_plot_donor_ibd.py")
_plot = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_plot)
STYLE, fit = _plot.STYLE, _plot.fit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pop", default="MXL")
    parser.add_argument("--reps", type=int, default=1000)
    args = parser.parse_args()

    data = {m: {"x": [], "original": [], "fixed": []} for m in MODELS}
    slopes = {m: {"original": [], "fixed": []} for m in MODELS}
    for i in range(1, args.reps + 1):
        print(f"replicate {i}/{args.reps}", flush=True)
        for model, (folder, tracts) in MODELS.items():
            ts = tskit.load(SIM_DIR / folder / "tree_seqs" / f"{args.pop}_rep{i}.ts.gz")
            adm = pd.read_csv(SIM_DIR / folder / "admix_csvs" / tracts.format(pop=args.pop, i=i))
            props = admix_props(adm, ts, args.pop)
            x = np.array([props[k] for k in props])
            data[model]["x"] += list(x)
            for key, fn in [("original", ibd_per_individual_original),
                            ("fixed", ibd_per_individual)]:
                ibd = fn(ts, args.pop)
                y = np.array([ibd[k] for k in props])
                data[model][key] += list(y)
                slopes[model][key].append(np.polyfit(x, y, 1)[0])

    rows = []
    for m in MODELS:
        row = {"model": m}
        for key in ("original", "fixed"):
            row[f"mean_ibd_cM_{key}"] = np.mean(data[m][key])
            row[f"rep_slope_mean_{key}"] = np.mean(slopes[m][key])
            row[f"rep_slope_sd_{key}"] = np.std(slopes[m][key], ddof=1)
        rows.append(row)
    table = pd.DataFrame(rows)
    print("\n" + table.round(2).to_string(index=False))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE_DIR / f"{args.pop}_IBD_original_vs_fixed.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, key, title in [(axes[0], "original", "Original (min_span 1 Mb on unjoined pieces)"),
                           (axes[1], "fixed", "Fixed (joined tracts, 2 cM cutoff)")]:
        for m in MODELS:
            xp, yp, lo, hi, _ = fit(np.array(data[m]["x"]), np.array(data[m][key]))
            ax.plot(xp, yp, color=STYLE[m]["color"], linestyle=STYLE[m]["ls"],
                    linewidth=2, label=m)
            ax.fill_between(xp, lo, hi, color=STYLE[m]["color"], alpha=0.25)
        ax.set_title(title)
        ax.set_xlabel("European admixture proportion", fontsize=12)
    axes[0].set_ylabel("IBD shared with other individuals (cM)", fontsize=12)
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle(f"{args.pop}, {args.reps} replicates")
    fig.tight_layout()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / f"{args.pop}_IBD_original_vs_fixed.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
