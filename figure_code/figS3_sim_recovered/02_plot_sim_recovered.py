#!/usr/bin/env python3
"""Figure S3: simulated archaic sequence recovered from European (solid) and
Indigenous American (dashed) tracts as haplotypes are added.

Same as the notebook (Chapter_1.ipynb cells 31-35): per replicate, add
haplotypes 1..100 in order and record the merged archaic sequence on each
background; plot the mean over replicates with a 95% confidence interval
(mean +- 1.96 SE).

Usage:
    python figure_code/figS3_sim_recovered/02_plot_sim_recovered.py              # 1000 replicates, standard model
    python figure_code/figS3_sim_recovered/02_plot_sim_recovered.py --reps 100

Output: $FIGURE_DIR/simulated_recovered_amr_eur_yay.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import FIGURE_DIR, SIM_DIR  # noqa: E402

POPS = ["PUR", "PEL", "MXL", "CLM"]
COLORS = {"CLM": "#A3D9FF", "MXL": "#DA3E52", "PEL": "#96E6B3", "PUR": "#7E6B8F"}


def merged_bp(df):
    if df is None or df.empty:
        return 0
    df = df.rename(columns={"start": "Start", "end": "End"}).copy()
    df["Chromosome"] = "chr1"
    m = pr.PyRanges(df[["Chromosome", "Start", "End"]]).merge().df
    return (m["End"] - m["Start"]).sum()


def curve_slow(df, n=100):
    """The notebook's version: merged sequence of hap_id 1..i, for i = 1..n."""
    return [merged_bp(df[df["hap_id"].isin(range(1, i + 1))]) if len(df) else 0 for i in range(1, n + 1)]


def curve(df, n=100):
    """Same result as curve_slow, in one pass: for each stretch of sequence, find
    the lowest hap_id (from 1 to n) covering it; the merged sequence of
    haplotypes 1..i is everything whose lowest covering hap_id is <= i."""
    if df is None or df.empty:
        return [0] * n
    df = df.rename(columns={"start": "Start", "end": "End"})
    df = df[(df["hap_id"] >= 1) & (df["hap_id"] <= n)]
    if df.empty:
        return [0] * n
    edges = np.unique(np.concatenate([df["Start"].to_numpy(), df["End"].to_numpy()]))
    first = np.full(len(edges) - 1, n + 1)
    lo = np.searchsorted(edges, df["Start"].to_numpy())
    hi = np.searchsorted(edges, df["End"].to_numpy())
    for h, a, b in zip(df["hap_id"].to_numpy(), lo, hi):
        np.minimum(first[a:b], h, out=first[a:b])
    width = np.diff(edges)
    by_hap = np.bincount(first, weights=width, minlength=n + 2)[1:n + 1]
    return list(np.cumsum(by_hap))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", default="no_bottleneck_donor")
    p.add_argument("--reps", type=int, default=1000)
    args = p.parse_args()
    d = SIM_DIR / args.scenario / "admix_csvs"

    plt.figure(figsize=(8, 6))
    for background, ls, suffix in (("eur", "-", " (EUR)"), ("amr", "--", " (AMR)")):
        for pop in POPS:
            ys = []
            for rep in range(1, args.reps + 1):
                path = d / f"{pop}_{background}_arc_ancestry_rep{rep}.csv"
                if not path.exists():
                    sys.exit(f"Missing {path}\nRun 01_archaic_by_ancestry.py {args.scenario} {pop} 1-{args.reps}")
                try:
                    df = pd.read_csv(path)
                except pd.errors.EmptyDataError:
                    df = pd.DataFrame()
                ys.append(curve(df))
            ys = np.array(ys, dtype=float)
            mean = ys.mean(0)
            se = ys.std(0) / np.sqrt(len(ys))
            x = np.arange(1, ys.shape[1] + 1)
            plt.plot(x, mean, color=COLORS[pop], linewidth=2, linestyle=ls, label=f"{pop}{suffix}")
            plt.fill_between(x, mean - 1.96 * se, mean + 1.96 * se, color=COLORS[pop], alpha=0.2)
            print(f"{pop} {background}: done", flush=True)

    plt.xlabel("Number of Individuals")
    plt.ylabel("Recovered Archaic Sequence (bp)")
    plt.ylim(0, 2.1e7)
    plt.legend()
    plt.tight_layout()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "simulated_recovered_amr_eur_yay.png"
    plt.savefig(out, dpi=300)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
