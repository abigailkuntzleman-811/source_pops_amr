#!/usr/bin/env python3
"""Figure S5 statistics: each simulated admixed population vs simulated IBS,
one value per simulation replicate.

Pairs within a replicate share simulated individuals, so they are not
independent; replicates are. Each replicate of a population's simulation
(standard model) contains both that admixed population and IBS, so the
comparison is paired by replicate:

    difference_r = statistic(admixed sample, replicate r) - statistic(IBS sample, replicate r)

for statistic = mean IBD per haplotype pair including zero pairs (mean_all) or
among pairs with any IBD (mean_nonzero). Reported: means over replicates, the
mean difference with a 95% interval (mean +- 1.96 SE over replicates), a
two-sided Wilcoxon signed-rank p value on the replicate differences, and Holm
correction over the four populations (separately for each statistic).

Needs $TABLE_DIR/IBD_boxplots_sims_standard_replicates.csv from
03_plot_sim_ibd_per_pair.py.

Usage:
    python figure_code/fig6/04_sim_ibd_stats.py

Output: $TABLE_DIR/IBD_boxplots_sims_standard_stats.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import TABLE_DIR  # noqa: E402

ADMIXED = ["CLM", "MXL", "PEL", "PUR"]
STATS = ["mean_all", "mean_nonzero"]


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, k in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[k]))
        adj[k] = running
    return adj


def main():
    path = TABLE_DIR / "IBD_boxplots_sims_standard_replicates.csv"
    if not path.exists():
        sys.exit(f"Missing {path}\nRun: python figure_code/fig6/03_plot_sim_ibd_per_pair.py --reps 1000")
    reps = pd.read_csv(path)

    rows = []
    for stat in STATS:
        block = []
        for pop in ADMIXED:
            sim = reps[reps["simulation"] == pop]
            w = sim.pivot(index="rep", columns="sample", values=stat)[[pop, "IBS"]].dropna()
            d = (w[pop] - w["IBS"]).to_numpy()
            se = d.std(ddof=1) / np.sqrt(len(d))
            block.append({
                "statistic": stat, "population": pop, "n_replicates": len(d),
                "value": w[pop].mean(), "value_sd_over_reps": w[pop].std(ddof=1),
                "IBS_value": w["IBS"].mean(),
                "difference_vs_IBS": d.mean(),
                "diff_ci_low": d.mean() - 1.96 * se, "diff_ci_high": d.mean() + 1.96 * se,
                "p": wilcoxon(d).pvalue if np.any(d != 0) else 1.0,
            })
        for r, adj in zip(block, holm([r["p"] for r in block])):
            r["p_holm"] = adj
        rows += block

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "IBD_boxplots_sims_standard_stats.csv", index=False)
    print(out.round(4).to_string(index=False))
    print(f"\nwrote {TABLE_DIR / 'IBD_boxplots_sims_standard_stats.csv'}")


if __name__ == "__main__":
    main()
