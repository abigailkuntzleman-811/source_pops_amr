#!/usr/bin/env python3
"""Figure S5: distribution of IBD per pair in the standard-model simulations,
drawn to match the real-data boxplot (Figure S4).

For each admixed population, every pair of simulated (haploid) individuals in
every replicate contributes its total IBD in cM. IBS is taken from the same
simulations (each population's model includes an IBS sample) and pooled.
As in S4, pairs with no IBD are left out of the boxes and counted in the table.

IBD uses sim_ibd.ibd_pairs: pieces joined into tracts before a 2 cM cutoff
(the hap-ibd minimum), recent shared ancestors only.

Usage:
    python figure_code/fig6/03_plot_sim_ibd_per_pair.py              # 1000 replicates
    python figure_code/fig6/03_plot_sim_ibd_per_pair.py --reps 5     # quick look
    python figure_code/fig6/03_plot_sim_ibd_per_pair.py --replot     # redraw only, from the saved pairs

Output:
    $FIGURE_DIR/IBD_boxplots_sims_standard.png
    $TABLE_DIR/IBD_boxplots_sims_standard.csv
    $TABLE_DIR/IBD_boxplots_sims_standard_pairs.npz   every pair's IBD, for --replot
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURE_DIR, SIM_DIR, TABLE_DIR  # noqa: E402

STANDARD = "no_bottleneck_donor"
ADMIXED = ["CLM", "MXL", "PEL", "PUR"]
# same order and colours as Figure S4
LABELS = ["IBS", "CLM", "MXL", "PEL", "PUR"]
COLORS = ["#1f4e79", "#A3D9FF", "#DA3E52", "#96E6B3", "#7E6B8F"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--replot", action="store_true",
                        help="redraw from the saved pairs (IBD_boxplots_sims_standard_pairs.npz) "
                             "instead of reading the tree sequences again")
    args = parser.parse_args()
    saved = TABLE_DIR / "IBD_boxplots_sims_standard_pairs.npz"

    if args.replot:
        if not saved.exists():
            sys.exit(f"Missing {saved}; run once without --replot first.")
        d = np.load(saved)
        plot({lab: d[lab] for lab in LABELS})
        return

    import tskit
    from sim_ibd import ibd_pairs

    tree = lambda pop, i: SIM_DIR / STANDARD / "tree_seqs" / f"{pop}_rep{i}.ts.gz"
    missing = [tree(p, i) for p in ADMIXED for i in range(1, args.reps + 1)
               if not tree(p, i).exists()]
    if missing:
        sys.exit(f"{len(missing)} tree sequence(s) missing for {args.reps} replicates, e.g.:\n  "
                 + "\n  ".join(map(str, missing[:5]))
                 + "\nRun fewer replicates with --reps, or finish the simulations.")

    pairs = {lab: [] for lab in LABELS}
    rep_means = {lab: [] for lab in LABELS}
    for pop in ADMIXED:
        for i in range(1, args.reps + 1):
            if i % 50 == 0 or i == 1:
                print(f"{pop}: replicate {i}/{args.reps}", flush=True)
            ts = tskit.load(tree(pop, i))
            for lab in (pop, "IBS"):
                v = ibd_pairs(ts, lab)
                pairs[lab] += v
                rep_means[lab].append(np.mean(v))

    # table with the same columns as Figure S4, plus the replicate-level spread
    rows = []
    for lab in LABELS:
        v = np.array(pairs[lab])
        rows.append({
            "Population": lab,
            "Total pairs": len(v),
            "Pairs with 0 IBD": int((v == 0).sum()),
            "Mean IBD across all pairs": v.mean(),
            "Mean IBD excluding 0-IBD pairs": v[v > 0].mean() if (v > 0).any() else np.nan,
            "SD of per-replicate means": np.std(rep_means[lab], ddof=1) if len(rep_means[lab]) > 1 else np.nan,
        })
    table = pd.DataFrame(rows)
    print("\n" + table.round(4).to_string(index=False))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE_DIR / "IBD_boxplots_sims_standard.csv", index=False)
    np.savez_compressed(saved, **{lab: np.array(pairs[lab]) for lab in LABELS})
    plot(pairs)


def plot(pairs):
    # boxplot styled like Figure S4: zero pairs excluded, log y axis
    fig, ax = plt.subplots(figsize=(8, 6))
    box = ax.boxplot(
        [[x for x in pairs[lab] if x > 0] for lab in LABELS],
        patch_artist=True,
        showfliers=False,  # with this many pairs, outlier circles merge into solid bars
    )
    ax.set_xticks(range(1, len(LABELS) + 1))
    ax.set_xticklabels(LABELS)
    for patch, color in zip(box["boxes"], COLORS):
        patch.set_facecolor(color)
    ax.set_yscale("log")
    # plain numbers (2, 5, 10, 20, ...) rather than powers of ten
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel("Total IBD (cM) per Pair")

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "IBD_boxplots_sims_standard.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
