#!/usr/bin/env python3
"""Figure S9: European admixture proportion against within-population IBD,
single-donor model vs. standard model, one panel per admixed population.

Same layout as plot_all_pops in Chapter_1.ipynb: each panel shows a least
squares line with its 95% confidence band for the donor (purple) and standard
(yellow) models. IBD is computed with sim_ibd.ibd_per_individual (pieces
joined before a 2 cM cutoff, in cM), replacing get_ibd_dict.

Two tests of whether the slope differs between the models are reported:
  pooled_interaction_p   OLS interaction term over all individuals pooled
                         across replicates, as in the notebook. Treats every
                         individual as independent, so it overstates certainty.
  replicate_paired_p     paired t-test on the per-replicate slope difference
                         (donor - standard), treating each replicate as one
                         observation. The more conservative of the two.

Usage:
    python figure_code/fig6/02_plot_donor_vs_standard_all_pops.py              # 1000 replicates
    python figure_code/fig6/02_plot_donor_vs_standard_all_pops.py --reps 5     # quick look
    python figure_code/fig6/02_plot_donor_vs_standard_all_pops.py --replot     # redraw only, from the saved values

Output:
    $FIGURE_DIR/sims_IBD_admix_prop.png
    $TABLE_DIR/sims_IBD_admix_prop_stats.csv
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import tskit
from matplotlib.lines import Line2D
from scipy.stats import ttest_rel

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURE_DIR, SIM_DIR, TABLE_DIR  # noqa: E402
from sim_ibd import admix_props, ibd_per_individual  # noqa: E402

POPS = ["CLM", "MXL", "PEL", "PUR"]

# model: (folder under SIM_DIR, admixture tract file pattern, colour), from the notebook
MODELS = {
    "Donor":    ("eur_donor",           "{pop}_10tracts_rep{i}_stitched.csv.gz", "#7E6B8F"),
    "Standard": ("no_bottleneck_donor", "{pop}_7tracts_rep{i}_stitched.csv.gz",  "#DCC48E"),
}


def paths(model, pop, i):
    folder, tracts, _ = MODELS[model]
    return (SIM_DIR / folder / "tree_seqs" / f"{pop}_rep{i}.ts.gz",
            SIM_DIR / folder / "admix_csvs" / tracts.format(pop=pop, i=i))


def fit(x, y):
    model = sm.OLS(y, sm.add_constant(x)).fit()
    xp = np.linspace(x.min(), x.max(), 200)
    pred = model.get_prediction(sm.add_constant(xp))
    lo, hi = pred.conf_int().T
    return xp, pred.predicted_mean, lo, hi, model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--replot", action="store_true",
                        help="redraw from the saved per-sample values (sims_IBD_admix_prop_data.npz) "
                             "instead of reading the tree sequences again")
    args = parser.parse_args()
    saved = TABLE_DIR / "sims_IBD_admix_prop_data.npz"

    if args.replot:
        if not saved.exists():
            sys.exit(f"Missing {saved}; run once without --replot first.")
        d = np.load(saved)
        plot({pop: {m: fit(d[f"{pop}_{k}_x"], d[f"{pop}_{k}_y"]) for k, m in enumerate(MODELS)}
              for pop in POPS})
        return

    # every model/population must have every replicate
    missing = [p for pop in POPS for m in MODELS for i in range(1, args.reps + 1)
               for p in paths(m, pop, i) if not p.exists()]
    if missing:
        sys.exit(f"{len(missing)} input file(s) missing for {args.reps} replicates, e.g.:\n  "
                 + "\n  ".join(map(str, missing[:5]))
                 + "\nRun fewer replicates with --reps, or finish the simulations.")

    results, rows, keep = {}, [], {}
    for pop in POPS:
        data = {m: {"x": [], "y": []} for m in MODELS}
        rep_slopes = {m: [] for m in MODELS}
        for i in range(1, args.reps + 1):
            if i % 50 == 0 or i == 1:
                print(f"{pop}: replicate {i}/{args.reps}", flush=True)
            for m in MODELS:
                ts_path, tract_path = paths(m, pop, i)
                ts = tskit.load(ts_path)
                props = admix_props(pd.read_csv(tract_path), ts, pop)
                ibd = ibd_per_individual(ts, pop)
                x = np.array([props[k] for k in props])
                y = np.array([ibd[k] for k in props])
                data[m]["x"] += list(x)
                data[m]["y"] += list(y)
                rep_slopes[m].append(np.polyfit(x, y, 1)[0])

        fits = {m: fit(np.array(d["x"]), np.array(d["y"])) for m, d in data.items()}
        for k, m in enumerate(MODELS):
            keep[f"{pop}_{k}_x"] = np.array(data[m]["x"])
            keep[f"{pop}_{k}_y"] = np.array(data[m]["y"])

        # pooled interaction test, as in the notebook
        df = pd.concat([pd.DataFrame({"IBD": data[m]["y"], "Prop": data[m]["x"],
                                      "Pop": float(m == "Standard")}) for m in MODELS])
        df["Interaction"] = df["Prop"] * df["Pop"]
        interact = sm.OLS(df["IBD"], sm.add_constant(df[["Prop", "Pop", "Interaction"]])).fit()

        # replicate-level test: paired difference in slope per replicate
        diff = np.array(rep_slopes["Donor"]) - np.array(rep_slopes["Standard"])
        paired_p = ttest_rel(rep_slopes["Donor"], rep_slopes["Standard"]).pvalue if args.reps > 1 else np.nan

        results[pop] = fits
        row = {"pop": pop}
        for m in MODELS:
            row[f"{m.lower()}_pooled_slope"] = fits[m][4].params[1]
            row[f"{m.lower()}_pooled_R2"] = fits[m][4].rsquared
            row[f"{m.lower()}_pooled_p"] = fits[m][4].pvalues[1]
            row[f"{m.lower()}_rep_slope_mean"] = np.mean(rep_slopes[m])
            row[f"{m.lower()}_rep_slope_sd"] = np.std(rep_slopes[m], ddof=1) if args.reps > 1 else np.nan
        row["pooled_interaction_p"] = interact.pvalues["Interaction"]
        row["slope_diff_rep_mean"] = diff.mean()
        row["replicate_paired_p"] = paired_p
        rows.append(row)

        print(f"{pop}:")
        for m in MODELS:
            print(f"  {m}: R² = {fits[m][4].rsquared:.3f}, slope = {fits[m][4].params[1]:.2f} cM, "
                  f"p = {fits[m][4].pvalues[1]:.3e}")
        print(f"  slope difference: pooled interaction p = {interact.pvalues['Interaction']:.3e}, "
              f"replicate-paired p = {paired_p:.3e}\n")

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TABLE_DIR / "sims_IBD_admix_prop_stats.csv", index=False)
    np.savez_compressed(saved, **keep)
    plot(results)


def plot(results):
    # shared y range across panels: 0 to the highest band + 5%
    highest = max(f[3].max() for fits in results.values() for f in fits.values())
    lowest = min(f[2].min() for fits in results.values() for f in fits.values())
    pad = 0.05 * (highest - min(lowest, 0))

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, pop in zip(axes.flatten(), POPS):
        for m in MODELS:
            xp, yp, lo, hi, _ = results[pop][m]
            ax.plot(xp, yp, color=MODELS[m][2], linewidth=2)
            ax.fill_between(xp, lo, hi, color=MODELS[m][2], alpha=0.25)
        ax.set_title(pop)
        ax.set_ylim(min(0, lowest - pad), highest + pad)

    fig.supxlabel("European admixture proportion", fontsize=16)
    fig.supylabel("Total IBD shared with other haplotypes (cM, per 50 Mb)", fontsize=16)
    handles = [Line2D([0], [0], color=MODELS[m][2], lw=2, label=m) for m in MODELS]
    # inside the upper-left panel, in its upper-left corner
    axes[0, 0].legend(handles=handles, loc="upper left", fontsize=14,
                      frameon=True, facecolor="white", framealpha=0.9)
    plt.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "sims_IBD_admix_prop.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
