#!/usr/bin/env python3
"""Figure 6: European admixture proportion against within-population IBD,
under the standard model and models with one, two or six European donor
pulses. Built up one model at a time (one PNG per step), as in the
progressive version in Chapter_1.ipynb.

A least squares line is fitted separately in each replicate; each plotted line
is the mean of those per-replicate lines and its band the 95% confidence
interval across replicates (mean +- t * SD / sqrt(replicates)), as in Figure S10. IBD is computed with sim_ibd.ibd_per_individual
(pieces joined before a 2 cM cutoff, in cM); see sim_ibd.py.

Usage:
    python figure_code/fig6/01_plot_donor_ibd.py                 # MXL, 1000 replicates
    python figure_code/fig6/01_plot_donor_ibd.py --reps 5        # quick look
    python figure_code/fig6/01_plot_donor_ibd.py --ymax 400      # fix the upper y limit (cM)
    python figure_code/fig6/01_plot_donor_ibd.py --replot        # redraw only, from the saved values

Output:
    $FIGURE_DIR/<pop>_IBD_progressive_<POP>_step<N>_<model>.png   (step4_plus_six is Figure 6)
    $TABLE_DIR/<POP>_IBD_progressive_slopes.csv   slope per model, mean and SD over replicates
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
from scipy.stats import t as t_dist

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURE_DIR, SIM_DIR, TABLE_DIR  # noqa: E402
from sim_ibd import MODELS, admix_props, ibd_per_individual  # noqa: E402

STYLE = {
    "Standard":                  dict(color="#0072B2", ls="-",          tag="standard"),
    "Single 10% European donor": dict(color="#CC79A7", ls=":",          tag="plus_single"),
    "Two 10% European donors":   dict(color="#009E73", ls="-.",         tag="plus_double"),
    "Six 10% European donors":   dict(color="#D55E00", ls=(0, (5, 1)),  tag="plus_six"),
}


def fit(x, y):
    model = sm.OLS(y, sm.add_constant(x)).fit()
    xp = np.linspace(x.min(), x.max(), 200)
    pred = model.get_prediction(sm.add_constant(xp))
    lo, hi = pred.conf_int().T
    return xp, pred.predicted_mean, lo, hi, model


def replicate_band(coefs, xmin, xmax):
    """Mean of the per-replicate fitted lines and its 95% CI across replicates.
    coefs: one (slope, intercept) row per replicate."""
    coefs = np.asarray(coefs)
    xp = np.linspace(xmin, xmax, 200)
    lines = coefs[:, [1]] + coefs[:, [0]] * xp
    mean = lines.mean(0)
    half = t_dist.ppf(0.975, len(lines) - 1) * lines.std(0, ddof=1) / np.sqrt(len(lines))
    return xp, mean, mean - half, mean + half, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pop", default="MXL")
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--ymax", type=float, default=None,
                        help="upper y-axis limit in cM (default: fit all four "
                             "confidence bands, with a little headroom)")
    parser.add_argument("--ymin", type=float, default=1.0,
                        help="lower y-axis limit in cM (default 1)")
    parser.add_argument("--replot", action="store_true",
                        help="redraw from the saved per-individual values (<POP>_IBD_progressive_data.npz) "
                             "instead of reading the tree sequences again")
    args = parser.parse_args()
    saved = TABLE_DIR / f"{args.pop}_IBD_progressive_data.npz"

    if args.replot:
        if not saved.exists():
            sys.exit(f"Missing {saved}; run once without --replot first.")
        d = np.load(saved)
        if "0_coef" not in d.files:
            sys.exit(f"{saved} is from before the per-replicate band; run once without --replot.")
        plot(args, {m: replicate_band(d[f"{k}_coef"], d[f"{k}_x"].min(), d[f"{k}_x"].max())
                    for k, m in enumerate(MODELS)})
        return

    # Every model must have every replicate, or the models aren't comparable.
    # Check before starting rather than failing hours in.
    missing = []
    for model, (folder, tracts) in MODELS.items():
        for i in range(1, args.reps + 1):
            for path in (SIM_DIR / folder / "tree_seqs" / f"{args.pop}_rep{i}.ts.gz",
                         SIM_DIR / folder / "admix_csvs" / tracts.format(pop=args.pop, i=i)):
                if not path.exists():
                    missing.append(path)
    if missing:
        sys.exit(f"{len(missing)} input file(s) missing for {args.reps} replicates, e.g.:\n  "
                 + "\n  ".join(map(str, missing[:5]))
                 + "\nRun fewer replicates with --reps, or finish the simulations.")

    data = {m: {"x": [], "y": []} for m in MODELS}
    slopes = {m: [] for m in MODELS}
    coefs = {m: [] for m in MODELS}
    for i in range(1, args.reps + 1):
        print(f"replicate {i}/{args.reps}", flush=True)
        for model, (folder, tracts) in MODELS.items():
            ts = tskit.load(SIM_DIR / folder / "tree_seqs" / f"{args.pop}_rep{i}.ts.gz")
            adm = pd.read_csv(SIM_DIR / folder / "admix_csvs" / tracts.format(pop=args.pop, i=i))
            props = admix_props(adm, ts, args.pop)
            ibd = ibd_per_individual(ts, args.pop)
            x = np.array([props[k] for k in props])
            y = np.array([ibd[k] for k in props])
            data[model]["x"] += list(x)
            data[model]["y"] += list(y)
            coef = np.polyfit(x, y, 1)  # slope, intercept
            coefs[model].append(coef)
            slopes[model].append(coef[0])

    fits = {m: fit(np.array(d["x"]), np.array(d["y"])) for m, d in data.items()}
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    keep = {f"{k}_{c}": np.array(data[m][c]) for k, m in enumerate(MODELS) for c in "xy"}
    keep.update({f"{k}_coef": np.array(coefs[m]) for k, m in enumerate(MODELS)})
    np.savez_compressed(saved, **keep)

    # slopes: pooled fit, and per-replicate mean +/- SD (the honest uncertainty)
    rows = [{"model": m, "pooled_slope": fits[m][4].params[1],
             "pooled_R2": fits[m][4].rsquared,
             "rep_slope_mean": np.mean(slopes[m]),
             "rep_slope_sd": np.std(slopes[m], ddof=1),
             "rep_slope_se": np.std(slopes[m], ddof=1) / np.sqrt(len(slopes[m])),
             "mean_ibd_cM": np.mean(data[m]["y"])} for m in MODELS]
    table = pd.DataFrame(rows)
    print("\n" + table.to_string(index=False))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE_DIR / f"{args.pop}_IBD_progressive_slopes.csv", index=False)
    # the figure shows the per-replicate band; the pooled fits stay in the table
    plot(args, {m: replicate_band(coefs[m], min(data[m]["x"]), max(data[m]["x"])) for m in MODELS})


def plot(args, fits):
    # One y range for every step, so the axes don't jump as models are added:
    # from 0 (or the lowest band, if below 0) to the highest band plus 5%.
    lowest = min(fits[m][2].min() for m in MODELS)
    highest = max(fits[m][3].max() for m in MODELS)
    pad = 0.05 * (highest - min(lowest, 0))
    ymin = args.ymin
    if lowest < ymin:
        print(f"\nnote: the confidence bands go down to {lowest:.2f} cM, below --ymin {ymin}; "
              f"parts of the plot will be cut off", file=sys.stderr)
    ymax = args.ymax if args.ymax is not None else highest + pad
    if highest > ymax:
        print(f"\nnote: the confidence bands reach {highest:.1f} cM, above --ymax {ymax}; "
              f"parts of the plot will be cut off", file=sys.stderr)

    xmin = min(fits[m][0].min() for m in MODELS)
    xmax = max(fits[m][0].max() for m in MODELS)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 7.5))
    for n, current in enumerate(MODELS, start=1):
        ax.clear()
        for m in list(MODELS)[:n]:
            xp, yp, lo, hi, _ = fits[m]
            ax.plot(xp, yp, color=STYLE[m]["color"], linestyle=STYLE[m]["ls"],
                    linewidth=2, label=m)
            ax.fill_between(xp, lo, hi, color=STYLE[m]["color"], alpha=0.25)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel("European admixture proportion", fontsize=18)
        ax.set_ylabel("Total IBD shared with other haplotypes (cM)", fontsize=18)
        ax.tick_params(labelsize=14)
        # inside, upper left (white box so it stays readable over a band)
        ax.legend(loc="upper left", fontsize=14, frameon=True, facecolor="white", framealpha=0.9)
        fig.tight_layout()
        # same names as the figure in Overleaf (e.g. mxl_IBD_progressive_MXL_step4_plus_six.png)
        out = FIGURE_DIR / f"{args.pop.lower()}_IBD_progressive_{args.pop}_step{n}_{STYLE[current]['tag']}.png"
        fig.savefig(out, dpi=300)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
