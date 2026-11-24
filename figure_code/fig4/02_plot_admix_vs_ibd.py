#!/usr/bin/env python3
"""Figure 4: continental admixture proportion against per-individual IBD.

A 2x2 panel, one per admixed American population. Each point is an individual:
their proportion of a given continental ancestry (x) against the average cM
they share with the rest of their population (y). Each panel carries a least
squares line with a 95% confidence band, and is titled with Spearman's rho and
its p value, highlighted when significant.

Reproduces admix_prop_vs_ibd in Chapter_1.ipynb (cell 53).

Prerequisites:
    python table1/00_ancestry_proportions.py          # needs LAI_DIR
    python figure_code/fig4/01_individual_ibd.py      # needs IBD output

Usage:
    python figure_code/fig4/02_plot_admix_vs_ibd.py
    python figure_code/fig4/02_plot_admix_vs_ibd.py --ancestry AMR
    python figure_code/fig4/02_plot_admix_vs_ibd.py --alpha 0.4      # more see-through points

Output:
    results/figures/<ANCESTRY>_prop_admix_by_ibd.png
    results/tables/<ANCESTRY>_prop_admix_by_ibd_stats.csv
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr
from scipy.stats import t as t_dist

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))

from config import DATA_PROCESSED, FIGURE_DIR, IBD_DIR, TABLE_DIR  # noqa: E402

# Share one definition of who counts as related with the other IBD scripts.
sys.path.insert(0, str(_root / "figure_code" / "fig3"))
_load_related_samples = __import__("01_within_population_ibd").load_related_samples

_related_cache = None


def load_related_samples() -> set:
    """Related sample IDs, read once and reused across panels."""
    global _related_cache
    if _related_cache is None:
        _related_cache = _load_related_samples()
    return _related_cache

# Panel order as in the paper. Titles use the population codes.
POP_LABELS = {
    "CLM": "CLM",
    "MXL": "MXL",
    "PEL": "PEL",
    "PUR": "PUR",
}

COLORS = {
    "CLM": "#A3D9FF",
    "MXL": "#DA3E52",
    "PEL": "#96E6B3",
    "PUR": "#7E6B8F",
}

AXIS_LABEL = {
    "EUR": "European admixture proportion",
    "AMR": "Indigenous American admixture proportion",
    "AFR": "African admixture proportion",
}


def load_panel_data(pop: str, ancestry: str, unit: str) -> pd.DataFrame:
    """Ancestry proportion and average IBD for every individual of one pop."""
    props_path = (DATA_PROCESSED / "ancestry_proportions" / f"{pop}_{unit}.csv")
    if not props_path.exists():
        sys.exit(
            f"Missing {props_path}\n"
            f"Run: python table1/00_ancestry_proportions.py --unit {unit}"
        )

    props = pd.read_csv(props_path, index_col=0)
    column = f"{ancestry}_prop"
    if column not in props.columns:
        sys.exit(f"{props_path.name} has no {column} column. "
                 f"Available: {sorted(props.columns)}")

    ibd_path = IBD_DIR / "individual_ibd.csv"
    if not ibd_path.exists():
        sys.exit(
            f"Missing {ibd_path}\n"
            "Run: python figure_code/fig4/01_individual_ibd.py"
        )
    ibd = pd.read_csv(ibd_path)
    ibd = ibd[ibd["pop"] == pop].set_index("individual")["avg_inpop_ibd"]

    # Drop related individuals explicitly rather than relying on their absence
    # from the IBD file. They have ancestry calls but were excluded upstream,
    # so a join leaves them missing; filling that with zero would plot them as
    # sharing no IBD with anyone, which is not a measurement.
    related = load_related_samples()
    excluded = sorted(set(props.index) & related)
    if excluded:
        print(f"{pop}: excluding {len(excluded)} related individual(s): "
              f"{', '.join(excluded)}", file=sys.stderr)
        props = props.loc[~props.index.isin(related)]

    df = pd.DataFrame({
        "prop": props[column],
        "avg_inpop_ibd": ibd.reindex(props.index),
    })

    # Anyone still missing has ancestry calls but no IBD record at all, which
    # should not happen once related individuals are removed.
    missing = df[df["avg_inpop_ibd"].isna()]
    if not missing.empty:
        print(f"{pop}: warning, {len(missing)} individual(s) have ancestry "
              f"calls but no IBD value and are dropped: "
              f"{', '.join(missing.index[:5])}", file=sys.stderr)

    return df.dropna()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ancestry", default="EUR",
                        help="continental ancestry on the x axis")
    parser.add_argument("--unit", default="individual",
                        choices=["individual", "haplotype"],
                        help="unit used by 00_ancestry_proportions.py")
    parser.add_argument("--ylim", nargs=2, type=float, default=(-5, 50),
                        help="y axis limits (default: -5 50)")
    parser.add_argument("--alpha", type=float, default=0.6,
                        help="point opacity, 0 (invisible) to 1 (solid); default 0.6")
    args = parser.parse_args()

    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    stats = []

    for ax, (pop, label) in zip(axes.flatten(), POP_LABELS.items()):
        df = load_panel_data(pop, args.ancestry, args.unit)
        x = df["prop"].to_numpy()
        y = df["avg_inpop_ibd"].to_numpy()

        if len(x) < 3:
            ax.set_title(f"{label}\n(insufficient data)")
            continue

        ax.scatter(x, y, s=18, alpha=args.alpha, color=COLORS[pop])

        # regression, band and title exactly as admix_prop_vs_ibd in
        # Chapter_1.ipynb (cell 62)
        result = linregress(x, y)
        m, b = result.slope, result.intercept
        rho, p_spearman = spearmanr(x, y)

        # --- 95% confidence band for the regression line ---
        n = len(x)
        x_mean = x.mean()
        Sxx = np.sum((x - x_mean) ** 2)
        residual_std = np.sqrt(np.sum((y - (m * x + b)) ** 2) / (n - 2))
        t_crit = t_dist.ppf(0.975, df=n - 2)

        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = m * x_line + b
        se_line = residual_std * np.sqrt(1/n + (x_line - x_mean)**2 / Sxx)

        ax.fill_between(x_line, y_line - t_crit * se_line, y_line + t_crit * se_line,
                        color=COLORS[pop], alpha=0.2)
        ax.plot(x_line, y_line, color=COLORS[pop], linewidth=2)

        # title: red if significant
        significant = p_spearman < 0.05
        p_str = f"{p_spearman:.2e}" if significant else f"{p_spearman:.2f}"
        ax.set_title(f"{label}\n$\\rho$ = {rho:.2f}, p = {p_str}",
                     color="red" if significant else "black")

        ax.set_xlim(0, 1)
        ax.set_ylim(*args.ylim)

        print(f"{label}   R² = {result.rvalue ** 2:.2f}   "
              f"p = {result.pvalue:.3e}   (n={n})")
        print(f"{label}   Spearman rho = {rho:.2f}   p = {p_spearman:.3e}")

        stats.append({
            "pop": pop,
            "label": label,
            "n": n,
            "slope": result.slope,
            "intercept": result.intercept,
            "r_squared": result.rvalue ** 2,
            "linregress_p": result.pvalue,
            "spearman_rho": rho,
            "spearman_p": p_spearman,
        })

    fig.supxlabel(AXIS_LABEL.get(args.ancestry,
                                 f"{args.ancestry} admixture proportion"))
    fig.supylabel(r"$\overline{\mathrm{IBD}}_i$  (cM)")
    fig.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / f"{args.ancestry}_prop_admix_by_ibd.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"\nwrote {out}")

    if stats:
        TABLE_DIR.mkdir(parents=True, exist_ok=True)
        stats_out = TABLE_DIR / f"{args.ancestry}_prop_admix_by_ibd_stats.csv"
        pd.DataFrame(stats).to_csv(stats_out, index=False)
        print(f"wrote {stats_out}")


if __name__ == "__main__":
    main()
