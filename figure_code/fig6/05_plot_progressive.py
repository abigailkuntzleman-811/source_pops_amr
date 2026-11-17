#!/usr/bin/env python3
"""Figure 6: IBD against European admixture proportion across scenarios.

One least squares fit per simulated scenario, each with a 95% confidence band,
drawn as a progressive reveal: one PNG per step, adding one scenario at a time,
with axes fixed across frames so nothing jumps.

Reproduces plot_IBD_eur_corr_CI_double_donor_recent from Chapter_1.ipynb,
restricted to the four scenarios that appear in the figure.

Prerequisite: figure_code/fig6/04_collect_stats.py

Usage:
    python figure_code/fig6/05_plot_progressive.py
    python figure_code/fig6/05_plot_progressive.py --final-only

Output:
    results/figures/mxl_IBD_progressive_<POP>_step<N>_<scenario>.png
    results/tables/mxl_IBD_progressive_fits.csv
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib.ticker import FuncFormatter

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import FIGURE_DIR, SIM_DIR, TABLE_DIR  # noqa: E402
from scenarios import PLOT_ORDER, SCENARIOS  # noqa: E402


def sci_formatter(val, _pos):
    scaled = val / 1e6
    if float(scaled).is_integer():
        return rf"{int(scaled)} × $10^{{6}}$"
    return rf"{scaled:.1f} × $10^{{6}}$"


def fit(x: np.ndarray, y: np.ndarray):
    """OLS fit with a 95% confidence band over the observed x range."""
    model = sm.OLS(y, sm.add_constant(x)).fit()
    x_pred = np.linspace(x.min(), x.max(), 200)
    prediction = model.get_prediction(sm.add_constant(x_pred))
    lo, hi = prediction.conf_int().T
    return model, x_pred, prediction.predicted_mean, lo, hi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pop", default="MXL")
    parser.add_argument("--stats", type=Path, default=None)
    parser.add_argument("--final-only", action="store_true",
                        help="write only the last frame, with all scenarios")
    args = parser.parse_args()

    stats_path = args.stats or (SIM_DIR / "simulation_stats.csv")
    if not stats_path.exists():
        sys.exit(
            f"Missing {stats_path}\n"
            "Run: python figure_code/fig6/04_collect_stats.py"
        )

    data = pd.read_csv(stats_path)

    layers, fits = [], []
    for scenario in PLOT_ORDER:
        subset = data[data["scenario"] == scenario]
        if subset.empty:
            print(f"warning: no data for {scenario}, omitted", file=sys.stderr)
            continue

        config = SCENARIOS[scenario]
        model, x_pred, y_pred, lo, hi = fit(subset["eur_prop"].to_numpy(),
                                            subset["ibd_bp"].to_numpy())

        layers.append(dict(scenario=scenario, x=x_pred, y=y_pred, lo=lo, hi=hi,
                           color=config["color"], ls=config["linestyle"],
                           label=config["label"]))

        print(f"{config['label']}: R² = {model.rsquared:.4f}, "
              f"slope = {model.params[1]:.3e}, p = {model.pvalues[1]:.3e}, "
              f"n = {len(subset):,}")

        fits.append({
            "scenario": scenario,
            "label": config["label"],
            "n": len(subset),
            "slope": model.params[1],
            "intercept": model.params[0],
            "r_squared": model.rsquared,
            "slope_p": model.pvalues[1],
        })

    if not layers:
        sys.exit("Nothing to plot.")

    # Fixed limits across frames so the axes do not move during the reveal.
    xmin = min(L["x"].min() for L in layers)
    xmax = max(L["x"].max() for L in layers)
    ymin = min(L["lo"].min() for L in layers)
    ymax = max(L["hi"].max() for L in layers)
    ypad = 0.05 * (ymax - ymin)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    steps = [len(layers)] if args.final_only else range(1, len(layers) + 1)

    fig, ax = plt.subplots(figsize=(7, 5))
    for n in steps:
        ax.clear()
        for layer in layers[:n]:
            ax.plot(layer["x"], layer["y"], color=layer["color"],
                    linewidth=2, linestyle=layer["ls"], label=layer["label"])
            ax.fill_between(layer["x"], layer["lo"], layer["hi"],
                            color=layer["color"], alpha=0.25)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ax.yaxis.set_major_formatter(FuncFormatter(sci_formatter))
        ax.set_xlabel("European admixture proportion", fontsize=14)
        ax.set_ylabel("IBD shared with other individuals", fontsize=14)
        ax.legend(loc="upper left")
        fig.tight_layout()

        out = (FIGURE_DIR /
               f"mxl_IBD_progressive_{args.pop}_step{n}_{layers[n-1]['scenario']}.png")
        fig.savefig(out, dpi=200)
        print(f"wrote {out.name}")

    plt.close(fig)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    fits_out = TABLE_DIR / "mxl_IBD_progressive_fits.csv"
    pd.DataFrame(fits).to_csv(fits_out, index=False)
    print(f"wrote {fits_out}")


if __name__ == "__main__":
    main()
