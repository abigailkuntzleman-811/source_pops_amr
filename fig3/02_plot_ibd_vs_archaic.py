#!/usr/bin/env python3
"""Figure 3: average shared archaic ancestry versus average shared IBD.

One point per population: mean shared archaic sequence between haplotype pairs
within the population (x) against mean shared IBD between haplotype pairs
within the population (y). Points are coloured and shaped by superpopulation.
Reproduces Chapter_1.ipynb cell 35.

Both axes are within-population averages, so this is the diagonal of the
Figure 1 heatmap plotted against the IBD equivalent.

Prerequisites:
    bash   fig1/03_pairwise_sharing.sh all        # archaic sharing
    python fig3/01_within_population_ibd.py        # IBD sharing

Usage:
    python fig3/02_plot_ibd_vs_archaic.py
    python fig3/02_plot_ibd_vs_archaic.py --archaic neanderthal

Output:
    results/figures/archaic_share_vs_ibd.png
    results/tables/archaic_share_vs_ibd.csv       the plotted values
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import FIGURE_DIR, IBD_DIR, PAIRWISE_DIR, TABLE_DIR  # noqa: E402

SUPERPOP = {
    "MXL": "AMR", "PEL": "AMR", "PUR": "AMR", "CLM": "AMR",
    "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "CDX": "EAS", "KHV": "EAS",
    "CEU": "EUR", "TSI": "EUR", "GBR": "EUR", "FIN": "EUR", "IBS": "EUR",
    "GIH": "SAS", "PJL": "SAS", "BEB": "SAS", "STU": "SAS", "ITU": "SAS",
}

STYLE = {
    "AMR": {"color": "#C45425", "marker": "o"},   # burnt sienna
    "EAS": {"color": "#2D6A4F", "marker": "s"},   # deep forest green
    "EUR": {"color": "#1F4E8C", "marker": "^"},   # deep slate blue
    "SAS": {"color": "#6A0DAD", "marker": "D"},   # deep violet
}

SUB_POPS = [
    "CHB", "JPT", "CHS", "CDX", "KHV",
    "CEU", "TSI", "GBR", "FIN", "IBS",
    "GIH", "PJL", "BEB", "STU", "ITU",
    "MXL", "PEL", "PUR", "CLM",
]


def mean_archaic_sharing(pop: str, archaic: str) -> float:
    """Mean shared archaic bp over within-population haplotype pairs.

    Same calculation as get_averages in the notebook: only lines where both
    IDs carry the sorted- prefix are counted.
    """
    path = PAIRWISE_DIR / archaic / f"all_{pop}.txt"
    if not path.exists():
        return float("nan")

    total, n = 0, 0
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            id1, id2, value = parts
            if "sorted" in id1 and "sorted" in id2:
                total += int(value)
                n += 1

    return total / n if n else float("nan")


def load_ibd_means() -> pd.Series:
    """Mean within-population IBD per haplotype pair, in cM."""
    path = IBD_DIR / "within_population_summary.csv"
    if not path.exists():
        sys.exit(
            f"Missing {path}\nRun: python fig3/01_within_population_ibd.py"
        )
    return pd.read_csv(path).set_index("pop")["mean_cM"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archaic", default="all",
                        choices=["all", "neanderthal", "denisovan"])
    parser.add_argument("--annotate", nargs="*", default=None,
                        help="populations to label; default is FIN plus all AMR")
    args = parser.parse_args()

    ibd_means = load_ibd_means()

    records = []
    for pop in SUB_POPS:
        records.append({
            "pop": pop,
            "superpop": SUPERPOP[pop],
            "archaic_bp": mean_archaic_sharing(pop, args.archaic),
            "ibd_cM": ibd_means.get(pop, float("nan")),
        })

    df = pd.DataFrame(records)

    missing = df[df[["archaic_bp", "ibd_cM"]].isna().any(axis=1)]["pop"].tolist()
    if missing:
        print(f"warning: no data for {', '.join(missing)}; omitted",
              file=sys.stderr)
    df = df.dropna(subset=["archaic_bp", "ibd_cM"])
    if df.empty:
        sys.exit("Nothing to plot.")

    if args.annotate is None:
        to_label = set(df.loc[df["superpop"] == "AMR", "pop"]) | {"FIN"}
    else:
        to_label = set(args.annotate)

    fig, ax = plt.subplots(figsize=(8, 6))

    for superpop, style in STYLE.items():
        group = df[df["superpop"] == superpop]
        if group.empty:
            continue
        ax.scatter(group["archaic_bp"], group["ibd_cM"],
                   color=style["color"], marker=style["marker"],
                   label=superpop, s=60)

    for _, row in df.iterrows():
        if row["pop"] in to_label:
            ax.annotate(row["pop"], (row["archaic_bp"], row["ibd_cM"]),
                        textcoords="offset points", xytext=(3, 7),
                        ha="center", fontsize=10)

    ax.set_xlabel("Average archaic sharing (bp)")
    ax.set_ylabel("Average shared IBD (cM)")
    ax.set_xlim(df["archaic_bp"].min() - 1_000_000,
                df["archaic_bp"].max() + 1_000_000)
    ax.set_ylim(df["ibd_cM"].min() - 1, df["ibd_cM"].max() + 1)
    ax.legend(title="Superpopulation", frameon=False)

    fig.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "archaic_share_vs_ibd.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"wrote {out}")

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    values_out = TABLE_DIR / "archaic_share_vs_ibd.csv"
    df.to_csv(values_out, index=False)
    print(f"wrote {values_out}")


if __name__ == "__main__":
    main()
