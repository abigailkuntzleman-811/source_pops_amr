#!/usr/bin/env python3
"""Figure 1: heatmap of pairwise shared archaic ancestry between populations.

Reads the pairwise sharing files written by 03_pairwise_sharing.sh, builds a
symmetric population-by-population matrix of mean shared bp per haplotype pair,
and draws the lower-triangle heatmap used in the paper. Diagonal entries are
within-population sharing.

Usage:
    python fig1/04_heatmap.py neanderthal
    python fig1/04_heatmap.py denisovan
    python fig1/04_heatmap.py all

Output:
    results/figures/<archaic>_shared_heatmap.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import FIGURE_DIR, PAIRWISE_DIR  # noqa: E402

# Display order for the heatmap axes: East Asian, European, South Asian,
# then admixed American.
SUB_POPS = [
    "CHB", "JPT", "CHS", "CDX", "KHV",
    "CEU", "TSI", "GBR", "FIN", "IBS",
    "GIH", "PJL", "BEB", "STU", "ITU",
    "MXL", "PEL", "PUR", "CLM",
]


def get_averages(path: Path):
    """Mean shared bp over haplotype pairs in one comparison file.

    Returns (n_pairs, mean_all, n_zero_pairs, mean_given_sharing).
    Only lines where both IDs carry the "sorted" prefix are counted, matching
    the original analysis.
    """
    values = []
    no_share = 0

    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            id1, id2, n = parts[0], parts[1], int(parts[2])
            if "sorted" in id1 and "sorted" in id2:
                values.append(n)
                if n == 0:
                    no_share += 1

    if not values:
        return 0, 0.0, 0, 0.0

    nonzero = [v for v in values if v > 0]
    mean_cond = sum(nonzero) / len(nonzero) if nonzero else 0.0
    return len(values), sum(values) / len(values), no_share, mean_cond


def build_matrix(archaic: str) -> np.ndarray:
    """Symmetric matrix of mean pairwise sharing, in SUB_POPS order."""
    directory = PAIRWISE_DIR / archaic
    if not directory.is_dir():
        sys.exit(
            f"Missing input directory: {directory}\n"
            f"Run: bash fig1/03_pairwise_sharing.sh {archaic}"
        )

    n = len(SUB_POPS)
    matrix = np.zeros((n, n))
    missing = []

    # Off-diagonal: between-population comparisons.
    for i in range(n):
        for j in range(i + 1, n):
            # The filename uses whichever population came first in the pairwise
            # run, so try both orderings.
            candidates = [
                directory / f"{SUB_POPS[i]}_{{{SUB_POPS[j]}}}.txt",
                directory / f"{SUB_POPS[j]}_{{{SUB_POPS[i]}}}.txt",
            ]
            path = next((p for p in candidates if p.is_file()), None)
            if path is None:
                missing.append(f"{SUB_POPS[i]}-{SUB_POPS[j]}")
                continue

            _, mean_all, _, _ = get_averages(path)
            matrix[i, j] = matrix[j, i] = mean_all

    # Diagonal: within-population comparisons.
    for i, pop in enumerate(SUB_POPS):
        path = directory / f"all_{pop}.txt"
        if path.is_file():
            _, mean_all, _, _ = get_averages(path)
            matrix[i, i] = mean_all
        else:
            missing.append(f"{pop}-{pop}")

    if missing:
        print(f"warning: {len(missing)} comparisons missing, left at zero: "
              f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}",
              file=sys.stderr)

    return matrix


def plot_comparisons(matrix: np.ndarray, archaic: str, out_path: Path) -> None:
    """Lower-triangle heatmap, matching the figure in the paper."""
    fig, ax = plt.subplots()

    # Mask the upper triangle; the matrix is symmetric.
    lower_mask = np.tri(matrix.shape[0], matrix.shape[1], k=-1)
    masked = np.ma.array(matrix, mask=lower_mask.T)

    im = ax.imshow(masked, cmap="YlGn", vmin=np.nanmin(matrix), vmax=np.nanmax(matrix))

    ax.set_xticks(range(len(SUB_POPS)), labels=SUB_POPS,
                  rotation=45, ha="right", rotation_mode="anchor", fontsize=12)
    ax.set_yticks(range(len(SUB_POPS)), labels=SUB_POPS, fontsize=12)

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x / 1e6:.1f}×10⁶")
    )

    plt.title(archaic.capitalize(), fontsize=16)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archaic", choices=["all", "neanderthal", "denisovan"])
    args = parser.parse_args()

    matrix = build_matrix(args.archaic)

    out_path = FIGURE_DIR / f"{args.archaic}_shared_heatmap.png"
    plot_comparisons(matrix, args.archaic, out_path)
    print(f"wrote {out_path}")

    # Save the underlying values so the figure can be checked without rerunning.
    csv_path = FIGURE_DIR / f"{args.archaic}_shared_heatmap_values.csv"
    np.savetxt(csv_path, matrix, delimiter=",",
               header=",".join(SUB_POPS), comments="")
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
