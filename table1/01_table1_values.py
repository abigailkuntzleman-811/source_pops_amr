#!/usr/bin/env python3
"""Table 1: European admixture proportion and shared archaic ancestry.

Two columns, from two sources:

  Shared archaic ancestry (bp)
      Mean +/- sample SD of pairwise shared bp over all haplotype pairs within
      a population, read from the within-population files written by
      fig1/03_pairwise_sharing.sh. Depends only on the public Skov (2024)
      calls, so this column is fully reproducible from this repository.

  European admixture proportion
      Mean +/- sample SD of the EUR ancestry proportion, read from the output
      of 00_ancestry_proportions.py. That step needs local ancestry inference
      files, which are NOT publicly available. Without them this column is
      left blank and the rest of the table is still produced.

Prerequisites:

    bash   fig1/03_pairwise_sharing.sh neanderthal      # sharing values
    python table1/00_ancestry_proportions.py            # admixture column

Usage:
    python table1/01_table1_values.py
    python table1/01_table1_values.py --archaic denisovan
    python table1/01_table1_values.py --unit haplotype
    python table1/01_table1_values.py --exclude-zeros

The --unit flag must match whichever unit 00_ancestry_proportions.py was run
with. The default, "individual", reproduces the published Table 1 exactly;
"haplotype" gives the same means and very slightly different SDs.

Output:
    results/tables/table1_<archaic>_<unit>.csv
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import (AMR_POPS, DATA_PROCESSED, PAIRWISE_DIR,  # noqa: E402
                    TABLE_DIR)

POP_NAMES = {
    "PEL": "Peruvians (PEL)",
    "MXL": "Mexicans (MXL)",
    "CLM": "Colombians (CLM)",
    "PUR": "Puerto Ricans (PUR)",
}



# get_distributions is copied from Chapter_1.ipynb (cell 17), with only the
# input directory taken from config. Note that a line with two fields (a pair
# with no third column) counts as 0, as in the original.
def get_distributions(archaic_type):
    input_dir = str(PAIRWISE_DIR / archaic_type)
    if not os.path.isdir(input_dir):
        sys.exit(
            f"Missing input directory: {input_dir}\n"
            f"Run: bash fig1/03_pairwise_sharing.sh {archaic_type}"
        )

    all_values = defaultdict(list)

    for fname in os.listdir(input_dir):
        filepath = os.path.join(input_dir, fname)

        # Match format like CLM_{ITU}.txt
        match_pairwise = re.match(r"^([A-Z]+)_\{([A-Z]+)\}\.txt$", fname)
        if match_pairwise:
            subpop1, subpop2 = match_pairwise.groups()
            pair_key = f"{subpop1}-{subpop2}"

            with open(filepath) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    value = int(parts[2]) if len(parts) == 3 else 0
                    all_values[pair_key].append(value)
            continue

        # Match format like all_MXL.txt
        match_all = re.match(r"^all_([A-Z]+)\.txt$", fname)
        if match_all:
            subpop = match_all.group(1)
            pair_key = f"{subpop}-{subpop}"

            with open(filepath) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    value = int(parts[2]) if len(parts) == 3 else 0
                    all_values[pair_key].append(value)

    # Make long-form DataFrame
    rows = []
    for pair, values in all_values.items():
        for v in values:
            rows.append((pair, v))

    df = pd.DataFrame(rows, columns=["subpop_pair", "value"])
    return df


def summarize_pairwise_sharing(df: pd.DataFrame, subpop_pair: str,
                               exclude_zeros: bool = False):
    """Mean and sample SD of shared bp for one population pair.

    Same calculation as summarize_pairwise_sharing in Chapter_1.ipynb (cell 16).
    """
    vals = df.loc[df["subpop_pair"] == subpop_pair, "value"]
    if exclude_zeros:
        vals = vals[vals > 0]
    if len(vals) == 0:
        return np.nan, np.nan, 0
    return float(np.mean(vals)), float(np.std(vals, ddof=1)), int(len(vals))


def summarize_admixture(pop: str, ancestry: str = "EUR",
                        unit: str = "individual"):
    """Mean and sample SD of one ancestry proportion.

    Reads the per-unit proportions written by 00_ancestry_proportions.py.
    Returns NaN if they have not been generated, so the rest of the table
    still builds from the public data alone.
    """
    path = (DATA_PROCESSED / "ancestry_proportions" / f"{pop}_{unit}.csv")
    if not path.exists():
        return float("nan"), float("nan"), 0

    props = pd.read_csv(path, index_col=0)
    column = f"{ancestry}_prop"
    if column not in props.columns:
        return float("nan"), float("nan"), 0

    values = props[column]
    return float(values.mean()), float(values.std(ddof=1)), int(len(values))


def fmt_bp(mean: float, sd: float) -> str:
    """Format as in the paper: 5.03 x 10^6 +/- 2.35 x 10^6."""
    if np.isnan(mean):
        return ""
    return f"{mean / 1e6:.2f} × 10⁶ ± {sd / 1e6:.2f} × 10⁶"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archaic", default="neanderthal",
                        choices=["all", "neanderthal", "denisovan"])
    parser.add_argument("--ancestry", default="EUR",
                        help="continental ancestry for the admixture column")
    parser.add_argument("--unit", default="individual",
                        choices=["individual", "haplotype"],
                        help="compute admixture proportions per individual "
                             "(default) or per haplotype")
    parser.add_argument("--exclude-zeros", action="store_true",
                        help="drop haplotype pairs that share nothing")
    args = parser.parse_args()

    shares = get_distributions(args.archaic)

    prop_dir = DATA_PROCESSED / "ancestry_proportions"
    if not any(prop_dir.glob(f"*_{args.unit}.csv")):
        print(f"note: no ancestry proportions found in {prop_dir};\n"
              f"      the {args.ancestry} admixture column will be blank.\n"
              f"      Run: python table1/00_ancestry_proportions.py "
              f"--unit {args.unit}\n",
              file=sys.stderr)

    rows = []
    for pop in AMR_POPS:
        mean_bp, sd_bp, n_pairs = summarize_pairwise_sharing(
            shares, f"{pop}-{pop}", args.exclude_zeros)
        mean_adm, sd_adm, n_units = summarize_admixture(
            pop, args.ancestry, args.unit)

        rows.append({
            "Population": POP_NAMES.get(pop, pop),
            f"{args.ancestry} admixture proportion":
                "" if np.isnan(mean_adm) else f"{mean_adm:.2f} ± {sd_adm:.2f}",
            f"Shared {args.archaic} ancestry (bp)": fmt_bp(mean_bp, sd_bp),
            f"n {args.unit}s": n_units,
            "n haplotype pairs": n_pairs,
            "mean_bp": "" if np.isnan(mean_bp) else round(mean_bp),
            "sd_bp": "" if np.isnan(sd_bp) else round(sd_bp),
        })

    table = pd.DataFrame(rows)

    # Order as in the paper: descending shared archaic ancestry.
    sortable = pd.to_numeric(table["mean_bp"], errors="coerce")
    table = table.loc[sortable.sort_values(ascending=False).index]

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    out = TABLE_DIR / f"table1_{args.archaic}_{args.unit}.csv"
    table.to_csv(out, index=False)

    display = table.drop(columns=["mean_bp", "sd_bp"])
    print(display.to_string(index=False))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
