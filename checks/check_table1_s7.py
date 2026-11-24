#!/usr/bin/env python3
"""Check Table 1 and Figure S7. Quick; fine on the login node.

Table 1 (shared archaic column)
    Runs get_distributions and summarize_pairwise_sharing from Chapter_1.ipynb
    (cells 16-17, unchanged apart from the input folder) and compares them with
    results/tables/table1_neanderthal_individual.csv. Also counts lines in the
    pairwise files that have only two fields: the notebook counts those as 0,
    Figure 1's get_averages does not, so if any exist the two disagree.

Table 1 (admixture column) is checked by check_fig3_fig4.py (same proportions).

Figure S7
    The code that made the original <POP>_<ANC>_IBD_ratios.csv files is not in
    the notebook, so the comparison is against those files: number of pairs,
    share with zero IBD on shared ancestry, and mean ratio, old vs new. The new
    run uses masked IBD and the KING list, so expect close, not identical.

    python checks/check_table1_s7.py
    python checks/check_table1_s7.py --old-ratios /users/akuntzle/data/akuntzle/IBD_ancestry_tracts
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_root))
from config import AMR_POPS, PAIRWISE_DIR, TABLE_DIR  # noqa: E402


# ------------------------------------------------- Chapter_1.ipynb cell 16
def summarize_pairwise_sharing(df, subpop_pair, exclude_zeros=False):
    vals = df[df['subpop_pair'] == subpop_pair]['value']

    if exclude_zeros:
        vals = vals[vals > 0]

    mean_val = np.mean(vals)
    std_val = np.std(vals, ddof=1)  # sample standard deviation

    # Convert to millions
    mean_million = mean_val / 1e6
    std_million = std_val / 1e6

    print(f"{subpop_pair}: {mean_million:.2f} ×10⁶ ± {std_million:.2f} ×10⁶")

    return mean_val, std_val


# ------------------------------------------------- Chapter_1.ipynb cell 17
def get_distributions(archaic_type):
    input_dir = str(PAIRWISE_DIR / archaic_type)

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
# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-ratios", default="/users/akuntzle/data/akuntzle/IBD_ancestry_tracts")
    args = ap.parse_args()
    problems = []

    print("=== Table 1: shared Neanderthal ancestry")
    nean_shares = get_distributions('neanderthal')
    repo_path = TABLE_DIR / "table1_neanderthal_individual.csv"
    repo = pd.read_csv(repo_path) if repo_path.exists() else None
    for pop in AMR_POPS:
        mean_val, std_val = summarize_pairwise_sharing(nean_shares, f'{pop}-{pop}')
        if repo is None:
            continue
        r = repo[repo["Population"].str.contains(f"({pop})", regex=False)].iloc[0]
        ok = round(mean_val) == r["mean_bp"] and round(std_val) == r["sd_bp"]
        print(f"  {'OK  ' if ok else 'DIFF'} {pop}: notebook {mean_val:,.0f} ± {std_val:,.0f} "
              f"vs repo {r['mean_bp']:,.0f} ± {r['sd_bp']:,.0f}")
        if not ok:
            problems.append(f"Table 1 {pop}")
    if repo is None:
        print(f"  (no {repo_path.name} yet; run python table1/01_table1_values.py)")

    counts = defaultdict(int)
    for path in (PAIRWISE_DIR / "neanderthal").glob("*.txt"):
        with open(path) as f:
            for line in f:
                counts[len(line.strip().split("\t"))] += 1
    print(f"  lines by number of fields in the pairwise files: {dict(counts)}")
    if counts.get(2):
        print("  NOTE: two-field lines exist; Table 1 counts them as 0, Figure 1's "
              "get_averages reuses the previous line's value for them")

    print("\n=== Figure S7: old ratio files vs new run")
    new_dir = TABLE_DIR
    old_dir = Path(args.old_ratios)
    for pop in AMR_POPS:
        new_path = new_dir / f"ancestry_aware_ibd_{pop}.csv.gz"
        new = pd.read_csv(new_path) if new_path.exists() else None
        for anc in ["EUR", "AMR", "AFR"]:
            line = f"  {pop} {anc}:"
            old_path = old_dir / f"{pop}_{anc}_IBD_ratios.csv"
            if old_path.exists():
                o = pd.read_csv(old_path)["ratio"]
                line += (f" old n={len(o):,} zero={100 * (o == 0).mean():.1f}% "
                         f"mean={o.mean():.4f} |")
            else:
                line += f" (no {old_path}) |"
            if new is not None:
                n = new[f"{anc}_ratio"].dropna()
                line += (f" new n={len(n):,} zero={100 * (n == 0).mean():.1f}% "
                         f"mean={n.mean():.4f}")
            print(line)

    print("\nTable 1 identical" if not problems else f"\nDIFFERENCES: {', '.join(problems)}")


if __name__ == "__main__":
    main()
