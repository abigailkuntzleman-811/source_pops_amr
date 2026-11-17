#!/usr/bin/env python3
"""Average within-population IBD per individual.

For each individual, sums the cM they share with every other member of their
population across the autosomes, then divides by (n - 1), the number of
possible within-population partners. This is the IBD_bar_i of Figure 4.

Note this is a different quantity from fig3/01_within_population_ibd.py, which
averages over haplotype *pairs*. Here the unit is the individual, matching the
ancestry proportions from table1/00_ancestry_proportions.py so the two can be
joined.

Related individuals are excluded, and are also left out of the (n - 1)
denominator. Set RELATED_EXCLUDE to a file of sample IDs; otherwise the
samples with a parent recorded in the 1000 Genomes pedigree are used.

Prerequisite: hapIBD_calls/01_call_ibd.sh for all 22 autosomes.

Usage:
    python figure_code/fig4/01_individual_ibd.py
    python figure_code/fig4/01_individual_ibd.py --pops MXL PEL PUR CLM

Output:
    data/processed/ibd/individual_ibd.csv    one row per individual
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))

from config import AMR_POPS, IBD_DIR  # noqa: E402

# Reuse the metadata helpers so both figures agree on populations and
# which individuals count as related.
sys.path.insert(0, str(_root / "figure_code" / "fig3"))
_fig3 = __import__("01_within_population_ibd")
load_population_map = _fig3.load_population_map
load_related_samples = _fig3.load_related_samples
read_ibd_segments = _fig3.read_ibd_segments


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pops", nargs="+", default=AMR_POPS)
    args = parser.parse_args()

    pop_of = load_population_map()
    related = load_related_samples()

    print("reading IBD segments")
    ibd = read_ibd_segments()
    print(f"  {len(ibd):,} segments total")

    ibd = ibd[~ibd["id1"].isin(related) & ~ibd["id2"].isin(related)]
    ibd["pop1"] = ibd["id1"].map(pop_of)
    ibd["pop2"] = ibd["id2"].map(pop_of)
    ibd = ibd[(ibd["pop1"] == ibd["pop2"]) & ibd["pop1"].isin(args.pops)]
    ibd = ibd[ibd["id1"] != ibd["id2"]]
    print(f"  {len(ibd):,} within-population segments after exclusions")

    rows = []
    for pop in args.pops:
        samples = sorted(set(pop_of[pop_of == pop].index) - related)
        if len(samples) < 2:
            print(f"{pop}: fewer than two unrelated samples, skipping",
                  file=sys.stderr)
            continue

        sub = ibd[ibd["pop1"] == pop]

        # An individual appears as id1 in some segments and id2 in others;
        # their total is the sum of both.
        as_id1 = sub.groupby("id1")["cM"].sum()
        as_id2 = sub.groupby("id2")["cM"].sum()
        totals = as_id1.add(as_id2, fill_value=0)

        # Average over the n - 1 possible partners, counting individuals with
        # no detected IBD as zero.
        avg = totals.reindex(samples).fillna(0) / (len(samples) - 1)

        rows.append(pd.DataFrame({
            "individual": avg.index,
            "pop": pop,
            "avg_inpop_ibd": avg.to_numpy(),
        }))

        print(f"{pop}: {len(samples)} individuals, "
              f"mean {avg.mean():.4f} cM, max {avg.max():.4f} cM")

    if not rows:
        sys.exit("No populations produced results.")

    IBD_DIR.mkdir(parents=True, exist_ok=True)
    out = IBD_DIR / "individual_ibd.csv"
    pd.concat(rows, ignore_index=True).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
