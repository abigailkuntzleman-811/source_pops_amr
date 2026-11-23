#!/usr/bin/env python3
"""UpSet plot: archaic sequence shared and unique across the admixed populations.

Same as the original scripts/shared_unique_between_admix_pops.py. Uses the
downsampled sets from figure_code/figS2_hmmix_recovered/01_recovered_from_hmmix.py
(every population downsampled to the smallest sample size). The denominator is
all sequence present in at least one population; the 15 regions (unique to one
population, shared by each combination of 2, 3, or all 4) partition it.

Usage:
    python figure_code/upset/01_upset_shared_unique.py all
    python figure_code/upset/01_upset_shared_unique.py neanderthal

Output:
    $TABLE_DIR/upset_regions_<type>.csv   bp and percent per region
    $FIGURE_DIR/upset_<type>.png
"""

import argparse
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
from upsetplot import UpSet

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import FIGURE_DIR, RECOVERED_DIR, TABLE_DIR  # noqa: E402

POPS = ["CLM", "PUR", "PEL", "MXL"]


def load(pop, archaic):
    paths = sorted(RECOVERED_DIR.glob(f"{pop}_n*_recovered_{archaic}.bed"))
    if len(paths) != 1:
        sys.exit(f"Expected one {pop}_n*_recovered_{archaic}.bed in {RECOVERED_DIR}, found {len(paths)}\n"
                 f"Run: python figure_code/figS2_hmmix_recovered/01_recovered_from_hmmix.py {archaic}")
    return pr.read_bed(str(paths[0])).merge()


def bp(gr):
    return int(gr.lengths().sum()) if len(gr) > 0 else 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("archaic", choices=["neanderthal", "denisovan", "all"])
    args = p.parse_args()

    beds = {pop: load(pop, args.archaic) for pop in POPS}
    total = bp(pr.concat(list(beds.values())).merge())
    regions = {}
    for pattern in itertools.product([True, False], repeat=len(POPS)):
        inside = [pop for pop, k in zip(POPS, pattern) if k]
        outside = [pop for pop, k in zip(POPS, pattern) if not k]
        if not inside:
            continue
        region = beds[inside[0]]
        for pop in inside[1:]:
            region = region.intersect(beds[pop])
        for pop in outside:
            region = region.subtract(beds[pop])
        regions[pattern] = bp(region)
    assert sum(regions.values()) == total, (sum(regions.values()), total)

    sizes = pd.Series(regions)
    sizes.index.names = POPS
    table = sizes.rename("bp").reset_index()
    table["percent"] = 100 * table["bp"] / total
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE_DIR / f"upset_regions_{args.archaic}.csv", index=False)
    print(f"{args.archaic}: union = {total:,} bp")
    print(table.to_string(index=False))

    UpSet(sizes, sort_by="cardinality", show_percentages=True).plot()
    plt.suptitle(f"{args.archaic.capitalize()} sequence shared across populations")
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / f"upset_{args.archaic}.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
