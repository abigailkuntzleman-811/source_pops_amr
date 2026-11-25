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
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
from matplotlib.ticker import FuncFormatter
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


def sci_label(x, _pos):
    """50000000 -> 5 x 10^7 (same style as the other figures)."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = round(x / 10 ** e, 2)
    if m == 10:
        m, e = 1, e + 1
    return rf"${m:g} \times 10^{{{e}}}$"


# "all" and "neanderthal" are drawn on the same axis ranges so they can be
# compared by eye; "denisovan" keeps its own (its bars are ~25x smaller).
SHARED_SCALE = ("all", "neanderthal")


def axis_maxima(table):
    """Largest intersection and largest population total (bp) in a regions table."""
    totals = [table.loc[table[pop], "bp"].sum() for pop in POPS]
    return table["bp"].max(), max(totals)


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

    axes = UpSet(sizes, sort_by="cardinality", show_percentages=True).plot()
    # bp axes as 5 x 10^7 rather than an offset "1e8"
    axes["intersections"].yaxis.set_major_formatter(FuncFormatter(sci_label))
    axes["totals"].xaxis.set_major_formatter(FuncFormatter(sci_label))
    axes["totals"].tick_params(axis="x", labelrotation=45)
    # percentage labels above the bars: vertical so neighbouring labels don't overlap
    ax_i = axes["intersections"]
    ax_t = axes["totals"]
    if args.archaic in SHARED_SCALE:
        maxima = [axis_maxima(table)]
        for other in SHARED_SCALE:
            path = TABLE_DIR / f"upset_regions_{other}.csv"
            if other != args.archaic and path.exists():
                maxima.append(axis_maxima(pd.read_csv(path)))
            elif other != args.archaic:
                print(f"note: {path.name} not found yet; run {other} too, then rerun "
                      f"{args.archaic} so both share the same axes")
        top_bar = max(m[0] for m in maxima)
        top_total = max(m[1] for m in maxima)
        ax_i.set_ylim(0, top_bar * 1.07)
        lo, hi = ax_t.get_xlim()
        ax_t.set_xlim((top_total * 1.05, 0) if lo > hi else (0, top_total * 1.05))
    ax_i.set_ylim(0, ax_i.get_ylim()[1] * 1.15)
    gap = 0.01 * ax_i.get_ylim()[1]
    for t in ax_i.texts:
        x, y = t.get_position()
        t.set_position((x, y + gap))
        t.set(rotation=90, rotation_mode="anchor", ha="left", va="center", fontsize=9)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / f"upset_{args.archaic}.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
