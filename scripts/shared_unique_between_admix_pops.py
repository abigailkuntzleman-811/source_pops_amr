## using pyranges, get all subset intersections between neanderthal, denisovan, and all
# from /users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop
# (named, i.e., MXL_n97_recovered_all.bed)
# I want to get those segments and count them, then do a ven diagram with % of total sequence present across
#everyone

#so denominator will be union of all segments across all pops, then numberator is each of the sections
#i.e. unique to each of the four, shared between all combos of 2 for each of the four, then 3, then all


#output to /users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/pop_intersections

import glob
import itertools
import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
from upsetplot import UpSet

DIR = '/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop'
OUT_DIR = f'{DIR}/pop_intersections'
os.makedirs(OUT_DIR, exist_ok=True)

POPS = ['CLM', 'PUR', 'PEL', 'MXL']
parser = argparse.ArgumentParser()
parser.add_argument('archaic', choices=['neanderthal', 'denisovan', 'all'])
args = parser.parse_args()
ARCHAIC_TYPES = [args.archaic]


def load(pop, archaic):
    # exactly one file per pop/type, e.g. MXL_n97_recovered_all.bed
    path, = glob.glob(f'{DIR}/{pop}_n*_recovered_{archaic}.bed')
    return pr.read_bed(path).merge()


def bp(gr):
    return int(gr.lengths().sum()) if len(gr) > 0 else 0


for archaic in ARCHAIC_TYPES:
    beds = {pop: load(pop, archaic) for pop in POPS}

    # denominator: all sequence present in at least one population
    total = bp(pr.concat(list(beds.values())).merge())

    # numerators: the 15 regions, e.g. (True, False, True, False) = in CLM and PEL only
    regions = {}
    for pattern in itertools.product([True, False], repeat=len(POPS)):
        inside = [pop for pop, k in zip(POPS, pattern) if k]
        outside = [pop for pop, k in zip(POPS, pattern) if not k]
        if not inside:
            continue

        region = beds[inside[0]]
        for pop in inside[1:]:
            region = region.intersect(beds[pop])   # keep what's shared
        for pop in outside:
            region = region.subtract(beds[pop])    # drop what any other pop has

        regions[pattern] = bp(region)

    # sanity check: the 15 regions don't overlap and together make up the union
    assert sum(regions.values()) == total, (sum(regions.values()), total)

    # one row per region, indexed by True/False membership in each pop
    sizes = pd.Series(regions)
    sizes.index.names = POPS

    # save the numbers
    table = sizes.rename('bp').reset_index()
    table['percent'] = 100 * table['bp'] / total
    table.to_csv(f'{OUT_DIR}/regions_{archaic}.csv', index=False)
    print(f'\n{archaic}: union = {total:,} bp')
    print(table.to_string(index=False))

    # UpSet plot; percentages are of the union (the sum of all 15 regions)
    UpSet(sizes, sort_by='cardinality', show_percentages=True).plot()
    plt.suptitle(f'{archaic.capitalize()} sequence shared across populations')
    plt.savefig(f'{OUT_DIR}/upset_{archaic}.png', dpi=300, bbox_inches='tight')
    plt.close()