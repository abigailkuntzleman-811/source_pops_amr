#!/usr/bin/env python3
"""Check Figure 2 against the original scripts/get_pop_nean_recovered.py.

Runs the original get_pop_archaic_recovered, copied unchanged except for the
two file paths (Jaz's old chr<N>_<IND>_<HAP>.bed folders -> the current
ARC_BED_DIR/arc_chr<N>_<IND>_<HAP>.bed and LAI_BED_DIR/lai_chr<N>_<IND>_<HAP>.bed),
with individuals in <POP>_lai.bed order as in the original. Compares it with
the curve written by figure_code/fig2/01_recovered_by_ancestry.py at a few
points along the curve (the original recomputes every point from scratch, so
the whole curve would take days).

    sbatch -J check_fig2 run.script checks/check_fig2.py
    sbatch -J check_fig2 run.script checks/check_fig2.py MXL EUR    # one panel
"""

import sys
from pathlib import Path

import pandas as pd
import pyranges as pr

_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_root))
from config import ARC_BED_DIR, LAI_BED_DIR, RECOVERED_DIR  # noqa: E402

ARC = str(ARC_BED_DIR / "arc_chr{i}_{ind}_{hap}.bed")
LAI = str(LAI_BED_DIR / "lai_chr{i}_{ind}_{hap}.bed")
# haplotype 1 and 2 as they appear in the file names (set in main: "1"/"2" or "hap1"/"hap2")
SUFFIX = {1: "1", 2: "2"}


# ---------------------------------------- from scripts/get_pop_nean_recovered.py
def get_pop_archaic_recovered(ind_ids, lai_pop):
    recovered_seq = 0
    #for each chromosomes/haplotype, I want to merge the sequence over all the individuals
    for i in range(1, 23):
        overlaps = []
        for hap in [1,2]:
            for ind in ind_ids:
                chr_bed = pd.read_csv(
                    ARC.format(i=i, ind=ind, hap=SUFFIX[hap]),
                    sep='\t'
                )
                # only introgressed regions
                intro_bed = chr_bed[chr_bed['ancestry'] == 'introgressed']


                eur_bed = pd.read_csv(LAI.format(i=i, ind=ind, hap=SUFFIX[hap]),
                                      sep = '\t'
                                     )
                #continental regions
                eur_bed = eur_bed[eur_bed['ancestry'] == lai_pop]

                intro_bed = intro_bed.rename(columns={
                    "chrom": "Chromosome",
                    "start": "Start",
                    "end": "End"
                })

                eur_bed = eur_bed.rename(columns={
                    "chrom": "Chromosome",
                    "start": "Start",
                    "end": "End"
                })
                intro_pr = pr.PyRanges(intro_bed)
                eur_pr   = pr.PyRanges(eur_bed)

                overlaps.append(intro_pr.intersect(eur_pr))
        overlap = pr.concat(overlaps)
        #now I merge all the sequence that is archaic and european in the population
        overlap = overlap.merge()
        if len(overlap) > 0:  # only addition: an empty chromosome adds 0 instead of crashing
            recovered_seq += (overlap.End - overlap.Start).sum()

    return recovered_seq
# --------------------------------------------------------------------------------


def repo_curve(pop, anc):
    path = RECOVERED_DIR / f"{pop}_only_{anc}_archaic_recovered.txt"
    lines = dict(line.split(": ", 1) for line in path.read_text().splitlines())
    x = [int(v) for v in lines["x"].split(", ")]
    y = [int(v) for v in lines["y"].split(", ")]
    return dict(zip(x, y))


def main():
    combos = ([tuple(sys.argv[1:3])] if len(sys.argv) >= 3 else
              [(p, a) for p in ("CLM", "MXL", "PEL", "PUR") for a in ("EUR", "AMR")])
    problems = []
    first = pd.read_csv(LAI_BED_DIR / f"{combos[0][0]}_lai.bed", sep='\t')['name'].iloc[0]
    for cand in (("1", "2"), ("hap1", "hap2"), ("0", "1")):
        if Path(LAI.format(i=1, ind=first, hap=cand[0])).exists():
            SUFFIX.update({1: cand[0], 2: cand[1]})
            break
    else:
        sys.exit(f"No {LAI.format(i=1, ind=first, hap='*')} file found")
    print(f"haplotype file suffixes: {SUFFIX[1]}, {SUFFIX[2]}")
    for pop, anc in combos:
        df = pd.read_csv(LAI_BED_DIR / f"{pop}_lai.bed", sep='\t')
        ids = df['name'].unique()
        repo = repo_curve(pop, anc)
        n = len(ids)
        print(f"{pop} {anc}: {n} individuals in {pop}_lai.bed, {len(repo)} points in the repo curve", flush=True)
        if n != len(repo):
            problems.append(f"{pop} {anc} length")
        for k in sorted({k for k in (1, 2, 5, 10, n // 2, n) if 1 <= k <= n}):
            original = int(get_pop_archaic_recovered(ids[:k], anc))
            ours = repo.get(k)
            ok = original == ours
            if not ok:
                problems.append(f"{pop} {anc} x={k}")
            print(f"  {'OK  ' if ok else 'DIFF'} x={k}: original {original:,} vs repo {ours if ours is None else f'{ours:,}'} bp",
                  flush=True)
    print("\nALL IDENTICAL" if not problems else f"\n{len(problems)} DIFFERENCES: {', '.join(problems)}")


if __name__ == "__main__":
    main()
