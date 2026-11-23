#!/usr/bin/env python3
"""Archaic sequence shared vs unique between ancestry backgrounds within a population.

Reviewer 2: "within a single population, are the same archaic haplotypes
observed across different ancestry backgrounds, or do the archaic haplotypes
tend to share the same source ancestry?"

Archaic (introgressed) sequence is split by the local ancestry of the
haplotype carrying it, as in Figure 2.

Steps, for one admixed population:
  1. Usable genome: positions where at least --k (default 20) haplotypes have
     European ancestry and at least --k have Indigenous American ancestry.
     Everything else is ignored.
  2. At each usable position, pick --k European and --k Indigenous haplotypes
     at random. The position is
        shared    if at least one picked European haplotype and at least one
                  picked Indigenous haplotype carry archaic sequence there
        EUR only  if only European ones do
        AMR only  if only Indigenous ones do
     Picking the same number from each background means a background with
     more haplotypes does not find more archaic sequence just by being
     sampled more.
  3. Of all archaic sequence (shared + EUR only + AMR only, in bp), the
     percentage in each class. Repeated --reps times (default 100) and
     averaged.
  4. Null: the same, but the EUR/AMR labels are first shuffled among the
     haplotypes at each position, i.e. what sharing would look like if
     archaic sequence did not depend on the ancestry of the tract it sits on.
     Repeated --shuffles times (default 1000).
  5. p-value: the fraction of shuffles with a shared percentage at or below
     the observed average. A small p means archaic sequence is shared between
     backgrounds less than expected: it tends to stay on one background.

"Shared" means archaic sequence at the same position on both backgrounds;
it does not check that the archaic haplotypes are identical.

Inputs (as in Figure 2; set in config.local.sh):
    $ARC_BED_DIR/arc_chr<N>_<IND>_<HAP>.bed   ancestry = introgressed / no_introgressed
    $LAI_BED_DIR/lai_chr<N>_<IND>_<HAP>.bed   ancestry = continental label (EUR, AMR, ...)
    $POP_PANEL                                1000 Genomes SampleID -> Population
All individuals in the population with LAI files are used, as in Figure 2.

Usage:
    python figure_code/archaic_by_background/01_shared_unique_by_background.py PEL
As a Slurm array over the four populations:
    sbatch --array=0-3 figure_code/archaic_by_background/01_shared_unique_by_background.script

Output ($TABLE_DIR):
    archaic_by_background_<POP>.csv   one row: usable genome, observed and
                                      null percentages, p-value
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import ARC_BED_DIR, CHROMOSOMES, LAI_BED_DIR, POP_PANEL, TABLE_DIR  # noqa: E402

BACKGROUNDS = ("EUR", "AMR")


def population_individuals(pop):
    if not POP_PANEL.exists():
        sys.exit(f"Missing input: {POP_PANEL}\nSet POP_PANEL in config.local.sh.")
    panel = pd.read_csv(POP_PANEL, sep=" ")
    ids = sorted(panel.loc[panel["Population"] == pop, "SampleID"].unique())
    if not ids:
        sys.exit(f"No individuals with Population == {pop} in {POP_PANEL}")
    return ids


def read_chromosome(chrom, inds):
    """(LAI tracts, archaic segments) for every haplotype with an LAI file on this chromosome."""
    haps = []
    for ind in inds:
        for lai_file in sorted(LAI_BED_DIR.glob(f"lai_{chrom}_{ind}_*.bed")):
            hap = lai_file.stem.rsplit("_", 1)[1]
            arc_file = ARC_BED_DIR / f"arc_{chrom}_{ind}_{hap}.bed"
            if not arc_file.exists():
                sys.exit(f"{lai_file.name} has no matching archaic file {arc_file}")
            lai = pd.read_csv(lai_file, sep="\t")
            arc = pd.read_csv(arc_file, sep="\t")
            haps.append((lai[["start", "end", "ancestry"]],
                         arc.loc[arc["ancestry"] == "introgressed", ["start", "end"]]))
    return haps


def counts_per_interval(haps):
    """Split the chromosome at every breakpoint of every haplotype; per interval return
    width (bp), and the number of EUR haplotypes, archaic EUR haplotypes, AMR
    haplotypes and archaic AMR haplotypes."""
    edges = np.unique(np.concatenate(
        [np.concatenate([l["start"], l["end"], a["start"], a["end"]]) for l, a in haps]).astype(np.int64))
    G = len(edges) - 1
    nE, aE, nA, aA = (np.zeros(G, np.int32) for _ in range(4))
    for l, a in haps:
        is_arc = np.zeros(G + 1, np.int32)
        np.add.at(is_arc, np.searchsorted(edges, a["start"]), 1)
        np.add.at(is_arc, np.searchsorted(edges, a["end"]), -1)
        is_arc = np.cumsum(is_arc)[:-1] > 0
        for bg, n, x in (("EUR", nE, aE), ("AMR", nA, aA)):
            on = np.zeros(G + 1, np.int32)
            t = l[l["ancestry"] == bg]
            np.add.at(on, np.searchsorted(edges, t["start"]), 1)
            np.add.at(on, np.searchsorted(edges, t["end"]), -1)
            on = np.cumsum(on)[:-1] > 0
            n += on
            x += on & is_arc
    return np.diff(edges), nE, aE, nA, aA


def percentages(width, eur_has, amr_has):
    """Percent of archaic bp that is shared, EUR only and AMR only."""
    shared = width[eur_has & amr_has].sum()
    eur = width[eur_has & ~amr_has].sum()
    amr = width[~eur_has & amr_has].sum()
    total = shared + eur + amr
    return 100 * np.array([shared, eur, amr]) / total if total else np.full(3, np.nan)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pop", help="admixed population, e.g. PEL")
    p.add_argument("--k", type=int, default=20, help="haplotypes picked per background (default 20)")
    p.add_argument("--reps", type=int, default=100, help="random picks to average (default 100)")
    p.add_argument("--shuffles", type=int, default=1000, help="shuffled-label repeats (default 1000)")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    k = args.k

    inds = population_individuals(args.pop)
    parts, called_bp, n_haps = [], 0, 0
    for chrom in CHROMOSOMES:
        haps = read_chromosome(chrom, inds)
        if not haps:
            continue
        n_haps = max(n_haps, len(haps))
        width, nE, aE, nA, aA = counts_per_interval(haps)
        called_bp += width[(nE > 0) & (nA > 0)].sum()
        use = (nE >= k) & (nA >= k)
        parts.append(np.column_stack([width, nE, aE, nA, aA])[use])
        print(f"{args.pop} {chrom}: {len(haps)} haplotypes, {width[use].sum() / 1e6:.1f} Mb usable", flush=True)
    if not parts:
        sys.exit(f"No lai_chr*_<IND>_*.bed files for {args.pop} in {LAI_BED_DIR}")
    width, nE, aE, nA, aA = np.concatenate(parts).T.astype(np.int64)

    # Picking k haplotypes at random from n, of which a are archaic: the number
    # of archaic ones picked is a hypergeometric draw.
    def pick(n, a):
        return rng.hypergeometric(a, n - a, k) > 0

    observed = np.array([percentages(width, pick(nE, aE), pick(nA, aA)) for _ in range(args.reps)])

    # Shuffled labels: pool the EUR and AMR haplotypes, deal k to "EUR", then k
    # of the rest to "AMR".
    N, a = nE + nA, aE + aA
    null = []
    for _ in range(args.shuffles):
        xE = rng.hypergeometric(a, N - a, k)
        xA = rng.hypergeometric(a - xE, (N - k) - (a - xE), k)
        null.append(percentages(width, xE > 0, xA > 0))
    null = np.array(null)

    obs = np.nanmean(observed, 0)
    row = {
        "pop": args.pop, "k": k, "haplotypes": n_haps,
        "usable_Mb": width.sum() / 1e6,
        "usable_pct_of_genome_with_both": 100 * width.sum() / called_bp if called_bp else np.nan,
        "archaic_Mb_mean": np.mean([width[pick(nE, aE) | pick(nA, aA)].sum() for _ in range(10)]) / 1e6,
        "shared_pct": obs[0], "eur_only_pct": obs[1], "amr_only_pct": obs[2],
        "null_shared_pct": np.nanmean(null[:, 0]),
        "null_eur_only_pct": np.nanmean(null[:, 1]), "null_amr_only_pct": np.nanmean(null[:, 2]),
        "p_shared_not_above_null": float((null[:, 0] <= obs[0]).mean()),
        "reps": args.reps, "shuffles": args.shuffles,
    }
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / f"archaic_by_background_{args.pop}.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    print()
    for key, v in row.items():
        print(f"  {key:32s} {v:.3f}" if isinstance(v, float) else f"  {key:32s} {v}")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
