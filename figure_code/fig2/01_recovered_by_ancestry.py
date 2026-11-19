#!/usr/bin/env python3
"""Archaic sequence recovered within a given continental ancestry background.

Intersects introgressed segments with local-ancestry tracts, so the archaic
sequence is partitioned by the continental ancestry of the haplotype carrying
it, then merges across individuals to give a cumulative recovery curve.

Inputs (joint-caller BED files, NOT public; set in config.local.sh):

    $ARC_BED_DIR/arc_chr<N>_<IND>_<HAP>.bed   ancestry = introgressed / no_introgressed
    $LAI_BED_DIR/lai_chr<N>_<IND>_<HAP>.bed   ancestry = continental label (EUR, AMR, ...)
    $POP_PANEL                                1000 Genomes SampleID -> Population

Usage:
    python fig2/01_recovered_by_ancestry.py <AMR_POP> <CONTINENTAL_ANCESTRY>

Example:
    python fig2/01_recovered_by_ancestry.py MXL EUR

As a Slurm array job covering all eight population/ancestry combinations:
    sbatch --array=0-7 fig2/01_recovered_by_ancestry.script
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import pyranges as pr

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import (ARC_BED_DIR, CHROMOSOMES, LAI_BED_DIR,  # noqa: E402
                    POP_PANEL, RECOVERED_DIR)

RENAME = {"chrom": "Chromosome", "start": "Start", "end": "End"}


def population_individuals(amr_pop: str) -> list:
    """Individual IDs for a population, from the 1000 Genomes sample panel."""
    if not POP_PANEL.exists():
        sys.exit(f"Missing input: {POP_PANEL}\nSet POP_PANEL in config.local.sh.")
    panel = pd.read_csv(POP_PANEL, sep=" ")
    ids = sorted(panel.loc[panel["Population"] == amr_pop, "SampleID"].unique())
    if not ids:
        sys.exit(f"No individuals with Population == {amr_pop} in {POP_PANEL}")
    return ids


def read_tracts(path: Path, ancestry: str) -> pd.DataFrame:
    """Rows of a per-haplotype BED file with the given ancestry label."""
    bed = pd.read_csv(path, sep="\t")
    return bed[bed["ancestry"] == ancestry].rename(columns=RENAME)


def archaic_in_ancestry(ind: str, continental_ancestry: str, labels_seen: set):
    """One individual's archaic sequence on tracts of the given ancestry.

    Returns a PyRanges over all autosomes and both haplotypes (unmerged), or
    None if the individual has no files. Haplotype suffixes are taken from the
    LAI filenames, and each must have a matching archaic file.
    """
    pieces = []
    n_pairs = 0
    for chrom in CHROMOSOMES:
        for lai_file in sorted(LAI_BED_DIR.glob(f"lai_{chrom}_{ind}_*.bed")):
            hap = lai_file.stem.rsplit("_", 1)[1]
            arc_file = ARC_BED_DIR / f"arc_{chrom}_{ind}_{hap}.bed"
            if not arc_file.exists():
                sys.exit(f"{lai_file.name} has no matching archaic file {arc_file}")
            n_pairs += 1

            anc = pd.read_csv(lai_file, sep="\t")
            labels_seen.update(anc["ancestry"].astype(str).unique())
            anc = anc[anc["ancestry"] == continental_ancestry].rename(columns=RENAME)
            arc = read_tracts(arc_file, "introgressed")
            if arc.empty or anc.empty:
                continue

            pieces.append(pr.PyRanges(arc[["Chromosome", "Start", "End"]])
                          .intersect(pr.PyRanges(anc[["Chromosome", "Start", "End"]])))

    if n_pairs == 0:
        return None
    pieces = [p for p in pieces if len(p) > 0]
    return pr.concat(pieces) if pieces else pr.PyRanges()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("amr_pop", help="admixed American population, e.g. MXL")
    parser.add_argument("continental_ancestry", help="ancestry label, e.g. EUR or AMR")
    args = parser.parse_args()

    RECOVERED_DIR.mkdir(parents=True, exist_ok=True)

    panel_ids = population_individuals(args.amr_pop)

    # Each individual's files are read once; the curve is then built by adding
    # one individual at a time to a running merged set.
    per_ind, labels_seen, no_files = [], set(), []
    for ind in panel_ids:
        ranges = archaic_in_ancestry(ind, args.continental_ancestry, labels_seen)
        if ranges is None:
            no_files.append(ind)
        else:
            per_ind.append((ind, ranges))

    print(f"{args.amr_pop}: {len(per_ind)} of {len(panel_ids)} panel individuals have files")
    if no_files:
        print(f"  no files for: {', '.join(no_files[:10])}"
              f"{' ...' if len(no_files) > 10 else ''}")
    if not per_ind:
        sys.exit(f"No lai_chr*_<IND>_*.bed files for {args.amr_pop} in {LAI_BED_DIR}")
    if args.continental_ancestry not in labels_seen:
        sys.exit(f"Ancestry label {args.continental_ancestry!r} never appears in the "
                 f"LAI files; labels found: {sorted(labels_seen)}")

    counts, recovered = [], []
    running = pr.PyRanges()
    for i, (ind, ranges) in enumerate(per_ind, start=1):
        if len(ranges) > 0:
            running = pr.concat([running, ranges]).merge() if len(running) > 0 else ranges.merge()
        counts.append(i)
        recovered.append(int(running.lengths().sum()) if len(running) > 0 else 0)
        print(f"  {i}/{len(per_ind)} {ind}: {recovered[-1]:,} bp", flush=True)

    out = (RECOVERED_DIR /
           f"{args.amr_pop}_only_{args.continental_ancestry}_archaic_recovered.txt")
    with open(out, "w") as f:
        f.write("x: " + ", ".join(map(str, counts)) + "\n")
        f.write("y: " + ", ".join(map(str, recovered)) + "\n")

    print(f"wrote {out}")


if __name__ == "__main__":
    main()
