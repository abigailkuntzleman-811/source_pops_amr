#!/usr/bin/env python3
"""Archaic sequence recovered within a given continental ancestry background.

Intersects introgressed segments with local-ancestry tracts, so the archaic
sequence is partitioned by the continental ancestry of the haplotype carrying
it, then merges across individuals to give a cumulative recovery curve.

Requires the local ancestry inference BED files, which are NOT public.
Set LAI_DIR in config.py (or export it) to point at them. Expected layout:

    $LAI_DIR/1kg_lai_global/<POP>.bedlist
    $LAI_DIR/1kg_bed_files/<POP>_lai.bed
    $LAI_DIR/1kg_bed_files/chr<N>_<IND>_<HAP>.bed       continental ancestry
    $LAI_DIR/1kg_arc_bed_files/chr<N>_<IND>_<HAP>.bed   archaic ancestry

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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CHROMOSOMES, LAI_DIR, RECOVERED_DIR  # noqa: E402

RENAME = {"chrom": "Chromosome", "start": "Start", "end": "End"}


def population_individuals(amr_pop: str) -> list:
    """Read the individual IDs for a population from its bedlist."""
    bedlist = LAI_DIR / "1kg_lai_global" / f"{amr_pop}.bedlist"
    if not bedlist.exists():
        sys.exit(
            f"Missing input: {bedlist}\n"
            "These local ancestry files are not public; set LAI_DIR in config.py."
        )
    ids = set()
    n_lines = 0
    with open(bedlist) as f:
        for line in f:
            if not line.strip():
                continue
            n_lines += 1
            parts = line.strip().split("_")
            if len(parts) > 3:
                ids.add(parts[3])

    # The sample ID is the fourth underscore-separated field of the filename.
    # If that assumption breaks, many lines collapse onto a few IDs and the
    # analysis silently runs on the wrong sample set, so check rather than
    # trust it.
    if not ids:
        sys.exit(f"No sample IDs parsed from {bedlist} ({n_lines} lines read).")
    if n_lines and len(ids) < max(2, n_lines // 100):
        sys.exit(
            f"Only {len(ids)} sample ID(s) parsed from {n_lines} lines in "
            f"{bedlist}: {sorted(ids)[:5]}\n"
            "The filename format has probably changed; this script expects the "
            "sample ID in the fourth underscore-separated field."
        )

    return sorted(ids)


def recovered_in_ancestry(individuals, continental_ancestry: str) -> int:
    """Unique bp that is both archaic and of the given continental ancestry."""
    total = 0

    for chrom in CHROMOSOMES:
        overlaps = []
        for hap in (1, 2):
            for ind in individuals:
                arc_file = LAI_DIR / "1kg_arc_bed_files" / f"{chrom}_{ind}_{hap}.bed"
                anc_file = LAI_DIR / "1kg_bed_files" / f"{chrom}_{ind}_{hap}.bed"
                if not (arc_file.exists() and anc_file.exists()):
                    continue

                arc = pd.read_csv(arc_file, sep="\t")
                arc = arc[arc["ancestry"] == "introgressed"].rename(columns=RENAME)

                anc = pd.read_csv(anc_file, sep="\t")
                anc = anc[anc["ancestry"] == continental_ancestry].rename(columns=RENAME)

                if arc.empty or anc.empty:
                    continue

                overlaps.append(pr.PyRanges(arc).intersect(pr.PyRanges(anc)))

        if not overlaps:
            continue

        merged = pr.concat(overlaps).merge()
        total += int((merged.End - merged.Start).sum())

    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("amr_pop", help="admixed American population, e.g. MXL")
    parser.add_argument("continental_ancestry", help="ancestry label, e.g. EUR or AMR")
    args = parser.parse_args()

    RECOVERED_DIR.mkdir(parents=True, exist_ok=True)

    individuals = population_individuals(args.amr_pop)
    print(f"{args.amr_pop}: {len(individuals)} individuals")

    counts, recovered = [], []
    for i in range(1, len(individuals) + 1):
        counts.append(i)
        recovered.append(recovered_in_ancestry(individuals[:i], args.continental_ancestry))
        print(f"  {i}/{len(individuals)}: {recovered[-1]:,} bp", flush=True)

    out = (RECOVERED_DIR /
           f"{args.amr_pop}_only_{args.continental_ancestry}_archaic_recovered.txt")
    with open(out, "w") as f:
        f.write("x: " + ", ".join(map(str, counts)) + "\n")
        f.write("y: " + ", ".join(map(str, recovered)) + "\n")

    print(f"wrote {out}")


if __name__ == "__main__":
    main()
