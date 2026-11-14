#!/usr/bin/env python3
"""Split filtered hmmix segments into per-individual, per-haplotype BED files.

Input is a filtered CSV from 01_filter_hmmix_calls.py. Output is one sorted BED
file per haplotype, grouped into a directory per population:

    <out_root>/<POP>/sorted-<IND>-<HAP>.bed

Files are written in sorted order because `bedops --intersect`, used by
03_pairwise_sharing.sh, requires sorted input. The "sorted-" prefix is also
what the downstream averaging step keys on.

This replaces the original organize_hmmix.py plus the sort-bed loop in
organize_hmmix.script. That version wrote an unsorted .bed and a sorted- copy
side by side, so the pairwise step intersected every pair four times and then
discarded the three involving an unsorted file. Writing only sorted files gives
identical averages for a quarter of the compute.

Handles both column conventions seen in this project: the uppercase schema of
tgp_hmmix_haplotype_tracts.csv (IND/CHR/START/END/POP/TRACT_ID) and the
lowercase schema of Skov's hg38_1000g_segments.txt (name/chrom/start/end/pop).

Usage:
    python fig1/02_segments_to_bed.py <POP> <filtered_csv> [out_root]

Example:
    python fig1/02_segments_to_bed.py MXL \
        data/processed/hmmix_filtered/hmmix_neanderthal_filtered.csv \
        data/processed/hmmix_by_pop/neanderthal
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import BY_POP_DIR  # noqa: E402


def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a frame with chrom / start / end / name / haplotype / pop."""
    cols = set(df.columns)

    if {"CHR", "START", "END"} <= cols:
        out = pd.DataFrame({
            "chrom": "chr" + df["CHR"].astype(str).str.removeprefix("chr"),
            "start": df["START"],
            "end": df["END"],
            "pop": df["POP"],
        })
        if "TRACT_ID" in cols:
            # TRACT_ID already looks like "NA19648-hap1"
            out["hap_id"] = df["TRACT_ID"]
        else:
            out["hap_id"] = df["IND"].astype(str) + "-hap" + df["HAP_ID"].astype(str)

    elif {"chrom", "start", "end"} <= cols:
        out = pd.DataFrame({
            "chrom": df["chrom"].astype(str),
            "start": df["start"],
            "end": df["end"],
            "pop": df["pop"],
        })
        out["hap_id"] = df["name"].astype(str) + "-" + df["haplotype"].astype(str)

    else:
        raise ValueError(
            "Could not find coordinate columns. Expected CHR/START/END or "
            f"chrom/start/end; got: {sorted(cols)}"
        )

    return out


def sort_key(df: pd.DataFrame) -> pd.DataFrame:
    """Sort as `sort-bed` does: lexicographic chromosome, then numeric start."""
    return df.sort_values(["chrom", "start", "end"], kind="mergesort")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pop", help="population code, e.g. MXL")
    parser.add_argument("segments", type=Path, help="filtered hmmix CSV")
    parser.add_argument(
        "out_root",
        type=Path,
        nargs="?",
        default=BY_POP_DIR,
        help=f"output root (default: {BY_POP_DIR})",
    )
    args = parser.parse_args()

    if not args.segments.exists():
        sys.exit(f"Missing input: {args.segments}\nRun 01_filter_hmmix_calls.py first.")

    df = normalise_columns(pd.read_csv(args.segments))

    pop_data = df[df["pop"] == args.pop]
    if pop_data.empty:
        print(f"No segments for population '{args.pop}'. Nothing to do.")
        return

    pop_dir = args.out_root / args.pop
    pop_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    for hap_id, group in pop_data.groupby("hap_id"):
        # Drop zero-length intervals; bedops rejects them.
        group = group[group["start"] != group["end"]]
        if group.empty:
            continue

        bed = sort_key(group)[["chrom", "start", "end"]]
        bed.to_csv(pop_dir / f"sorted-{hap_id}.bed", sep="\t",
                   index=False, header=False)
        written += 1

    print(f"{args.pop}: wrote {written} haplotype BED files to {pop_dir}")


if __name__ == "__main__":
    main()
