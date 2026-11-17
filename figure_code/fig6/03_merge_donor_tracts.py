#!/usr/bin/env python3
"""Merge per-donor ancestry tracts into one file per replicate.

Scenarios with more than one European donor population produce one stitched
tract file per donor. Figure 6 treats all of them as European ancestry, so
they are combined here: tracts are pooled per haplotype and any overlapping or
abutting intervals merged, so shared material is not counted twice.

NOTE: this reproduces a step whose original script was not available. Before
relying on it, compare its output against the existing merged files; see the
--compare option.

Usage:
    python figure_code/fig6/03_merge_donor_tracts.py <scenario> <rep> \
        --sources 8 9 --out-name "{pop}_merged_89tracts_rep{rep}.csv.gz"

    python figure_code/fig6/03_merge_donor_tracts.py mxl_donor6 1 \
        --sources 10 11 12 13 14 15 \
        --out-name "{pop}_totalIBSmerged_rep{rep}.csv.gz"

    # check against a file produced by the original pipeline
    python figure_code/fig6/03_merge_donor_tracts.py mxl_donor6 1 \
        --sources 10 11 12 13 14 15 --compare

Output:
    $SIM_DIR/<scenario>/admix_csvs/<out-name>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import SIM_DIR  # noqa: E402
from scenarios import SCENARIOS  # noqa: E402


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping or abutting tracts within each haplotype."""
    rows = []
    for hap_id, group in df.groupby("hap_id"):
        group = group.sort_values("start")
        start = end = None

        for row in group.itertuples(index=False):
            if start is None:
                start, end = row.start, row.end
            elif row.start <= end:
                # overlapping or abutting: extend the open tract
                end = max(end, row.end)
            else:
                rows.append((hap_id, start, end, end - start))
                start, end = row.start, row.end

        if start is not None:
            rows.append((hap_id, start, end, end - start))

    # Column names match the original merge step's output, so files produced
    # either way are interchangeable.
    out = pd.DataFrame(rows, columns=["hap_id", "Start", "End", "length"])
    out.insert(0, "Chromosome", "chr1")
    return out[["Chromosome", "Start", "End", "hap_id", "length"]]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=sorted(SCENARIOS))
    parser.add_argument("rep", type=int)
    parser.add_argument("--pop", default="MXL")
    parser.add_argument("--sources", nargs="+", type=int, required=True,
                        help="population ids of the donor populations")
    parser.add_argument("--out-name", default=None,
                        help="output filename pattern; defaults to the "
                             "scenario's configured tracts filename")
    parser.add_argument("--compare", action="store_true",
                        help="compare against an existing file of that name "
                             "instead of writing")
    args = parser.parse_args()

    admix_dir = SIM_DIR / args.scenario / "admix_csvs"
    pattern = args.out_name or SCENARIOS[args.scenario]["tracts"]
    out_path = admix_dir / pattern.format(pop=args.pop, rep=args.rep)

    frames = []
    for source_id in args.sources:
        path = admix_dir / f"{args.pop}_{source_id}tracts_rep{args.rep}_stitched.csv.gz"
        if not path.exists():
            sys.exit(
                f"Missing {path}\n"
                f"Run: python figure_code/fig6/02_extract_tracts.py "
                f"{args.scenario} {args.rep} --source <name of population "
                f"{source_id}>"
            )
        df = pd.read_csv(path)
        df["hap_id"] = df["hap_id"].astype(int)
        frames.append(df)
        print(f"  source {source_id}: {len(df):,} tracts")

    merged = merge_intervals(pd.concat(frames, ignore_index=True))
    total_bp = merged["length"].sum()
    print(f"merged: {len(merged):,} tracts, {total_bp:,.0f} bp")

    if args.compare:
        if not out_path.exists():
            sys.exit(f"Nothing to compare against: {out_path} does not exist")
        existing = pd.read_csv(out_path)
        existing["hap_id"] = existing["hap_id"].astype(int)
        print(f"existing: {len(existing):,} tracts, "
              f"{existing['length'].sum():,.0f} bp")

        same_n = len(merged) == len(existing)
        same_bp = abs(total_bp - existing["length"].sum()) < 1
        print(f"\nsame tract count: {same_n}")
        print(f"same total bp:    {same_bp}")
        if not (same_n and same_bp):
            print("\nThe reconstruction differs from the original pipeline; "
                  "do not use it without working out why.")
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False, compression="gzip")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
