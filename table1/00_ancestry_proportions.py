#!/usr/bin/env python3
"""Continental ancestry proportions from local ancestry inference calls.

For each admixed American population, sums local ancestry tract lengths per
ancestry and divides by the total, giving the proportion of each individual's
(or haplotype's) genome assigned to each continental ancestry.

Two units are available, because they give different standard deviations:

    individual   group by sample, summing both haplotypes. One row per person.
    haplotype    group by sample and haplotype. One row per haplotype, twice
                 as many rows, and a wider spread.

Means are nearly identical either way; only the spread differs. Use whichever
matches the published table, and say which in the methods.

INPUT: joint-caller local ancestry BED files. These are not publicly
available. Set LAI_BED_DIR in config.local.sh to wherever you have them.
Expected, one file per population:

    $LAI_BED_DIR/<POP>_lai.bed

with columns chrom, start, end, ancestry, name, haplotype (pop optional).
Every individual in the file is used, related individuals included.

Usage:
    python table1/00_ancestry_proportions.py
    python table1/00_ancestry_proportions.py --unit haplotype

Output:
    data/processed/ancestry_proportions/<POP>_<unit>.csv   per-unit proportions
    data/processed/ancestry_proportions/summary_<unit>.csv mean and SD per pop
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import AMR_POPS, DATA_PROCESSED, LAI_BED_DIR  # noqa: E402

OUT_DIR = DATA_PROCESSED / "ancestry_proportions"


def lai_path(pop: str) -> Path:
    return LAI_BED_DIR / f"{pop}_lai.bed"


def ancestry_proportions(pop: str, unit: str) -> pd.DataFrame:
    """Proportion of each continental ancestry, per individual or haplotype."""
    path = lai_path(pop)
    if not path.exists():
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t")

    missing = {"name", "haplotype", "start", "end", "ancestry"} - set(df.columns)
    if missing:
        sys.exit(f"{path} is missing expected columns: {sorted(missing)}")

    df["seg_len"] = df["end"] - df["start"]

    if unit == "individual":
        key = df["name"].astype(str)
    else:
        key = df["name"].astype(str) + "-" + df["haplotype"].astype(str)
    df = df.assign(**{unit: key})

    # bp per ancestry, then divide by that unit's total bp
    wide = df.groupby([unit, "ancestry"])["seg_len"].sum().unstack(fill_value=0)
    totals = wide.sum(axis=1)

    wide = wide.loc[totals > 0]
    totals = totals[totals > 0]

    props = wide.div(totals, axis=0).add_suffix("_prop")
    out = wide.add_suffix("_bp").join(props)
    out.index.name = unit
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit", default="individual",
                        choices=["individual", "haplotype"])
    args = parser.parse_args()

    if not LAI_BED_DIR.is_dir():
        sys.exit(
            f"Local ancestry files not found under {LAI_BED_DIR}.\n"
            "These data are not public. Set LAI_BED_DIR in config.local.sh or export it:\n"
            "    export LAI_BED_DIR=/path/to/autosome_lai_bedfiles"
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    summary = []
    for pop in AMR_POPS:
        props = ancestry_proportions(pop, args.unit)
        if props.empty:
            print(f"{pop}: no data at {lai_path(pop)}, skipping", file=sys.stderr)
            continue

        out = OUT_DIR / f"{pop}_{args.unit}.csv"
        props.to_csv(out)

        row = {"pop": pop, "n": len(props)}
        for column in sorted(c for c in props.columns if c.endswith("_prop")):
            ancestry = column.removesuffix("_prop")
            row[f"{ancestry}_mean"] = props[column].mean()
            row[f"{ancestry}_sd"] = props[column].std(ddof=1)
        summary.append(row)

        shown = ", ".join(
            f"{c.removesuffix('_prop')} {props[c].mean():.2f}±{props[c].std(ddof=1):.2f}"
            for c in sorted(props.columns) if c.endswith("_prop")
        )
        print(f"{pop}: n={len(props)} {args.unit}s   {shown}")

    if not summary:
        sys.exit("No populations had usable data.")

    summary_path = OUT_DIR / f"summary_{args.unit}.csv"
    pd.DataFrame(summary).to_csv(summary_path, index=False)
    print(f"\nwrote per-{args.unit} proportions and {summary_path.name} to {OUT_DIR}")


if __name__ == "__main__":
    main()
