#!/usr/bin/env python3
"""Merge the per-donor tract files of a multi-donor scenario into one file.

Same as the original notebook step (Chapter_1.ipynb cells 86-87): for each
haplotype, pool its tracts from every donor and merge overlapping or touching
intervals, so all donors count together as European ancestry.

    mxl_double_donor_recent: donors 8, 9     -> MXL_merged_89tracts_rep<rep>.csv.gz
    mxl_donor6:              donors 8 to 13  -> MXL_totalIBSmerged_rep<rep>.csv.gz

Usage:
    python figure_code/simulations/03_merge_donor_tracts.py <scenario> <POP> <rep>
    python figure_code/simulations/03_merge_donor_tracts.py mxl_donor6 MXL 1
    python figure_code/simulations/03_merge_donor_tracts.py mxl_donor6 MXL 1-1000   # a range

Output columns as before: Chromosome, Start, End, hap_id, length.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import pyranges as pr
import tskit

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import SIM_DIR  # noqa: E402
from scenarios import SCENARIOS, tract_file  # noqa: E402


def merge_one(scenario, pop, rep, source_ids):
    d = SIM_DIR / scenario / "admix_csvs"
    dfs = []
    for sid in source_ids:
        path = d / tract_file(pop, sid, rep)
        if not path.exists():
            sys.exit(f"Missing {path}\nRun 02_extract_tracts.py {scenario} {pop} {rep} first.")
        df = pd.read_csv(path)
        df["Chromosome"] = "chr1"
        dfs.append(df.rename(columns={"start": "Start", "end": "End"}))

    result = []
    for hap in set().union(*(set(df["hap_id"]) for df in dfs)):
        combined = pd.concat([df[df["hap_id"] == hap] for df in dfs], ignore_index=True)
        if combined.empty:
            continue
        merged = pr.PyRanges(combined[["Chromosome", "Start", "End"]]).merge().df
        merged["hap_id"] = hap
        merged["length"] = merged["End"] - merged["Start"]
        result.append(merged)
    final = pd.concat(result, ignore_index=True) if result else \
        pd.DataFrame(columns=["Chromosome", "Start", "End", "hap_id", "length"])
    out = d / SCENARIOS[scenario]["merged"].format(pop=pop, rep=rep)
    final.to_csv(out, index=False, compression="gzip")
    return out, len(final)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scenario", choices=[s for s, c in SCENARIOS.items() if "merged" in c])
    p.add_argument("pop")
    p.add_argument("reps", help="one replicate (5) or a range (1-1000)")
    args = p.parse_args()
    lo, _, hi = args.reps.partition("-")
    reps = range(int(lo), int(hi or lo) + 1)

    # donor ids from the model, via the first replicate's tree sequence
    ts = tskit.load(SIM_DIR / args.scenario / "tree_seqs" / f"{args.pop}_rep{reps[0]}.ts.gz")
    names = {p.metadata.get("name"): p.id for p in ts.populations()}
    ids = [names[s] for s in SCENARIOS[args.scenario]["sources"]]
    for rep in reps:
        out, n = merge_one(args.scenario, args.pop, rep, ids)
        print(f"rep {rep}: {n:,} merged tracts -> {out.name}", flush=True)


if __name__ == "__main__":
    main()
