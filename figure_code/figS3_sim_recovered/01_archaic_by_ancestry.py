#!/usr/bin/env python3
"""Simulated archaic (Neanderthal) sequence on European and on Indigenous
American tracts, per haplotype, for one replicate (Figure S3).

Same as the original scripts/get_admixed_introed_real.py: for each haplotype,
intersect its Neanderthal tracts with its European tracts and with its
Indigenous American (MXB) tracts.

Needs the tract files from figure_code/simulations/02_extract_tracts.py for the
European source and for NEA and MXB:
    python figure_code/simulations/02_extract_tracts.py no_bottleneck_donor MXL 1
    python figure_code/simulations/02_extract_tracts.py no_bottleneck_donor MXL 1 --sources NEA MXB

Usage:
    python figure_code/figS3_sim_recovered/01_archaic_by_ancestry.py <scenario> <POP> <reps>
    python figure_code/figS3_sim_recovered/01_archaic_by_ancestry.py no_bottleneck_donor MXL 1-100

Output ($SIM_DIR/<scenario>/admix_csvs/):
    <POP>_eur_arc_ancestry_rep<rep>.csv, <POP>_amr_arc_ancestry_rep<rep>.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import pyranges as pr
import tskit

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(_root / "figure_code" / "simulations"))
from config import SIM_DIR  # noqa: E402
from scenarios import SCENARIOS, tract_file  # noqa: E402


def read_tracts(path):
    df = pd.read_csv(path)
    df["hap_id"] = df["hap_id"].astype(int)
    df["Chromosome"] = "chr1"
    return df.rename(columns={"start": "Start", "end": "End"})


def arc_on(arc_df, anc_df):
    out = []
    for hap in set(arc_df.hap_id) & set(anc_df.hap_id):
        a = pr.PyRanges(arc_df.loc[arc_df.hap_id == hap, ["Chromosome", "Start", "End"]])
        b = pr.PyRanges(anc_df.loc[anc_df.hap_id == hap, ["Chromosome", "Start", "End"]])
        overlap = a.intersect(b)
        if len(overlap) > 0:
            df = overlap.df.copy()
            df["hap_id"] = hap
            out.append(df)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scenario", choices=sorted(SCENARIOS))
    p.add_argument("pop")
    p.add_argument("reps", help="one replicate (5) or a range (1-100)")
    args = p.parse_args()
    lo, _, hi = args.reps.partition("-")
    d = SIM_DIR / args.scenario / "admix_csvs"

    for rep in range(int(lo), int(hi or lo) + 1):
        ts = tskit.load(SIM_DIR / args.scenario / "tree_seqs" / f"{args.pop}_rep{rep}.ts.gz")
        ids = {q.metadata.get("name"): q.id for q in ts.populations()}
        eur_id = ids[SCENARIOS[args.scenario]["sources"][0]]
        arc = read_tracts(d / tract_file(args.pop, ids["NEA"], rep))
        eur = read_tracts(d / tract_file(args.pop, eur_id, rep))
        mxb = read_tracts(d / tract_file(args.pop, ids["MXB"], rep))
        arc_on(arc, eur).to_csv(d / f"{args.pop}_eur_arc_ancestry_rep{rep}.csv", index=False)
        arc_on(arc, mxb).to_csv(d / f"{args.pop}_amr_arc_ancestry_rep{rep}.csv", index=False)
        print(f"{args.pop} rep {rep} done", flush=True)


if __name__ == "__main__":
    main()
