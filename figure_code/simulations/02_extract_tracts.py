#!/usr/bin/env python3
"""Extract European ancestry tracts from one simulated tree sequence.

Same algorithm as the original scripts/call_intro.py, with paths from
config.py and the source population looked up by name:
  1. from the recorded migrations, find tree intervals carrying material that
     migrated from the source population into the admixed population's samples;
  2. stitch consecutive intervals of each haplotype into tracts.

As in the original, the stitching does not merge overlapping intervals, so a
haplotype can have overlapping tracts. The analyses merge them before summing
(sim_ibd.covered_bp, and 03_merge_donor_tracts.py for multi-donor scenarios).

Usage:
    python figure_code/simulations/02_extract_tracts.py <scenario> <POP> <rep>
    python figure_code/simulations/02_extract_tracts.py eur_donor MXL 1
For multi-donor scenarios every donor is extracted, one file each.
Other sources by name, e.g. Neanderthal and Indigenous American tracts for
Figure S3 (NEA is id 2, MXB id 9 in the standard models):
    python figure_code/simulations/02_extract_tracts.py no_bottleneck_donor MXL 1 --sources NEA MXB
As a Slurm array over replicates:
    sbatch --array=1-1000 figure_code/simulations/02_extract_tracts.script <scenario> <POP>

Input:  $SIM_DIR/<scenario>/tree_seqs/<POP>_rep<rep>.ts.gz
Output: $SIM_DIR/<scenario>/admix_csvs/<POP>_<source id>tracts_rep<rep>_stitched.csv.gz
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import tskit

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import SIM_DIR  # noqa: E402
from scenarios import SCENARIOS, tract_file  # noqa: E402


def population_id(ts, name):
    for p in ts.populations():
        if p.metadata.get("name") == name:
            return p.id
    sys.exit(f"No population {name!r} in this tree sequence; found "
             f"{[p.metadata.get('name') for p in ts.populations()]}")


def get_intro(dest, target_populations, ts):
    """Tree intervals of target haplotypes covered by migrations into `dest`
    (backwards in time: lineages moving into the source population)."""
    mig = ts.tables.migrations
    mig = mig[mig.dest == dest]
    targets = sorted(h for pop in target_populations for h in ts.samples(population=pop))
    nodes, lefts, rights, times = mig.node, mig.left, mig.right, mig.time
    rec = {k: [] for k in ("hap_id", "t_l", "t_r", "i_l", "i_r", "o_len", "o_perct", "time")}
    for tree in ts.trees(tracked_samples=targets):
        t_l, t_r = tree.interval
        for j, node in enumerate(nodes):
            if not (t_l < rights[j] and lefts[j] < t_r):
                continue
            if tree.num_tracked_samples(node) == 0:
                continue
            leaves = set(tree.leaves(node))
            o_len = min(t_r, rights[j]) - max(t_l, lefts[j])
            for target in targets:
                if target in leaves:
                    rec["hap_id"].append(target)
                    rec["t_l"].append(t_l)
                    rec["t_r"].append(t_r)
                    rec["i_l"].append(lefts[j])
                    rec["i_r"].append(rights[j])
                    rec["o_len"].append(o_len)
                    rec["o_perct"].append(o_len / (t_r - t_l))
                    rec["time"].append(times[j])
    return pd.DataFrame(rec)


def extract_introgressed_tracts(intro_df):
    """Join each haplotype's consecutive intervals (end == next start) into tracts."""
    out = {k: [] for k in ("hap_id", "start", "end", "length", "time")}
    for hap_id in intro_df.hap_id.unique():
        hap_df = intro_df[intro_df.hap_id == hap_id]
        start = end = time = None
        for row in hap_df.itertuples(index=False):
            left = row.t_l if row.o_perct == 1 else max(row.t_l, row.i_l)
            right = row.t_r if row.o_perct == 1 else min(row.t_r, row.i_r)
            time = row.time
            if start is None and end is None:
                start, end = left, right
            elif end == left:
                end = right
            else:
                for k, v in zip(("hap_id", "time", "start", "end", "length"),
                                (hap_id, time, start, end, end - start)):
                    out[k].append(v)
                start, end = left, right
        if start is not None and end is not None:
            for k, v in zip(("hap_id", "time", "start", "end", "length"),
                            (hap_id, time, start, end, end - start)):
                out[k].append(v)
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scenario", choices=sorted(SCENARIOS))
    p.add_argument("pop", help="admixed population, e.g. MXL")
    p.add_argument("rep", type=int)
    p.add_argument("--sources", nargs="+", default=None,
                   help="source populations by name (default: the scenario's European source(s))")
    p.add_argument("--skip-existing", action="store_true",
                   help="leave tract files that already exist untouched")
    args = p.parse_args()

    ts_path = SIM_DIR / args.scenario / "tree_seqs" / f"{args.pop}_rep{args.rep}.ts.gz"
    if not ts_path.exists():
        sys.exit(f"Missing {ts_path}")
    ts = tskit.load(ts_path)
    out_dir = SIM_DIR / args.scenario / "admix_csvs"
    out_dir.mkdir(parents=True, exist_ok=True)
    target = population_id(ts, args.pop)
    for source in args.sources or SCENARIOS[args.scenario]["sources"]:
        sid = population_id(ts, source)
        out = out_dir / tract_file(args.pop, sid, args.rep)
        if args.skip_existing and out.exists():
            print(f"{source} (id {sid}): {out.name} exists, skipped")
            continue
        tracts = extract_introgressed_tracts(get_intro(sid, [target], ts))
        tracts.to_csv(out, index=False, compression="gzip")
        print(f"{source} (id {sid}): {len(tracts):,} tracts -> {out}")


if __name__ == "__main__":
    main()
