#!/usr/bin/env python3
"""Extract contiguous ancestry tracts from a simulated tree sequence.

Reads the recorded migrations, finds the tree intervals carrying material that
migrated from the source population into the target haplotypes, then stitches
adjacent intervals into contiguous tracts.

Usage:
    python figure_code/fig6/02_extract_tracts.py <scenario> <rep> [--pop MXL]

As a Slurm array over replicates:
    sbatch --array=1-100 figure_code/fig6/02_extract_tracts.script eur_donor

Output:
    $SIM_DIR/<scenario>/admix_csvs/<POP>_<SOURCE_ID>tracts_rep<REP>_stitched.csv.gz

The number in that filename is the source population's id in the demographic
model, not a count. Scenarios whose ancestry comes from several donor
populations need those per-source files merged afterwards; this script does
not do that.
"""

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
import tskit

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import SIM_DIR  # noqa: E402
from scenarios import SCENARIOS  # noqa: E402


def population_id(ts: tskit.TreeSequence, name: str) -> int:
    """Numeric id of a named population in the simulated demography."""
    for pop in ts.populations():
        if pop.metadata.get("name") == name:
            return pop.id
    names = sorted(p.metadata.get("name", "?") for p in ts.populations())
    sys.exit(f"No population named {name!r}. Available: {names}")


def get_intro(ts: tskit.TreeSequence, source_id: int, target_id: int) -> pd.DataFrame:
    """Tree intervals carrying material migrated from source into target."""
    mig = ts.tables.migrations
    mig = mig[mig.dest == source_id]

    target_haplotypes = sorted(ts.samples(population=target_id))

    nodes, lefts, rights, times = mig.node, mig.left, mig.right, mig.time

    records = {k: [] for k in
               ("hap_id", "t_l", "t_r", "i_l", "i_r", "o_len", "o_perct", "time")}

    for tree in ts.trees(tracked_samples=target_haplotypes):
        t_l, t_r = tree.interval
        for j, node in enumerate(nodes):
            # Does this migrating segment overlap the current tree?
            if not (t_l < rights[j] and lefts[j] < t_r):
                continue
            if tree.num_tracked_samples(node) == 0:
                continue

            leaves = set(tree.get_leaves(node))
            overlap = min(t_r, rights[j]) - max(t_l, lefts[j])

            for target in target_haplotypes:
                if target not in leaves:
                    continue
                records["hap_id"].append(target)
                records["t_l"].append(t_l)
                records["t_r"].append(t_r)
                records["i_l"].append(lefts[j])
                records["i_r"].append(rights[j])
                records["o_len"].append(overlap)
                records["o_perct"].append(overlap / (t_r - t_l))
                records["time"].append(times[j])

    return pd.DataFrame(records)


def extract_tracts(intro_df: pd.DataFrame, merge: bool = True) -> pd.DataFrame:
    """Stitch intervals into contiguous tracts per haplotype.

    With merge=True (the default) the intervals are sorted by position and
    overlapping ones are combined. Without that, a haplotype covered by
    several migration records yields overlapping tracts whose lengths are
    then counted twice: summed tract length can exceed the sequence length,
    giving ancestry proportions above 1. The original implementation did not
    sort or merge; merge=False reproduces it for comparison.
    """
    tracts = {k: [] for k in ("hap_id", "start", "end", "length", "time")}

    for hap_id in intro_df.hap_id.unique():
        hap_df = intro_df[intro_df.hap_id == hap_id]
        if merge:
            hap_df = hap_df.assign(
                _l=hap_df.apply(
                    lambda r: r.t_l if r.o_perct == 1 else max(r.t_l, r.i_l),
                    axis=1),
            ).sort_values("_l")
        start = end = time = None

        for row in hap_df.itertuples(index=False):
            left = row.t_l if row.o_perct == 1 else max(row.t_l, row.i_l)
            right = row.t_r if row.o_perct == 1 else min(row.t_r, row.i_r)
            time = row.time

            if start is None:
                start, end = left, right
            elif (left <= end) if merge else (end == left):
                # contiguous with, or overlapping, the tract being built
                end = max(end, right)
            else:
                tracts["hap_id"].append(hap_id)
                tracts["time"].append(time)
                tracts["start"].append(start)
                tracts["end"].append(end)
                tracts["length"].append(end - start)
                start, end = left, right

        if start is not None:
            tracts["hap_id"].append(hap_id)
            tracts["time"].append(time)
            tracts["start"].append(start)
            tracts["end"].append(end)
            tracts["length"].append(end - start)

    return pd.DataFrame(tracts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=sorted(SCENARIOS))
    parser.add_argument("rep", type=int)
    parser.add_argument("--pop", default="MXL")
    parser.add_argument("--no-merge", action="store_true",
                        help="reproduce the original behaviour: do not sort "
                             "or merge overlapping intervals. Produces "
                             "ancestry proportions above 1 for some "
                             "haplotypes.")
    parser.add_argument("--source", default=None,
                        help="source population name; defaults to the "
                             "scenario's admix_pop")
    args = parser.parse_args()

    config = SCENARIOS[args.scenario]
    source_name = args.source or config.get("admix_pop")
    if source_name is None:
        sys.exit(
            f"Scenario {args.scenario!r} draws ancestry from several donor "
            "populations and has no single source. Pass --source explicitly, "
            "once per donor, and merge the results."
        )

    ts_path = SIM_DIR / args.scenario / "tree_seqs" / f"{args.pop}_rep{args.rep}.ts.gz"
    if not ts_path.exists():
        sys.exit(f"Missing {ts_path}\nRun figure_code/fig6/01_simulate.py first.")

    with gzip.open(ts_path, "rb") as f:
        ts = tskit.load(f)

    source_id = population_id(ts, source_name)
    target_id = population_id(ts, args.pop)

    out_dir = SIM_DIR / args.scenario / "admix_csvs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{args.pop}_{source_id}tracts_rep{args.rep}_stitched.csv.gz"

    if out_path.exists():
        print(f"{out_path.name} exists, nothing to do")
        return

    print(f"extracting {source_name} (id {source_id}) tracts in {args.pop}")
    intro_df = get_intro(ts, source_id, target_id)
    print(f"  {len(intro_df):,} overlapping tree intervals")

    tract_df = extract_tracts(intro_df, merge=not args.no_merge)
    tract_df.to_csv(out_path, index=False, compression="gzip")

    print(f"  {len(tract_df):,} tracts, {tract_df['length'].sum():,.0f} bp total")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
