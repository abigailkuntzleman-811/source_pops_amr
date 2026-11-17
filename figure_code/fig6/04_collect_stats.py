#!/usr/bin/env python3
"""Collect European ancestry proportion and within-population IBD per individual.

For every replicate of every scenario, loads the tree sequence, computes how
much IBD each individual shares with the rest of the simulated population, and
reads their European ancestry proportion from the stitched tract file. Writes
one tidy table, so the plotting step does not have to reload hundreds of tree
sequences.

This is the expensive step: each replicate loads a 50 Mb tree sequence and
runs an IBD scan.

Prerequisites: 01_simulate.py, 02_extract_tracts.py, and 03_merge_donor_tracts.py
for scenarios with several donors.

Usage:
    python figure_code/fig6/04_collect_stats.py
    python figure_code/fig6/04_collect_stats.py --scenarios eur_donor --reps 10

Output:
    $SIM_DIR/simulation_stats.csv    scenario, rep, individual, eur_prop, ibd_bp
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

from config import SIM_DIR, SIM_REPS  # noqa: E402
from scenarios import PLOT_ORDER, SCENARIOS  # noqa: E402

# Minimum IBD segment length, as in the original analysis.
MIN_SPAN = 1_000_000


def population_id(ts: tskit.TreeSequence, name: str) -> int:
    for pop in ts.populations():
        if pop.metadata.get("name") == name:
            return pop.id
    names = sorted(p.metadata.get("name", "?") for p in ts.populations())
    sys.exit(f"No population named {name!r}. Available: {names}")


def ibd_per_individual(ts: tskit.TreeSequence, target_pop: str) -> dict:
    """Total IBD each individual shares with others in their population."""
    target_id = population_id(ts, target_pop)

    node_to_ind = {n: ind.id for ind in ts.individuals() for n in ind.nodes}

    target_nodes = [
        n for ind in ts.individuals()
        if ts.node(ind.nodes[0]).population == target_id
        for n in ind.nodes
    ]

    totals = {
        ind.id: 0.0 for ind in ts.individuals()
        if ts.node(ind.nodes[0]).population == target_id
    }

    segments = ts.ibd_segments(within=target_nodes, min_span=MIN_SPAN,
                               store_segments=True)

    for (node1, node2), segment_list in segments.items():
        ind1, ind2 = node_to_ind.get(node1), node_to_ind.get(node2)
        for seg in segment_list:
            span = seg.right - seg.left
            if ind1 is not None:
                totals[ind1] += span
            if ind2 is not None:
                totals[ind2] += span

    return totals


def admixture_proportions(tract_path: Path, individuals, seq_length: float) -> dict:
    """European ancestry proportion per individual, from stitched tracts.

    The denominator is the simulated sequence length, taken from the tree
    sequence rather than assumed, so a run at a different length still gives
    proportions.
    """
    df = pd.read_csv(tract_path)
    df["hap_id"] = df["hap_id"].astype(int)
    by_hap = df.groupby("hap_id")["length"].sum()
    # Samples are haploid, so hap_id indexes individuals directly.
    return {ind: by_hap.get(ind, 0.0) / seq_length for ind in individuals}


def process_replicate(scenario: str, rep: int, pop: str):
    config = SCENARIOS[scenario]

    ts_path = SIM_DIR / scenario / "tree_seqs" / f"{pop}_rep{rep}.ts.gz"
    tract_path = (SIM_DIR / scenario / "admix_csvs" /
                  config["tracts"].format(pop=pop, rep=rep))

    if not ts_path.exists() or not tract_path.exists():
        return None

    with gzip.open(ts_path, "rb") as f:
        ts = tskit.load(f)

    ibd = ibd_per_individual(ts, pop)
    props = admixture_proportions(tract_path, ibd.keys(), ts.sequence_length)

    return pd.DataFrame({
        "scenario": scenario,
        "rep": rep,
        "individual": list(ibd.keys()),
        "eur_prop": [props[i] for i in ibd],
        "ibd_bp": [ibd[i] for i in ibd],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenarios", nargs="+", default=PLOT_ORDER,
                        choices=sorted(SCENARIOS))
    parser.add_argument("--pop", default="MXL")
    parser.add_argument("--reps", type=int, default=SIM_REPS)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    frames, missing = [], []
    for scenario in args.scenarios:
        n_done = 0
        for rep in range(1, args.reps + 1):
            result = process_replicate(scenario, rep, args.pop)
            if result is None:
                missing.append(f"{scenario}/rep{rep}")
                continue
            frames.append(result)
            n_done += 1
            if n_done % 10 == 0:
                print(f"  {scenario}: {n_done} replicates", flush=True)
        print(f"{scenario}: {n_done} replicates collected")

    if missing:
        print(f"\nwarning: {len(missing)} replicate(s) missing inputs, e.g. "
              f"{', '.join(missing[:5])}", file=sys.stderr)
    if not frames:
        sys.exit("Nothing collected. Run 01_simulate.py and 02_extract_tracts.py.")

    out = args.out or (SIM_DIR / "simulation_stats.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
