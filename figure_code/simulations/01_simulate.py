#!/usr/bin/env python3
"""Simulate one replicate of one scenario with msprime.

Same simulation as the original script: the scenario's demes model
(data/demographic_models/), 100 haploid samples from the admixed population
(plus 100 IBS haploid samples with --with-ibs), 50 Mb, recombination rate
1e-8 per bp per generation, migrations recorded (02_extract_tracts.py reads
them to find European tracts).

The original runs did not set a random seed. --seed makes a run repeatable;
without it, each run differs, as before.

Usage:
    python figure_code/simulations/01_simulate.py <scenario> <POP> <rep> [--with-ibs] [--seed N]
    python figure_code/simulations/01_simulate.py no_bottleneck_donor MXL 1 --with-ibs

Scenarios (see scenarios.py): no_bottleneck_donor (standard), eur_donor
(single donor) for CLM, MXL, PEL, PUR; mxl_double_donor_recent and
mxl_donor6 for MXL. Figure S5 needs the IBS samples in the standard runs.

Output: $SIM_DIR/<scenario>/tree_seqs/<POP>_rep<rep>.ts.gz
"""

import argparse
import gzip
import sys
from pathlib import Path

import demes
import msprime

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import SIM_DIR  # noqa: E402
from scenarios import SCENARIOS  # noqa: E402

MODEL_DIR = _root / "data" / "demographic_models"
N_HAPLOID = 100
SEQUENCE_LENGTH = 50_000_000
RECOMBINATION_RATE = 1e-8


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scenario", choices=sorted(SCENARIOS))
    p.add_argument("pop", choices=["CLM", "MXL", "PEL", "PUR"])
    p.add_argument("rep", type=int)
    p.add_argument("--with-ibs", action="store_true", help="also sample 100 IBS haplotypes")
    p.add_argument("--seed", type=int, default=None, help="random seed (default: none, as originally)")
    p.add_argument("--sequence-length", type=float, default=SEQUENCE_LENGTH,
                   help="bp (default 50 Mb; smaller only for testing)")
    args = p.parse_args()

    model = MODEL_DIR / SCENARIOS[args.scenario]["demes"].format(pop_lower=args.pop.lower())
    if not model.exists():
        sys.exit(f"No model {model.name} for {args.pop} in scenario {args.scenario}")
    demography = msprime.Demography.from_demes(demes.load(model))

    samples = [msprime.SampleSet(N_HAPLOID, ploidy=1, population=args.pop)]
    if args.with_ibs:
        samples.append(msprime.SampleSet(N_HAPLOID, ploidy=1, population="IBS"))
    ts = msprime.sim_ancestry(samples=samples, demography=demography,
                              sequence_length=args.sequence_length,
                              recombination_rate=RECOMBINATION_RATE,
                              record_migrations=True, random_seed=args.seed)

    out = SIM_DIR / args.scenario / "tree_seqs" / f"{args.pop}_rep{args.rep}.ts.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wb") as f:
        ts.dump(f)
    print(f"{args.scenario} {args.pop} rep {args.rep}: {ts.num_trees:,} trees -> {out}")


if __name__ == "__main__":
    main()
