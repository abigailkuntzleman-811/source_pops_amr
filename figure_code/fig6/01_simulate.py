#!/usr/bin/env python3
"""Simulate one replicate of one demographic scenario.

Loads a demes model, simulates 100 haploid samples from the admixed American
population over 50 Mb with migrations recorded, and writes a gzipped tree
sequence. Migrations are recorded because 02_extract_tracts.py reads them to
work out which segments came from which source population.

Usage:
    python figure_code/fig6/01_simulate.py <scenario> <rep> [--pop MXL]

Example:
    python figure_code/fig6/01_simulate.py eur_donor 1

As a Slurm array over replicates:
    sbatch --array=1-100 figure_code/fig6/01_simulate.script eur_donor

Output:
    $SIM_DIR/<scenario>/tree_seqs/<POP>_rep<REP>.ts.gz
"""

import argparse
import gzip
import sys
from pathlib import Path

import demes
import msprime

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import DEMES_DIR, SIM_DIR  # noqa: E402
from scenarios import SCENARIOS  # noqa: E402

SEQUENCE_LENGTH = 50_000_000
RECOMBINATION_RATE = 1e-8
N_SAMPLES = 100


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=sorted(SCENARIOS))
    parser.add_argument("rep", type=int, help="replicate number")
    parser.add_argument("--pop", default="MXL",
                        help="admixed population to sample (default: MXL)")
    parser.add_argument("--with-ibs", action="store_true",
                        help="also sample 100 haplotypes from IBS")
    parser.add_argument("--sequence-length", type=float, default=SEQUENCE_LENGTH,
                        help="simulated sequence length in bp "
                             "(default: 50 Mb, as used for the paper; "
                             "lower values are for testing only)")
    parser.add_argument("--seed", type=int, default=None,
                        help="random seed; default derives one from the "
                             "scenario and replicate so runs are repeatable")
    args = parser.parse_args()

    config = SCENARIOS[args.scenario]
    model_path = DEMES_DIR / config["demes"]
    if not model_path.exists():
        sys.exit(f"Missing demographic model: {model_path}")

    out_dir = SIM_DIR / args.scenario / "tree_seqs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{args.pop}_rep{args.rep}.ts.gz"

    if out_path.exists():
        print(f"{out_path.name} exists, nothing to do")
        return

    # A deterministic seed, so a replicate can be regenerated identically.
    seed = args.seed if args.seed is not None else (
        abs(hash((args.scenario, args.pop, args.rep))) % (2 ** 31 - 1)
    )

    demography = msprime.Demography.from_demes(demes.load(model_path))

    samples = [msprime.SampleSet(N_SAMPLES, ploidy=1, population=args.pop)]
    if args.with_ibs:
        samples.append(msprime.SampleSet(N_SAMPLES, ploidy=1, population="IBS"))

    print(f"simulating {args.scenario} {args.pop} rep {args.rep} (seed {seed})")
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=args.sequence_length,
        recombination_rate=RECOMBINATION_RATE,
        record_migrations=True,
        random_seed=seed,
    )

    with gzip.open(out_path, "wb") as f:
        ts.dump(f)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
