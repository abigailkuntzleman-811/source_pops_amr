#!/usr/bin/env python3
"""Archaic sequence recovered per admixed population from the hmmix calls (Figure S2),
and the downsampled sets used for the UpSet plot.

Same as the original scripts/get_archaic_recovered_from_hmmix.py:
  - recovery curve: add individuals one at a time (in file order) and record
    the merged archaic sequence (bp) across both haplotypes and the autosomes;
  - downsampling: draw the smallest population's sample size from each
    population (numpy default_rng(42), populations in the order CLM, PUR, PEL,
    MXL), and write their merged archaic sequence as a BED file.

Usage:
    python figure_code/figS2_hmmix_recovered/01_recovered_from_hmmix.py all
    python figure_code/figS2_hmmix_recovered/01_recovered_from_hmmix.py neanderthal
    python figure_code/figS2_hmmix_recovered/01_recovered_from_hmmix.py denisovan

Input:  data/processed/hmmix_filtered/hmmix[_<type>]_filtered.csv  (figure_code/fig1/01_filter_hmmix_calls.py)
Output ($RECOVERED_DIR):
    <POP>_total_<type>_recovered_from_hmmix_pyranges.txt   curve (x: individuals, y: bp)
    <POP>_n<N>_recovered_<type>.bed                         downsampled merged sequence
    <POP>_n<N>_recovered_<type>_samples.txt                 the sampled individuals
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import DATA_PROCESSED, RECOVERED_DIR  # noqa: E402

FILES = {"all": "hmmix_filtered.csv",
         "neanderthal": "hmmix_neanderthal_filtered.csv",
         "denisovan": "hmmix_denisovan_filtered.csv"}
POPS = ["CLM", "PUR", "PEL", "MXL"]     # order matters for the seeded downsampling


def recovered(pop_df, inds):
    chroms = [f"chr{i}" for i in range(1, 23)]
    df = pop_df[pop_df["name"].isin(inds) & pop_df["chrom"].isin(chroms)
                & pop_df["haplotype"].isin(["hap1", "hap2"])]
    df = df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    return pr.PyRanges(df[["Chromosome", "Start", "End"]]).merge()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("archaic", choices=sorted(FILES))
    args = p.parse_args()

    path = DATA_PROCESSED / "hmmix_filtered" / FILES[args.archaic]
    if not path.exists():
        sys.exit(f"Missing {path}\nRun: python figure_code/fig1/01_filter_hmmix_calls.py")
    hmmix = pd.read_csv(path)
    rng = np.random.default_rng(42)

    by_pop = {pop: hmmix[hmmix["pop"] == pop] for pop in POPS}
    inds = {pop: by_pop[pop]["name"].unique() for pop in POPS}
    n_min = min(len(v) for v in inds.values())
    RECOVERED_DIR.mkdir(parents=True, exist_ok=True)

    for pop in POPS:
        x, y = [], []
        for i in range(1, len(inds[pop]) + 1):
            x.append(i)
            y.append(recovered(by_pop[pop], inds[pop][:i]).lengths().sum())
        with open(RECOVERED_DIR / f"{pop}_total_{args.archaic}_recovered_from_hmmix_pyranges.txt", "w") as f:
            f.write("x: " + ", ".join(map(str, x)) + "\n")
            f.write("y: " + ", ".join(map(str, y)) + "\n")

        sampled = rng.choice(np.sort(inds[pop]), size=n_min, replace=False)
        merged = recovered(by_pop[pop], sampled)
        pd.DataFrame({"chrom": merged.Chromosome, "start": merged.Start, "end": merged.End}).to_csv(
            RECOVERED_DIR / f"{pop}_n{n_min}_recovered_{args.archaic}.bed", sep="\t", header=False, index=False)
        (RECOVERED_DIR / f"{pop}_n{n_min}_recovered_{args.archaic}_samples.txt").write_text(
            "\n".join(sorted(sampled)) + "\n")
        print(f"{pop}: {len(inds[pop])} individuals, {y[-1]:,} bp; downsampled to {n_min}", flush=True)


if __name__ == "__main__":
    main()
