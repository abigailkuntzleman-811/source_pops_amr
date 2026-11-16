#!/usr/bin/env python3
"""Within-population IBD sharing per haplotype pair.

Sums hap-ibd segment lengths (cM) across autosomes for every pair of
haplotypes belonging to two individuals of the same population, then averages
over all such pairs.

Pairs sharing no IBD count as zero, which matters: restricting to pairs with
detected IBD would bias the average upward and differently in each population.
Every possible within-population haplotype pair is enumerated first, so the
denominator is the number of possible pairs rather than the number observed.

Related individuals are excluded. Set RELATED_EXCLUDE to a file of sample IDs
(one per line, or PLINK-style FID/IID pairs as KING writes). If unset, the
related samples named in the 1000 Genomes pedigree are used instead.

Prerequisite: hapIBD_calls/01_call_ibd.sh for all 22 autosomes.

Usage:
    python fig3/01_within_population_ibd.py
    python fig3/01_within_population_ibd.py --pops MXL PEL PUR CLM

Output:
    data/processed/ibd/within_population_pairs.csv.gz   every pair and its cM
    data/processed/ibd/within_population_summary.csv    mean, SD and n per pop
"""

import argparse
import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents
             if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import (DATA_RAW, IBD_DIR, RELATED_EXCLUDE,  # noqa: E402
                    CHROMOSOMES)

IBD_COLUMNS = ["id1", "hap1", "id2", "hap2", "chrom", "start_bp", "end_bp", "cM"]

DEFAULT_POPS = [
    "CHB", "JPT", "CHS", "CDX", "KHV",
    "CEU", "TSI", "GBR", "FIN", "IBS",
    "GIH", "PJL", "BEB", "STU", "ITU",
    "MXL", "PEL", "PUR", "CLM",
]


def load_population_map() -> pd.Series:
    """Sample ID -> population, from the 3,202-sample pedigree.

    The pedigree covers the whole 30x release. The phase 3 panel is only a
    fallback: it lists the 2,504 unrelated samples alone, so any population
    that gained related individuals in the 30x release would be averaged over
    the wrong sample set. CEU is the clearest case, being the CEPH collection:
    118 samples in the pedigree against 90 that survive a panel-based lookup.
    """
    meta = DATA_RAW / "metadata"
    ped = meta / "20130606_g1k_3202_samples_ped_population.txt"

    if ped.exists():
        df = pd.read_csv(ped, sep=r"\s+")
        if {"SampleID", "Population"} <= set(df.columns):
            return df.set_index("SampleID")["Population"]
        print(f"warning: {ped.name} lacks SampleID/Population columns; "
              "falling back to the phase 3 panel", file=sys.stderr)

    panel = meta / "integrated_call_samples_v3.20130502.ALL.panel"
    if panel.exists():
        print("note: using the phase 3 panel for population labels. It covers "
              "only the 2,504 unrelated samples, so related individuals in the "
              "IBD output will be dropped.", file=sys.stderr)
        return pd.read_csv(panel, sep=r"\s+").set_index("sample")["pop"]

    sys.exit(
        f"No sample metadata found in {meta}\n"
        "Fetch the 3,202-sample pedigree:\n"
        "  mkdir -p data/raw/metadata && wget -P data/raw/metadata \\\n"
        "    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000G_2504_high_coverage/working/"
        "20130606_g1k_3202_samples_ped_population.txt"
    )


def load_related_samples() -> set:
    """Sample IDs to exclude as related."""
    if RELATED_EXCLUDE is not None:
        if not RELATED_EXCLUDE.exists():
            sys.exit(f"RELATED_EXCLUDE points at a missing file: {RELATED_EXCLUDE}")
        ids = set()
        for line in RELATED_EXCLUDE.read_text().splitlines():
            fields = line.split()
            if fields:
                # one ID per line, or PLINK-style FID IID
                ids.add(fields[-1])
        print(f"excluding {len(ids):,} samples listed in {RELATED_EXCLUDE.name}")
        return ids

    ped = DATA_RAW / "metadata" / "20130606_g1k_3202_samples_ped_population.txt"
    if not ped.exists():
        sys.exit(
            f"Missing {ped}\n"
            "Run fig1/00_download_data.sh, or set RELATED_EXCLUDE."
        )
    df = pd.read_csv(ped, sep=r"\s+")
    # A sample with a recorded parent is a child in a trio or duo.
    has_parent = (df["fatherID"].astype(str) != "0") | \
                 (df["motherID"].astype(str) != "0")
    ids = set(df.loc[has_parent, "SampleID"])
    print(f"excluding {len(ids):,} samples with a parent in the 1kGP pedigree")
    return ids


def read_ibd_segments() -> pd.DataFrame:
    """All autosomal hap-ibd segments, concatenated."""
    frames = []
    for chrom in CHROMOSOMES:
        n = chrom.removeprefix("chr")
        path = IBD_DIR / f"kgp_chr{n}_segments.ibd.gz"
        if not path.exists():
            sys.exit(
                f"Missing {path}\nRun: sbatch --array=1-22 hapIBD_calls/01_call_ibd.script"
            )
        frames.append(pd.read_csv(path, sep="\t", header=None, names=IBD_COLUMNS))
        print(f"  chr{n}: {len(frames[-1]):,} segments", flush=True)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pops", nargs="+", default=DEFAULT_POPS)
    args = parser.parse_args()

    pop_of = load_population_map()
    related = load_related_samples()

    print("reading IBD segments")
    ibd = read_ibd_segments()
    print(f"  {len(ibd):,} segments total")

    # Drop related individuals and anything outside the requested populations.
    ibd = ibd[~ibd["id1"].isin(related) & ~ibd["id2"].isin(related)]
    ibd["pop1"] = ibd["id1"].map(pop_of)
    ibd["pop2"] = ibd["id2"].map(pop_of)
    ibd = ibd[(ibd["pop1"] == ibd["pop2"]) & ibd["pop1"].isin(args.pops)]
    print(f"  {len(ibd):,} within-population segments after exclusions")

    # Haplotype pair key, order-independent so id order never matters.
    h1 = ibd["id1"] + "-" + ibd["hap1"].astype(str)
    h2 = ibd["id2"] + "-" + ibd["hap2"].astype(str)
    lo = np.minimum(h1, h2)
    hi = np.maximum(h1, h2)
    ibd = ibd.assign(pair=lo + "_" + hi)

    observed = ibd.groupby(["pop1", "pair"])["cM"].sum()

    rows, summary = [], []
    for pop in args.pops:
        samples = sorted(set(pop_of[pop_of == pop].index) - related)
        if len(samples) < 2:
            print(f"{pop}: fewer than two unrelated samples, skipping",
                  file=sys.stderr)
            continue

        # Every possible within-population haplotype pair, so pairs with no
        # detected IBD are counted as zero rather than dropped.
        pairs = []
        for a, b in combinations(samples, 2):
            for ha, hb in product(("1", "2"), repeat=2):
                x, y = f"{a}-{ha}", f"{b}-{hb}"
                pairs.append(f"{min(x, y)}_{max(x, y)}")

        values = pd.Series(0.0, index=pd.Index(pairs, name="pair"))
        if pop in observed.index.get_level_values(0):
            found = observed.loc[pop]
            values.update(found.reindex(values.index).dropna())

        rows.append(pd.DataFrame({"pop": pop, "pair": values.index,
                                  "cM": values.to_numpy()}))
        summary.append({
            "pop": pop,
            "n_individuals": len(samples),
            "n_pairs": len(values),
            "mean_cM": float(values.mean()),
            "sd_cM": float(values.std(ddof=1)),
            "frac_pairs_zero": float((values == 0).mean()),
        })
        print(f"{pop}: {len(samples)} individuals, {len(values):,} haplotype "
              f"pairs, mean {values.mean():.4f} cM")

    if not rows:
        sys.exit("No populations produced results.")

    IBD_DIR.mkdir(parents=True, exist_ok=True)

    pairs_out = IBD_DIR / "within_population_pairs.csv.gz"
    pd.concat(rows, ignore_index=True).to_csv(pairs_out, index=False)

    summary_out = IBD_DIR / "within_population_summary.csv"
    pd.DataFrame(summary).to_csv(summary_out, index=False)

    print(f"\nwrote {pairs_out}")
    print(f"wrote {summary_out}")


if __name__ == "__main__":
    main()
