#!/usr/bin/env python3
"""Find IBD pileups in the hap-ibd calls and measure how much they contribute.

A pileup is a short region where far more pairs appear to share an IBD
segment than anywhere else, usually an artifact of segmental duplications or
genetic map errors near centromeres. This script only reports; it changes
nothing.

1. For every 1 Mb bin of each autosome, counts the hap-ibd segments (all
   samples, all pairs) that overlap it.
2. Flags bins whose count is more than --fold times the genome-wide median
   (median over bins with any segments), and merges neighbouring flagged bins
   into regions.
3. For the admixed populations, using unrelated individuals only, reports how
   much within-population IBD comes from segments that overlap a flagged
   region, and the mean per-individual IBD (as in Figure 4) with and without
   those segments.

Usage:
    python hapIBD_calls/02_find_pileups.py
    python hapIBD_calls/02_find_pileups.py --fold 5

Output (in IBD_DIR):
    ibd_pileup_mask.bed            flagged regions: chrom, start, end, peak count, fold
    ibd_pileup_bins.csv.gz         segment count for every 1 Mb bin
    ibd_pileup_impact.csv          per-population effect on within-population IBD
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import AMR_POPS, CHROMOSOMES, IBD_DIR  # noqa: E402

sys.path.insert(0, str(_root / "figure_code" / "fig3"))
_fig3 = __import__("01_within_population_ibd")
IBD_COLUMNS = _fig3.IBD_COLUMNS
load_population_map = _fig3.load_population_map
load_related_samples = _fig3.load_related_samples

BIN = 1_000_000


def bin_counts(seg: pd.DataFrame) -> np.ndarray:
    """Number of segments overlapping each 1 Mb bin."""
    first = (seg["start_bp"].to_numpy() // BIN).astype(int)
    last = ((seg["end_bp"].to_numpy() - 1) // BIN).astype(int)
    diff = np.zeros(last.max() + 2, dtype=np.int64)
    np.add.at(diff, first, 1)
    np.add.at(diff, last + 1, -1)
    return np.cumsum(diff)[:-1]


def merge_flagged(chrom: str, counts: np.ndarray, flagged: np.ndarray, median: float):
    """Runs of consecutive flagged bins, as BED rows."""
    rows, start = [], None
    for b, f in enumerate(list(flagged) + [False]):
        if f and start is None:
            start = b
        elif not f and start is not None:
            peak = int(counts[start:b].max())
            rows.append((chrom, start * BIN, b * BIN, peak, round(peak / median, 1)))
            start = None
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fold", type=float, default=3.0,
                        help="flag bins above this multiple of the median count (default 3)")
    parser.add_argument("--pops", nargs="+", default=AMR_POPS)
    args = parser.parse_args()

    pop_of = load_population_map()
    related = load_related_samples()
    unrelated = {pop: sorted(set(pop_of[pop_of == pop].index) - related) for pop in args.pops}
    keep_ids = set().union(*unrelated.values())

    # Pass 1: coverage per bin, plus the within-population segments we need later.
    counts, within = {}, []
    for chrom in CHROMOSOMES:
        n = chrom.removeprefix("chr")
        path = IBD_DIR / f"kgp_chr{n}_segments.ibd.gz"
        if not path.exists():
            sys.exit(f"Missing {path}")
        seg = pd.read_csv(path, sep="\t", header=None, names=IBD_COLUMNS)
        counts[chrom] = bin_counts(seg)
        print(f"{chrom}: {len(seg):,} segments, busiest bin "
              f"{int(counts[chrom].argmax())} Mb with {int(counts[chrom].max()):,}", flush=True)

        w = seg[seg["id1"].isin(keep_ids) & seg["id2"].isin(keep_ids)
                & (seg["id1"] != seg["id2"])].copy()
        w["pop1"] = w["id1"].map(pop_of)
        w = w[w["pop1"] == w["id2"].map(pop_of)]
        w["chrom_name"] = chrom
        within.append(w)

    # Genome-wide median over bins that have any segments (skips centromere gaps).
    all_counts = np.concatenate(list(counts.values()))
    median = float(np.median(all_counts[all_counts > 0]))
    threshold = args.fold * median
    print(f"\nmedian segments per 1 Mb bin: {median:,.0f}; flagging bins above {threshold:,.0f}")

    mask_rows, bin_rows = [], []
    for chrom, c in counts.items():
        flagged = c > threshold
        mask_rows += merge_flagged(chrom, c, flagged, median)
        bin_rows.append(pd.DataFrame({"chrom": chrom, "bin_start": np.arange(len(c)) * BIN,
                                      "segments": c, "flagged": flagged}))
    mask = pd.DataFrame(mask_rows, columns=["chrom", "start", "end", "peak_segments", "fold_over_median"])

    IBD_DIR.mkdir(parents=True, exist_ok=True)
    mask.to_csv(IBD_DIR / "ibd_pileup_mask.bed", sep="\t", index=False, header=False)
    pd.concat(bin_rows).to_csv(IBD_DIR / "ibd_pileup_bins.csv.gz", index=False)

    print(f"\n{len(mask)} flagged region(s), {(mask['end'] - mask['start']).sum() / 1e6:.0f} Mb in total:")
    print(mask.to_string(index=False) if len(mask) else "  none")

    # Pass 2: how much within-population IBD overlaps a flagged region.
    within = pd.concat(within, ignore_index=True)
    hits = np.zeros(len(within), dtype=bool)
    for chrom, start, end, *_ in mask_rows:
        on = within["chrom_name"].to_numpy() == chrom
        hits |= on & (within["start_bp"].to_numpy() < end) & (within["end_bp"].to_numpy() > start)
    within["in_pileup"] = hits

    rows = []
    for pop in args.pops:
        sub = within[within["pop1"] == pop]
        n = len(unrelated[pop])
        if n < 2:
            continue
        total, masked = sub["cM"].sum(), sub.loc[sub["in_pileup"], "cM"].sum()
        # mean over individuals of Figure 4's IBD_bar_i = 2 * total / (n (n - 1))
        per_ind = lambda cm: 2 * cm / (n * (n - 1))
        rows.append({
            "pop": pop, "n_unrelated": n,
            "within_segments": len(sub),
            "segments_in_pileups": int(sub["in_pileup"].sum()),
            "cM_total": round(total, 1),
            "cM_in_pileups": round(masked, 1),
            "pct_cM_in_pileups": round(100 * masked / total, 1) if total else np.nan,
            "mean_IBD_bar_i_now": round(per_ind(total), 2),
            "mean_IBD_bar_i_masked": round(per_ind(total - masked), 2),
        })
    impact = pd.DataFrame(rows)
    impact.to_csv(IBD_DIR / "ibd_pileup_impact.csv", index=False)
    print("\nwithin-population IBD, unrelated individuals:")
    print(impact.to_string(index=False))
    print(f"\nwrote {IBD_DIR / 'ibd_pileup_mask.bed'}, ibd_pileup_bins.csv.gz, ibd_pileup_impact.csv")


if __name__ == "__main__":
    main()
