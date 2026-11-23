#!/usr/bin/env python3
"""Table of the IBD pileup regions masked by 03_mask_pileups.py.

For each masked region: location, size, peak segment coverage and how many
times the genome-wide median that is; where the region sits relative to the
genetic map (inside it, before its first position, after its last); the map's
recombination rate in the region (blank outside the map); the rate hap-ibd
actually used there, as the median cM per Mb of the raw segments lying entirely
inside the region; and the chromosome's average rate for comparison.

Outside the map's range hap-ibd extrapolates, so the map itself says nothing
about those regions and only the hap-ibd rate is meaningful.

Usage:
    python hapIBD_calls/04_masked_regions_table.py
    python hapIBD_calls/04_masked_regions_table.py --map-template /path/plink.chr{CHROM}.GRCh38.map

Output (in IBD_DIR):
    ibd_peaks_masked_table.csv   the table
    ibd_peaks_masked_table.tex   the same rows as a LaTeX tabular (no caption)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import GENETIC_MAP_DIR, IBD_DIR  # noqa: E402

sys.path.insert(0, str(_root / "figure_code" / "fig3"))
IBD_COLUMNS = __import__("01_within_population_ibd").IBD_COLUMNS


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--map-template", default=str(GENETIC_MAP_DIR / "plink.chr{CHROM}.GRCh38.map"),
                   help="genetic map path with {CHROM} (default: GENETIC_MAP_DIR/plink.chr{CHROM}.GRCh38.map)")
    args = p.parse_args()

    bed = IBD_DIR / "ibd_peaks_masked.bed"
    if not bed.exists():
        sys.exit(f"Missing {bed}\nRun: python hapIBD_calls/03_mask_pileups.py")
    peaks = pd.read_csv(bed, sep="\t", header=None,
                        names=["chrom", "start", "end", "peak_segments", "fold_over_median"])

    rows, maps, segs = [], {}, {}
    for r in peaks.itertuples(index=False):
        n = r.chrom.removeprefix("chr")
        if r.chrom not in maps:
            path = Path(args.map_template.replace("{CHROM}", n))
            if not path.exists():
                sys.exit(f"Missing genetic map {path}; pass --map-template")
            m = pd.read_csv(path, sep=r"\s+", header=None, usecols=[2, 3], names=["cM", "bp"]).sort_values("bp")
            maps[r.chrom] = (m["bp"].to_numpy(float), m["cM"].to_numpy(float))
        bp, cm = maps[r.chrom]
        if r.chrom not in segs:
            segs[r.chrom] = pd.read_csv(IBD_DIR / f"kgp_chr{n}_segments.ibd.gz", sep="\t", header=None,
                                        names=IBD_COLUMNS, usecols=["start_bp", "end_bp", "cM"])
        size_mb = (r.end - r.start) / 1e6
        chrom_rate = (cm.max() - cm.min()) / ((bp.max() - bp.min()) / 1e6)
        if r.end <= bp[0]:
            status, region_rate = "before map start", np.nan
        elif r.start >= bp[-1]:
            status, region_rate = "after map end", np.nan
        elif r.start < bp[0] or r.end > bp[-1]:
            status, region_rate = "partly outside map", np.nan
        else:
            status = "inside map"
            region_rate = (np.interp(r.end, bp, cm) - np.interp(r.start, bp, cm)) / size_mb
        s = segs[r.chrom]
        inside = s[(s["start_bp"] >= r.start) & (s["end_bp"] <= r.end)]
        hapibd_rate = float((inside["cM"] / ((inside["end_bp"] - inside["start_bp"]) / 1e6)).median()) \
            if len(inside) else np.nan
        rows.append({
            "chrom": r.chrom,
            "start_Mb": round(r.start / 1e6, 1),
            "end_Mb": round(r.end / 1e6, 1),
            "size_Mb": round(size_mb, 1),
            "peak_segments": int(r.peak_segments),
            "fold_over_median": round(float(r.fold_over_median), 1),
            "map_status": status,
            "map_cM_per_Mb": round(region_rate, 2),
            "hapibd_cM_per_Mb": round(hapibd_rate, 2),
            "chrom_cM_per_Mb": round(chrom_rate, 2),
            "hapibd_rate_ratio": round(hapibd_rate / chrom_rate, 2),
        })

    table = pd.DataFrame(rows)
    table.to_csv(IBD_DIR / "ibd_peaks_masked_table.csv", index=False)

    fmt = lambda x: "--" if pd.isna(x) else f"{x:.2f}"
    lines = [r"\begin{tabular}{lrrrrrlrrr}", r"\toprule",
             r"Chr & Start (Mb) & End (Mb) & Size (Mb) & Peak segments & Fold over median & Map & "
             r"cM/Mb (map) & cM/Mb (hap-ibd) & cM/Mb (chromosome) \\", r"\midrule"]
    for r in table.itertuples(index=False):
        lines.append(f"{r.chrom.removeprefix('chr')} & {r.start_Mb:.1f} & {r.end_Mb:.1f} & {r.size_Mb:.1f} & "
                     f"{r.peak_segments:,} & {r.fold_over_median:,.1f} & {r.map_status} & "
                     f"{fmt(r.map_cM_per_Mb)} & {fmt(r.hapibd_cM_per_Mb)} & {r.chrom_cM_per_Mb:.2f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (IBD_DIR / "ibd_peaks_masked_table.tex").write_text("\n".join(lines) + "\n")

    print(table.to_string(index=False))
    print(f"\ntotal masked: {table['size_Mb'].sum():.1f} Mb in {len(table)} regions")
    print(f"wrote ibd_peaks_masked_table.csv and .tex to {IBD_DIR}")


if __name__ == "__main__":
    main()
