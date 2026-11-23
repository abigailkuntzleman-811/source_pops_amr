#!/usr/bin/env python3
"""Mask IBD pileups in the hap-ibd calls.

Segment handling follows Guo et al. (2024):
    Guo B, et al. (2024) Strong positive selection biases identity-by-descent-
    based inferences of recent demography and population structure in
    Plasmodium falciparum. Nature Communications 15:2499.
    doi:10.1038/s41467-024-46659-0

Steps
  1. IBD coverage: for each autosome, the number of segments (all samples, all
     pairs) overlapping each --bin bp window.
  2. Peaks: bins with coverage more than --fold times the genome-wide median
     (median over bins with any segments, so centromere gaps are ignored).
     Neighbouring flagged bins are merged into one region.

     This differs from Guo et al., who flagged candidates above the trimmed
     mean + 2 trimmed SD and extended them to the chromosome median, then kept
     only candidates with a signal of positive selection. Without that second
     filter, the 2-SD rule flags ordinary variation in IBD rate: in testing
     with a single planted pileup it masked about a quarter of the genome.
     The pileups here are technical artifacts tens to thousands of times
     above background, so a fold-over-median threshold isolates them.
  3. Segments (as in Guo et al.): a segment entirely inside a peak is removed.
     A segment partly overlapping a peak has the overlapping part cut out, and
     each remaining piece is kept if it is at least --min-cm long. Piece
     lengths in cM come from the genetic map used to call IBD; if it cannot be
     found, they are scaled from the segment's own cM per bp, and the script
     says so.

Masked segments are written next to the originals, in the same hap-ibd format:
    $IBD_DIR/kgp_chr<N>_segments.masked.ibd.gz
The figure scripts read these by default (see IBD_MASKED in config.py).

Usage:
    python hapIBD_calls/03_mask_pileups.py
    python hapIBD_calls/03_mask_pileups.py --map-template /path/plink.chr{CHROM}.GRCh38.map

Also writes:
    $IBD_DIR/ibd_peaks_masked.bed           masked regions (chrom, start, end, peak, fold)
    $IBD_DIR/ibd_peaks_masking_summary.csv  per chromosome: segments and cM before/after
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import CHROMOSOMES, GENETIC_MAP_DIR, IBD_DIR  # noqa: E402

sys.path.insert(0, str(_root / "figure_code" / "fig3"))
IBD_COLUMNS = __import__("01_within_population_ibd").IBD_COLUMNS


def raw_path(n):
    return IBD_DIR / f"kgp_chr{n}_segments.ibd.gz"


def masked_path(n):
    return IBD_DIR / f"kgp_chr{n}_segments.masked.ibd.gz"


# ---------------------------------------------------------------------------
# coverage and peaks
# ---------------------------------------------------------------------------
def coverage(seg: pd.DataFrame, bin_bp: int) -> np.ndarray:
    first = (seg["start_bp"].to_numpy() // bin_bp).astype(int)
    last = ((seg["end_bp"].to_numpy() - 1) // bin_bp).astype(int)
    diff = np.zeros(last.max() + 2, dtype=np.int64)
    np.add.at(diff, first, 1)
    np.add.at(diff, last + 1, -1)
    return np.cumsum(diff)[:-1]


def find_peaks(cov: np.ndarray, threshold: float):
    """Runs of consecutive bins with coverage above threshold, as [start, end) bin indices."""
    peaks, start = [], None
    for b, above in enumerate(list(cov > threshold) + [False]):
        if above and start is None:
            start = b
        elif not above and start is not None:
            peaks.append((start, b))
            start = None
    return peaks


# ---------------------------------------------------------------------------
# genetic map
# ---------------------------------------------------------------------------
def map_from_log(n):
    """The map file hap-ibd recorded in its log for this chromosome, if any."""
    log = IBD_DIR / f"kgp_chr{n}_segments.log"
    if not log.exists():
        return None
    m = re.search(r"\bmap\s*[=:]\s*(\S+)", log.read_text(errors="ignore"))
    return Path(m.group(1)) if m else None


def load_map(n, template):
    """(bp, cM) arrays for one chromosome, or None if no map can be found."""
    candidates = []
    if template:
        candidates.append(Path(template.replace("{CHROM}", str(n))))
    candidates.append(map_from_log(n))
    candidates.append(GENETIC_MAP_DIR / f"plink.chr{n}.GRCh38.map")
    for path in candidates:
        if path is not None and path.exists():
            m = pd.read_csv(path, sep=r"\s+", header=None, usecols=[2, 3], names=["cM", "bp"])
            m = m.sort_values("bp")
            return path, m["bp"].to_numpy(float), m["cM"].to_numpy(float)
    return None


# ---------------------------------------------------------------------------
# masking
# ---------------------------------------------------------------------------
def mask_segments(seg, peaks_bp, gmap, min_cm):
    """Drop segments inside peaks and trim those crossing them."""
    if not peaks_bp:
        return seg, 0, 0
    s, e = seg["start_bp"].to_numpy(), seg["end_bp"].to_numpy()
    touches = np.zeros(len(seg), bool)
    for ps, pe in peaks_bp:
        touches |= (s < pe) & (e > ps)

    untouched = seg[~touches]
    pieces = []
    n_dropped = n_trimmed = 0
    for row in seg[touches].itertuples(index=False):
        keep = [(row.start_bp, row.end_bp)]
        for ps, pe in peaks_bp:                   # subtract each peak
            nxt = []
            for a, b in keep:
                if b <= ps or a >= pe:
                    nxt.append((a, b))
                    continue
                if a < ps:
                    nxt.append((a, ps))
                if b > pe:
                    nxt.append((pe, b))
            keep = nxt
        kept_any = False
        for a, b in keep:
            if gmap is not None:
                cm = float(np.interp(b, gmap[0], gmap[1]) - np.interp(a, gmap[0], gmap[1]))
            else:
                cm = row.cM * (b - a) / (row.end_bp - row.start_bp)
            if cm >= min_cm:
                pieces.append(row._replace(start_bp=a, end_bp=b, cM=round(cm, 3)))
                kept_any = True
        if kept_any:
            n_trimmed += 1
        else:
            n_dropped += 1
    trimmed = pd.DataFrame(pieces, columns=seg.columns)
    out = pd.concat([untouched, trimmed], ignore_index=True) if len(trimmed) else untouched
    return out.sort_values(["start_bp", "end_bp"]), n_dropped, n_trimmed


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", type=int, default=100_000, help="coverage bin size in bp (default 100 kb)")
    p.add_argument("--fold", type=float, default=3.0,
                   help="flag bins above this multiple of the genome-wide median coverage (default 3)")
    p.add_argument("--min-cm", type=float, default=2.0, help="minimum length of a trimmed piece (default 2 cM)")
    p.add_argument("--map-template", default=None,
                   help="genetic map path with {CHROM}; default: the map in each hap-ibd log, "
                        "then GENETIC_MAP_DIR/plink.chr{CHROM}.GRCh38.map")
    args = p.parse_args()

    # pass 1: coverage for every chromosome
    cov = {}
    for chrom in CHROMOSOMES:
        n = chrom.removeprefix("chr")
        if not raw_path(n).exists():
            sys.exit(f"Missing {raw_path(n)}")
        seg = pd.read_csv(raw_path(n), sep="\t", header=None, names=IBD_COLUMNS,
                          usecols=["start_bp", "end_bp"])
        cov[chrom] = coverage(seg, args.bin)
        print(f"read {chrom}: {len(seg):,} segments", flush=True)

    nonzero = np.concatenate([c[c > 0] for c in cov.values()]).astype(float)
    median = float(np.median(nonzero))
    threshold = args.fold * median
    print(f"\nmedian coverage per {args.bin // 1000} kb bin: {median:,.0f}; "
          f"flagging bins above {threshold:,.0f} ({args.fold:g}x)")

    # pass 2: peaks and masking, chromosome by chromosome
    beds, rows, map_note = [], [], set()
    for chrom in CHROMOSOMES:
        n = chrom.removeprefix("chr")
        c = cov[chrom]
        runs = find_peaks(c, threshold)
        peaks = [(a * args.bin, b * args.bin) for a, b in runs]
        beds += [(chrom, a, b, int(c[i:j].max()), round(c[i:j].max() / median, 1))
                 for (a, b), (i, j) in zip(peaks, runs)]

        seg = pd.read_csv(raw_path(n), sep="\t", header=None, names=IBD_COLUMNS)
        found = load_map(n, args.map_template) if peaks else None
        gmap = (found[1], found[2]) if found else None
        if peaks:
            map_note.add(str(found[0].parent) if found else "none")
        out, n_drop, n_trim = mask_segments(seg, peaks, gmap, args.min_cm)
        out.to_csv(masked_path(n), sep="\t", header=False, index=False, compression="gzip")

        rows.append({"chrom": chrom, "peaks": len(peaks),
                     "masked_Mb": sum(b - a for a, b in peaks) / 1e6,
                     "segments_before": len(seg), "segments_after": len(out),
                     "segments_removed": n_drop, "segments_trimmed": n_trim,
                     "cM_before": round(seg["cM"].sum(), 1), "cM_after": round(out["cM"].sum(), 1)})
        flag = (f"  {len(peaks)} peak(s), {rows[-1]['masked_Mb']:.1f} Mb; "
                f"removed {n_drop:,}, trimmed {n_trim:,}") if peaks else ""
        print(f"{chrom}: {len(seg):,} -> {len(out):,} segments{flag}", flush=True)

    pd.DataFrame(beds).to_csv(IBD_DIR / "ibd_peaks_masked.bed", sep="\t", header=False, index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(IBD_DIR / "ibd_peaks_masking_summary.csv", index=False)

    print(f"\n{len(beds)} masked region(s), {summary['masked_Mb'].sum():.1f} Mb in total:")
    for b in beds:
        print(f"  {b[0]}:{b[1] / 1e6:.1f}-{b[2] / 1e6:.1f} Mb   peak {b[3]:,} segments ({b[4]}x median)")
    if "none" in map_note:
        print("\nnote: no genetic map found for some chromosomes with peaks; trimmed piece "
              "lengths there are scaled from each segment's cM per bp. Pass --map-template "
              "to use the map instead.")
    else:
        print(f"\ntrimmed lengths from the genetic map in: {', '.join(sorted(map_note)) or '(no peaks)'}")
    print(f"wrote kgp_chr*_segments.masked.ibd.gz, ibd_peaks_masked.bed and "
          f"ibd_peaks_masking_summary.csv to {IBD_DIR}")


if __name__ == "__main__":
    main()
