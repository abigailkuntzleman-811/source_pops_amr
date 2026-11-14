#!/usr/bin/env python3
"""Filter the published hmmix archaic introgression calls.

Reads the segment file from Skov (2024), doi:10.5281/zenodo.14136628, retains
segments whose mean posterior probability meets the threshold, and writes three
files: Neanderthal segments, Denisovan segments, and all retained segments.

The threshold (default 0.80) follows the recommendation in the Zenodo record.
Change it in config.py or by setting HMMIX_MIN_PROB.

Usage:
    python fig1/01_filter_hmmix_calls.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import HMMIX_DIR, HMMIX_MIN_PROB, DATA_PROCESSED  # noqa: E402

SEGMENTS_FILE = HMMIX_DIR / "hg38_1000g_segments.txt"
OUT_DIR = DATA_PROCESSED / "hmmix_filtered"

# ND_type value in the source file -> output filename
ARCHAIC_TYPES = {
    "Neanderthal": "hmmix_neanderthal_filtered.csv",
    "Denisova": "hmmix_denisovan_filtered.csv",
}


def main() -> None:
    if not SEGMENTS_FILE.exists():
        sys.exit(
            f"Missing input: {SEGMENTS_FILE}\n"
            "Run fig1/00_download_data.sh first."
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"reading {SEGMENTS_FILE}")
    hmmix = pd.read_csv(SEGMENTS_FILE, sep=r"\s+")
    print(f"  {len(hmmix):,} segments, {hmmix['name'].nunique():,} individuals")

    passing = hmmix[hmmix["mean_prob"] >= HMMIX_MIN_PROB]
    print(f"  {len(passing):,} segments with mean_prob >= {HMMIX_MIN_PROB}")

    # Per-archaic-source files
    for nd_type, filename in ARCHAIC_TYPES.items():
        subset = passing[passing["ND_type"] == nd_type]
        out = OUT_DIR / filename
        subset.to_csv(out, index=False)
        print(f"  wrote {len(subset):>9,} {nd_type:<12} -> {out.name}")

    # Combined file, all archaic sources
    out = OUT_DIR / "hmmix_filtered.csv"
    passing.to_csv(out, index=False)
    print(f"  wrote {len(passing):>9,} {'all':<12} -> {out.name}")

    # Record the sample set actually analysed, for the methods section.
    samples = sorted(passing["name"].unique())
    sample_file = OUT_DIR / "analysed_samples.txt"
    sample_file.write_text("\n".join(samples) + "\n")
    print(f"  wrote {len(samples):>9,} {'samples':<12} -> {sample_file.name}")


if __name__ == "__main__":
    main()
