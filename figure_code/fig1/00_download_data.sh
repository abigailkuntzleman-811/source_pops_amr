#!/usr/bin/env bash
#
# Downloads the public input data for the Figure 1 / Table 1 pipeline.
# Safe to re-run: complete files are skipped, partial downloads resumed.
#
# Usage:
#     bash fig1/00_download_data.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../../config.sh"

mkdir -p "$HMMIX_DIR"

# ---------------------------------------------------------------------------
# Archaic introgression calls (hmmix), GRCh38.
#
# Skov, L. (2024). Archaic introgression for HGDP and 1000genomes in hg38
# (v1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14136628
# Licence: CC BY 4.0
#
# These calls are based on the 1000 Genomes 30x high-coverage callset on
# GRCh38 (3,202 samples; Byrska-Bishop et al. 2022, Cell 185:3426-3440.e19).
# ---------------------------------------------------------------------------

ZENODO="https://zenodo.org/records/14136628/files"

for f in hg38_1000g_segments.txt hg38_1000g_SNPS.txt; do
  if [[ -s "$HMMIX_DIR/$f" ]]; then
    echo "skip  $f"
  else
    echo "fetch $f"
    wget -c --show-progress -O "$HMMIX_DIR/$f" "$ZENODO/${f}?download=1"
  fi
done

# ---------------------------------------------------------------------------
# 1000 Genomes sample metadata (small; useful for population labels).
# ---------------------------------------------------------------------------

META_DIR="$DATA_RAW/metadata"
mkdir -p "$META_DIR"

IGSR="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"

[[ -s "$META_DIR/integrated_call_samples_v3.20130502.ALL.panel" ]] || \
  wget -c --show-progress -P "$META_DIR" \
    "$IGSR/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

[[ -s "$META_DIR/20130606_g1k_3202_samples_ped_population.txt" ]] || \
  wget -c --show-progress -P "$META_DIR" \
    "$IGSR/data_collections/1000G_2504_high_coverage/working/20130606_g1k_3202_samples_ped_population.txt"

# ---------------------------------------------------------------------------
# Checksums, so others can confirm identical inputs.
# ---------------------------------------------------------------------------

( cd "$DATA_RAW" && find . -type f ! -name 'checksums.sha256' -print0 \
    | sort -z | xargs -0 sha256sum > checksums.sha256 )

echo
echo "done. verify later with:  (cd $DATA_RAW && sha256sum -c checksums.sha256)"
echo
echo "NOTE: local ancestry inference files are not public and are not"
echo "downloaded here. Set LAI_DIR in config.sh if you have them."
