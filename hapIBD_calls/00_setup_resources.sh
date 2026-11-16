#!/usr/bin/env bash
#
# Downloads the tools needed for IBD calling.
#
# Usage:
#     bash hapIBD_calls/00_setup_resources.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

# ---------------------------------------------------------------------------
# hap-ibd
#
# Zhou, Y., Browning, S.R., Browning, B.L. (2020). A fast and simple method
# for detecting identity by descent segments in large-scale data.
# American Journal of Human Genetics, 106(4), 426-437.
# https://doi.org/10.1016/j.ajhg.2020.02.010
# ---------------------------------------------------------------------------

mkdir -p "$(dirname "$HAPIBD_JAR")"

if [[ -s "$HAPIBD_JAR" ]]; then
  echo "skip  hap-ibd.jar already present"
else
  echo "fetch hap-ibd.jar"
  wget -c --show-progress -O "$HAPIBD_JAR" \
    "https://faculty.washington.edu/browning/hap-ibd.jar"
fi

java -jar "$HAPIBD_JAR" 2>&1 | head -2 || {
  echo "ERROR: could not run hap-ibd. Java 1.8 or later is required." >&2
  exit 1
}

# ---------------------------------------------------------------------------
# Genetic maps
#
# hap-ibd needs PLINK-format genetic maps in cM for GRCh38, one per
# chromosome. These are distributed from the Beagle genetic maps page:
#
#     https://faculty.washington.edu/browning/beagle/beagle.html
#
# Download the PLINK-format GRCh38 maps, unpack them into
# $GENETIC_MAP_DIR, and make sure the filenames match
# $GENETIC_MAP_TEMPLATE (default: plink.chr{CHROM}.GRCh38.map).
#
# They are not downloaded automatically here because the archive name has
# changed between releases; pinning a URL that later moves is worse than
# saying plainly where to get them.
# ---------------------------------------------------------------------------

mkdir -p "$GENETIC_MAP_DIR"

missing=0
for CHROM in {1..22}; do
  map_name="$(printf "%s" "$GENETIC_MAP_TEMPLATE" | sed "s/{CHROM}/$CHROM/g")"
  map="$GENETIC_MAP_DIR/$map_name"
  [[ -s "$map" ]] || missing=$((missing + 1))
done

if (( missing > 0 )); then
  echo
  echo "NOTE: $missing of 22 genetic maps are missing from $GENETIC_MAP_DIR"
  echo "      Expected filenames like: $(printf "%s" "$GENETIC_MAP_TEMPLATE" | sed "s/{CHROM}/1/g")"
  echo "      See the comment block in this script for where to get them."
else
  echo "ok    all 22 genetic maps present"
fi
