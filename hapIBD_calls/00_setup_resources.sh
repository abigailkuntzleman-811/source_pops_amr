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
# hap-ibd needs PLINK-format genetic maps in cM, with chromosome identifiers
# matching the VCF. These are the HapMap maps distributed alongside Beagle;
# the GRCh38 archive is linked from the hap-ibd README.
# ---------------------------------------------------------------------------

mkdir -p "$GENETIC_MAP_DIR"

MAP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
MAP_ZIP="$GENETIC_MAP_DIR/$(basename "$MAP_URL")"

first_map="$(printf "%s" "$GENETIC_MAP_TEMPLATE" | sed "s/{CHROM}/1/g")"
if [[ -s "$GENETIC_MAP_DIR/$first_map" ]]; then
  echo "skip  genetic maps already present"
else
  echo "fetch $(basename "$MAP_URL")"
  wget -c --show-progress -O "$MAP_ZIP" "$MAP_URL"
  command -v unzip >/dev/null \
    || { echo "unzip not found; unpack $MAP_ZIP manually" >&2; exit 1; }
  unzip -o -d "$GENETIC_MAP_DIR" "$MAP_ZIP"
  rm -f "$MAP_ZIP"

  # The archive unpacks into chr_in_chrom_field/ and no_chr_in_chrom_field/.
  # 01_call_ibd.sh adds a chr prefix when the map lacks one, so either works;
  # lift one set up to $GENETIC_MAP_DIR so the default path finds them.
  for sub in no_chr_in_chrom_field chr_in_chrom_field; do
    if [[ -d "$GENETIC_MAP_DIR/$sub" ]]; then
      echo "unpacked into $sub/; linking those maps into $GENETIC_MAP_DIR"
      for f in "$GENETIC_MAP_DIR/$sub"/*.map; do
        [[ -e "$f" ]] && ln -sf "$sub/$(basename "$f")" \
                               "$GENETIC_MAP_DIR/$(basename "$f")"
      done
      break
    fi
  done
fi

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
