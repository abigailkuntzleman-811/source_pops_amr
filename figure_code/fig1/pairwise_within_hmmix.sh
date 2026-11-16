#!/usr/bin/env bash
#
# Within-population pairwise sharing.
#
# Intersects every unordered pair of distinct haplotypes within one population
# and records the total shared bp. Writes one line per pair to
#
#     <output_dir>/all_<pop>.txt           columns: id1, id2, shared_bp
#
# These within-population values are the diagonal of the Figure 1 heatmap and
# the "shared Neanderthal ancestry" column of Table 1.
#
# Usage: pairwise_within_hmmix.sh <pop> <input_dir> <output_dir>

set -uo pipefail

pop="${1:?}"
input_dir="${2:?}"
output_dir="${3:?}"

dir="$input_dir/$pop"
mkdir -p "$output_dir"

individuals_file="$output_dir/all_${pop}.txt"
: > "$individuals_file"

total_intersection=0
total_combinations=0

mapfile -t bed_files < <(find "$dir" -name "sorted-*.bed" | sort)

n=${#bed_files[@]}
for ((a = 0; a < n; a++)); do
  for ((b = a + 1; b < n; b++)); do
    bed1="${bed_files[$a]}"
    bed2="${bed_files[$b]}"

    id1=$(basename "$bed1" .bed)
    id2=$(basename "$bed2" .bed)

    intersect_bp=$(bedops --intersect "$bed1" "$bed2" \
                   | awk '{sum += ($3 - $2)} END {print sum + 0}')

    printf '%s\t%s\t%s\n' "$id1" "$id2" "$intersect_bp" >> "$individuals_file"

    total_intersection=$((total_intersection + intersect_bp))
    total_combinations=$((total_combinations + 1))
  done
done

if ((total_combinations > 0)); then
  avg=$(echo "scale=6; $total_intersection / $total_combinations" | bc)
  echo -e "${pop}\t${avg}" > "$output_dir/${pop}_pairwise_average_sharing.txt"
  echo "$pop within: $total_combinations pairs, mean ${avg} bp"
else
  echo "$pop within: no BED files found" >&2
fi
