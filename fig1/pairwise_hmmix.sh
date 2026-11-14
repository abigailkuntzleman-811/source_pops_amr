#!/usr/bin/env bash
#
# Between-population pairwise sharing.
#
# For every haplotype in pop1 and every haplotype in pop2, intersects their
# archaic segments and records the total shared bp. Writes one line per pair to
#
#     <output_dir>/<pop1>_{<pop2>}.txt      columns: id1, id2, shared_bp
#
# The brace in the filename is intentional: the downstream parsing in
# 05_table1_values.py matches that pattern.
#
# Usage: pairwise_hmmix.sh <pop1> <pop2> <input_dir> <output_dir>

set -uo pipefail

pop1="${1:?}"
pop2="${2:?}"
input_dir="${3:?}"
output_dir="${4:?}"

dir1="$input_dir/$pop1"
dir2="$input_dir/$pop2"

mkdir -p "$output_dir"

individuals_file="$output_dir/${pop1}_{$pop2}.txt"
: > "$individuals_file"

total_intersection=0
total_combinations=0

# Only sorted- files: bedops requires sorted input, and the averaging step
# counts only pairs where both filenames carry the sorted- prefix.
mapfile -t bed_files1 < <(find "$dir1" -name "sorted-*.bed" | sort)
mapfile -t bed_files2 < <(find "$dir2" -name "sorted-*.bed" | sort)

for bed1 in "${bed_files1[@]}"; do
  for bed2 in "${bed_files2[@]}"; do
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
  echo -e "${pop1}\t${pop2}\t${avg}" \
    > "$output_dir/${pop1}_${pop2}_pairwise_average_sharing.txt"
  echo "$pop1 vs $pop2: $total_combinations pairs, mean ${avg} bp"
else
  echo "$pop1 vs $pop2: no BED files found" >&2
fi
