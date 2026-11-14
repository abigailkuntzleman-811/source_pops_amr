#!/usr/bin/env bash
#
# Runs all between- and within-population pairwise intersections of archaic
# segments, for one archaic source. Produces the sharing values that become
# Figure 1 and the shared-ancestry column of Table 1.
#
# Replaces pairwise_all.script. Two changes from the original:
#   * paths come from config.sh instead of being hardcoded
#   * the output directory passed to the helper scripts was $rootoutdir, which
#     was never assigned -- the helpers received an empty string. It now passes
#     the directory the script actually defines.
#
# Usage:
#     bash fig1/03_pairwise_sharing.sh <archaic>        # all | neanderthal | denisovan
#
# As a Slurm array job (190 tasks = 171 between + 19 within):
#     sbatch --array=1-190 fig1/03_pairwise_sharing.sh neanderthal
#
#SBATCH -J pairwise_intersections
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -o logs/pairwise_intersections_%A_%a.out
#SBATCH -e logs/pairwise_intersections_%A_%a.err

set -euo pipefail

# Locate the repository root. When submitted with sbatch, BASH_SOURCE points
# at a copy of this script in the Slurm spool directory, so the path relative
# to it does not exist; fall back to the submission directory.
find_repo_root() {
  local candidate
  for candidate in "$(dirname "${BASH_SOURCE[0]}")/.." \
                   "${SLURM_SUBMIT_DIR:-}" "${SLURM_SUBMIT_DIR:-}/.." \
                   "$PWD" "$PWD/.."; do
    if [[ -n "$candidate" && -f "$candidate/config.sh" ]]; then
      (cd "$candidate" && pwd)
      return 0
    fi
  done
  echo "ERROR: could not find config.sh. Submit from the repository root," >&2
  echo "       or set REPO_ROOT to the repository path." >&2
  exit 1
}
REPO_ROOT="${REPO_ROOT:-$(find_repo_root)}"
source "$REPO_ROOT/config.sh"

ARCHAIC="${1:?usage: 03_pairwise_sharing.sh <all|neanderthal|denisovan>}"

INPUT_DIR="$BY_POP_DIR/$ARCHAIC"
OUTPUT_DIR="$PAIRWISE_DIR/$ARCHAIC"
mkdir -p "$OUTPUT_DIR"

command -v bedops >/dev/null 2>&1 || { echo "bedops not found (try: ml bedops)" >&2; exit 1; }

POPULATIONS=(BEB CDX CEU CHB CHS CLM FIN GBR GIH IBS ITU JPT KHV MXL PEL PJL PUR STU TSI)

# Between-population pairs, in the same order as the original script.
pairwise_combos=()
within_combos=()
for ((i = 0; i < ${#POPULATIONS[@]}; i++)); do
  for ((j = i + 1; j < ${#POPULATIONS[@]}; j++)); do
    pairwise_combos+=("${POPULATIONS[$i]} ${POPULATIONS[$j]}")
  done
  within_combos+=("${POPULATIONS[$i]}")
done

num_between=${#pairwise_combos[@]}   # 171

run_task() {
  local task_id="$1"
  if (( task_id <= num_between )); then
    local pair="${pairwise_combos[$((task_id - 1))]}"
    echo "BETWEEN: $pair"
    bash "$REPO_ROOT/fig1/pairwise_hmmix.sh" $pair "$INPUT_DIR" "$OUTPUT_DIR"
  else
    local pop="${within_combos[$((task_id - num_between - 1))]}"
    echo "WITHIN: $pop"
    bash "$REPO_ROOT/fig1/pairwise_within_hmmix.sh" "$pop" "$INPUT_DIR" "$OUTPUT_DIR"
  fi
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  run_task "$SLURM_ARRAY_TASK_ID"
else
  # Run everything serially. Slow; prefer the array job on a cluster.
  total=$((num_between + ${#within_combos[@]}))
  echo "running $total comparisons serially"
  for ((t = 1; t <= total; t++)); do
    run_task "$t"
  done
fi
