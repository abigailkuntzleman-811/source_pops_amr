#!/usr/bin/env bash
#
# Driver for 02_segments_to_bed.py: builds per-haplotype BED files for every
# population, for each archaic source.
#
# Replaces organize_hmmix.script. Paths come from config.sh.
#
# Usage:
#     bash fig1/02_segments_to_bed.sh
#
# On a Slurm cluster, submit with the header below or wrap in your own job file.
#SBATCH -J segments_to_bed
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -o logs/segments_to_bed.out
#SBATCH -e logs/segments_to_bed.err

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

FILTERED="$DATA_PROCESSED/hmmix_filtered"

# archaic label -> filtered CSV produced by 01_filter_hmmix_calls.py
declare -A INPUTS=(
  [all]="$FILTERED/hmmix_filtered.csv"
  [neanderthal]="$FILTERED/hmmix_neanderthal_filtered.csv"
  [denisovan]="$FILTERED/hmmix_denisovan_filtered.csv"
)

SUB_POPS=(CHB JPT CHS CDX KHV
          CEU TSI GBR FIN IBS
          GIH PJL BEB STU ITU
          MXL PEL PUR CLM)

for archaic in "${!INPUTS[@]}"; do
  input="${INPUTS[$archaic]}"
  if [[ ! -f "$input" ]]; then
    echo "missing $input -- run 01_filter_hmmix_calls.py first" >&2
    exit 1
  fi

  out_root="$BY_POP_DIR/$archaic"

  for pop in "${SUB_POPS[@]}"; do
    # Start clean so a re-run never mixes old and new files.
    rm -rf "${out_root:?}/${pop}"
    python3 "$REPO_ROOT/fig1/02_segments_to_bed.py" "$pop" "$input" "$out_root"
  done
done

echo "done. per-population BED files under $BY_POP_DIR"
