#!/usr/bin/env bash
#
# Central path configuration for this project.
#
# Everything downstream reads these variables, so the pipeline can be run
# anywhere by editing this one file (or by exporting the variables before
# calling the scripts). No analysis script should contain an absolute path.
#
# Usage in a shell script:   source "$(dirname "$0")/../config.sh"
# Python scripts read the same variables via config.py.

# Repository root, resolved relative to this file.
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT

# Machine-specific settings (paths, virtual environment) live in
# config.local.sh, which is not tracked by git. Copy config.local.sh.example
# and edit it. Sourced first so its values take precedence over the defaults
# below, which all use ${VAR:-default}.
if [[ -f "$PROJECT_ROOT/config.local.sh" ]]; then
  source "$PROJECT_ROOT/config.local.sh"
fi

# --- Inputs -----------------------------------------------------------------
# Downloaded public data (see code/00_download_data.sh).
export DATA_RAW="${DATA_RAW:-$PROJECT_ROOT/data/raw}"

# Skov (2024) hmmix archaic introgression calls.
export HMMIX_DIR="${HMMIX_DIR:-$DATA_RAW/skov_archaic}"

# Local ancestry inference BED files.
# NOTE: not publicly available. Set this to wherever you have them.
export LAI_DIR="${LAI_DIR:-$DATA_RAW/local_ancestry}"

# --- Outputs ----------------------------------------------------------------
export DATA_PROCESSED="${DATA_PROCESSED:-$PROJECT_ROOT/data/processed}"
export RESULTS_DIR="${RESULTS_DIR:-$PROJECT_ROOT/results}"
export FIGURE_DIR="${FIGURE_DIR:-$RESULTS_DIR/figures}"
export TABLE_DIR="${TABLE_DIR:-$RESULTS_DIR/tables}"

# Per-population hmmix BED files, input to the pairwise sharing scripts.
export BY_POP_DIR="${BY_POP_DIR:-$DATA_PROCESSED/hmmix_by_pop}"
export PAIRWISE_DIR="${PAIRWISE_DIR:-$DATA_PROCESSED/pairwise_comparisons}"

# --- Analysis parameters ----------------------------------------------------
# Minimum hmmix mean posterior probability for a segment to be retained.
export HMMIX_MIN_PROB="${HMMIX_MIN_PROB:-0.80}"

mkdir -p "$DATA_RAW" "$DATA_PROCESSED" "$FIGURE_DIR" "$TABLE_DIR"
