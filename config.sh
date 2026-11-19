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

# Joint-caller per-haplotype archaic and local ancestry BED files (Figure 2).
# NOTE: not publicly available. Set these to wherever you have them.
export ARC_BED_DIR="${ARC_BED_DIR:-$DATA_RAW/autosome_arc_bedfiles}"
export LAI_BED_DIR="${LAI_BED_DIR:-$DATA_RAW/autosome_lai_bedfiles}"

# 1000 Genomes sample panel, used to assign individuals to populations.
export POP_PANEL="${POP_PANEL:-$DATA_RAW/20130606_g1k_3202_samples_ped_population.txt}"

# msprime simulations of the admixture models (Figure 6).
export SIM_DIR="${SIM_DIR:-$DATA_RAW/simulations_donor_eur}"

# --- Outputs ----------------------------------------------------------------
export DATA_PROCESSED="${DATA_PROCESSED:-$PROJECT_ROOT/data/processed}"
export RESULTS_DIR="${RESULTS_DIR:-$PROJECT_ROOT/results}"
export FIGURE_DIR="${FIGURE_DIR:-$RESULTS_DIR/figures}"
export TABLE_DIR="${TABLE_DIR:-$RESULTS_DIR/tables}"

# Per-population hmmix BED files, input to the pairwise sharing scripts.
export BY_POP_DIR="${BY_POP_DIR:-$DATA_PROCESSED/hmmix_by_pop}"
export PAIRWISE_DIR="${PAIRWISE_DIR:-$DATA_PROCESSED/pairwise_comparisons}"

# --- IBD ---
# Downloaded 1000 Genomes VCFs (reused by the IBD pipeline).
export VCF_DIR="${VCF_DIR:-$DATA_RAW/1kgp_30x_grch38}"

# hap-ibd executable. ibd/00_setup_resources.sh downloads it here.
export HAPIBD_JAR="${HAPIBD_JAR:-$DATA_RAW/tools/hap-ibd.jar}"

# PLINK-format genetic maps in cM for GRCh38, one per chromosome.
# Available from the Beagle genetic maps page; see ibd/00_setup_resources.sh.
export GENETIC_MAP_DIR="${GENETIC_MAP_DIR:-$DATA_RAW/genetic_maps}"
# Assigned with an if rather than ${VAR:-default}: the default contains
# {CHROM}, and that inner brace closes the expansion early, so the variable
# would be exported already mangled as plink.chr{CHROM.GRCh38.map}.
if [[ -z "${GENETIC_MAP_TEMPLATE:-}" ]]; then
  GENETIC_MAP_TEMPLATE='plink.chr{CHROM}.GRCh38.map'
fi
export GENETIC_MAP_TEMPLATE

# hap-ibd output: kgp_chr<N>_segments.ibd.gz
export IBD_DIR="${IBD_DIR:-$DATA_PROCESSED/ibd}"

# Optional file of sample IDs to exclude as related (one per line, or a
# PLINK-style FID/IID pair per line, as KING writes). Leave unset to fall
# back to the 698 related samples named in the 1000 Genomes pedigree.
export RELATED_EXCLUDE="${RELATED_EXCLUDE:-}"

# Minimum MAF for a variant to enter IBD calling.
export IBD_MIN_MAF="${IBD_MIN_MAF:-0.05}"

# --- Analysis parameters ----------------------------------------------------
# Minimum hmmix mean posterior probability for a segment to be retained.
export HMMIX_MIN_PROB="${HMMIX_MIN_PROB:-0.80}"

mkdir -p "$DATA_RAW" "$DATA_PROCESSED" "$FIGURE_DIR" "$TABLE_DIR"
