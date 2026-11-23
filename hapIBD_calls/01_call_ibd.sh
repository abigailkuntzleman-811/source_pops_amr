#!/usr/bin/env bash
#
# Calls IBD segments on one chromosome of the 1000 Genomes 30x phased panel.
#
# Filters the panel to common biallelic SNPs, then runs hap-ibd. Each step
# checks for its output first, so an interrupted run resumes rather than
# repeating work.
#
#   1    download the phased panel VCF, if not already present
#   2    keep variants with MAF > $IBD_MIN_MAF
#   3    drop positions appearing more than once
#   4    keep SNPs only
#   5    run hap-ibd
#
# Usage:
#     bash hapIBD_calls/01_call_ibd.sh <CHROM>
#
# As a Slurm array over the autosomes:
#     sbatch --array=1-22 hapIBD_calls/01_call_ibd.script
#
# Output:
#     $IBD_DIR/kgp_chr<CHROM>_segments.ibd.gz   IBD segments
#     $IBD_DIR/kgp_chr<CHROM>_segments.hbd.gz   homozygosity-by-descent
#     $IBD_DIR/kgp_chr<CHROM>_segments.log

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

CHROM="${1:?usage: 01_call_ibd.sh <CHROM>}"

command -v bcftools >/dev/null || { echo "bcftools not found (try: module load bcftools)" >&2; exit 1; }
command -v java >/dev/null     || { echo "java not found" >&2; exit 1; }
[[ -s "$HAPIBD_JAR" ]]         || { echo "hap-ibd.jar not at $HAPIBD_JAR; run hapIBD_calls/00_setup_resources.sh" >&2; exit 1; }

mkdir -p "$VCF_DIR" "$IBD_DIR"

WORK="$DATA_PROCESSED/ibd_tmp/chr${CHROM}"
mkdir -p "$WORK"

# ---------------------------------------------------------------------------
# 1. Phased panel VCF
# ---------------------------------------------------------------------------
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
PANEL="1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

echo "STEP 1: phased panel"
if [[ -s "$VCF_DIR/$PANEL" ]]; then
  echo "  $PANEL exists"
else
  wget -c --show-progress -P "$VCF_DIR" "$BASE/$PANEL"
  wget -c --show-progress -P "$VCF_DIR" "$BASE/${PANEL}.tbi"
fi
curfile="$VCF_DIR/$PANEL"

# ---------------------------------------------------------------------------
# 2. Common variants only
# ---------------------------------------------------------------------------
echo "STEP 2: MAF > $IBD_MIN_MAF"
target="$WORK/maf_chr${CHROM}.vcf.gz"
if [[ -s "$target" ]]; then
  echo "  $(basename "$target") exists"
else
  bcftools view -i "MAF>${IBD_MIN_MAF}" "$curfile" -O z -o "$target"
  bcftools index -f "$target"
fi
curfile="$target"

# ---------------------------------------------------------------------------
# 3. Drop duplicated positions
#
# hap-ibd expects one record per site; duplicated positions produce
# inconsistent haplotype indices.
# ---------------------------------------------------------------------------
echo "STEP 3: remove duplicate positions"
target="$WORK/maf_nodup_chr${CHROM}.vcf.gz"
if [[ -s "$target" ]]; then
  echo "  $(basename "$target") exists"
else
  bcftools view -T <(bcftools query -f "%CHROM %POS\n" "$curfile" \
                     | sort -n | uniq -c \
                     | awk '{ if ($1 == 1) print $2"\t"$3 }') \
                "$curfile" -O z -o "$target"
  bcftools index -f "$target"
fi
curfile="$target"

# ---------------------------------------------------------------------------
# 4. SNPs only
# ---------------------------------------------------------------------------
echo "STEP 4: SNPs only"
target="$WORK/maf_nodup_snps_chr${CHROM}.vcf.gz"
if [[ -s "$target" ]]; then
  echo "  $(basename "$target") exists"
else
  bcftools view -v snps "$curfile" -O z -o "$target"
  bcftools index -f "$target"
fi
curfile="$target"

# ---------------------------------------------------------------------------
# 5. hap-ibd
#
# The VCF uses "chr1" style contig names while the PLINK maps use "1", so a
# chr-prefixed copy of the map is made here. The original is left untouched.
# ---------------------------------------------------------------------------
echo "STEP 5: hap-ibd"

# Substitute with sed rather than bash parameter expansion: the brace
# escaping in ${VAR/\{CHROM\}/...} is not portable across bash versions
# and silently produced a mangled filename on some systems.
MAP_NAME="$(printf "%s" "$GENETIC_MAP_TEMPLATE" | sed "s/{CHROM}/$CHROM/g")"
MAP_SRC="$GENETIC_MAP_DIR/$MAP_NAME"
[[ -s "$MAP_SRC" ]] || { echo "genetic map not found: $MAP_SRC" >&2; exit 1; }

MAP="$WORK/$(basename "$MAP_SRC" .map).chrprefix.map"
if [[ ! -s "$MAP" ]]; then
  if head -1 "$MAP_SRC" | grep -q '^chr'; then
    cp "$MAP_SRC" "$MAP"
  else
    sed "s/^/chr/" "$MAP_SRC" > "$MAP"
  fi
fi

OUT_PREFIX="$IBD_DIR/kgp_chr${CHROM}_segments"
if [[ -s "${OUT_PREFIX}.ibd.gz" ]]; then
  echo "  ${OUT_PREFIX}.ibd.gz exists, nothing to do"
else
  java -Xmx"${IBD_JAVA_MEM:-16}"g -jar "$HAPIBD_JAR" \
       gt="$curfile" \
       map="$MAP" \
       out="$OUT_PREFIX"
fi

# ---------------------------------------------------------------------------
# Keep the filtered VCF (input to hap-ibd), which KING also uses
# (hapIBD_calls/05_king_unrelated.script), under the name of the originals.
# ---------------------------------------------------------------------------
mkdir -p "$FILTERED_VCF_DIR"
KEPT="$FILTERED_VCF_DIR/filtered_1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
if [[ ! -s "$KEPT" ]]; then
  cp "$curfile" "$KEPT"
  bcftools index -f "$KEPT"
fi

# ---------------------------------------------------------------------------
# Clean up intermediates. Set KEEP_TMP=1 to retain them for debugging.
# ---------------------------------------------------------------------------
if [[ "${KEEP_TMP:-0}" != "1" ]]; then
  rm -rf "$WORK"
fi

echo "done: ${OUT_PREFIX}.ibd.gz"
