#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# 03_variant_calling.sh — Variant calling with GATK HaplotypeCaller
#
# Calls SNPs per individual with haploid ploidy (--ploidy 1), appropriate for
# non-recombining markers (MSY scaffold and mtDNA). Applies hard filters with
# GATK VariantFiltration: DP < 7 (low coverage), DP > 30 (excess coverage),
# QUAL < 28 (low base quality).
#
# Tools: GATK v4.1.5.0, SAMtools v1.10, GNU Parallel
# HPC modules: GATK/4.1.5.0-GCCcore-9.3.0-Java-11, SAMtools/1.10-GCC-8.3.0
#
# Usage: bash 03_variant_calling.sh
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
REFERENCE="path/to/reference.fasta"      # MSY: NC_087473.1 | mtDNA: NC_009849.1
BAM_DIR="path/to/bams_with_rg"           # output of 02_add_readgroups.sh
VCF_DIR="path/to/vcf_output"
NUM_THREADS=80
THREADS_PER_JOB=20
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# module load SAMtools/1.10-GCC-8.3.0
# module load GATK/4.1.5.0-GCCcore-9.3.0-Java-11

mkdir -p "$VCF_DIR"

# Index reference and create sequence dictionary if needed
[ ! -f "${REFERENCE}.fai" ] && samtools faidx "$REFERENCE"
DICT="${REFERENCE%.fasta}.dict"
[ ! -f "$DICT" ] && gatk CreateSequenceDictionary -R "$REFERENCE"

process_bam() {
    local bam=$1
    local base=$(basename "$bam" .bam)
    local vcf_raw="${VCF_DIR}/${base}_gatk_variantes.vcf"
    local vcf_filtered="${VCF_DIR}/${base}_gatk_variantes_filtered.vcf"

    if [ -f "$vcf_filtered" ]; then
        echo "  [SKIP] $base"
        return
    fi

    echo "  HaplotypeCaller: $base..."
    gatk HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$bam" \
        -O "$vcf_raw" \
        --ploidy 1 \
        --native-pair-hmm-threads "$THREADS_PER_JOB"

    echo "  VariantFiltration: $base..."
    gatk VariantFiltration \
        -R "$REFERENCE" \
        -V "$vcf_raw" \
        -O "$vcf_filtered" \
        --filter-name "LowDP"   --filter-expression "DP < 7" \
        --filter-name "HighDP"  --filter-expression "DP > 30" \
        --filter-name "LowQual" --filter-expression "QUAL < 28"

    echo "  Done: $base"
}
export -f process_bam
export REFERENCE VCF_DIR THREADS_PER_JOB

find "$BAM_DIR" -name "*_rg.bam" | \
    parallel -j $((NUM_THREADS / THREADS_PER_JOB)) process_bam {}

echo "Done. Filtered VCFs in: $VCF_DIR"
