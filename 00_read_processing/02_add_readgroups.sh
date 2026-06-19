#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# 02_add_readgroups.sh — Add read group tags to BAM files
#
# Adds @RG tags (ID, LB, PL, PU, SM) required by GATK HaplotypeCaller.
# Skips files that already contain read groups.
#
# Tools: SAMtools v1.10, GNU Parallel
# HPC modules: SAMtools/1.10-GCC-8.3.0
#
# Usage: bash 02_add_readgroups.sh
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
BAM_DIR="path/to/sorted_bams"            # output of 01_mapping.sh
OUTPUT_DIR="path/to/bams_with_rg"
NUM_THREADS=60
THREADS_PER_JOB=4
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# module load SAMtools/1.10-GCC-8.3.0

mkdir -p "$OUTPUT_DIR"

process_bam() {
    local bam=$1
    local base=$(basename "$bam" .bam)
    local out="${OUTPUT_DIR}/${base}_rg.bam"

    if [ -f "$out" ]; then
        echo "  [SKIP] $base — already processed"
        return
    fi

    if samtools view -H "$bam" | grep -q '@RG'; then
        echo "  [COPY] $base — read groups already present"
        cp "$bam" "$out"
    else
        echo "  [RG]   $base — adding read groups..."
        samtools addreplacerg \
            -@ "$THREADS_PER_JOB" \
            -r "ID:${base}" \
            -r "LB:lib1" \
            -r "PL:ILLUMINA" \
            -r "PU:unit1" \
            -r "SM:${base}" \
            -o "$out" \
            "$bam"
    fi

    samtools index -@ "$THREADS_PER_JOB" "$out"
    echo "  Done: $base"
}
export -f process_bam
export OUTPUT_DIR THREADS_PER_JOB

find "$BAM_DIR" -name "*_sorted.bam" | \
    parallel -j $((NUM_THREADS / THREADS_PER_JOB)) process_bam {}

echo "Done. BAMs with read groups in: $OUTPUT_DIR"
