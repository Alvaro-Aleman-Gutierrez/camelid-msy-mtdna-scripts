#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# 05_merge_vcf.sh — Merge per-individual VCFs into a population-level VCF
#
# Produces the combined VCF used for FST and PCA analyses.
# Requires all individual VCFs to be bgzip-compressed and tabix-indexed.
#
# Tools: BCFtools v1.10.2
# HPC modules: BCFtools/1.10.2-GCC-8.3.0
#
# Usage: bash 05_merge_vcf.sh
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
VCF_DIR="path/to/vcf_compressed"         # directory with *.vcf.gz + *.tbi files
OUTPUT_DIR="path/to/merged_vcf"
MERGED_OUT="${OUTPUT_DIR}/combined.vcf.gz"
NUM_THREADS=20
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# module load BCFtools/1.10.2-GCC-8.3.0

mkdir -p "$OUTPUT_DIR"

bcftools merge \
    --threads "$NUM_THREADS" \
    "$VCF_DIR"/*.vcf.gz \
    -o "$MERGED_OUT"

echo "Merged VCF written: $MERGED_OUT"
