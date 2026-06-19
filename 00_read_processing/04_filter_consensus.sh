#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# 04_filter_consensus.sh — SNP filtering and consensus FASTA generation
#
# For each individual VCF:
#   1. Retain only PASS variants (bcftools view -f PASS,.)
#   2. Keep only biallelic SNPs (bcftools view -v snps)
#   3. Compress and index with bgzip/tabix
#   4. Generate per-individual consensus FASTA (bcftools consensus)
#   5. Concatenate all individuals into a multi-FASTA
#
# Tools: BCFtools v1.10.2, SAMtools v1.10
# HPC modules: BCFtools/1.10.2-GCC-8.3.0, SAMtools/1.10-GCC-8.3.0
#
# Usage: bash 04_filter_consensus.sh
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
VCF_DIR="path/to/vcf_filtered"           # output of 03_variant_calling.sh
OUTPUT_DIR="path/to/fasta_consensus"
REFERENCE="path/to/reference.fasta"      # same reference used in mapping
MULTIFASTA_OUT="${OUTPUT_DIR}/csa_male.fasta"   # final concatenated FASTA
NUM_THREADS=60
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# module load BCFtools/1.10.2-GCC-8.3.0
# module load SAMtools/1.10-GCC-8.3.0

mkdir -p "$OUTPUT_DIR"
> "$MULTIFASTA_OUT"

for VCF in "$VCF_DIR"/*_gatk_variantes_filtered.vcf; do
    BASE=$(basename "$VCF" .vcf)
    SEQ_NAME=$(echo "$BASE" | cut -d'_' -f1)   # use sample ID as FASTA header
    FASTA_OUT="${OUTPUT_DIR}/${BASE}_consensus.fasta"

    echo "Processing $SEQ_NAME..."

    # 1. Retain PASS variants
    PASS_VCF="${OUTPUT_DIR}/${BASE}_pass.vcf"
    bcftools view -f 'PASS,.' "$VCF" -O v --threads "$NUM_THREADS" -o "$PASS_VCF"

    # 2. Keep SNPs only
    SNP_VCF="${OUTPUT_DIR}/${BASE}_snps.vcf"
    bcftools view -v snps "$PASS_VCF" -O v -o "$SNP_VCF"

    # 3. Compress and index
    bgzip -f "$SNP_VCF"
    tabix -p vcf "${SNP_VCF}.gz"

    # 4. Consensus FASTA
    bcftools consensus -f "$REFERENCE" "${SNP_VCF}.gz" > "$FASTA_OUT"

    # 5. Rename header to sample ID
    sed -i "s/^>.*$/>$SEQ_NAME/" "$FASTA_OUT"

    cat "$FASTA_OUT" >> "$MULTIFASTA_OUT"
    echo "  Done: $SEQ_NAME"
done

echo "Multi-FASTA written: $MULTIFASTA_OUT"
echo "Add outgroup sequence manually before alignment:"
echo "  cat outgroup.fasta >> $MULTIFASTA_OUT"
