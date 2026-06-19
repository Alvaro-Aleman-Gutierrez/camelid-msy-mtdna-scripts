#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# 06_align_mafft.sh — Multiple sequence alignment with MAFFT
#
# Aligns the multi-FASTA (44 ingroup + 1 outgroup) produced by
# 04_filter_consensus.sh. Uses --auto mode (MAFFT selects the best strategy)
# with --inputorder to preserve the original sequence order.
#
# Tools: MAFFT v7.427
# HPC modules: MAFFT/7.427-foss-2018b-with-extensions
#
# Usage: bash 06_align_mafft.sh
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
INPUT_FASTA="path/to/csa_male_w_outgroup.fasta"   # multi-FASTA with outgroup
OUTPUT_ALN="path/to/csa_male_w_outgroup_aln.fasta"
THREADS=60
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# module load MAFFT/7.427-foss-2018b-with-extensions
# conda activate mafft

echo "Starting MAFFT alignment..."
echo "  Input : $INPUT_FASTA"
echo "  Output: $OUTPUT_ALN"
echo "  Threads: $THREADS"

mafft --auto --inputorder --thread "$THREADS" "$INPUT_FASTA" > "$OUTPUT_ALN"

echo "Alignment complete: $OUTPUT_ALN"
