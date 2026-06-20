#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# IQ-TREE 3 — Maximum-likelihood phylogeny of the complete mitochondrial genome
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Input  : Multiple-sequence alignment of complete mtDNA (16,643 bp)
#          (45 sequences including C. dromedarius outgroup)
# Output : ML tree + consensus tree + bootstrap support
# Model  : TPM2u+F+R3 (selected by ModelFinder under BIC)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
INPUT_FASTA="path/to/csa_male_mt_w_outgroup.fasta"    # 45-seq mtDNA alignment
RESULTS_DIR="results/iqtree_mtdna"
BASE_NAME="csa_male_mtdna26"
OUTGROUP="NC_009849.1_C_dromedarius_mt"   # C. dromedarius mtDNA accession

BOOTSTRAP_REPLICAS=1000
RUNS=10
THREADS=30
MEM=80GB
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "$RESULTS_DIR"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: FASTA not found: $INPUT_FASTA"
    exit 1
fi

echo "Input : $INPUT_FASTA  ($(du -sh "$INPUT_FASTA" | cut -f1))"
echo "Output: $RESULTS_DIR/$BASE_NAME.*"

if ! command -v iqtree3 &>/dev/null; then
    echo "ERROR: iqtree3 not found. Install with: conda install -c bioconda iqtree"
    exit 1
fi
echo "Using: $(which iqtree3)  $(iqtree3 --version 2>&1 | head -1)"

echo "Starting IQ-TREE 3 — $(date)"

iqtree3 \
    -s "$INPUT_FASTA" \
    -m MFP \
    -o "$OUTGROUP" \
    -B "$BOOTSTRAP_REPLICAS" \
    --alrt "$BOOTSTRAP_REPLICAS" \
    -T "$THREADS" \
    --runs "$RUNS" \
    --mem "$MEM" \
    -pre "${RESULTS_DIR}/${BASE_NAME}" \
    -bnni \
    --redo

echo "IQ-TREE finished — $(date)"
echo "Consensus tree: ${RESULTS_DIR}/${BASE_NAME}.contree"
ls -lh "${RESULTS_DIR}/${BASE_NAME}".*
