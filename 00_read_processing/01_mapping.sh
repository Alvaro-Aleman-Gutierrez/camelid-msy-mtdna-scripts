#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# 01_mapping.sh — Read mapping with BWA-MEM v0.7.17
#
# Maps paired-end FASTQ reads to a reference genome (MSY scaffold or mtDNA).
# Output: sorted, indexed BAM files per individual.
#
# Tools: BWA v0.7.17, SAMtools v1.10, GNU Parallel
# HPC modules: BWA/0.7.17-foss-2018b, SAMtools/1.10-GCC-8.3.0
#
# Usage: bash 01_mapping.sh
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
READS_DIR="path/to/fastq_reads"          # directory with *_R1.fastq.gz files
REFERENCE="path/to/reference.fasta"      # MSY: NC_087473.1 | mtDNA: NC_009849.1
OUTPUT_DIR="path/to/bam_output"
NCORES=80                                # total CPU cores available
CORES_PER_TASK=20                        # threads per BWA-MEM job
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# Load modules (HPC) or activate conda environment
# module load BWA/0.7.17-foss-2018b
# module load SAMtools/1.10-GCC-8.3.0
# conda activate bwa_env

mkdir -p "$OUTPUT_DIR"

# Index reference if not already indexed
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Indexing reference..."
    bwa index "$REFERENCE"
fi

# ── Chunk 1: Mapping ──────────────────────────────────────────────────────────
map_reads() {
    local read1=$1
    local read2="${read1/_R1.fastq.gz/_R2.fastq.gz}"
    local base=$(basename "$read1" _R1.fastq.gz)
    local sam_file="${OUTPUT_DIR}/${base}.sam"

    if [ -f "$sam_file" ]; then
        echo "  [SKIP] $base — SAM already exists"
        return
    fi
    echo "  Mapping $base..."
    zcat "$read1" | bwa mem -t "$CORES_PER_TASK" "$REFERENCE" - "$read2" > "$sam_file"
}
export -f map_reads
export OUTPUT_DIR REFERENCE CORES_PER_TASK

echo "Step 1: Mapping reads..."
find "$READS_DIR" -name "*_R1.fastq.gz" | \
    parallel -j $((NCORES / CORES_PER_TASK)) map_reads {}

# ── Chunk 2: SAM → BAM ───────────────────────────────────────────────────────
convert_sam_to_bam() {
    local sam=$1
    local bam="${sam%.sam}.bam"
    echo "  Converting $(basename $sam)..."
    samtools view -Sb "$sam" > "$bam"
}
export -f convert_sam_to_bam

echo "Step 2: Converting SAM to BAM..."
find "$OUTPUT_DIR" -name "*.sam" | parallel -j $((NCORES / CORES_PER_TASK)) convert_sam_to_bam {}

# ── Chunk 3: Sort and index BAM ───────────────────────────────────────────────
sort_bam() {
    local bam=$1
    local base=$(basename "$bam" .bam)
    local sorted="${OUTPUT_DIR}/${base}_sorted.bam"
    echo "  Sorting $(basename $bam)..."
    samtools sort -@ "$CORES_PER_TASK" -o "$sorted" "$bam"
    samtools index "$sorted"
}
export -f sort_bam

echo "Step 3: Sorting and indexing BAM files..."
find "$OUTPUT_DIR" -name "*.bam" ! -name "*sorted*" | \
    parallel -j $((NCORES / CORES_PER_TASK)) sort_bam {}

echo "Done. Sorted BAM files in: $OUTPUT_DIR"
