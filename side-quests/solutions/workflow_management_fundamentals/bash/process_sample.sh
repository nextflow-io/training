#!/bin/bash
# process_sample.sh - Process a single RNA-seq sample
#
# This is the foundation script that learners build in Section 1.2.
# It runs the complete analysis pipeline for ONE sample.
#
# Usage: ./process_sample.sh <sample_id> <fastq_r1_url> <fastq_r2_url>

set -e

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <sample_id> <fastq_r1_url> <fastq_r2_url>"
    echo "Example: $0 WT_REP1 https://...R1.fastq.gz https://...R2.fastq.gz"
    exit 1
fi

SAMPLE_ID=$1
FASTQ_R1_URL=$2
FASTQ_R2_URL=$3

echo "Processing sample: $SAMPLE_ID"

# Create output directories
mkdir -p data/fastq
mkdir -p results/fastqc
mkdir -p results/fastp
mkdir -p results/salmon

# Step 1: Download FASTQ files (if not already present)
echo "  Downloading FASTQ files..."
if [ ! -f "data/fastq/${SAMPLE_ID}_R1.fastq.gz" ]; then
    curl -sL "$FASTQ_R1_URL" -o "data/fastq/${SAMPLE_ID}_R1.fastq.gz"
fi
if [ ! -f "data/fastq/${SAMPLE_ID}_R2.fastq.gz" ]; then
    curl -sL "$FASTQ_R2_URL" -o "data/fastq/${SAMPLE_ID}_R2.fastq.gz"
fi

# Step 2: Run FastQC (quality control)
echo "  Running FastQC..."
fastqc -q -o results/fastqc \
    "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
    "data/fastq/${SAMPLE_ID}_R2.fastq.gz"

# Step 3: Run fastp (adapter trimming)
echo "  Running fastp..."
fastp \
    -i "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
    -I "data/fastq/${SAMPLE_ID}_R2.fastq.gz" \
    -o "results/fastp/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
    -O "results/fastp/${SAMPLE_ID}_trimmed_R2.fastq.gz" \
    -j "results/fastp/${SAMPLE_ID}.fastp.json" \
    -h "results/fastp/${SAMPLE_ID}.fastp.html" \
    2>/dev/null

# Step 4: Download salmon index (if not present)
if [ ! -d "data/salmon_index/salmon" ]; then
    echo "  Downloading salmon index..."
    mkdir -p data/salmon_index
    curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz \
        -o data/salmon_index/salmon.tar.gz
    tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
    rm data/salmon_index/salmon.tar.gz
fi

# Step 5: Run Salmon (transcript quantification)
echo "  Running Salmon..."
salmon quant \
    --index data/salmon_index/salmon \
    --libType A \
    --mates1 "results/fastp/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
    --mates2 "results/fastp/${SAMPLE_ID}_trimmed_R2.fastq.gz" \
    --output "results/salmon/${SAMPLE_ID}" \
    --threads 2 \
    --quiet

echo "Completed: $SAMPLE_ID"
