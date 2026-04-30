#!/bin/bash
# Process a single RNA-seq sample
#
# Usage: ./process_sample.sh <sample_id> <fastq_r1_url> <fastq_r2_url>

set -e  # Exit on error

SAMPLE_ID=$1
FASTQ_R1_URL=$2
FASTQ_R2_URL=$3

echo "Processing sample: $SAMPLE_ID"

# Create output directories
mkdir -p data/fastq results/fastqc results/fastp results/salmon data/salmon_index

# Step 1: Download FASTQ files
# TODO: Add curl commands to download R1 and R2 files

# Step 2: Run FastQC
# TODO: Add fastqc command

# Step 3: Run fastp
# TODO: Add fastp command

# Step 4: Download salmon index (if needed)
# TODO: Add conditional download of salmon index

# Step 5: Run Salmon quantification
# TODO: Add salmon quant command

echo "Completed: $SAMPLE_ID"
