#!/bin/bash
set -e  # Exit on error

echo "Processing sample_01..."

# Create output directory
mkdir -p results/fastqc

# Run FastQC (mock - just echoes for learning purposes)
echo "Running FastQC on sample_01_R1.fastq.gz..."
sleep 1
echo "Running FastQC on sample_01_R2.fastq.gz..."
sleep 1

echo "Done!"
