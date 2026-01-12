#!/bin/bash
# Process a single RNA-seq sample

set -e  # Exit on error

SAMPLE_ID=$1
FASTQ_R1_URL=$2
FASTQ_R2_URL=$3

echo "Processing sample: $SAMPLE_ID"

# Create output directories
mkdir -p data/fastq results/fastqc results/fastp results/salmon

# Download FASTQ files
echo "  Downloading FASTQ files..."
curl -sL "$FASTQ_R1_URL" -o "data/fastq/${SAMPLE_ID}_R1.fastq.gz"
curl -sL "$FASTQ_R2_URL" -o "data/fastq/${SAMPLE_ID}_R2.fastq.gz"

# Run FastQC
echo "  Running FastQC..."
fastqc -q -o results/fastqc \
    "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
    "data/fastq/${SAMPLE_ID}_R2.fastq.gz"

# Run fastp
echo "  Running fastp..."
fastp \
    -i "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
    -I "data/fastq/${SAMPLE_ID}_R2.fastq.gz" \
    -o "results/fastp/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
    -O "results/fastp/${SAMPLE_ID}_trimmed_R2.fastq.gz" \
    -j "results/fastp/${SAMPLE_ID}.fastp.json" \
    -h "results/fastp/${SAMPLE_ID}.fastp.html" \
    2>/dev/null

# Download salmon index (if not present)
if [ ! -d "data/salmon_index/salmon" ]; then
    echo "  Downloading salmon index..."
    mkdir -p data/salmon_index
    curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz \
        -o data/salmon_index/salmon.tar.gz
    tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
    rm data/salmon_index/salmon.tar.gz
fi

# Run Salmon
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
