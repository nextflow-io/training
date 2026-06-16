#!/bin/bash
# Process all RNA-seq samples sequentially
#
# Usage: ./pipeline_sequential.sh [samples.csv]

set -e

SAMPLES_FILE=${1:-data/samples.csv}

echo "=========================================="
echo "RNA-seq Pipeline (Sequential)"
echo "=========================================="

mkdir -p data/fastq results/fastqc results/fastp results/salmon data/salmon_index

# Download salmon index once (shared by all samples)
if [ ! -d "data/salmon_index/salmon" ]; then
    echo "Downloading salmon index..."
    curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz \
        -o data/salmon_index/salmon.tar.gz
    tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
    rm data/salmon_index/salmon.tar.gz
fi

# Process each sample from the CSV (skip header line)
tail -n +2 "$SAMPLES_FILE" | while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
    echo ""
    echo "Processing: $sample_id"

    # TODO: Add the processing steps for each sample
    # - Download FASTQ files
    # - Run FastQC
    # - Run fastp
    # - Run Salmon

    echo "  Done: $sample_id"
done

echo ""
echo "Pipeline complete!"
