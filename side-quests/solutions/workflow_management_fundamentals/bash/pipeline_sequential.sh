#!/bin/bash
# pipeline_sequential.sh - Process all RNA-seq samples sequentially
#
# This is the loop-based script that learners build in Section 1.3.
# It processes each sample one-by-one, waiting for each to complete.
#
# Usage: ./pipeline_sequential.sh <samples.csv>

set -e

SAMPLES_FILE=${1:-data/samples.csv}

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file not found: $SAMPLES_FILE"
    exit 1
fi

echo "=========================================="
echo "RNA-seq Pipeline (Sequential)"
echo "=========================================="
echo ""

# Create output directories
mkdir -p data/fastq results/fastqc results/fastp results/salmon data/salmon_index

# Download salmon index once (shared by all samples)
if [ ! -d "data/salmon_index/salmon" ]; then
    echo "Downloading salmon index..."
    curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz \
        -o data/salmon_index/salmon.tar.gz
    tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
    rm data/salmon_index/salmon.tar.gz
    echo ""
fi

# Read samples from CSV (skip header)
START_TIME=$(date +%s)

tail -n +2 "$SAMPLES_FILE" | while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
    echo "Processing: $sample_id"

    # Download FASTQ files
    echo "  Downloading FASTQ files..."
    if [ ! -f "data/fastq/${sample_id}_R1.fastq.gz" ]; then
        curl -sL "$fastq_r1" -o "data/fastq/${sample_id}_R1.fastq.gz"
    fi
    if [ ! -f "data/fastq/${sample_id}_R2.fastq.gz" ]; then
        curl -sL "$fastq_r2" -o "data/fastq/${sample_id}_R2.fastq.gz"
    fi

    # FastQC
    echo "  Running FastQC..."
    fastqc -q -o results/fastqc \
        "data/fastq/${sample_id}_R1.fastq.gz" \
        "data/fastq/${sample_id}_R2.fastq.gz"

    # fastp
    echo "  Running fastp..."
    fastp \
        -i "data/fastq/${sample_id}_R1.fastq.gz" \
        -I "data/fastq/${sample_id}_R2.fastq.gz" \
        -o "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        -O "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${sample_id}.fastp.json" \
        -h "results/fastp/${sample_id}.fastp.html" \
        2>/dev/null

    # Salmon
    echo "  Running Salmon..."
    salmon quant \
        --index data/salmon_index/salmon \
        --libType A \
        --mates1 "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${sample_id}" \
        --threads 2 \
        --quiet

    echo "  Done: $sample_id"
    echo ""
done

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Run MultiQC to aggregate reports
echo "Running MultiQC..."
multiqc results/ -o results/ --quiet --force
echo ""

echo "=========================================="
echo "Pipeline complete!"
echo "Total time: ${DURATION}s"
echo "=========================================="
echo ""
echo "Notice how each sample waited for the previous one to finish."
echo "Your CPUs were mostly idle during this time."
