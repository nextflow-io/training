#!/bin/bash
# Process all RNA-seq samples in PARALLEL
#
# Usage: ./pipeline_parallel.sh [samples.csv]

set -e

SAMPLES_FILE=${1:-data/samples.csv}

echo "=========================================="
echo "RNA-seq Pipeline (Parallel)"
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

# Function to process one sample (all steps for a single sample)
process_sample() {
    local sample_id=$1
    local fastq_r1=$2
    local fastq_r2=$3

    echo "[${sample_id}] Starting..."

    # Download
    curl -sL "$fastq_r1" -o "data/fastq/${sample_id}_R1.fastq.gz"
    curl -sL "$fastq_r2" -o "data/fastq/${sample_id}_R2.fastq.gz"

    # FastQC
    fastqc -q -o results/fastqc "data/fastq/${sample_id}_R1.fastq.gz" "data/fastq/${sample_id}_R2.fastq.gz"

    # fastp
    fastp -i "data/fastq/${sample_id}_R1.fastq.gz" -I "data/fastq/${sample_id}_R2.fastq.gz" \
        -o "results/fastp/${sample_id}_trimmed_R1.fastq.gz" -O "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${sample_id}.fastp.json" -h "results/fastp/${sample_id}.fastp.html" 2>/dev/null

    # Salmon
    salmon quant --index data/salmon_index/salmon --libType A \
        --mates1 "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${sample_id}" --threads 2 --quiet

    echo "[${sample_id}] Complete!"
}

export -f process_sample

echo "Launching samples..."

# Process each sample
# TODO: Modify this loop to run samples in PARALLEL
# Hint: Add & after the function call to run in background
while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
    process_sample "$sample_id" "$fastq_r1" "$fastq_r2"
done < <(tail -n +2 "$SAMPLES_FILE")

# TODO: Add wait command to wait for all background jobs

echo ""
echo "Pipeline complete!"
