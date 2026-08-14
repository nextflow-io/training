#!/bin/bash
# pipeline_parallel.sh - Process all RNA-seq samples in parallel
#
# This is the parallelized script that learners build in Section 1.4.
# It uses & to background jobs and wait to synchronize.
#
# Usage: ./pipeline_parallel.sh <samples.csv>
#
# WARNING: This script launches ALL samples simultaneously!
# With many samples, this could exhaust memory or overwhelm the system.
# A production script would need job limiting - see the discussion in the tutorial.

set -e

SAMPLES_FILE=${1:-data/samples.csv}

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file not found: $SAMPLES_FILE"
    exit 1
fi

echo "=========================================="
echo "RNA-seq Pipeline (Parallel)"
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

# Function to process a single sample
process_sample() {
    local sample_id=$1
    local fastq_r1=$2
    local fastq_r2=$3

    echo "[${sample_id}] Starting..."

    # Download FASTQ files
    if [ ! -f "data/fastq/${sample_id}_R1.fastq.gz" ]; then
        curl -sL "$fastq_r1" -o "data/fastq/${sample_id}_R1.fastq.gz"
    fi
    if [ ! -f "data/fastq/${sample_id}_R2.fastq.gz" ]; then
        curl -sL "$fastq_r2" -o "data/fastq/${sample_id}_R2.fastq.gz"
    fi

    # FastQC
    fastqc -q -o results/fastqc \
        "data/fastq/${sample_id}_R1.fastq.gz" \
        "data/fastq/${sample_id}_R2.fastq.gz"

    # fastp
    fastp \
        -i "data/fastq/${sample_id}_R1.fastq.gz" \
        -I "data/fastq/${sample_id}_R2.fastq.gz" \
        -o "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        -O "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${sample_id}.fastp.json" \
        -h "results/fastp/${sample_id}.fastp.html" \
        2>/dev/null

    # Salmon
    salmon quant \
        --index data/salmon_index/salmon \
        --libType A \
        --mates1 "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${sample_id}" \
        --threads 2 \
        --quiet

    echo "[${sample_id}] Complete!"
}

# Export function so subshells can use it
export -f process_sample

START_TIME=$(date +%s)

# Launch all samples in parallel using & (backgrounding)
echo "Launching samples in parallel..."
PIDS=()

while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
    process_sample "$sample_id" "$fastq_r1" "$fastq_r2" &
    PIDS+=($!)
done < <(tail -n +2 "$SAMPLES_FILE")

# Wait for all background jobs to complete
wait

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Run MultiQC to aggregate reports
echo ""
echo "Running MultiQC..."
multiqc results/ -o results/ --quiet --force
echo ""

echo "=========================================="
echo "Pipeline complete!"
echo "Total time: ${DURATION}s"
echo "=========================================="
echo ""
echo "All samples ran in parallel - notice the speedup!"
echo ""
echo "BUT... what if you had 50 samples? Or 500?"
echo "They would ALL start at once, potentially:"
echo "  - Exhausting memory (each Salmon uses ~4GB)"
echo "  - Overwhelming disk I/O"
echo "  - Crashing your system"
echo ""
echo "To fix this, you'd need to add job limiting logic..."
echo "...which is exactly what workflow managers do automatically."
