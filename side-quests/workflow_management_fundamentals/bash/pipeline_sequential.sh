#!/bin/bash
# Process all RNA-seq samples sequentially

set -e

SAMPLES_FILE=${1:-data/samples.csv}

echo "=========================================="
echo "RNA-seq Pipeline (Sequential)"
echo "=========================================="

mkdir -p data/fastq results/fastqc results/fastp results/salmon data/salmon_index

# Download salmon index once
if [ ! -d "data/salmon_index/salmon" ]; then
    echo "Downloading salmon index..."
    curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz \
        -o data/salmon_index/salmon.tar.gz
    tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
    rm data/salmon_index/salmon.tar.gz
fi

# Process each sample
tail -n +2 "$SAMPLES_FILE" | while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
    echo ""
    echo "Processing: $sample_id"

    # Download
    curl -sL "$fastq_r1" -o "data/fastq/${sample_id}_R1.fastq.gz"
    curl -sL "$fastq_r2" -o "data/fastq/${sample_id}_R2.fastq.gz"

    # FastQC
    echo "  Running FastQC..."
    fastqc -q -o results/fastqc "data/fastq/${sample_id}_R1.fastq.gz" "data/fastq/${sample_id}_R2.fastq.gz"

    # fastp
    echo "  Running fastp..."
    fastp -i "data/fastq/${sample_id}_R1.fastq.gz" -I "data/fastq/${sample_id}_R2.fastq.gz" \
        -o "results/fastp/${sample_id}_trimmed_R1.fastq.gz" -O "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${sample_id}.fastp.json" -h "results/fastp/${sample_id}.fastp.html" 2>/dev/null

    # Salmon
    echo "  Running Salmon..."
    salmon quant --index data/salmon_index/salmon --libType A \
        --mates1 "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${sample_id}" --threads 2 --quiet

    echo "  Done: $sample_id"
done

# Aggregate reports
echo ""
echo "Running MultiQC..."
multiqc results/ -o results/ --quiet --force

echo ""
echo "Pipeline complete!"
