#!/bin/bash
set -e  # Exit on error

echo "Starting bacterial genome analysis pipeline"
echo "==========================================="

# Create output directories
mkdir -p results/fastqc
mkdir -p results/trimmed
mkdir -p results/assemblies
mkdir -p results/quast

# Read the sample CSV and process each sample
tail -n +2 data/samples.csv | while IFS=',' read -r sample_id organism read1 read2; do

    echo ""
    echo "Processing $sample_id ($organism)..."
    echo "-----------------------------------"

    # Step 1: Quality control with FastQC
    echo "Running FastQC..."
    fastqc -q -o results/fastqc $read1 $read2

    # Step 2: Trim adapters and filter with fastp
    echo "Running fastp..."
    fastp \
        -i $read1 \
        -I $read2 \
        -o results/trimmed/${sample_id}_R1.fastq.gz \
        -O results/trimmed/${sample_id}_R2.fastq.gz \
        --json results/trimmed/${sample_id}.json \
        --html results/trimmed/${sample_id}.html \
        --thread 4

    # Step 3: Genome assembly with SPAdes
    echo "Running SPAdes assembly..."
    spades.py \
        -1 results/trimmed/${sample_id}_R1.fastq.gz \
        -2 results/trimmed/${sample_id}_R2.fastq.gz \
        -o results/assemblies/${sample_id} \
        --threads 8 \
        --memory 16

    # Step 4: Assembly quality assessment with QUAST
    echo "Running QUAST..."
    quast.py \
        results/assemblies/${sample_id}/contigs.fasta \
        -o results/quast/${sample_id} \
        --threads 4

    echo "Completed $sample_id"
done

echo ""
echo "Pipeline complete!"
echo "Results available in results/"
