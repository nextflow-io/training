#!/bin/bash
# Bacterial Genome Analysis Pipeline - Bash Version
# Demonstrates the limitations that workflow managers solve

set -e

# Parse arguments
FAIL_AT=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --fail-at)
            FAIL_AT="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

echo "Starting bacterial genome analysis pipeline"
echo "==========================================="
echo ""

# Create output directories
mkdir -p results/fastqc
mkdir -p results/trimmed
mkdir -p results/assemblies
mkdir -p results/quast

# Process each sample SEQUENTIALLY (this is the problem!)
tail -n +2 data/samples.csv | while IFS=',' read -r sample_id organism read1 read2; do

    echo "Processing $sample_id ($organism)..."

    # Step 1: FastQC
    echo -n "  Running FastQC... "
    sleep 2  # Simulate processing time
    # In real life: fastqc -q -o results/fastqc $read1 $read2
    touch "results/fastqc/${sample_id}_R1_fastqc.html"
    touch "results/fastqc/${sample_id}_R2_fastqc.html"
    echo "done (2s)"

    # Step 2: fastp
    echo -n "  Running fastp... "
    sleep 3  # Simulate processing time
    # In real life: fastp -i $read1 -I $read2 -o trimmed_R1.fq.gz -O trimmed_R2.fq.gz
    cp "$read1" "results/trimmed/${sample_id}_trimmed_R1.fastq.gz"
    cp "$read2" "results/trimmed/${sample_id}_trimmed_R2.fastq.gz"
    echo "done (3s)"

    # Step 3: SPAdes (check for simulated failure)
    echo -n "  Running SPAdes... "
    if [[ "$FAIL_AT" == "${sample_id}_spades" ]]; then
        echo "ERROR: Out of memory"
        echo ""
        echo "Pipeline failed! All progress for remaining samples is lost."
        echo "To retry, you must either:"
        echo "  1. Re-run the entire pipeline (wasting completed work)"
        echo "  2. Manually edit the script to skip completed samples"
        exit 1
    fi
    sleep 5  # Simulate processing time (SPAdes is slow!)
    # In real life: spades.py -1 trimmed_R1.fq.gz -2 trimmed_R2.fq.gz -o assembly/
    mkdir -p "results/assemblies/${sample_id}"
    echo ">contig_1" > "results/assemblies/${sample_id}/contigs.fasta"
    echo "ATCGATCG" >> "results/assemblies/${sample_id}/contigs.fasta"
    echo "done (5s)"

    # Step 4: QUAST
    echo -n "  Running QUAST... "
    sleep 1  # Simulate processing time
    # In real life: quast.py assembly/contigs.fasta -o quast_report/
    mkdir -p "results/quast/${sample_id}"
    echo "Assembly: ${sample_id}" > "results/quast/${sample_id}/report.tsv"
    echo "done (1s)"

    echo "Completed $sample_id"
    echo ""

done

echo "Pipeline complete!"
echo ""
echo "Total time: ~33 seconds for 3 samples (11 seconds each, sequential)"
echo ""
echo "Notice: Each sample waited for the previous one to completely finish."
echo "With 50 samples, this would take 9+ minutes of unnecessary waiting."
