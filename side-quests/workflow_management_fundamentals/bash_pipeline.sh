#!/bin/bash
# RNA-seq Analysis Pipeline - Bash Version
#
# This script represents the "old way" of doing bioinformatics:
# - Assumes all tools are installed in your environment
# - Processes samples sequentially
# - No built-in resume, parallelization, or provenance
#
# If tools aren't installed, this script will fail - that's the point!
# Compare this to the Nextflow version which handles everything automatically.

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

echo "Starting RNA-seq analysis pipeline (bash version)"
echo "=================================================="
echo ""

# Check for required tools
echo "Checking for required tools..."
for tool in fastqc fastp salmon multiqc curl; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool is not installed!"
        echo ""
        echo "This is exactly the problem workflow managers solve."
        echo "With Nextflow, each process runs in its own container"
        echo "with the exact tool version specified - no installation needed."
        echo ""
        echo "Try the Nextflow version instead:"
        echo "  cd nextflow_pipeline && nextflow run main.nf"
        exit 1
    fi
done
echo "All tools found."
echo ""

# Create output directories
mkdir -p data/fastq
mkdir -p results/fastqc
mkdir -p results/fastp
mkdir -p results/salmon
mkdir -p data/salmon_index

# Download salmon index if not present
if [ ! -d "data/salmon_index/salmon" ]; then
    echo "Downloading salmon index..."
    curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz -o data/salmon_index/salmon.tar.gz
    tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
    rm data/salmon_index/salmon.tar.gz
    echo "Salmon index ready."
    echo ""
fi

# Sample information (hardcoded - another limitation vs workflow managers)
declare -a SAMPLES=("WT_REP1" "WT_REP2" "RAP1_IAA_30M_REP1")
declare -A FASTQ_R1=(
    ["WT_REP1"]="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz"
    ["WT_REP2"]="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz"
    ["RAP1_IAA_30M_REP1"]="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz"
)
declare -A FASTQ_R2=(
    ["WT_REP1"]="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz"
    ["WT_REP2"]="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz"
    ["RAP1_IAA_30M_REP1"]="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz"
)

echo "Samples: ${SAMPLES[*]}"
echo ""

# Process each sample SEQUENTIALLY (this is the problem!)
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."

    # Download FASTQ files if not present
    if [ ! -f "data/fastq/${sample}_R1.fastq.gz" ]; then
        echo "  Downloading FASTQ files..."
        curl -sL "${FASTQ_R1[$sample]}" -o "data/fastq/${sample}_R1.fastq.gz"
        curl -sL "${FASTQ_R2[$sample]}" -o "data/fastq/${sample}_R2.fastq.gz"
    fi

    # Step 1: FastQC (quality control)
    echo -n "  Running FastQC... "
    fastqc -q -o results/fastqc "data/fastq/${sample}_R1.fastq.gz" "data/fastq/${sample}_R2.fastq.gz"
    echo "done"

    # Step 2: fastp (adapter trimming)
    echo -n "  Running fastp... "
    if [[ "$FAIL_AT" == "${sample}_fastp" ]]; then
        echo "ERROR: Simulated failure (disk full)"
        echo ""
        echo "Pipeline failed! All progress for remaining samples is lost."
        echo "To retry, you must either:"
        echo "  1. Re-run the entire pipeline (wasting completed work)"
        echo "  2. Manually edit the script to skip completed samples"
        exit 1
    fi
    fastp \
        -i "data/fastq/${sample}_R1.fastq.gz" \
        -I "data/fastq/${sample}_R2.fastq.gz" \
        -o "results/fastp/${sample}_trimmed_R1.fastq.gz" \
        -O "results/fastp/${sample}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${sample}.fastp.json" \
        -h "results/fastp/${sample}.fastp.html" \
        --quiet
    echo "done"

    # Step 3: Salmon (transcript quantification)
    echo -n "  Running Salmon... "
    if [[ "$FAIL_AT" == "${sample}_salmon" ]]; then
        echo "ERROR: Simulated failure (out of memory)"
        echo ""
        echo "Pipeline failed! All progress for remaining samples is lost."
        echo "To retry, you must either:"
        echo "  1. Re-run the entire pipeline (wasting completed work)"
        echo "  2. Manually edit the script to skip completed samples"
        exit 1
    fi
    salmon quant \
        --index data/salmon_index/salmon \
        --libType A \
        --mates1 "results/fastp/${sample}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${sample}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${sample}" \
        --threads 1 \
        --quiet
    echo "done"

    echo "Completed $sample"
    echo ""
done

# Step 4: MultiQC (aggregate reports)
echo -n "Running MultiQC... "
multiqc results/ -o results/ --quiet --force
echo "done"
echo ""

echo "Pipeline complete!"
echo ""
echo "Notice: Each sample was processed sequentially."
echo "With workflow management, all samples could run in parallel."
echo ""
echo "Try the Nextflow version to see the difference:"
echo "  cd nextflow_pipeline && nextflow run main.nf"
