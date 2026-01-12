# Workflow Management Fundamentals

If you're coming from a background of writing shell scripts for data processing, workflow management tools might seem like unnecessary complexity.
Why learn a whole new framework when your bash script does the job?

This side quest answers that question through direct experience.
You'll build a real analysis pipeline in bash, hit its limitations firsthand, then rebuild it in Nextflow to see what workflow management actually provides.

**What you'll build:**

An RNA-seq analysis pipeline that runs:

1. Quality control (FastQC)
2. Adapter trimming (fastp)
3. Transcript quantification (Salmon)
4. Report aggregation (MultiQC)

**How this tutorial works:**

| Part | What You'll Do | What You'll Learn |
|------|----------------|-------------------|
| **Part 1: Bash** | Build a pipeline from scratch | Where bash scripts break down |
| **Part 2: Nextflow** | Rebuild the same pipeline | How workflow managers solve each problem |

By the end, you'll understand *why* workflow managers exist - not because someone told you, but because you experienced the problems yourself.

---

## 1. Building a Bash Pipeline

In this part, you'll build an RNA-seq analysis pipeline using bash.
Along the way, you'll encounter real limitations that plague shell-based pipelines.

---

### 1.1. Setup

#### The Scenario

You're a bioinformatician analyzing RNA-seq data from a yeast gene expression study.
You have 3 samples today. Your PI mentions that 50 more samples are coming next week.

#### Navigate to the Working Directory

```bash
cd side-quests/workflow_management_fundamentals
```

#### Explore the Starting Point

```bash
ls -la
```

```console title="Output"
.
├── bash/           # You'll build scripts here
├── data/
│   └── samples.csv # Sample metadata
└── nextflow/       # For Part 2
```

#### Examine the Sample Data

```bash
cat data/samples.csv
```

```console title="Output"
sample,fastq_1,fastq_2
WT_REP1,https://raw.githubusercontent.com/.../SRR6357070_1.fastq.gz,https://...
WT_REP2,https://raw.githubusercontent.com/.../SRR6357072_1.fastq.gz,https://...
RAP1_IAA_30M_REP1,https://raw.githubusercontent.com/.../SRR6357076_1.fastq.gz,https://...
```

Each row is a sample with URLs to paired-end FASTQ files.

---

### 1.2. Installing Your Tools

Before you can run any analysis, you need to install the bioinformatics tools.

#### Create a Conda Environment

=== "Mamba (Recommended)"

    ```bash
    mamba create -n rnaseq-bash fastqc fastp salmon multiqc -c bioconda -c conda-forge -y
    ```

=== "Conda"

    ```bash
    conda create -n rnaseq-bash fastqc fastp salmon multiqc -c bioconda -c conda-forge -y
    ```

Watch the output. Notice:

- Dependencies being resolved (sometimes this takes minutes)
- Packages being downloaded
- Hope that there are no conflicts...

#### Activate the Environment

```bash
mamba activate rnaseq-bash
```

!!! warning "Remember this step"

    Every time you open a new terminal, you must activate this environment again.
    Forget, and your script fails with "command not found".

#### Verify Installation

```bash
fastqc --version
fastp --version
salmon --version
multiqc --version
```

#### Reflection: The Installation Tax

You just spent time on:

- Choosing a package manager (conda vs mamba vs pip?)
- Finding the right channel (bioconda)
- Waiting for dependency resolution
- Hoping nothing conflicts

This is before writing a single line of analysis code.
**Keep this in mind for later.**

---

### 1.3. Processing Your First Sample

Let's start simple: process just ONE sample.

#### Create Your First Script

Create a new file `bash/process_sample.sh`:

```bash title="bash/process_sample.sh" linenums="1"
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
```

#### Make It Executable and Run

```bash
chmod +x bash/process_sample.sh
```

Now run it for one sample (extract values from samples.csv):

```bash
./bash/process_sample.sh WT_REP1 \
    "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz" \
    "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz"
```

```console title="Output"
Processing sample: WT_REP1
  Downloading FASTQ files...
  Running FastQC...
  Running fastp...
  Downloading salmon index...
  Running Salmon...
Completed: WT_REP1
```

#### Check the Results

```bash
ls results/
```

```console title="Output"
fastqc/  fastp/  salmon/
```

It worked! But you have 3 samples (and 50 more coming).
Running this command manually for each one isn't practical.

---

### 1.4. Processing Multiple Samples

Let's wrap the logic in a loop to process all samples.

#### Create a Sequential Pipeline

Create `bash/pipeline_sequential.sh`:

```bash title="bash/pipeline_sequential.sh" linenums="1"
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
```

#### Run and Time It

```bash
chmod +x bash/pipeline_sequential.sh
time ./bash/pipeline_sequential.sh
```

Watch the output. Each sample waits for the previous one to finish completely:

```
Processing: WT_REP1
  Running FastQC...
  Running fastp...
  Running Salmon...
  Done: WT_REP1

Processing: WT_REP2      <- Starts only after WT_REP1 is completely done
  Running FastQC...
  ...
```

#### The Problem: Sequential Execution

```
Timeline (sequential):
WT_REP1:          [====FastQC====][====fastp====][======Salmon======]
WT_REP2:                                                              [====FastQC====][====fastp====][======Salmon======]
RAP1_IAA_30M_REP1:                                                                                                        [====FastQC====]...
```

These samples are **completely independent**.
Why is WT_REP2 waiting for WT_REP1's Salmon to finish before it can start FastQC?

With 50 samples, you're looking at **50× the runtime**.
Your CPUs sit idle while samples wait in line.

---

### 1.5. Adding Parallelization

Let's fix the sequential bottleneck by running samples in parallel.

#### Create a Parallel Pipeline

Create `bash/pipeline_parallel.sh`:

```bash title="bash/pipeline_parallel.sh" linenums="1"
#!/bin/bash
# Process all RNA-seq samples in PARALLEL

set -e

SAMPLES_FILE=${1:-data/samples.csv}

echo "=========================================="
echo "RNA-seq Pipeline (Parallel)"
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

# Function to process one sample
process_sample() {
    local sample_id=$1
    local fastq_r1=$2
    local fastq_r2=$3

    echo "[${sample_id}] Starting..."

    curl -sL "$fastq_r1" -o "data/fastq/${sample_id}_R1.fastq.gz"
    curl -sL "$fastq_r2" -o "data/fastq/${sample_id}_R2.fastq.gz"

    fastqc -q -o results/fastqc "data/fastq/${sample_id}_R1.fastq.gz" "data/fastq/${sample_id}_R2.fastq.gz"

    fastp -i "data/fastq/${sample_id}_R1.fastq.gz" -I "data/fastq/${sample_id}_R2.fastq.gz" \
        -o "results/fastp/${sample_id}_trimmed_R1.fastq.gz" -O "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${sample_id}.fastp.json" -h "results/fastp/${sample_id}.fastp.html" 2>/dev/null

    salmon quant --index data/salmon_index/salmon --libType A \
        --mates1 "results/fastp/${sample_id}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${sample_id}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${sample_id}" --threads 2 --quiet

    echo "[${sample_id}] Complete!"
}

export -f process_sample

echo "Launching all samples in parallel..."

# Launch each sample in background with &
while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
    process_sample "$sample_id" "$fastq_r1" "$fastq_r2" &
done < <(tail -n +2 "$SAMPLES_FILE")

# Wait for all background jobs
wait

echo ""
echo "Running MultiQC..."
multiqc results/ -o results/ --quiet --force

echo ""
echo "Pipeline complete!"
```

#### Run and Compare

```bash
chmod +x bash/pipeline_parallel.sh
time ./bash/pipeline_parallel.sh
```

Notice the interleaved output - all samples running at once:

```
[WT_REP1] Starting...
[WT_REP2] Starting...
[RAP1_IAA_30M_REP1] Starting...
[WT_REP1] Complete!
[RAP1_IAA_30M_REP1] Complete!
[WT_REP2] Complete!
```

The timeline now looks like:

```
Timeline (parallel):
WT_REP1:          [====FastQC====][====fastp====][======Salmon======]
WT_REP2:          [====FastQC====][====fastp====][======Salmon======]
RAP1_IAA_30M_REP1:[====FastQC====][====fastp====][======Salmon======]
                                                                      [MultiQC]
```

Much faster! But...

#### The Hidden Problem

What happens with 50 samples? Or 500?

```bash
# This would launch 500 Salmon jobs simultaneously!
# Each Salmon uses ~4GB RAM
# That's 500 × 4GB = 2TB RAM required!
```

Your machine would crash from memory exhaustion.

??? question "How would you fix this?"

    To limit concurrent jobs, you'd need to:

    - Track how many jobs are running
    - Wait when at the limit
    - Start new jobs as others complete

    This requires:

    - PID tracking arrays
    - A counting loop
    - Sleep/polling logic
    - Handling jobs that finish out of order

    It's doable, but it's 20-30 lines of code that has nothing to do with your science.
    **This is infrastructure code masquerading as analysis.**

---

### 1.6. When Things Go Wrong

Let's simulate what happens when a job fails mid-pipeline.

#### Simulate a Failure

Imagine Salmon crashes on WT_REP2 (out of memory, network timeout, etc.):

```bash
# Re-run but imagine WT_REP2's Salmon step failed
# What happens?
```

With our parallel script, if one sample fails:

- Other samples might complete
- Some might fail
- You don't know which finished successfully
- No checkpoint of progress

#### Your Options After Failure

1. **Re-run everything** - WT_REP1 runs again unnecessarily
2. **Manually track progress** - Comment out completed samples, hope you don't make mistakes
3. **Build checkpoint logic** - Track state in files, check before each step...

??? question "What would checkpointing require?"

    To implement proper checkpointing, you'd need:

    ```bash
    # Track completed steps
    STATE_FILE=".pipeline_state"

    process_sample() {
        local sample_id=$1

        # Check if FastQC already done
        if ! grep -q "${sample_id}:fastqc:done" "$STATE_FILE" 2>/dev/null; then
            fastqc ...
            echo "${sample_id}:fastqc:done" >> "$STATE_FILE"
        fi

        # Check if fastp already done
        if ! grep -q "${sample_id}:fastp:done" "$STATE_FILE" 2>/dev/null; then
            fastp ...
            echo "${sample_id}:fastp:done" >> "$STATE_FILE"
        fi

        # ... and so on for every step
    }
    ```

    This adds:

    - State tracking for every step of every sample
    - File locking for concurrent writes
    - Logic to determine what's "safe to skip"
    - Handling of partial failures within a sample

    **More infrastructure code. More maintenance burden. More chances for bugs.**

#### Reflection: The Failure Tax

Failures are inevitable in computational pipelines:

- Network timeouts
- Out of memory
- Disk full
- Corrupted files
- Tool bugs

Every minute spent re-running already-completed work is wasted.
Every bug in your checkpoint logic could corrupt your results.

---

### 1.7. Sharing With a Colleague

Your colleague wants to run your pipeline. You send them the script.

#### They Try to Run It

```console
$ bash pipeline_parallel.sh
pipeline_parallel.sh: line 34: fastqc: command not found
```

They don't have FastQC installed.

#### After Installing

```console
$ salmon --version
salmon 1.4.0
```

But you used salmon 1.10.3. Different versions can produce different results.

#### The Reproducibility Problem

To share your analysis reproducibly, you'd need to:

- Document every tool and version
- Document the installation steps
- Hope their system resolves conda packages the same way
- Test on their machine
- Repeat when anything updates

??? question "How bad is this really?"

    In a 2016 Nature survey, **70% of researchers** reported failing to reproduce another scientist's experiments.
    Software environments are a major contributor.

    "It worked on my machine" is not science.

---

### 1.8. Scaling to the Cluster

Your PI's 50 samples arrive. Your laptop can't handle it. Time for the cluster.

#### The Cluster Uses SLURM

Your laptop script won't work. The cluster needs:

```bash
#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00

salmon quant ...
```

#### Your Collaborator's Cluster Uses PBS

Different syntax entirely:

```bash
#!/bin/bash
#PBS -N salmon
#PBS -l ncpus=4
#PBS -l mem=8gb
#PBS -l walltime=2:00:00

salmon quant ...
```

#### The Portability Problem

You now need to maintain:

- A laptop version
- A SLURM version
- A PBS version
- Maybe an AWS version

**Same analysis logic. Multiple scripts. Diverging codebases.**

---

### 1.9. Organizing Results

Look at your results directory after running the pipeline:

```bash
ls -la results/
```

```console title="Output"
results/
├── fastqc/
│   ├── WT_REP1_R1_fastqc.html
│   ├── WT_REP1_R1_fastqc.zip
│   ├── WT_REP1_R2_fastqc.html
│   ├── ... (many more files)
├── fastp/
│   ├── WT_REP1.fastp.html
│   ├── WT_REP1.fastp.json
│   ├── WT_REP1_trimmed_R1.fastq.gz
│   ├── ... (mixed outputs and intermediates)
├── salmon/
│   ├── WT_REP1/
│   │   ├── aux_info/
│   │   ├── quant.sf
│   │   └── ... (nested structure)
```

#### The Organization Problem

Your results are scattered across directories you manually created with `mkdir -p`.
But consider:

- What if you want QC reports separate from intermediate files?
- What if you want to keep trimmed FASTQs but not raw downloads?
- What about organizing by sample vs. by tool?
- Who cleans up the `work/` intermediates vs. final outputs?

In your bash script, every output location is hardcoded:

```bash
fastqc -o results/fastqc ...
fastp -o results/fastp/${sample}_trimmed_R1.fastq.gz ...
salmon quant --output results/salmon/${sample} ...
```

#### Changing Organization Requires Code Changes

Want to reorganize? You need to:

- Edit every output path in the script
- Update any downstream code that reads those paths
- Hope you don't break anything
- Re-run everything (since you don't have resume)

**Output organization is tangled with analysis logic.**

---

### 1.10. Part 1 Summary: The True Cost

Let's inventory what we built and what it cost:

| Component | Lines of Code | Purpose |
|-----------|---------------|---------|
| Sample processing logic | ~30 | Actual science |
| Sequential loop | ~10 | Infrastructure |
| Parallelization | ~20 | Infrastructure |
| Checkpointing (sketched) | ~40+ | Infrastructure |
| Job limiting (not built) | ~30+ | Infrastructure |
| Multi-platform support | ~50+ per platform | Infrastructure |
| Results organization | ~15 | Infrastructure |

**Most of your code would be infrastructure, not science.**

### What We Achieved

- [x] Process multiple samples
- [x] Basic parallelization

### What We're Missing

- [ ] Smart job scheduling (fit jobs to available resources)
- [ ] Resume after failure (without re-running everything)
- [ ] Automatic retries (retry just the failed task, not everything)
- [ ] Reproducible environments (same tools, same versions, everywhere)
- [ ] Transparent remote file handling (no manual downloads)
- [ ] Cluster/cloud execution (without rewriting everything)
- [ ] Audit trail / provenance (what ran, when, with what inputs)
- [ ] Clean separation of outputs from intermediates
- [ ] Simple workflow logic (separate "what" from "how")

### The Fundamental Problem

Bash scripts mix **what** to compute with **how** to compute it:

- Analysis logic tangled with parallelization code
- Tool commands mixed with error handling
- Science buried under infrastructure

**There has to be a better way.**

---

## 2. Building a Nextflow Pipeline

Now you'll rebuild the same pipeline using Nextflow.
At each step, we'll call back to the bash equivalent and see what's different.

---

### 2.1. What is Nextflow?

Nextflow is a workflow manager. It handles the infrastructure so you can focus on the science:

- **You declare** what processes exist and what data they need
- **Nextflow figures out** parallelization, scheduling, error handling

#### Key Concepts

| Concept | What It Is |
|---------|------------|
| **Process** | A unit of work (like running FastQC on one sample) |
| **Channel** | A queue of data flowing between processes |
| **Workflow** | How processes connect together |

#### Separation of Concerns

In bash, everything is tangled together - tool commands, parallelization logic, error handling, file management.

In Nextflow, you separate:

- **Process definitions** (in module files): What each tool does, its container, resource needs
- **Workflow logic** (in main.nf): How processes connect - just a few lines showing data flow

This means your workflow logic becomes remarkably simple:

```groovy
FASTQC(samples)
FASTP(samples)
SALMON_QUANT(FASTP.out.reads, index)
```

Three lines. That's the entire pipeline flow. All the complexity of parallelization, containers, file staging, and error handling is declared in the process definitions and handled by Nextflow.

#### No Tool Installation Needed

Remember installing FastQC, fastp, Salmon, MultiQC with conda?

With Nextflow, each process declares its own container.
The tools install themselves, on-demand, with exact versions.

Let's see how.

---

### 2.2. Your First Process: FastQC

#### Create the FastQC Module

Create `nextflow/modules/fastqc.nf`:

```groovy title="nextflow/modules/fastqc.nf" linenums="1"
process FASTQC {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip

    script:
    """
    fastqc --quiet --threads 2 ${reads}
    """
}
```

Let's break this down:

| Line | What It Does |
|------|--------------|
| `tag "$meta.id"` | Label for logs (shows sample name) |
| `container '...'` | The Docker image with FastQC 0.12.1 |
| `publishDir` | Copy outputs to results folder |
| `input:` | What data this process needs |
| `output:` | What files it produces |
| `script:` | The actual command |

!!! tip "Remember installing FastQC with conda?"

    That entire conda environment setup is now **one line**: `container '...'`

    The version is locked. Your colleague gets the same version. Always.

!!! tip "Remember the results organization mess?"

    In bash, you had to manually create directories (`mkdir -p results/fastqc`) and carefully construct output paths.

    With `publishDir`, you declare **where results go** once, and Nextflow handles the rest.
    The organization is part of the process definition, not scattered through your script.

#### Update main.nf

Edit `nextflow/main.nf`:

```groovy title="nextflow/main.nf" linenums="1" hl_lines="8 17-23 27"
#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline - Nextflow Version
 */

// Include the process
include { FASTQC } from './modules/fastqc'

// Pipeline parameters
params.samples = '../data/samples.csv'
params.outdir = '../results'

// Main workflow
workflow {
    // Create channel from sample sheet
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // Run FastQC on all samples
    FASTQC(ch_samples)
}
```

#### Run It

```bash
cd nextflow
nextflow run main.nf
```

```console title="Output"
N E X T F L O W  ~  version 24.10.0

executor >  local (3)
[a1/b2c3d4] FASTQC (WT_REP1)           [100%] 3 of 3 ✔
[e5/f6g7h8] FASTQC (WT_REP2)           [100%] 3 of 3 ✔
[i9/j0k1l2] FASTQC (RAP1_IAA_30M_REP1) [100%] 3 of 3 ✔
```

Notice: **All 3 samples ran in parallel automatically.**

You wrote zero parallelization code. No `&`, no `wait`, no job tracking.
Nextflow saw that the samples were independent and parallelized them.

!!! tip "Remember all that curl downloading in bash?"

    Look at your `samples.csv` - those are **remote URLs**, not local files.

    Nextflow treats remote files like local ones. When a process needs a file, Nextflow:

    - Downloads it automatically (you may see "Staging foreign file..." messages)
    - Caches it so it's not re-downloaded on resume
    - Handles HTTP, HTTPS, FTP, S3, Azure Blob, Google Cloud Storage...

    **No curl commands. No download scripts. Just use the URL.**

---

### 2.3. Adding More Processes

#### Create the fastp Module

Create `nextflow/modules/fastp.nf`:

```groovy title="nextflow/modules/fastp.nf" linenums="1"
process FASTP {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.html"), emit: html

    script:
    def prefix = meta.id
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${prefix}_trimmed_R1.fastq.gz \\
        --out2 ${prefix}_trimmed_R2.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread 4
    """
}
```

!!! tip "Different container, no conflicts"

    fastp uses a completely different container than FastQC.
    They could even require incompatible Python versions - doesn't matter.
    Each process is isolated.

#### Update main.nf

```groovy title="nextflow/main.nf" linenums="1" hl_lines="9 29"
#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline - Nextflow Version
 */

// Include processes
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'

// Pipeline parameters
params.samples = '../data/samples.csv'
params.outdir = '../results'

// Main workflow
workflow {
    // Create channel from sample sheet
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // Run processes
    FASTQC(ch_samples)
    FASTP(ch_samples)
}
```

#### Run It

```bash
nextflow run main.nf
```

FastQC and FASTP run **in parallel** because they both only need the raw reads.

---

### 2.4. Connecting Processes

Now let's add Salmon, which needs fastp's trimmed reads.

#### Create the Salmon Module

Create `nextflow/modules/salmon.nf`:

```groovy title="nextflow/modules/salmon.nf" linenums="1"
process UNTAR {
    container 'ubuntu:22.04'

    input:
    path archive

    output:
    path "salmon", emit: index

    script:
    """
    tar -xzf $archive
    """
}

process SALMON_QUANT {
    tag "$meta.id"
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    """
    salmon quant \\
        --index $index \\
        --libType A \\
        --mates1 ${reads[0]} \\
        --mates2 ${reads[1]} \\
        --output ${meta.id} \\
        --threads $task.cpus
    """
}
```

Notice the **resource declarations**:

```groovy
cpus 4
memory '8.GB'
```

Each process can declare its own resource requirements. Nextflow tracks these and schedules jobs to fit available resources.

!!! tip "Remember the bash parallelization problem?"

    With 50 samples, all Salmon jobs would start at once and crash from OOM.

    Nextflow uses these declarations to **schedule intelligently**:

    - Machine has 16 CPUs? Run 4 Salmon jobs at a time (each needs 4)
    - Machine has 64GB RAM? Run 8 Salmon jobs at a time (each needs 8GB)
    - Different processes can have different requirements - FastQC might use 2 CPUs while Salmon uses 4

    You can also use `maxForks` to explicitly limit concurrent instances of a process, or configure executor-level limits in `nextflow.config`.

    **Declarative resource management. No semaphore code. No manual job limiting.**

#### Update main.nf

```groovy title="nextflow/main.nf" linenums="1" hl_lines="10 15 33-37"
#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline - Nextflow Version
 */

// Include processes
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { UNTAR; SALMON_QUANT } from './modules/salmon'

// Pipeline parameters
params.samples = '../data/samples.csv'
params.outdir = '../results'
params.salmon_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz'

// Main workflow
workflow {
    // Create channel from sample sheet
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // Run processes
    FASTQC(ch_samples)
    FASTP(ch_samples)

    // Salmon needs the index AND fastp's output
    ch_salmon_index = Channel.fromPath(params.salmon_index)
    UNTAR(ch_salmon_index)

    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())
}
```

#### Run and Watch

```bash
nextflow run main.nf
```

Nextflow automatically determines:

- FastQC and FASTP can run in parallel (both only need raw reads)
- SALMON_QUANT must wait for FASTP (needs trimmed reads)
- Each sample's SALMON_QUANT starts as soon as ITS FASTP finishes

```
Timeline (Nextflow automatic):
UNTAR:     [===]
WT_REP1:         [FastQC]    [FASTP]    [======SALMON======]
WT_REP2:         [FastQC]    [FASTP]    [======SALMON======]
RAP1:            [FastQC]    [FASTP]    [======SALMON======]
```

**You wrote zero scheduling code.**
Nextflow inferred the dependencies from the data flow.

---

### 2.5. The Magic of Resume

#### Run the Pipeline

```bash
nextflow run main.nf
```

#### Modify Something

Let's say you want to change a Salmon parameter. Edit the script, then:

```bash
nextflow run main.nf -resume
```

```console title="Output"
[a1/b2c3d4] UNTAR              [100%] 1 of 1, cached: 1 ✔
[e5/f6g7h8] FASTQC (WT_REP1)   [100%] 3 of 3, cached: 3 ✔
[i9/j0k1l2] FASTP (WT_REP1)    [100%] 3 of 3, cached: 3 ✔
[m3/n4o5p6] SALMON_QUANT (WT_REP1) [100%] 3 of 3 ✔  <- Re-ran
```

Only SALMON_QUANT re-ran. Everything else used cached results.

!!! tip "Remember the checkpointing nightmare?"

    You'd need 40+ lines of state tracking, file locking, and recovery logic.

    Nextflow gives you `-resume`. **One flag.**

#### How It Works

Nextflow creates a unique hash for each task based on:

- The script content
- Input file contents (not just names)
- Parameters

If nothing changed, the cached result is used.
This is **scientifically rigorous** - you know exactly what ran vs. what was reused.

#### Automatic Task Retries

What if a task fails due to a transient error - network glitch, temporary resource shortage, or a tool that occasionally crashes?

In bash, the whole pipeline fails and you start debugging.

In Nextflow, you can configure automatic retries per process:

```groovy
process SALMON_QUANT {
    errorStrategy 'retry'
    maxRetries 3

    // Increase memory on each retry
    memory { 8.GB * task.attempt }

    // ... rest of process
}
```

If Salmon fails (maybe it ran out of memory), Nextflow automatically:

1. Retries the **specific failed task** (not the whole pipeline)
2. Increases memory on retry (8GB → 16GB → 24GB)
3. Keeps successful tasks cached
4. Continues other independent tasks while retrying

!!! tip "In bash, a single failure stops everything"

    You'd need to implement retry logic, track which tasks failed, increase resources manually, and restart. With 50 samples, finding and re-running just the 2 that failed is painful.

    Nextflow handles this automatically, per-task, with configurable strategies.

---

### 2.6. Configuration Profiles

#### Update nextflow.config

```groovy title="nextflow/nextflow.config" linenums="1"
// Parameters
params {
    samples = '../data/samples.csv'
    salmon_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz'
    outdir = '../results'
}

// Automatic reports
timeline {
    enabled = true
    file = "${params.outdir}/timeline.html"
}
report {
    enabled = true
    file = "${params.outdir}/report.html"
}
trace {
    enabled = true
    file = "${params.outdir}/trace.txt"
}

// Execution profiles
profiles {
    standard {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    slurm {
        process.executor = 'slurm'
        process.queue = 'general'
        singularity.enabled = true
    }

    aws {
        process.executor = 'awsbatch'
        process.queue = 'genomics-queue'
        aws.region = 'us-east-1'
    }
}
```

#### Run Anywhere

```bash
# On your laptop
nextflow run main.nf -profile standard

# On SLURM cluster
nextflow run main.nf -profile slurm

# On AWS
nextflow run main.nf -profile aws
```

**Same workflow. Same results. Different infrastructure.**

Your processes don't change. The `container`, `cpus`, and `memory` declarations work everywhere. Nextflow handles:

- Submitting jobs to SLURM/PBS/LSF/SGE schedulers
- Spinning up cloud instances on AWS/Google/Azure
- Pulling containers via Docker, Singularity, or Podman
- Transferring input files to remote workers and results back

!!! tip "Imagine doing this in bash..."

    To run your bash pipeline on a SLURM cluster, you'd need to:

    - Rewrite each step as a separate SLURM job script
    - Add `sbatch` commands with `--dependency` flags for job ordering
    - Handle output file staging between jobs
    - Track job IDs and poll for completion
    - Parse scheduler output to detect failures
    - Implement retry logic for transient failures

    For AWS Batch? Start over with a completely different approach.

    **With Nextflow: change one flag.** The same workflow runs on your laptop, your institution's cluster, or a cloud with thousands of CPUs - no code changes.

---

### 2.7. Automatic Provenance

After running the pipeline, check the results directory:

```bash
ls ../results/
```

```console title="Output"
fastqc/          fastp/          salmon/
report.html      timeline.html   trace.txt
```

#### Execution Timeline

Open `timeline.html` in a browser. You'll see:

- When each task ran
- Which tasks ran in parallel
- Where the bottlenecks are

#### Detailed Trace

```bash
head ../results/trace.txt
```

Every task tracked:

- Duration and wall-clock time
- CPU usage percentage
- Memory consumption
- Exit status

!!! tip "Remember we never added logging to the bash script?"

    With Nextflow, comprehensive provenance tracking is built in.
    Zero extra code.

---

### 2.8. Part 2 Summary: What Nextflow Gave Us

### Side-by-Side Comparison

| Feature | Bash | Nextflow |
|---------|------|----------|
| Tool installation | Manual conda/mamba | `container` directive |
| Parallelization | `&` + `wait` + job limiting | Automatic from data flow |
| Failure recovery | Manual checkpoint logic | `-resume` flag |
| Task retries | Manual retry loops | `errorStrategy 'retry'` per process |
| Resource management | Manual, error-prone | Declarative per-process (`cpus`, `memory`, `maxForks`) |
| Remote files | Manual curl/wget scripts | Transparent (URLs work like local paths) |
| Results organization | Manual mkdir + path construction | `publishDir` directive |
| Cluster/cloud execution | Complete rewrite per platform | Same code, different profile |
| Workflow logic | Tangled with infrastructure | Clean process calls (3 lines for whole pipeline) |
| Provenance | Manual logging | Built-in |

### Code Comparison

| Component | Bash Lines | Nextflow Lines |
|-----------|-----------|----------------|
| Core analysis logic | ~30 | ~50 |
| Parallelization | ~20 | 0 |
| Job limiting | ~30 | 0 (declarative resources) |
| Resume/checkpoint | ~40 | 0 (`-resume` built-in) |
| Cluster/cloud support | ~100+/platform | ~10 (profiles) |
| Logging/provenance | ~30 | 0 (built-in) |

**Nextflow code is almost all science. Bash code is mostly infrastructure.**

Running on a cluster or cloud with bash would require completely separate scripts with scheduler-specific job submission, dependency tracking, and file staging - hundreds of lines of infrastructure code per platform. With Nextflow, you write your workflow once and add a few lines of configuration.

---

## When to Use What

### Use a Workflow Manager When:

- Processing more than a handful of samples
- Pipeline has multiple steps
- Results need to be reproducible
- Running on different machines/clusters
- Others need to run your pipeline
- You want to focus on science, not infrastructure

### Bash Scripts Are Fine For:

- Quick one-off explorations
- Single-sample, single-step operations
- Interactive work with immediate feedback
- Simple tasks where setup overhead isn't worth it

---

## Next Steps

This tutorial used Nextflow, but the concepts apply to any workflow manager:

- **Nextflow** - What we used here; strong in bioinformatics
- **Snakemake** - Python-based; popular in genomics
- **WDL** - Broad Institute standard; used in Terra/Cromwell
- **CWL** - Highly portable; vendor-neutral standard

**Continue learning:**

- [Hello Nextflow](../hello_nextflow/index.md) - Full Nextflow introduction
- [nf-core](https://nf-co.re/) - Production-ready community pipelines

---

## Quick Reference

### The Problems Solved

| Problem | Bash Reality | Workflow Solution |
|---------|-------------|-------------------|
| Sequential processing | Samples wait in line | Automatic parallelization |
| No crash recovery | Start over from scratch | Resume from failure point |
| Transient failures | Pipeline crashes, manual restart | Automatic per-task retries |
| Environment chaos | "Works on my machine" | Per-process containers |
| Scaling pain | Rewrite for more samples | Declarative resources (`cpus`, `memory`, `maxForks`) |
| Remote data access | Manual curl/wget everywhere | URLs just work (auto-cached) |
| Results scattered | Manual mkdir + paths everywhere | Declarative publishDir |
| Tangled code | Infrastructure mixed with science | Clean separation (3-line workflow) |
| No audit trail | What ran? Who knows? | Automatic provenance |
| Cluster/cloud execution | Complete rewrite per scheduler | Same workflow, `-profile slurm` or `-profile aws` |

### Nextflow Commands

```bash
# Run workflow
nextflow run main.nf

# Resume from cached results
nextflow run main.nf -resume

# Use specific profile
nextflow run main.nf -profile slurm

# Override parameters
nextflow run main.nf --samples data/more_samples.csv
```
