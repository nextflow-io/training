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

**How this tutorial works:**

| Part | What You'll Do | What You'll Learn |
|------|----------------|-------------------|
| **Part 1: Bash** | Build a pipeline incrementally | Where bash scripts break down |
| **Part 2: Nextflow** | Fill in process definitions | How workflow managers solve each problem |

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

#### Explore the Starter Files

```bash
ls -la bash/
```

```console title="Output"
bash/
├── process_sample.sh      # Single sample script (has TODOs)
├── pipeline_sequential.sh # Multi-sample loop (has TODOs)
└── pipeline_parallel.sh   # Parallel version (has TODOs)
```

You'll fill in these starter files as you progress through the tutorial.

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
    mamba create -n rnaseq-bash fastqc fastp salmon -c bioconda -c conda-forge -y
    ```

=== "Conda"

    ```bash
    conda create -n rnaseq-bash fastqc fastp salmon -c bioconda -c conda-forge -y
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

### 1.3. Building Your First Script

Open the starter file `bash/process_sample.sh`. It has the structure ready - you just need to add the tool commands.

```bash
cat bash/process_sample.sh
```

```bash title="bash/process_sample.sh"
#!/bin/bash
# Process a single RNA-seq sample
set -e

SAMPLE_ID=$1
FASTQ_R1_URL=$2
FASTQ_R2_URL=$3

echo "Processing sample: $SAMPLE_ID"

# Create output directories
mkdir -p data/fastq results/fastqc results/fastp results/salmon data/salmon_index

# Step 1: Download FASTQ files
# TODO: Add curl commands to download R1 and R2 files

# Step 2: Run FastQC
# TODO: Add fastqc command

# Step 3: Run fastp
# TODO: Add fastp command

# Step 4: Download salmon index (if needed)
# TODO: Add conditional download of salmon index

# Step 5: Run Salmon quantification
# TODO: Add salmon quant command

echo "Completed: $SAMPLE_ID"
```

#### Understanding the Starter Script

The script accepts three arguments: a sample ID and two URLs for paired-end FASTQ files. The `set -e` line means the script will stop immediately if any command fails.

The TODO comments mark where you'll add each analysis step:

1. **Download** - Get the raw sequencing data
2. **FastQC** - Check the quality of the raw reads
3. **fastp** - Trim adapters and filter low-quality reads
4. **Salmon** - Quantify transcript expression levels

This is a typical RNA-seq preprocessing workflow. Now let's fill in each step.

#### 1.3.1. Add the Download Step

First, we need to download the FASTQ files from their URLs. We'll use `curl` with `-sL` (silent mode, follow redirects) to fetch each file and save it with a consistent naming pattern.

=== "After"

    ```bash title="bash/process_sample.sh" linenums="17" hl_lines="2-4"
    # Step 1: Download FASTQ files
    echo "  Downloading FASTQ files..."
    curl -sL "$FASTQ_R1_URL" -o "data/fastq/${SAMPLE_ID}_R1.fastq.gz"
    curl -sL "$FASTQ_R2_URL" -o "data/fastq/${SAMPLE_ID}_R2.fastq.gz"
    ```

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="17"
    # Step 1: Download FASTQ files
    # TODO: Add curl commands to download R1 and R2 files
    ```

#### 1.3.2. Add FastQC

FastQC analyzes sequencing data and generates quality reports. We run it on both read files (`R1` and `R2`) and save the output to our results directory. The `-q` flag suppresses progress output.

=== "After"

    ```bash title="bash/process_sample.sh" linenums="22" hl_lines="3-5"
    # Step 2: Run FastQC
    echo "  Running FastQC..."
    fastqc -q -o results/fastqc \
        "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
        "data/fastq/${SAMPLE_ID}_R2.fastq.gz"
    ```

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="22"
    # Step 2: Run FastQC
    # TODO: Add fastqc command
    ```

#### 1.3.3. Add fastp

fastp trims adapter sequences and filters out low-quality reads. It takes paired input files (`-i`, `-I`) and produces trimmed output files (`-o`, `-O`), plus JSON and HTML reports for quality metrics.

=== "After"

    ```bash title="bash/process_sample.sh" linenums="28" hl_lines="3-10"
    # Step 3: Run fastp
    echo "  Running fastp..."
    fastp \
        -i "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
        -I "data/fastq/${SAMPLE_ID}_R2.fastq.gz" \
        -o "results/fastp/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
        -O "results/fastp/${SAMPLE_ID}_trimmed_R2.fastq.gz" \
        -j "results/fastp/${SAMPLE_ID}.fastp.json" \
        -h "results/fastp/${SAMPLE_ID}.fastp.html" \
        2>/dev/null
    ```

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="28"
    # Step 3: Run fastp
    # TODO: Add fastp command
    ```

#### 1.3.4. Add Salmon Index Download

Salmon needs a pre-built index of the reference transcriptome. We'll download a pre-built index (to save time) only if it doesn't already exist. This avoids re-downloading for every sample.

=== "After"

    ```bash title="bash/process_sample.sh" linenums="39" hl_lines="2-8"
    # Step 4: Download salmon index (if needed)
    if [ ! -d "data/salmon_index/salmon" ]; then
        echo "  Downloading salmon index..."
        curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz \
            -o data/salmon_index/salmon.tar.gz
        tar -xzf data/salmon_index/salmon.tar.gz -C data/salmon_index/
        rm data/salmon_index/salmon.tar.gz
    fi
    ```

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="39"
    # Step 4: Download salmon index (if needed)
    # TODO: Add conditional download of salmon index
    ```

#### 1.3.5. Add Salmon Quantification

Salmon quantifies transcript expression by pseudo-aligning reads to the index. It takes the **trimmed** reads from fastp (not the raw reads) and outputs expression counts per transcript.

=== "After"

    ```bash title="bash/process_sample.sh" linenums="48" hl_lines="2-10"
    # Step 5: Run Salmon quantification
    echo "  Running Salmon..."
    salmon quant \
        --index data/salmon_index/salmon \
        --libType A \
        --mates1 "results/fastp/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
        --mates2 "results/fastp/${SAMPLE_ID}_trimmed_R2.fastq.gz" \
        --output "results/salmon/${SAMPLE_ID}" \
        --threads 2 \
        --quiet
    ```

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="48"
    # Step 5: Run Salmon quantification
    # TODO: Add salmon quant command
    ```

#### Test Your Script

Make it executable and run on one sample:

```bash
chmod +x bash/process_sample.sh
./bash/process_sample.sh WT_REP1 \
    "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz" \
    "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz"
```

Check the results:

```bash
ls results/
```

It worked! But you have 3 samples (and 50 more coming).
Running this command manually for each one isn't practical.

---

### 1.4. Processing Multiple Samples

Now that we can process one sample, we need to handle all samples from our CSV file. Open `bash/pipeline_sequential.sh` - it already has the loop structure that reads the CSV and iterates over each sample.

The starter script handles:

- Reading the CSV file and skipping the header row
- Parsing each line into sample ID and FASTQ URLs
- Creating output directories and downloading the shared Salmon index

Your task is to add the actual processing commands inside the loop - the same steps you just wrote, but using the loop variables (`$sample_id`, `$fastq_r1`, `$fastq_r2`) instead of the script arguments.

#### Add Processing Logic to the Loop

=== "After"

    ```bash title="bash/pipeline_sequential.sh" linenums="26" hl_lines="5-24"
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
    ```

=== "Before"

    ```bash title="bash/pipeline_sequential.sh" linenums="26"
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
    ```

#### Run and Time It

```bash
chmod +x bash/pipeline_sequential.sh
time ./bash/pipeline_sequential.sh
```

Watch the output. Each sample waits for the previous one to finish:

```
Processing: WT_REP1
  Running FastQC...
  Running fastp...
  Running Salmon...
  Done: WT_REP1

Processing: WT_REP2      <- Starts only after WT_REP1 is completely done
```

#### The Problem: Sequential Execution

```
Timeline (sequential):
WT_REP1:          [====FastQC====][====fastp====][======Salmon======]
WT_REP2:                                                              [====FastQC====]...
```

These samples are **completely independent**.
Why is WT_REP2 waiting for WT_REP1's Salmon to finish before it can start FastQC?

With 50 samples, you're looking at **50× the runtime**.

---

### 1.5. Adding Parallelization

The sequential script works, but it's slow - each sample waits for the previous one to finish completely. Since samples are independent, they could run simultaneously.

Open `bash/pipeline_parallel.sh`. This version wraps the processing logic in a function called `process_sample()`. The starter script has everything except the parallelization itself.

Your task is to make the loop run samples in parallel by:

1. Adding `&` after each function call to run it in the background
2. Adding `wait` after the loop to pause until all background jobs complete

#### Add Background Execution

=== "After"

    ```bash title="bash/pipeline_parallel.sh" linenums="58" hl_lines="4 7"
    # Process each sample in parallel
    echo "Launching all samples in parallel..."
    while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
        process_sample "$sample_id" "$fastq_r1" "$fastq_r2" &
    done < <(tail -n +2 "$SAMPLES_FILE")

    wait
    ```

=== "Before"

    ```bash title="bash/pipeline_parallel.sh" linenums="58"
    # Process each sample
    # TODO: Modify this loop to run samples in PARALLEL
    # Hint: Add & after the function call to run in background
    while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
        process_sample "$sample_id" "$fastq_r1" "$fastq_r2"
    done < <(tail -n +2 "$SAMPLES_FILE")

    # TODO: Add wait command to wait for all background jobs
    ```

The key changes:

- **`&`** at the end of the function call runs it in the background
- **`wait`** pauses until all background jobs complete

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

Much faster!

#### The Hidden Problem

What happens with 50 samples? Or 500?

```bash
# This would launch 500 Salmon jobs simultaneously!
# Each Salmon uses ~4GB RAM
# That's 500 × 4GB = 2TB RAM required!
```

Your machine would crash from memory exhaustion.

!!! question "How would you limit concurrent jobs?"

    You'd need to track running jobs, wait when at the limit, and start new jobs as others complete. That's 20-30 lines of infrastructure code that has nothing to do with your science.

---

### 1.6. The Pain Points

We won't implement these, but consider what you'd need for a production pipeline:

#### Failure Recovery

If Salmon fails on sample 47 of 50, what happens?

- Everything stops
- No record of what completed
- You'd need checkpoint logic (40+ lines)

#### Reproducibility

When your colleague tries to run your script:

```
salmon: command not found
```

Or worse - they have salmon 1.4.0 but you used 1.10.3. Results differ silently.

#### Scaling to Cluster

Your SLURM cluster needs completely different scripts:

```bash
#SBATCH --job-name=salmon
#SBATCH --cpus-per-task=4
```

Your collaborator's PBS cluster needs yet another rewrite.

---

### 1.7. Part 1 Summary

**What we achieved:**

- [x] Process multiple samples
- [x] Basic parallelization

**What we're missing:**

- [ ] Smart job scheduling (no OOM crashes)
- [ ] Resume after failure
- [ ] Reproducible environments
- [ ] Cluster/cloud portability
- [ ] Audit trail

**The fundamental problem:** Bash scripts mix *what* to compute with *how* to compute it.

**There has to be a better way.**

---

## 2. Building a Nextflow Pipeline

Now you'll rebuild the same pipeline using Nextflow.
At each step, we'll call back to the bash equivalent.

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

#### No Tool Installation Needed

Remember installing FastQC, fastp, Salmon with conda?

With Nextflow, each process declares its own container.
The tools install themselves, on-demand, with exact versions.

---

### 2.2. Your First Process: FastQC

In Nextflow, each tool runs inside a **process** - a self-contained unit that declares its inputs, outputs, and command. You don't write loops or job management; you just describe what each process needs and produces.

Open `nextflow/modules/fastqc.nf` to see the process structure:

```bash
cat nextflow/modules/fastqc.nf
```

```groovy title="nextflow/modules/fastqc.nf"
process FASTQC {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    // TODO: Define input - a tuple with sample metadata and read files
    ???

    output:
    // TODO: Define outputs - HTML reports and ZIP files
    ???

    script:
    // TODO: Add the fastqc command
    """
    ???
    """
}
```

#### Understanding the Process Structure

The starter file provides the "infrastructure" parts you'd write the same way every time:

- **`tag`** - Labels each job with the sample ID (visible in logs)
- **`container`** - The Docker image containing FastQC (no manual installation!)
- **`publishDir`** - Where to copy final outputs

Your job is to define the **data flow** - what goes in, what comes out, and what command runs. These are the learning outcomes for this section.

#### 2.2.1. Fill in the Input

Nextflow processes receive data through **channels**. We need to declare what shape of data this process expects.

=== "After"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="6" hl_lines="2"
    input:
    tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="6"
    input:
    // TODO: Define input - a tuple with sample metadata and read files
    ???
    ```

This declaration tells Nextflow: "I expect a tuple (a grouped set of values) containing two things":

- **`val(meta)`** - A value (metadata map like `[id: "WT_REP1"]`) - stays in memory
- **`path(reads)`** - File path(s) to stage into the working directory

The `val()` vs `path()` distinction matters: `path()` tells Nextflow to make files available in the task's working directory, while `val()` passes data directly.

#### 2.2.2. Fill in the Output

Outputs declare what files the process produces. This is how Nextflow knows what to pass to downstream processes.

=== "After"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="10" hl_lines="2-3"
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    ```

=== "Before"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="10"
    output:
    // TODO: Define outputs - HTML reports and ZIP files
    ???
    ```

FastQC produces HTML reports and ZIP files. We capture both and pass along the metadata so downstream processes know which sample each file belongs to.

The **`emit: html`** and **`emit: zip`** labels create named output channels. Later, you can access these as `FASTQC.out.html` or `FASTQC.out.zip`.

#### 2.2.3. Fill in the Script

The script block contains the actual command - the same command you'd run in bash. Variables from the input block are available here.

=== "After"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="14" hl_lines="3"
    script:
    """
    fastqc --quiet --threads 2 ${reads}
    """
    ```

=== "Before"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="14"
    script:
    // TODO: Add the fastqc command
    """
    ???
    """
    ```

The `${reads}` variable expands to both read files (they were staged as a list). This is the same FastQC command from your bash script, just with Nextflow variable substitution.

!!! tip "Remember installing FastQC with conda?"

    That entire setup is now **one line**: `container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'`

    The version is locked. Your colleague gets the same version. Always.

#### Call FASTQC in main.nf

Open `nextflow/main.nf` and add the FASTQC call:

=== "After"

    ```groovy title="nextflow/main.nf" linenums="34" hl_lines="1"
    FASTQC(ch_samples)
    ```

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="34"
    // TODO: Call FASTQC process with ch_samples
    ```

#### Test It

```bash
cd nextflow
nextflow run main.nf
```

```console title="Output"
N E X T F L O W  ~  version 24.10.0

executor >  local (3)
[a1/b2c3d4] FASTQC (WT_REP1)           [100%] 3 of 3 ✔
```

**All 3 samples ran in parallel automatically.** No `&`, no `wait`, no job tracking.

---

### 2.3. Adding fastp

Now that you understand the process structure, let's add fastp. This process follows the same pattern but introduces one key difference: its output will be used by another process (Salmon needs the trimmed reads).

Open `nextflow/modules/fastp.nf`. The structure is identical to FASTQC - fill in the input, output, and script blocks.

#### Complete fastp.nf

Pay attention to the `emit: reads` on the trimmed reads output. This is how we'll connect fastp to Salmon.

=== "After"

    ```groovy title="nextflow/modules/fastp.nf" hl_lines="7 10-12 15 17-24"
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

=== "Before"

    ```groovy title="nextflow/modules/fastp.nf"
    process FASTP {
        tag "$meta.id"
        container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
        publishDir "${params.outdir}/fastp", mode: 'copy'

        input:
        ???

        output:
        ???

        script:
        """
        ???
        """
    }
    ```

The key line is `emit: reads` on the trimmed read files output. This creates a named output channel that Salmon will consume. When you later write `FASTP.out.reads`, you're accessing this specific output.

#### Call FASTP in main.nf

Now connect FASTP to the sample channel in the workflow:

=== "After"

    ```groovy title="nextflow/main.nf" linenums="36" hl_lines="1"
    FASTP(ch_samples)
    ```

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="36"
    // TODO: Call FASTP process with ch_samples
    ```

!!! tip "Different container, no conflicts"

    fastp uses a completely different container than FastQC.
    They could even require incompatible Python versions - doesn't matter.

---

### 2.4. Connecting to Salmon

This is where data flow gets interesting. Salmon needs **two** inputs:

1. **Trimmed reads** - Different for each sample (from FASTP)
2. **Reference index** - Same for all samples (from UNTAR)

This is a common pattern: per-sample data combined with a shared reference.

#### Complete salmon.nf

The UNTAR process (which extracts the pre-built Salmon index) is already complete. Your task is to fill in SALMON_QUANT, which has **two input declarations** - one for the reads tuple and one for the index path:

=== "After"

    ```groovy title="nextflow/modules/salmon.nf" linenums="17" hl_lines="10-11 14 18-24"
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

=== "Before"

    ```groovy title="nextflow/modules/salmon.nf" linenums="17"
    process SALMON_QUANT {
        tag "$meta.id"
        container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
        publishDir "${params.outdir}/salmon", mode: 'copy'

        cpus 4
        memory '8.GB'

        input:
        ???

        output:
        ???

        script:
        """
        ???
        """
    }
    ```

Notice two things in the script:

- **`$task.cpus`** - Nextflow makes declared resources available as variables. If you later change `cpus 4` to `cpus 8`, the script automatically uses 8 threads.
- **Two inputs** - The process receives the trimmed reads tuple and the index path separately.

#### Wire It Up in main.nf

Now comes the crucial part - connecting the processes. SALMON_QUANT needs two arguments:

=== "After"

    ```groovy title="nextflow/main.nf" linenums="38" hl_lines="1"
    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())
    ```

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="38"
    // TODO: Call SALMON_QUANT with FASTP output and the index
    // Hint: Use FASTP.out.reads for the trimmed reads
    // Hint: Use UNTAR.out.index.first() for the index
    ```

Understanding this line:

- **`FASTP.out.reads`** - The trimmed reads channel (one item per sample)
- **`UNTAR.out.index.first()`** - The `.first()` operator converts the channel to a single value that gets reused for every sample

Without `.first()`, Nextflow would try to pair each sample with a separate index item. Since there's only one index, we tell Nextflow to reuse it.

#### Run and Watch

Run the complete pipeline:

```bash
nextflow run main.nf
```

Watch what happens - Nextflow automatically determines the execution order from the data flow you defined:

- **FastQC and FASTP can run in parallel** - Both only need the raw reads, no dependency between them
- **SALMON_QUANT must wait for FASTP** - It needs the trimmed reads output
- **Each sample runs independently** - Sample 2's Salmon can start as soon as Sample 2's fastp finishes, even if Sample 1 is still running

!!! tip "Remember `&` and `wait`?"

    Gone. Nextflow infers parallelization from the data flow.

---

### 2.5. The Magic of Resume

Run the pipeline, then modify something and run with `-resume`:

```bash
nextflow run main.nf -resume
```

```console title="Output"
[a1/b2c3d4] UNTAR              [100%] 1 of 1, cached: 1 ✔
[e5/f6g7h8] FASTQC (WT_REP1)   [100%] 3 of 3, cached: 3 ✔
[i9/j0k1l2] FASTP (WT_REP1)    [100%] 3 of 3, cached: 3 ✔
[m3/n4o5p6] SALMON_QUANT       [100%] 3 of 3 ✔  <- Only this re-ran
```

!!! tip "Remember the checkpointing nightmare?"

    You'd need 40+ lines of state tracking.
    Nextflow gives you `-resume`. **One flag.**

---

### 2.6. Configuration Profiles

The `nextflow.config` file lets you run the same workflow anywhere:

```bash
# On your laptop
nextflow run main.nf -profile standard

# On SLURM cluster
nextflow run main.nf -profile slurm

# On AWS
nextflow run main.nf -profile aws
```

!!! tip "Remember rewriting scripts for SLURM vs PBS?"

    Same workflow. Same results. Different infrastructure.

---

### 2.7. Part 2 Summary

| Feature | Bash | Nextflow |
|---------|------|----------|
| Tool installation | Manual conda | `container` directive |
| Parallelization | `&` + `wait` + job limiting | Automatic from data flow |
| Failure recovery | Manual checkpoint logic | `-resume` flag |
| Resource management | Manual, error-prone | Declarative per-process |
| Cluster/cloud | Complete rewrite | Same code, different profile |

**Nextflow code is almost all science. Bash code is mostly infrastructure.**

---

## Quick Reference

### Nextflow Commands

```bash
# Run workflow
nextflow run main.nf

# Resume from cached results
nextflow run main.nf -resume

# Use specific profile
nextflow run main.nf -profile slurm
```

### Process Template

```groovy
process TOOL_NAME {
    tag "$meta.id"
    container 'quay.io/biocontainers/tool:version'
    publishDir "${params.outdir}/tool", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: results

    script:
    """
    tool command ${reads}
    """
}
```
