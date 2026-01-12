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

Now let's fill in each step.

#### 1.3.1. Add the Download Step

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="17"
    # Step 1: Download FASTQ files
    # TODO: Add curl commands to download R1 and R2 files
    ```

=== "After"

    ```bash title="bash/process_sample.sh" linenums="17" hl_lines="2-4"
    # Step 1: Download FASTQ files
    echo "  Downloading FASTQ files..."
    curl -sL "$FASTQ_R1_URL" -o "data/fastq/${SAMPLE_ID}_R1.fastq.gz"
    curl -sL "$FASTQ_R2_URL" -o "data/fastq/${SAMPLE_ID}_R2.fastq.gz"
    ```

#### 1.3.2. Add FastQC

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="22"
    # Step 2: Run FastQC
    # TODO: Add fastqc command
    ```

=== "After"

    ```bash title="bash/process_sample.sh" linenums="22" hl_lines="3-5"
    # Step 2: Run FastQC
    echo "  Running FastQC..."
    fastqc -q -o results/fastqc \
        "data/fastq/${SAMPLE_ID}_R1.fastq.gz" \
        "data/fastq/${SAMPLE_ID}_R2.fastq.gz"
    ```

#### 1.3.3. Add fastp

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="28"
    # Step 3: Run fastp
    # TODO: Add fastp command
    ```

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

#### 1.3.4. Add Salmon Index Download

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="39"
    # Step 4: Download salmon index (if needed)
    # TODO: Add conditional download of salmon index
    ```

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

#### 1.3.5. Add Salmon Quantification

=== "Before"

    ```bash title="bash/process_sample.sh" linenums="48"
    # Step 5: Run Salmon quantification
    # TODO: Add salmon quant command
    ```

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

Open `bash/pipeline_sequential.sh`. The loop structure is already there - you need to add the processing steps inside.

#### Add Processing Logic to the Loop

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

Open `bash/pipeline_parallel.sh`. The function is already there - you just need to add `&` and `wait` to run samples in parallel.

#### Add Background Execution

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

=== "After"

    ```bash title="bash/pipeline_parallel.sh" linenums="58" hl_lines="4 7"
    # Process each sample in parallel
    echo "Launching all samples in parallel..."
    while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
        process_sample "$sample_id" "$fastq_r1" "$fastq_r2" &
    done < <(tail -n +2 "$SAMPLES_FILE")

    wait
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

Open `nextflow/modules/fastqc.nf`. The process structure is there - you need to fill in the input, output, and script.

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

Notice that `tag`, `container`, and `publishDir` are already set. These aren't learning outcomes - you just need to define the data flow and command.

#### 2.2.1. Fill in the Input

=== "Before"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="6"
    input:
    // TODO: Define input - a tuple with sample metadata and read files
    ???
    ```

=== "After"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="6" hl_lines="2"
    input:
    tuple val(meta), path(reads)
    ```

This declares that FASTQC receives a tuple containing:

- `meta` - a map with sample info (like `[id: "WT_REP1"]`)
- `reads` - the FASTQ file paths

#### 2.2.2. Fill in the Output

=== "Before"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="10"
    output:
    // TODO: Define outputs - HTML reports and ZIP files
    ???
    ```

=== "After"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="10" hl_lines="2-3"
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    ```

FastQC produces HTML reports and ZIP files. We capture both and pass along the metadata.

#### 2.2.3. Fill in the Script

=== "Before"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="14"
    script:
    // TODO: Add the fastqc command
    """
    ???
    """
    ```

=== "After"

    ```groovy title="nextflow/modules/fastqc.nf" linenums="14" hl_lines="3"
    script:
    """
    fastqc --quiet --threads 2 ${reads}
    """
    ```

The `${reads}` variable expands to both read files.

!!! tip "Remember installing FastQC with conda?"

    That entire setup is now **one line**: `container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'`

    The version is locked. Your colleague gets the same version. Always.

#### Call FASTQC in main.nf

Open `nextflow/main.nf` and add the FASTQC call:

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="34"
    // TODO: Call FASTQC process with ch_samples
    ```

=== "After"

    ```groovy title="nextflow/main.nf" linenums="34" hl_lines="1"
    FASTQC(ch_samples)
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

Open `nextflow/modules/fastp.nf` and fill in the placeholders using the same pattern as FASTQC.

#### Complete fastp.nf

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

Note the `emit: reads` - this names the output channel so downstream processes can use it.

#### Call FASTP in main.nf

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="36"
    // TODO: Call FASTP process with ch_samples
    ```

=== "After"

    ```groovy title="nextflow/main.nf" linenums="36" hl_lines="1"
    FASTP(ch_samples)
    ```

!!! tip "Different container, no conflicts"

    fastp uses a completely different container than FastQC.
    They could even require incompatible Python versions - doesn't matter.

---

### 2.4. Connecting to Salmon

Now the interesting part - Salmon needs fastp's output (trimmed reads) plus a reference index.

#### Complete salmon.nf

The UNTAR process is already complete. Fill in SALMON_QUANT:

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

Notice `$task.cpus` - Nextflow makes the declared resources available to your script.

#### Wire It Up in main.nf

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="38"
    // TODO: Call SALMON_QUANT with FASTP output and the index
    // Hint: Use FASTP.out.reads for the trimmed reads
    // Hint: Use UNTAR.out.index.first() for the index
    ```

=== "After"

    ```groovy title="nextflow/main.nf" linenums="38" hl_lines="1"
    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())
    ```

The `.first()` converts the index channel to a value that's reused for all samples.

#### Run and Watch

```bash
nextflow run main.nf
```

Nextflow automatically determines:

- FastQC and FASTP can run in parallel (both only need raw reads)
- SALMON_QUANT must wait for FASTP (needs trimmed reads)
- Each sample's SALMON runs as soon as ITS FASTP finishes

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
