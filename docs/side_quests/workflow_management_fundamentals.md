# Workflow Management Fundamentals

If you're coming from a background of writing shell scripts, Python scripts, or other sequential data processing scripts, Nextflow might seem like unnecessary complexity. Why learn a whole new framework when your bash script does the job?

This side quest takes a different approach: **you'll build a bash script yourself**, experience its limitations firsthand, and then discover how Nextflow solves each problem as it emerges. By the end, you'll understand not just *what* workflow management does, but *why* it exists.

This tutorial uses a single exemplar throughout: a quality control and assembly pipeline for bacterial genome sequencing data.

**What you'll learn:**

- **Why parallelization matters:** Experience sequential processing delays, then see automatic parallelization
- **Why software isolation matters:** Hit version conflicts, then see container-based solutions
- **Why resume capability matters:** Lose hours of work to a failure, then see instant recovery
- **Why portability matters:** See environment-dependent code, then see write-once-run-anywhere
- **Why resource management matters:** Experience resource bottlenecks, then see declarative solutions

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Be comfortable with command-line basics (file navigation, running commands, writing simple bash scripts)
- Have basic familiarity with bioinformatics concepts (FASTQ files, quality control, assembly)
- Understand what containers are at a high level (isolated software environments)

You do NOT need:

- Prior Nextflow experience (we'll teach you as we go)
- Deep knowledge of Groovy or any programming language
- Experience with workflow management systems

### 0.2. The Scenario

You're analyzing bacterial genome sequences. Your workflow involves:

1. Quality control of raw sequencing reads (FastQC)
2. Adapter trimming and filtering (fastp)
3. Genome assembly (SPAdes)
4. Assembly quality assessment (QUAST)

You have 3 samples to process initially, but your PI just said there are 20 more coming next week...

### 0.3. Starting Point

Navigate to the project directory:

```bash title="Navigate to project directory"
cd side-quests/workflow_management_fundamentals
```

The directory contains sample data:

```console title="Directory contents"
> tree
.
├── data
│   ├── samples.csv
│   └── reads
│       ├── sample_01_R1.fastq.gz
│       ├── sample_01_R2.fastq.gz
│       ├── sample_02_R1.fastq.gz
│       ├── sample_02_R2.fastq.gz
│       ├── sample_03_R1.fastq.gz
│       └── sample_03_R2.fastq.gz
└── README.md

3 directories, 8 files
```

Our sample CSV contains metadata:

```csv title="data/samples.csv"
sample_id,organism,read1,read2
sample_01,E.coli,data/reads/sample_01_R1.fastq.gz,data/reads/sample_01_R2.fastq.gz
sample_02,S.aureus,data/reads/sample_02_R1.fastq.gz,data/reads/sample_02_R2.fastq.gz
sample_03,P.aeruginosa,data/reads/sample_03_R1.fastq.gz,data/reads/sample_03_R2.fastq.gz
```

!!! note "About running the commands"

    This tutorial is conceptual - we won't actually execute these computationally intensive commands. Focus on understanding the patterns and problems. The code is realistic and representative of what you'd write in practice.

---

## 1. Building the Bash Script

### 1.1. Start Simple: Process One Sample

Let's begin by processing just one sample. Create a file called `process_one.sh`:

=== "Your Task"

    Write a bash script that:

    1. Creates an output directory called `results/`
    2. Runs FastQC on `data/reads/sample_01_R1.fastq.gz` and `data/reads/sample_01_R2.fastq.gz`
    3. Outputs results to `results/fastqc/`

=== "Solution"

    ```bash title="process_one.sh" linenums="1"
    #!/bin/bash
    set -e  # Exit on error

    echo "Processing sample_01..."

    # Create output directory
    mkdir -p results/fastqc

    # Run FastQC
    fastqc -q -o results/fastqc \
        data/reads/sample_01_R1.fastq.gz \
        data/reads/sample_01_R2.fastq.gz

    echo "Done!"
    ```

This works! For a single sample, this is perfectly fine. But you have 3 samples...

### 1.2. Scale to Multiple Samples: Add a Loop

Now modify your script to process all samples. Read from `data/samples.csv` and loop through them.

=== "Your Task"

    Modify `process_one.sh` to:

    1. Read sample information from `data/samples.csv`
    2. Loop through each sample
    3. Run FastQC on each sample's reads

    Hint: Use `tail -n +2` to skip the CSV header and `while IFS=',' read` to parse CSV lines

=== "Solution"

    ```bash title="process_all.sh" linenums="1"
    #!/bin/bash
    set -e  # Exit on error

    echo "Processing all samples..."
    echo "========================"

    # Create output directory
    mkdir -p results/fastqc

    # Read CSV and process each sample
    tail -n +2 data/samples.csv | while IFS=',' read -r sample_id organism read1 read2; do
        echo ""
        echo "Processing $sample_id ($organism)..."

        fastqc -q -o results/fastqc $read1 $read2

        echo "Completed $sample_id"
    done

    echo ""
    echo "All samples processed!"
    ```

Run this script:

```bash title="Run the script"
chmod +x process_all.sh
./process_all.sh
```

**Expected output:**

```console title="Sequential processing"
Processing all samples...
========================

Processing sample_01 (E.coli)...
Completed sample_01

Processing sample_02 (S.aureus)...
Completed sample_02

Processing sample_03 (P.aeruginosa)...
Completed sample_03

All samples processed!
```

### 1.3. 🛑 Problem #1: It's Slow

Time your script:

```bash title="Time the script"
time ./process_all.sh
```

Notice something? The samples process **one at a time**. If each FastQC takes 5 minutes:

- Sample 1: 5 minutes
- Sample 2: 5 minutes (waits for sample 1)
- Sample 3: 5 minutes (waits for sample 2)
- **Total: 15 minutes**

But these samples are completely independent! They could run **simultaneously**.

**Think about it:** With 20 samples, you're looking at 100 minutes of sequential waiting. On a machine with 16 CPUs, you're using only 2 CPUs at a time. The other 14 sit idle.

You *could* parallelize this with `&` backgrounding or GNU parallel, but that adds complexity. Let's see this problem in action first, then we'll return to it.

### 1.4. Add More Steps: The Full Pipeline

Now let's add the remaining steps. Modify your script to:

1. Run FastQC
2. Trim adapters with fastp
3. Assemble with SPAdes
4. Assess quality with QUAST

=== "Your Task"

    Extend your script to include all four pipeline steps for each sample.

=== "Solution"

    ```bash title="process_complete.sh" linenums="1"
    #!/bin/bash
    set -e  # Exit on error

    echo "Starting bacterial genome analysis pipeline"
    echo "==========================================="

    # Create output directories
    mkdir -p results/fastqc
    mkdir -p results/trimmed
    mkdir -p results/assemblies
    mkdir -p results/quast

    # Read CSV and process each sample
    tail -n +2 data/samples.csv | while IFS=',' read -r sample_id organism read1 read2; do

        echo ""
        echo "Processing $sample_id ($organism)..."
        echo "-----------------------------------"

        # Step 1: Quality control
        echo "Running FastQC..."
        fastqc -q -o results/fastqc $read1 $read2

        # Step 2: Trim adapters
        echo "Running fastp..."
        fastp \
            -i $read1 \
            -I $read2 \
            -o results/trimmed/${sample_id}_R1.fastq.gz \
            -O results/trimmed/${sample_id}_R2.fastq.gz \
            --json results/trimmed/${sample_id}.json \
            --html results/trimmed/${sample_id}.html \
            --thread 4

        # Step 3: Genome assembly
        echo "Running SPAdes assembly..."
        spades.py \
            -1 results/trimmed/${sample_id}_R1.fastq.gz \
            -2 results/trimmed/${sample_id}_R2.fastq.gz \
            -o results/assemblies/${sample_id} \
            --threads 8 \
            --memory 16

        # Step 4: Assembly quality assessment
        echo "Running QUAST..."
        quast.py \
            results/assemblies/${sample_id}/contigs.fasta \
            -o results/quast/${sample_id} \
            --threads 4

        echo "Completed $sample_id"
    done

    echo ""
    echo "Pipeline complete!"
    ```

### 1.5. 🛑 Problem #2: Resource Conflicts

Look at the thread allocations:

- `fastp --thread 4`
- `spades.py --threads 8`
- `quast.py --threads 4`

These are **hardcoded**. What if:

- You run this on your laptop with only 4 cores? SPAdes will fail or thrash
- You run this on an HPC node with 64 cores? You're wasting 56 cores
- You run this on two different machines simultaneously? They'll compete for resources

**The problem:** Resource requirements are baked into your code, not adapted to the environment.

### 1.6. 🛑 Problem #3: The Crash

You start processing your samples. Everything is going well...

```console
Processing sample_01 (E.coli)...
-----------------------------------
Running FastQC...
Running fastp...
Running SPAdes assembly...
Running QUAST...
Completed sample_01

Processing sample_02 (S.aureus)...
-----------------------------------
Running FastQC...
Running fastp...
Running SPAdes assembly...
```

Then disaster strikes at 2 AM:

```console
ERROR: SPAdes terminated with an error
Out of memory: Killed process 12234 (spades.py)
```

**You've lost all the work for sample_02** and any samples that would have come after it. Your options:

1. **Re-run everything** - Waste all completed work (sample_01 runs again unnecessarily)
2. **Manually edit the CSV** - Remove completed samples, re-run script (error-prone)
3. **Build checkpoint logic** - Track what completed, skip finished work (hours of additional coding)

None of these are good options.

### 1.7. 🛑 Problem #4: Environment Dependencies

Your script assumes that `fastqc`, `fastp`, `spades.py`, and `quast.py` are:

- Installed on the system
- In your PATH
- The correct versions
- Compatible with each other

Try to run your script on a colleague's machine:

```console
$ ./process_complete.sh
Processing sample_01 (E.coli)...
Running FastQC...
bash: fastqc: command not found
```

Or worse - it runs with different versions:

```console
# Your machine: fastp version 0.23.4
# Colleague's machine: fastp version 0.20.0
# Different results, different outputs, not reproducible!
```

**The problem:** Your script is **environment-dependent**. It won't run the same way everywhere.

### 1.8. 🛑 Problem #5: No Visibility

Your script is running. Questions you can't easily answer:

- How much memory is SPAdes actually using?
- Which step is the bottleneck?
- How long did sample_02 take compared to sample_01?
- What were the exact parameters used for sample_03?
- If I need to optimize, where should I focus?

You'd need to manually add:

- Logging at every step
- Timestamp tracking
- Resource monitoring
- Parameter recording

That's a lot of boilerplate code that has nothing to do with your science.

### 1.9. Taking Stock

Let's summarize what we've experienced:

| Problem | Impact | Your Options |
|---------|--------|--------------|
| **Sequential processing** | 3× slower than necessary | Add complex parallelization logic |
| **Hardcoded resources** | Fails on some machines, wastes resources on others | Add environment detection, make configurable |
| **No resume capability** | Re-run everything after failures | Build checkpoint tracking system |
| **Environment dependencies** | Not reproducible, won't run elsewhere | Create conda environments, document meticulously |
| **No provenance** | Can't optimize, hard to debug | Add extensive logging code |

You *could* solve all of these. And after weeks of work, you'd have built... a workflow management system.

**This is why Nextflow exists.**

---

## 2. Enter Nextflow: Solving the Problems

### 2.1. The Nextflow Equivalent

Let's see how the same pipeline looks in Nextflow. We'll start with the main workflow file:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

// Include process definitions
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { SPADES } from './modules/spades'
include { QUAST } from './modules/quast'

// Pipeline parameters
params.samples = 'data/samples.csv'
params.outdir = 'results'

// Main workflow
workflow {

    // Read sample sheet and create channel
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id, organism: row.organism]
            def reads = [file(row.read1), file(row.read2)]
            return [meta, reads]
        }

    // Run the pipeline steps
    FASTQC(ch_samples)
    FASTP(ch_samples)
    SPADES(FASTP.out.reads)
    QUAST(SPADES.out.assembly)
}
```

And one of the process definitions:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    tag "$meta.id"
    container 'biocontainers/fastp:0.23.4'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_R{1,2}.fastq.gz"), emit: reads
    path "*.{json,html}", emit: reports

    script:
    def prefix = meta.id
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${prefix}_R1.fastq.gz \\
        -O ${prefix}_R2.fastq.gz \\
        --json ${prefix}.json \\
        --html ${prefix}.html \\
        --thread $task.cpus
    """
}
```

Let's run it:

```bash title="Run the Nextflow pipeline"
nextflow run main.nf
```

### 2.2. ✅ Solution #1: Automatic Parallelization

Watch the execution:

```console title="Nextflow automatically parallelizes"
N E X T F L O W  ~  version 24.10.0
Launching `main.nf` [friendly_darwin] DSL2 - revision: a1b2c3d4e5

executor >  local (12)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3 ✔
[8e/9f0a1b] FASTP (sample_01)      [100%] 3 of 3 ✔
[2c/3d4e5f] SPADES (sample_01)     [100%] 3 of 3 ✔
[7g/8h9i0j] QUAST (sample_01)      [100%] 3 of 3 ✔

Completed at: 11-Oct-2025 14:23:45
Duration    : 28m 15s
CPU hours   : 2.5
Succeeded   : 12
```

**Notice the difference:**

- **Bash script**: 90 minutes (sequential: 30 min × 3 samples)
- **Nextflow**: 28 minutes (parallel where possible)

**What happened?** Nextflow analyzed the dependencies:

1. All 3 FASTQC tasks are independent → **Run simultaneously**
2. All 3 FASTP tasks depend only on input data → **Run simultaneously**
3. Each SPADES depends on its FASTP output → **Run when ready**
4. Each QUAST depends on its SPADES output → **Run when ready**

You wrote **zero parallelization code**. Nextflow figured it out from the data dependencies.

### 2.3. ✅ Solution #2: Declarative Resource Requirements

Look at the process definition again:

```groovy title="Process with resource declarations"
process FASTP {
    cpus 4
    memory '8.GB'

    // ...
}
```

These aren't hardcoded values - they're **requirements**. Nextflow uses them differently based on where you run:

**On your laptop (4 cores, 16 GB RAM):**

- Schedules tasks so they don't exceed available resources
- Might run 1 FASTP task at a time (needs 4 CPUs)
- Queues remaining tasks until resources free up

**On HPC cluster (64 cores, 256 GB RAM):**

- Generates appropriate SLURM job submission script
- Requests 4 CPUs and 8 GB per task
- Can run 16 FASTP tasks in parallel

**On cloud (elastic scaling):**

- Provisions appropriately-sized instances
- Scales up when needed, down when done
- You pay only for what you use

**Same code. Different environments. Automatic adaptation.**

### 2.4. ✅ Solution #3: Built-in Resume Capability

Remember the crash scenario? Let's see it with Nextflow:

```console title="Pipeline fails on sample_02"
N E X T F L O W  ~  version 24.10.0

executor >  local (9)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3 ✔
[8e/9f0a1b] FASTP (sample_01)      [100%] 3 of 3 ✔
[2c/3d4e5f] SPADES (sample_01)     [100%] 2 of 3 ✔
[2c/a1b2c3] SPADES (sample_02)     [  0%] 0 of 1, failed: 1 ✘
[7g/8h9i0j] QUAST (sample_01)      [100%] 2 of 2 ✔

Error executing process > 'SPADES (sample_02)'
```

**Key observations:**

1. Samples 01 and 03 completed successfully
2. Only sample 02 failed
3. Downstream work (QUAST) continued for successful samples

**Fix the problem** (increase memory):

```groovy title="modules/spades.nf" hl_lines="3"
process SPADES {
    tag "$meta.id"
    memory '32.GB'  // Increased from 16.GB
    // ... rest of process
}
```

**Resume from where it failed:**

```bash title="Resume with one flag"
nextflow run main.nf -resume
```

```console title="Only failed work re-runs"
N E X T F L O W  ~  version 24.10.0

executor >  local (2)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3, cached: 3 ✔
[8e/9f0a1b] FASTP (sample_01)      [100%] 3 of 3, cached: 3 ✔
[2c/3d4e5f] SPADES (sample_01)     [100%] 3 of 3, cached: 2 ✔
[7g/8h9i0j] QUAST (sample_01)      [100%] 3 of 3, cached: 2 ✔

Completed at: 11-Oct-2025 14:45:12
Duration    : 32m 5s
```

**Only SPADES(sample_02) and QUAST(sample_02) re-ran.** Everything else used cached results.

**Zero checkpoint logic required.** Nextflow tracks this automatically using content-based hashing.

### 2.5. ✅ Solution #4: Container-Based Software Isolation

Look at the process definition:

```groovy title="Process with container specification"
process FASTP {
    container 'biocontainers/fastp:0.23.4'

    // ... rest of process
}
```

**What Nextflow does automatically:**

1. Downloads the container image (first run only)
2. Executes the task inside the container
3. Ensures exact version 0.23.4 is used
4. Prevents conflicts with other tools

**Running on a colleague's machine:**

```bash title="Works identically"
nextflow run main.nf
```

```console
Pulling biocontainers/fastp:0.23.4 ... done
Pulling biocontainers/spades:3.15.5 ... done
Pulling biocontainers/quast:5.2.0 ... done

[... pipeline runs identically ...]
```

**Zero installation steps.** No PATH configuration. No version conflicts. Perfect reproducibility.

### 2.6. ✅ Solution #5: Comprehensive Provenance

Generate execution reports:

```bash title="Run with reporting flags"
nextflow run main.nf -with-report -with-timeline -with-trace
```

**Execution Report (`report.html`):**

- Task-level resource usage (CPU, memory, time)
- Success/failure status for each task
- Resource efficiency metrics
- Total cost estimates (for cloud)

**Timeline (`timeline.html`):**

- Visual timeline of when each task ran
- Shows parallelization patterns
- Identifies bottlenecks
- Reveals idle time

**Trace file (`trace.txt`):**

```csv
task_id,hash,process,status,exit,duration,realtime,%cpu,%mem,peak_rss,peak_vmem
1,3a/4b5c6d,FASTQC (sample_01),COMPLETED,0,1m 23s,1m 18s,98.5%,12.3%,1.8 GB,2.1 GB
2,8e/9f0a1b,FASTP (sample_01),COMPLETED,0,3m 45s,3m 40s,390.2%,8.7%,2.4 GB,3.1 GB
3,2c/3d4e5f,SPADES (sample_01),COMPLETED,0,28m 12s,27m 51s,785.3%,45.2%,12.1 GB,15.8 GB
```

Perfect for:

- Identifying which samples need more resources
- Optimizing resource allocations
- Cost analysis
- Publication methods sections

**Zero logging code written.** It's built in.

---

## 3. Side-by-Side Comparison

### 3.1. The Same Task, Two Approaches

Let's directly compare what we built:

**Bash Script Approach:**

```bash
./process_complete.sh
```

- ❌ Sequential execution (90 minutes for 3 samples)
- ❌ Hardcoded resource requirements
- ❌ No resume capability
- ❌ Environment-dependent (requires manual software installation)
- ❌ No execution tracking
- ✅ Simple to understand initially
- ✅ No new tools to learn

**Nextflow Approach:**

```bash
nextflow run main.nf -resume
```

- ✅ Automatic parallelization (28 minutes for 3 samples)
- ✅ Adaptive resource management
- ✅ Built-in resume capability
- ✅ Container-based software isolation
- ✅ Comprehensive provenance tracking
- ⚠️ Initial learning curve
- ✅ Same code works everywhere

### 3.2. Real-World Impact: The Numbers

Let's quantify the differences:

**Time Investment:**

| Task | Bash Script | Nextflow |
|------|-------------|----------|
| Initial implementation | 30 minutes | 60 minutes |
| Add parallelization | +2 hours | Already included |
| Add resume logic | +3 hours | Already included |
| Add logging/tracking | +2 hours | Already included |
| Fix environment issues | +1 hour per machine | Works everywhere |
| **Total development** | **8.5 hours** | **1 hour** |

**Execution Time (3 samples):**

| Approach | Time | Speedup |
|----------|------|---------|
| Bash (sequential) | 90 minutes | 1× |
| Nextflow (parallel) | 28 minutes | 3.2× |

**Execution Time (20 samples):**

| Approach | Time | Speedup |
|----------|------|---------|
| Bash (sequential) | 10 hours | 1× |
| Nextflow (parallel) | 2 hours | 5× |

**After a failure at sample 15:**

| Approach | Time to Resume |
|----------|----------------|
| Bash | 7.5 hours (re-run samples 1-15) |
| Nextflow | 30 minutes (only sample 15) |

### 3.3. The Breaking Point

When should you use workflow management instead of bash scripts?

**Consider Nextflow when:**

✅ You have more than 3-5 samples to process
✅ Your pipeline has multiple steps with dependencies
✅ You need to share your work with collaborators
✅ You run analyses on different computers/clusters
✅ Your analyses take more than 30 minutes
✅ You need to re-run analyses with different parameters
✅ You care about reproducibility

**Bash scripts are fine when:**

✅ One-off exploratory analysis you'll never repeat
✅ Single sample, single step
✅ Interactive work where you adjust based on output
✅ Extremely simple operations (< 5 minutes total)

**The transition point:** When you find yourself thinking "I wish I could..." about your bash script - parallelization, resuming, tracking - that's when to switch to Nextflow.

---

## 4. Hands-On: Converting Your Script

Now it's your turn. Let's convert part of the bash script to Nextflow.

### Exercise 1: Create the FastQC Process

Take the FastQC step from the bash script and convert it to a Nextflow process.

=== "Your Task"

    Create `modules/fastqc.nf` with a process that:

    1. Takes sample metadata and read files as input
    2. Runs FastQC on the reads
    3. Outputs the HTML and ZIP reports
    4. Uses the `biocontainers/fastqc:0.12.1` container
    5. Uses 2 CPUs and 2 GB memory

    Use this template:

    ```groovy title="modules/fastqc.nf"
    process FASTQC {
        tag "$meta.id"
        container '???'
        publishDir "${params.outdir}/fastqc", mode: 'copy'

        cpus ???
        memory ???

        input:
        tuple val(meta), path(reads)

        output:
        // What outputs does FastQC produce?

        script:
        """
        fastqc -q -t $task.cpus ${reads}
        """
    }
    ```

=== "Solution"

    ```groovy title="modules/fastqc.nf"
    process FASTQC {
        tag "$meta.id"
        container 'biocontainers/fastqc:0.12.1'
        publishDir "${params.outdir}/fastqc", mode: 'copy'

        cpus 2
        memory '2.GB'

        input:
        tuple val(meta), path(reads)

        output:
        tuple val(meta), path("*.html"), emit: html
        tuple val(meta), path("*.zip"), emit: zip

        script:
        """
        fastqc -q -t $task.cpus ${reads}
        """
    }
    ```

    **Key elements:**

    - `tag "$meta.id"`: Labels output with sample ID
    - `container`: Specifies exact software version
    - `publishDir`: Where to save final results
    - `cpus`/`memory`: Resource requirements (not hardcoded thread counts!)
    - `input`: Takes metadata + read files as a tuple
    - `output`: Declares what files are produced, with named outputs (`emit`)
    - `$task.cpus`: Uses the declared CPU count (adapts to environment)

### Exercise 2: Create the Main Workflow

Now create a workflow that reads the CSV and runs FastQC on all samples.

=== "Your Task"

    Create `main.nf` that:

    1. Includes the FASTQC process
    2. Reads `data/samples.csv`
    3. Creates a channel with [meta, reads] tuples
    4. Runs FASTQC on all samples

    Use this template:

    ```groovy title="main.nf"
    #!/usr/bin/env nextflow

    include { FASTQC } from './modules/fastqc'

    params.samples = 'data/samples.csv'
    params.outdir = 'results'

    workflow {
        ch_samples = Channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .map { row ->
                // Create metadata map
                // Create reads list
                // Return tuple
            }

        FASTQC(ch_samples)
    }
    ```

=== "Solution"

    ```groovy title="main.nf"
    #!/usr/bin/env nextflow

    include { FASTQC } from './modules/fastqc'

    params.samples = 'data/samples.csv'
    params.outdir = 'results'

    workflow {
        ch_samples = Channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample_id, organism: row.organism]
                def reads = [file(row.read1), file(row.read2)]
                return [meta, reads]
            }

        FASTQC(ch_samples)
    }
    ```

    **How it works:**

    - `Channel.fromPath()`: Creates channel from file
    - `.splitCsv(header: true)`: Parses CSV, skips header
    - `.map { row -> ... }`: Transforms each row
    - `meta`: Metadata map (sample ID, organism, etc.)
    - `reads`: List of read file paths
    - `[meta, reads]`: Tuple that matches process input
    - `FASTQC(ch_samples)`: Runs process on all items in channel

### Exercise 3: Add the FASTP Process

Now add the trimming step.

=== "Your Task"

    1. Create `modules/fastp.nf` (similar to FASTQC but with different outputs)
    2. Add it to the workflow after FASTQC
    3. Make sure it receives the input reads (hint: use the original `ch_samples`)

=== "Solution"

    ```groovy title="modules/fastp.nf"
    process FASTP {
        tag "$meta.id"
        container 'biocontainers/fastp:0.23.4'
        publishDir "${params.outdir}/trimmed", mode: 'copy'

        cpus 4
        memory '8.GB'

        input:
        tuple val(meta), path(reads)

        output:
        tuple val(meta), path("*_R{1,2}.fastq.gz"), emit: reads
        path "*.{json,html}", emit: reports

        script:
        def prefix = meta.id
        """
        fastp \\
            -i ${reads[0]} \\
            -I ${reads[1]} \\
            -o ${prefix}_R1.fastq.gz \\
            -O ${prefix}_R2.fastq.gz \\
            --json ${prefix}.json \\
            --html ${prefix}.html \\
            --thread $task.cpus
        """
    }
    ```

    ```groovy title="main.nf" hl_lines="3 17"
    #!/usr/bin/env nextflow

    include { FASTQC } from './modules/fastqc'
    include { FASTP } from './modules/fastp'

    params.samples = 'data/samples.csv'
    params.outdir = 'results'

    workflow {
        ch_samples = Channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample_id, organism: row.organism]
                def reads = [file(row.read1), file(row.read2)]
                return [meta, reads]
            }

        FASTQC(ch_samples)
        FASTP(ch_samples)
    }
    ```

    **Note:** Both FASTQC and FASTP receive `ch_samples`. They'll run in parallel since neither depends on the other.

### Exercise 4: Chain Processes with Dependencies

Now add SPADES, which depends on FASTP output.

=== "Your Task"

    Create the SPADES process and connect it to receive trimmed reads from FASTP.

    Hint: FASTP outputs `emit: reads` - you can reference this as `FASTP.out.reads`

=== "Solution"

    ```groovy title="modules/spades.nf"
    process SPADES {
        tag "$meta.id"
        container 'biocontainers/spades:3.15.5'
        publishDir "${params.outdir}/assemblies", mode: 'copy'

        cpus 8
        memory '16.GB'
        time '6.h'

        input:
        tuple val(meta), path(reads)

        output:
        tuple val(meta), path("${meta.id}/contigs.fasta"), emit: assembly

        script:
        """
        spades.py \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${meta.id} \\
            --threads $task.cpus \\
            --memory ${task.memory.toGiga()}
        """
    }
    ```

    ```groovy title="main.nf" hl_lines="4 21"
    #!/usr/bin/env nextflow

    include { FASTQC } from './modules/fastqc'
    include { FASTP } from './modules/fastp'
    include { SPADES } from './modules/spades'

    params.samples = 'data/samples.csv'
    params.outdir = 'results'

    workflow {
        ch_samples = Channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample_id, organism: row.organism]
                def reads = [file(row.read1), file(row.read2)]
                return [meta, reads]
            }

        FASTQC(ch_samples)
        FASTP(ch_samples)
        SPADES(FASTP.out.reads)
    }
    ```

    **The dependency chain:**

    - `ch_samples` flows into FASTQC and FASTP (parallel)
    - `FASTP.out.reads` flows into SPADES (sequential - SPADES waits for FASTP)
    - Nextflow automatically handles the scheduling

---

## 5. Configuration: Flexibility Without Code Changes

### 5.1. Separating Configuration from Logic

One of Nextflow's powerful features: **configuration lives in separate files**.

Your workflow (`main.nf`) describes **what** to do.
Config files describe **how** to do it.

Create `nextflow.config`:

```groovy title="nextflow.config" linenums="1"
// Project defaults
params {
    samples = 'data/samples.csv'
    outdir = 'results'
}

// Process defaults
process {
    cpus = 2
    memory = '4.GB'
    time = '2.h'

    // Override for specific processes
    withName: 'SPADES' {
        cpus = 8
        memory = '16.GB'
        time = '6.h'
    }
}

// Container settings
docker.enabled = true

// Generate execution reports
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
```

Now you can run with just:

```bash
nextflow run main.nf
```

All settings come from the config file.

### 5.2. Environment-Specific Profiles

Add execution profiles for different environments:

```groovy title="nextflow.config (add to end)" linenums="35"
profiles {

    standard {
        process.executor = 'local'
        docker.enabled = true
    }

    cluster {
        process.executor = 'slurm'
        process.queue = 'general'
        process.clusterOptions = '--account=myproject'
        singularity.enabled = true
        singularity.cacheDir = '/shared/containers'
        docker.enabled = false
    }

    cloud {
        process.executor = 'awsbatch'
        process.queue = 'genomics-queue'
        workDir = 's3://my-bucket/work'
        aws.region = 'us-east-1'
    }

    test {
        params.samples = 'test/mini_samples.csv'
        process.cpus = 1
        process.memory = '2.GB'
    }
}
```

Switch environments with one flag:

```bash title="Different execution environments"
# Local execution
nextflow run main.nf -profile standard

# HPC cluster
nextflow run main.nf -profile cluster

# AWS cloud
nextflow run main.nf -profile cloud

# Quick test
nextflow run main.nf -profile test
```

**The workflow code never changes.**

---

## 6. Conclusion: The Paradigm Shift

### 6.1. What We've Learned

We started by building a bash script and experienced its limitations:

1. ⏱️ **Sequential processing** - Samples waited unnecessarily
2. 🖥️ **Hardcoded resources** - Failed on different machines
3. 💥 **No resume capability** - Lost hours of work to failures
4. 📦 **Environment dependencies** - "Works on my machine" problems
5. 📊 **No provenance** - Couldn't optimize or debug effectively

Then we saw how Nextflow solves each problem:

1. ✅ **Automatic parallelization** - 3× faster with zero code changes
2. ✅ **Adaptive resources** - Same code works everywhere
3. ✅ **Built-in caching** - Resume from failures instantly
4. ✅ **Container isolation** - Perfect reproducibility
5. ✅ **Comprehensive tracking** - Complete provenance built-in

### 6.2. When to Make the Switch

**You should use Nextflow when:**

- Processing 3+ samples in a multi-step pipeline
- Sharing analyses with collaborators
- Running on different computers/clusters
- Analyses take > 30 minutes
- You need reproducibility

**Bash scripts are fine for:**

- Quick one-off explorations
- Single-sample, single-step operations
- Interactive work with immediate feedback
- Very simple tasks (< 5 minutes)

### 6.3. The Investment

**Learning Nextflow:**

- 1-2 days to become productive
- 1 week to feel comfortable
- After 2-3 workflows, you'll prefer it by default

**Return on investment:**

- Hours saved per analysis (scaling with sample count)
- Near-perfect reproducibility
- Seamless collaboration
- Infinite code reuse

### 6.4. Next Steps

**Continue learning:**

1. [Groovy Essentials](./groovy_essentials.md) - Data manipulation and advanced patterns
2. [nf-core](./nf-core.md) - Community best practices and ready-made pipelines
3. [nf-test](./nf-test.md) - Testing your workflows
4. [Nextflow Patterns](https://nextflow-io.github.io/patterns/) - Common workflow patterns

**Join the community:**

- [Nextflow Slack](https://www.nextflow.io/slack-invite.html) - Active community support
- [nf-core Slack](https://nf-co.re/join) - Pipeline-specific help
- [Nextflow Training](https://training.nextflow.io/) - Official workshops

### 6.5. The Mindset Shift

Moving from scripts to workflows isn't just about using a new tool - it's a new way of thinking:

**Script thinking:**
> "How do I run these commands in sequence on my data?"

**Workflow thinking:**
> "What are the data dependencies? How will this scale? How can others reproduce this?"

Once you start thinking in terms of processes, channels, and data flow, you'll wonder how you ever managed with bash scripts.

Welcome to workflow management. Your future self will thank you.

---

## Quick Reference

### Nextflow vs Bash: Command Comparison

| Task | Bash | Nextflow |
|------|------|----------|
| Process one sample | `fastqc file.fq` | `FASTQC(ch_sample)` |
| Process multiple samples | `for f in *.fq; do fastqc $f; done` | `FASTQC(ch_samples)` (automatic) |
| Chain commands | `cmd1 | cmd2 | cmd3` | `CMD1(input); CMD2(CMD1.out); CMD3(CMD2.out)` |
| Parallel execution | `cmd &` or GNU parallel | Automatic based on dependencies |
| Resume after failure | Manual checkpoints | `-resume` flag |
| Specify resources | Hardcoded threads | `cpus = X, memory = 'Y.GB'` |
| Software environment | `module load` or PATH | `container = 'image:version'` |

### Common Nextflow Commands

```bash
# Run workflow
nextflow run main.nf

# Resume from cached results
nextflow run main.nf -resume

# Override parameters
nextflow run main.nf --samples my_data.csv --outdir my_results

# Use specific profile
nextflow run main.nf -profile cluster

# Generate reports
nextflow run main.nf -with-report -with-timeline -with-trace

# Clean up work directory
nextflow clean -f

# List previous runs
nextflow log

# View run details
nextflow log <run_name> -f script,workdir,hash,duration
```

### Process Template

```groovy
process PROCESS_NAME {
    tag "$meta.id"                    // Label for logs
    container 'image:version'         // Software environment
    publishDir "${params.outdir}/dir" // Where to save outputs

    cpus 4                            // CPU requirement
    memory '8.GB'                     // Memory requirement
    time '2.h'                        // Time limit

    input:
    tuple val(meta), path(files)      // Input declaration

    output:
    tuple val(meta), path("*.out"), emit: results

    script:
    """
    command --input $files --output output.out --threads $task.cpus
    """
}
```
