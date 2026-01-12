# Workflow Management Fundamentals

If you're coming from a background of writing shell scripts or Python scripts for data processing, formal workflow management tools might seem like unnecessary complexity.
Why learn a whole new framework when your bash script does the job?

This side quest answers that question through direct experience.
You'll hit real limitations with bash scripts, then see how workflow management tools solve each one.
We use Nextflow as our example, but the concepts apply to any workflow management system (Snakemake, WDL, CWL, etc.).

**What you'll experience:**

| Problem | What Goes Wrong | How Workflow Management Helps |
|---------|-----------------|-------------------------------|
| Sequential processing | Samples wait in line; CPUs sit idle | Automatic parallelization |
| Crash recovery | Failure means restart from scratch | Built-in caching and resume |
| Environment chaos | "Works on my machine" syndrome | Per-process containers |
| Scaling pain | 3 samples → 300 requires rewrite | Handle any scale, same code |
| No audit trail | What ran? When? With what? | Automatic provenance tracking |
| Portability lock-in | Tied to one machine/cluster | Write once, run anywhere |

---

## 0. Setup

### 0.1. The Scenario

You're analyzing RNA-seq data from a yeast gene expression study.
Your workflow involves:

1. Quality control (FastQC)
2. Adapter trimming (fastp)
3. Transcript quantification (Salmon)
4. Report aggregation (MultiQC)

You have 3 samples today. Your PI says 50 more are coming next week.

!!! note "Real data, real tools"

    This tutorial uses real RNA-seq data from the nf-core test-datasets repository.
    The FASTQ files are from a published yeast study (GSE110004).
    All processes run actual bioinformatics tools in containers.

### 0.2. Navigate to the Project Directory

```bash
cd side-quests/workflow_management_fundamentals
```

### 0.3. Explore the Starting Point

```bash
ls -la
```

```console title="Output"
.
├── data/
│   └── samples.csv
├── bash_pipeline.sh
└── nextflow_pipeline/
    ├── main.nf
    ├── modules/
    └── nextflow.config
```

We've provided both a bash script (`bash_pipeline.sh`) and a Nextflow workflow (`nextflow_pipeline/`).

### 0.4. Install Tools for the Bash Script

To run the bash script, you need to install the required bioinformatics tools:

=== "Conda/Mamba (Recommended)"

    ```bash
    mamba create -n rnaseq-bash fastqc fastp salmon multiqc -c bioconda -c conda-forge -y
    mamba activate rnaseq-bash
    ```

=== "Conda"

    ```bash
    conda create -n rnaseq-bash fastqc fastp salmon multiqc -c bioconda -c conda-forge -y
    conda activate rnaseq-bash
    ```

This installation step is **exactly the kind of friction** that workflow managers eliminate.
Notice you had to:

- Know which package manager to use
- Know the correct channel (`bioconda`)
- Hope there are no dependency conflicts
- Remember to activate the environment

Later, you'll see how Nextflow handles this automatically.

---

## 1. Problem: Sequential Processing

### 1.1. Run the Bash Pipeline

```bash
time bash bash_pipeline.sh
```

```console title="Output"
Starting RNA-seq analysis pipeline (bash version)
==================================================

Checking for required tools...
All tools found.

Samples: WT_REP1 WT_REP2 RAP1_IAA_30M_REP1

Processing WT_REP1...
  Downloading FASTQ files...
  Running FastQC... done
  Running fastp... done
  Running Salmon... done
Completed WT_REP1

Processing WT_REP2...
  Downloading FASTQ files...
  Running FastQC... done
  Running fastp... done
  Running Salmon... done
Completed WT_REP2

Processing RAP1_IAA_30M_REP1...
  Downloading FASTQ files...
  Running FastQC... done
  Running fastp... done
  Running Salmon... done
Completed RAP1_IAA_30M_REP1

Running MultiQC... done

Pipeline complete!

real    2m34.567s
```

Watch the output carefully: each sample waits for the previous one to finish completely.

### 1.2. The Problem

Look at what happened:

```
Timeline (bash - sequential):
WT_REP1:          [FastQC][fastp][====Salmon====]
WT_REP2:                                         [FastQC][fastp][====Salmon====]
RAP1_IAA_30M_REP1:                                                               [FastQC][fastp][====Salmon====]
                                                                                                                 [MultiQC]
```

These samples are **completely independent**.
Why is WT_REP2 waiting for WT_REP1 to finish Salmon before it can start FastQC?

!!! warning "The scaling nightmare"

    With 50 samples processed sequentially, you're looking at **50× the runtime**.
    On a machine with 16 CPUs, you're using 1-4 at a time. The other 12 sit idle.

You *could* add parallelization to your bash script with `&` and `wait`, but that adds complexity, error handling, and resource management headaches.

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you modify the bash script to run samples in parallel?
    - What happens if you try to run 50 Salmon jobs simultaneously on a laptop with 16GB RAM?
    - How would you limit concurrent jobs to match available resources?

### 1.3. Solution: Automatic Parallelization

Now run the same pipeline with Nextflow:

```bash
cd nextflow_pipeline
nextflow run main.nf
```

```console title="Output"
N E X T F L O W  ~  version 24.10.0

executor >  local (11)
[5a/8b2c3d] UNTAR              [100%] 1 of 1 ✔
[3a/4b5c6d] FASTQC (WT_REP1)   [100%] 3 of 3 ✔
[8e/9f0a1b] FASTP (WT_REP2)    [100%] 3 of 3 ✔
[2c/3d4e5f] SALMON_QUANT (RAP1_IAA_30M_REP1) [100%] 3 of 3 ✔
[7g/8h9i0j] MULTIQC            [100%] 1 of 1 ✔

Completed at: 2024-10-11 14:23:45
Duration    : 45s
Succeeded   : 11
```

All samples processed **in parallel** - significantly faster than sequential!

```
Timeline (Nextflow - parallel):
WT_REP1:          [FastQC][fastp][====Salmon====]
WT_REP2:          [FastQC][fastp][====Salmon====]
RAP1_IAA_30M_REP1:[FastQC][fastp][====Salmon====]
                                                 [MultiQC]
```

### 1.4. How It Works

Nextflow analyzed your data dependencies:

- All 3 FastQC tasks need only raw reads → **run in parallel**
- All 3 fastp tasks need only raw reads → **run in parallel**
- Each Salmon needs its sample's fastp output → **run when ready**
- MultiQC needs all reports → **wait for everything, then run**

**You wrote zero parallelization code.**
The workflow manager figured it out from the data flow.

### Takeaway

Sequential bash scripts waste resources.
Workflow managers parallelize automatically based on data dependencies.

---

## 2. Problem: No Resume After Failure

### 2.1. Simulate a Crash

Let's say Salmon runs out of memory on WT_REP2 (this happens with real large datasets).

```bash
cd ..
bash bash_pipeline.sh --fail-at WT_REP2_salmon
```

```console title="Output"
Processing WT_REP1...
  Running FastQC... done (5s)
  Running fastp... done (8s)
  Running Salmon... done (15s)
Completed WT_REP1

Processing WT_REP2...
  Running FastQC... done (5s)
  Running fastp... done (8s)
  Running Salmon... ERROR: Out of memory

Pipeline failed! All progress for remaining samples is lost.
```

### 2.2. The Problem

Your options with bash:

1. **Re-run everything** - WT_REP1 runs again unnecessarily
2. **Manually track progress** - Comment out completed samples, hope you don't make mistakes
3. **Build checkpoint logic** - Hours of additional coding for something that isn't your science

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you implement checkpointing in a bash script?
    - What if the failure happened on sample 47 of 50? How much time would you lose?
    - How would you know which tasks are safe to skip vs. which need re-running?

### 2.3. Solution: Built-in Resume

With Nextflow, fix the issue (increase memory) and resume:

```bash
cd nextflow_pipeline
nextflow run main.nf -resume
```

```console title="Output"
executor >  local (2)
[5a/8b2c3d] UNTAR              [100%] 1 of 1, cached: 1 ✔
[3a/4b5c6d] FASTQC (WT_REP1)   [100%] 3 of 3, cached: 3 ✔
[8e/9f0a1b] FASTP (WT_REP2)    [100%] 3 of 3, cached: 3 ✔
[2c/3d4e5f] SALMON_QUANT (WT_REP2) [100%] 3 of 3, cached: 2 ✔
[7g/8h9i0j] MULTIQC            [100%] 1 of 1 ✔
```

**Only the failed task and its downstream dependencies re-run.**
Everything else uses cached results.

### 2.4. How It Works

Nextflow creates a unique hash for each task based on:

- The script content
- The input files (content, not just names)
- The parameters

If nothing changed, it skips the task and uses cached outputs.
This isn't just convenient - it's **scientifically rigorous**.
You know exactly what ran vs. what was reused.

### Takeaway

Failures are inevitable in computational pipelines.
Workflow managers track completed work automatically, so you never lose progress.

---

## 3. Problem: "Works on My Machine"

### 3.1. The Environment Problem

Remember the setup step where you installed the tools?

```bash
mamba create -n rnaseq-bash fastqc fastp salmon multiqc -c bioconda -c conda-forge -y
mamba activate rnaseq-bash
```

This is **your responsibility** to manage. You must:

- Remember which environment to activate before running
- Document what you installed (and hope the documentation stays current)
- Reinstall everything when you move to a new machine
- Debug version conflicts when tools need incompatible dependencies

What happens when a colleague tries to run your script?

```console
$ bash bash_pipeline.sh
Checking for required tools...
ERROR: salmon is not installed!
```

Or worse - they have different versions:

```console
# Your machine
$ salmon --version
salmon 1.10.3

# Colleague's machine
$ salmon --version
salmon 1.4.0   # Different version = potentially different results!
```

### 3.2. The Reproducibility Crisis

This causes real problems:

- **Can't reproduce your own results** after system updates
- **Colleagues get different results** with different tool versions
- **Paper reviewers can't verify** your analysis
- **Future you** has no idea what versions you used
- **"It worked last month"** - but something changed and you don't know what

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you share your analysis with a collaborator so they get identical results?
    - What if your pipeline needs tools that require incompatible dependencies (e.g., Python 2.7 vs Python 3.11)?
    - How would you document *exactly* which software versions were used for a publication?

### 3.3. Solution: Declarative Software Environments

In Nextflow, each process declares its own software environment:

```groovy title="modules/fastqc.nf" linenums="1" hl_lines="3"
process FASTQC {
    tag "$meta.id"
    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    // ...
}
```

```groovy title="modules/salmon.nf" linenums="1" hl_lines="3"
process SALMON_QUANT {
    tag "$meta.id"
    container 'biocontainers/salmon:1.10.3--h6dccd9a_2'

    // ...
}
```

**The workflow manages the software, not you.**

- Each process specifies exactly what it needs
- Versions are locked and documented in the code itself
- No manual environment activation required
- No "did I install this?" questions

### 3.4. Flexibility in How Environments Are Provided

Nextflow can use containers (Docker, Singularity) OR conda - you choose via profiles:

```bash
# Use Docker containers
nextflow run main.nf -profile docker

# Use Singularity (common on HPC)
nextflow run main.nf -profile singularity

# Use Conda environments
nextflow run main.nf -profile conda
```

The workflow code stays the same. The software management is handled automatically.

### 3.5. Mix Incompatible Tools

Ever tried to install tools that need different Python versions in the same environment?
With workflow managers, each process is isolated:

```groovy
process TOOL_A {
    container 'tool-a:1.0'  // Needs Python 2.7
    // ...
}

process TOOL_B {
    container 'tool-b:2.0'  // Needs Python 3.11
    // ...
}
```

They run in separate environments. No conflicts. No workarounds.

### 3.6. Effortless Collaboration

Sharing your analysis becomes trivial:

```bash
# You send your colleague:
git clone https://github.com/your-lab/rnaseq-analysis
cd rnaseq-analysis
nextflow run main.nf -profile docker --samples their_data.csv
```

**Same workflow. Same software. Same results.**
No "can you send me your installation notes?" emails.

This enables:

- **Reproducible science** - Reviewers can verify your analysis
- **Collaborative research** - Teams get identical results
- **Shared pipelines** - Communities like [nf-core](https://nf-co.re/) publish production-ready workflows

### Takeaway

With bash scripts, **you** manage software environments - installation, activation, documentation, troubleshooting.
With workflow managers, the **workflow** manages environments declaratively - you just run it.

---

## 4. Problem: Scaling Requires Rewriting

### 4.1. The Scaling Problem

Your bash script works for 3 samples.
Now you have 50 samples.

For bash, you need to think about:

- How many parallel jobs can you run?
- How do you manage resource contention?
- What if different steps need different resources?
- How do you handle partial failures across 50 samples?

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you modify your bash script to handle 50 samples? 500 samples?
    - Salmon needs 8GB RAM, FastQC needs 2GB. How do you prevent memory exhaustion?
    - What changes are needed to move from your laptop to a cluster with 100 nodes?

### 4.2. Solution: Declarative Scaling

With Nextflow, going from 3 to 50 samples requires **zero code changes**:

```bash
# 3 samples
nextflow run main.nf --samples data/samples_3.csv

# 50 samples - same command, same code
nextflow run main.nf --samples data/samples_50.csv

# 500 samples - still the same
nextflow run main.nf --samples data/samples_500.csv
```

The workflow manager handles:

- **Scheduling** - Queue tasks, respect resource limits
- **Resource allocation** - Don't run 50 Salmon jobs simultaneously if you only have 64GB RAM
- **Progress tracking** - Know exactly which samples completed
- **Failure isolation** - One sample's failure doesn't affect others

### 4.3. Declarative Resources

Each process declares what it needs:

```groovy
process SALMON_QUANT {
    cpus 4
    memory '8.GB'

    // ...
}
```

Nextflow uses this to schedule intelligently:

- **On your laptop (16GB RAM)**: Run 2 Salmon jobs at a time
- **On a cluster (256GB RAM)**: Run 32 Salmon jobs in parallel
- **On cloud**: Provision appropriately-sized instances

**Same code. Automatic adaptation to available resources.**

### Takeaway

Bash scripts require rewriting to scale.
Workflow managers handle scaling automatically through declarative resource management.

---

## 5. Problem: No Audit Trail

### 5.1. The Traceability Problem

Your bash script ran. Now answer these questions:

- How long did Salmon take for WT_REP2?
- Which samples used more than 8GB memory?
- What were the exact commands run for each step?
- When did the pipeline run, and how long did each phase take?

With bash, you'd need to manually add logging, timing, and resource monitoring.
That's a lot of code that has nothing to do with your science.

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you write the "Methods" section of a paper based on your bash script?
    - If a reviewer asks "how long did quantification take for each sample?", can you answer?
    - How would you identify which step is causing out-of-memory errors?

### 5.2. Solution: Automatic Provenance

Nextflow generates comprehensive reports automatically:

```bash
nextflow run main.nf
ls results/
```

```console title="Output"
results/
├── fastqc/
├── fastp/
├── salmon/
├── multiqc_report.html
├── report.html      # Execution report
├── timeline.html    # Visual timeline
└── trace.txt        # Detailed metrics
```

### 5.3. Execution Timeline

Open `results/timeline.html` in your browser:

See exactly when each task ran, which ran in parallel, and where bottlenecks are.

### 5.4. Detailed Trace

```bash
head results/trace.txt
```

```console title="Output"
task_id  hash       name                   status     exit  duration  realtime  %cpu   peak_rss  peak_vmem
1        5a/8b2c3d  UNTAR                  COMPLETED  0     12s       10s       45%    128 MB    256 MB
2        3a/4b5c6d  FASTQC (WT_REP1)       COMPLETED  0     45s       42s       180%   512 MB    1.2 GB
3        8e/9f0a1b  FASTP (WT_REP1)        COMPLETED  0     1m 12s    1m 8s     350%   1.8 GB    2.4 GB
4        2c/3d4e5f  SALMON_QUANT (WT_REP1) COMPLETED  0     2m 45s    2m 38s    380%   4.2 GB    6.1 GB
```

Every task tracked: duration, CPU usage, memory consumption, exit status.
Perfect for:

- **Optimization** - Which step is the bottleneck?
- **Resource planning** - How much memory does Salmon actually need?
- **Publication** - Methods section writes itself
- **Debugging** - Why did sample_42 fail?

### Takeaway

Bash scripts provide no visibility into execution.
Workflow managers track everything automatically - timing, resources, status, and provenance.

---

## 6. Problem: Tied to One Environment

### 6.1. The Portability Problem

Your bash script runs on your laptop. Now you need to:

- Run on your institution's HPC cluster (uses SLURM)
- Run on a collaborator's cluster (uses PBS)
- Run on AWS for a large dataset

Each environment needs different:

- Job submission commands
- Resource request syntax
- File path conventions
- Container runtimes (Docker vs Singularity)

??? question "Think about it"

    Before seeing the solution, consider:

    - How many different versions of your script would you need to maintain?
    - What if the cluster uses Singularity but your laptop uses Docker?
    - How do you ensure identical results across different execution environments?

### 6.2. Solution: Configuration Profiles

With Nextflow, the **workflow never changes**. Only the configuration:

```groovy title="nextflow.config"
profiles {

    standard {
        process.executor = 'local'
        docker.enabled = true
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

### 6.3. Run Anywhere

```bash
# On your laptop
nextflow run main.nf -profile standard

# On SLURM cluster
nextflow run main.nf -profile slurm

# On AWS
nextflow run main.nf -profile aws
```

**Same workflow. Same results. Different infrastructure.**

### Takeaway

Bash scripts lock you into one environment.
Workflow managers separate logic from infrastructure through configuration profiles.

---

## 7. Summary: Why Workflow Management?

### 7.1. The Problems with Custom Scripts

| Problem | Impact | DIY Solution | Time Cost |
|---------|--------|--------------|-----------|
| Sequential processing | Wasted time and resources | Add parallelization | Days |
| No resume | Re-run everything after failures | Build checkpoint system | Days |
| Environment chaos | Not reproducible | Write installation docs | Hours per machine |
| Scaling pain | Rewrite for each scale | Add resource management | Days |
| No audit trail | Can't debug or optimize | Add comprehensive logging | Days |
| Portability | Tied to one system | Rewrite for each platform | Days per platform |

**Total: Weeks of engineering that isn't your science.**

### 7.2. What Workflow Managers Provide

All of this comes built-in:

- **Automatic parallelization** from data dependencies
- **Resume and caching** with content-based hashing
- **Container isolation** for reproducibility
- **Declarative scaling** for any sample count
- **Comprehensive provenance** for every execution
- **Portable execution** across any infrastructure

### 7.3. When to Use Workflow Management

**Use a workflow manager when:**

- Processing more than a handful of samples
- Pipeline has multiple steps
- Results need to be reproducible
- Running on different machines/clusters
- Analysis takes more than a few minutes
- Others need to run your pipeline

**Bash scripts are fine for:**

- Quick one-off explorations
- Single-sample, single-step operations
- Interactive work with immediate feedback

### 7.4. Next Steps

This tutorial used Nextflow, but the concepts apply broadly:

- **Nextflow** - What we used here; strong in bioinformatics
- **Snakemake** - Python-based; popular in genomics
- **WDL** - Broad Institute standard; used in Terra/Cromwell
- **CWL** - Highly portable; vendor-neutral standard

Choose based on your community and infrastructure. The principles are the same.

**Continue learning:**

- [Hello Nextflow](../hello_nextflow/index.md) - Full Nextflow introduction
- [Snakemake Tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [WDL Documentation](https://openwdl.org/)

---

## Quick Reference

### The Six Problems Solved by Workflow Management

| # | Problem | Bash Reality | Workflow Solution |
|---|---------|-------------|-------------------|
| 1 | Sequential processing | Samples wait in line | Automatic parallelization |
| 2 | No crash recovery | Start over from scratch | Resume from failure point |
| 3 | Environment chaos | "Works on my machine" | Per-process containers |
| 4 | Scaling pain | Rewrite for more samples | Declarative scaling |
| 5 | No audit trail | What ran? Who knows? | Automatic provenance |
| 6 | Portability lock-in | Tied to one machine | Write once, run anywhere |

### Nextflow Commands Used

```bash
# Run workflow
nextflow run main.nf

# Resume from cached results
nextflow run main.nf -resume

# Use specific profile
nextflow run main.nf -profile slurm

# Override parameters
nextflow run main.nf --samples my_data.csv
```
