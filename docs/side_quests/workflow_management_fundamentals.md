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

You're analyzing bacterial genome sequences.
Your workflow involves:

1. Quality control (FastQC)
2. Adapter trimming (fastp)
3. Genome assembly (SPAdes)
4. Assembly quality assessment (QUAST)

You have 3 samples today. Your PI says 50 more are coming next week.

### 0.2. Navigate to the Project Directory

```bash
cd side-quests/workflow_management_fundamentals
```

### 0.3. Explore the Starting Point

```bash
tree
```

```console title="Output"
.
├── data
│   ├── reads
│   │   ├── sample_01_R1.fastq.gz
│   │   ├── sample_01_R2.fastq.gz
│   │   ├── sample_02_R1.fastq.gz
│   │   ├── sample_02_R2.fastq.gz
│   │   ├── sample_03_R1.fastq.gz
│   │   └── sample_03_R2.fastq.gz
│   └── samples.csv
├── bash_pipeline.sh
├── nextflow_pipeline/
│   ├── main.nf
│   ├── modules/
│   └── nextflow.config
└── solutions/
```

We've provided both a bash script (`bash_pipeline.sh`) and a Nextflow workflow (`nextflow_pipeline/`) that do the same thing.
You'll run both and compare.

---

## 1. Problem: Sequential Processing

### 1.1. Run the Bash Pipeline

```bash
time bash bash_pipeline.sh
```

```console title="Output"
Starting bacterial genome analysis pipeline
===========================================

Processing sample_01 (E.coli)...
  Running FastQC... done (2s)
  Running fastp... done (3s)
  Running SPAdes... done (5s)
  Running QUAST... done (1s)
Completed sample_01

Processing sample_02 (S.aureus)...
  Running FastQC... done (2s)
  Running fastp... done (3s)
  Running SPAdes... done (5s)
  Running QUAST... done (1s)
Completed sample_02

Processing sample_03 (P.aeruginosa)...
  Running FastQC... done (2s)
  Running fastp... done (3s)
  Running SPAdes... done (5s)
  Running QUAST... done (1s)
Completed sample_03

Pipeline complete!

real    0m33.0s
```

**33 seconds** for 3 samples. Each sample waits for the previous one to finish completely.

### 1.2. The Problem

Look at what happened:

```
Timeline (bash):
sample_01: [====FastQC====][======fastp======][==========SPAdes==========][==QUAST==]
sample_02:                                                                            [====FastQC====][======fastp======][==========SPAdes==========][==QUAST==]
sample_03:                                                                                                                                                        [====FastQC====]...
           |-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
           0s      5s     10s     15s     20s     25s     30s     35s
```

These samples are **completely independent**.
Why is sample_02 waiting for sample_01 to finish QUAST before it can start FastQC?

!!! warning "The scaling nightmare"

    With 50 samples at 11 seconds each: **9+ minutes of sequential waiting**.
    On a machine with 16 CPUs, you're using 1-2 at a time. The other 14 sit idle.

You *could* add parallelization to your bash script with `&` and `wait`, but that adds complexity, error handling, and resource management headaches.

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you modify the bash script to run samples in parallel?
    - What happens if you try to run 50 SPAdes assemblies simultaneously on a laptop with 16GB RAM?
    - How would you limit concurrent jobs to match available resources?

### 1.3. Solution: Automatic Parallelization

Now run the same pipeline with Nextflow:

```bash
cd nextflow_pipeline
nextflow run main.nf
```

```console title="Output"
N E X T F L O W  ~  version 24.10.0

executor >  local (12)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3 ✔
[8e/9f0a1b] FASTP (sample_02)      [100%] 3 of 3 ✔
[2c/3d4e5f] SPADES (sample_03)     [100%] 3 of 3 ✔
[7g/8h9i0j] QUAST (sample_01)      [100%] 3 of 3 ✔

Completed at: 2024-10-11 14:23:45
Duration    : 12s
Succeeded   : 12
```

**12 seconds** vs 33 seconds. Same work, 3× faster.

```
Timeline (Nextflow):
sample_01: [==FastQC==][===fastp===][=====SPAdes=====][QUAST]
sample_02: [==FastQC==][===fastp===][=====SPAdes=====][QUAST]
sample_03: [==FastQC==][===fastp===][=====SPAdes=====][QUAST]
           |-------|-------|-------|-------|
           0s      5s     10s     15s
```

### 1.4. How It Works

Nextflow analyzed your data dependencies:

- All 3 FastQC tasks need only raw reads → **run in parallel**
- All 3 fastp tasks need only raw reads → **run in parallel**
- Each SPAdes needs its sample's fastp output → **run when ready**
- Each QUAST needs its sample's SPAdes output → **run when ready**

**You wrote zero parallelization code.**
The workflow manager figured it out from the data flow.

### Takeaway

Sequential bash scripts waste resources.
Workflow managers parallelize automatically based on data dependencies.

---

## 2. Problem: No Resume After Failure

### 2.1. Simulate a Crash

Let's say SPAdes runs out of memory on sample_02 (this happens constantly with real assemblies).

```bash
cd ..
bash bash_pipeline.sh --fail-at sample_02_spades
```

```console title="Output"
Processing sample_01 (E.coli)...
  Running FastQC... done
  Running fastp... done
  Running SPAdes... done
  Running QUAST... done
Completed sample_01

Processing sample_02 (S.aureus)...
  Running FastQC... done
  Running fastp... done
  Running SPAdes... ERROR: Out of memory
```

### 2.2. The Problem

Your options with bash:

1. **Re-run everything** - sample_01 runs again unnecessarily
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
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3, cached: 3 ✔
[8e/9f0a1b] FASTP (sample_02)      [100%] 3 of 3, cached: 3 ✔
[2c/3d4e5f] SPADES (sample_01)     [100%] 3 of 3, cached: 2 ✔
[7g/8h9i0j] QUAST (sample_01)      [100%] 3 of 3, cached: 2 ✔
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

Your bash script assumes tools are installed:

```bash
# From bash_pipeline.sh
fastqc -q -o results/fastqc $read1 $read2
fastp -i $read1 -I $read2 -o trimmed_R1.fq.gz -O trimmed_R2.fq.gz
spades.py -1 trimmed_R1.fq.gz -2 trimmed_R2.fq.gz -o assembly/
quast.py assembly/contigs.fasta -o quast_report/
```

Try running this on a colleague's machine:

```console
$ bash bash_pipeline.sh
Running FastQC...
bash: fastqc: command not found
```

Or worse - it runs but with different versions:

```console
# Your machine
$ fastp --version
fastp 0.23.4

# Colleague's machine
$ fastp --version
fastp 0.20.0   # Different version = potentially different results!
```

### 3.2. The Reproducibility Crisis

This causes real problems:

- **Can't reproduce your own results** after system updates
- **Colleagues get different results** with different tool versions
- **Paper reviewers can't verify** your analysis
- **Future you** has no idea what versions you used

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you share your analysis with a collaborator so they get identical results?
    - What if your pipeline needs tools that require incompatible dependencies (e.g., Python 2.7 vs Python 3.11)?
    - How would you document *exactly* which software versions were used for a publication?

### 3.3. Solution: Per-Process Containers

In Nextflow, each process declares its exact software environment:

```groovy title="modules/fastqc.nf" linenums="1" hl_lines="3"
process FASTQC {
    tag "$meta.id"
    container 'biocontainers/fastqc:0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    // ...
}
```

```groovy title="modules/fastp.nf" linenums="1" hl_lines="3"
process FASTP {
    tag "$meta.id"
    container 'biocontainers/fastp:0.23.4'

    // ...
}
```

**Each process gets its own isolated container** with the exact tool version specified.

### 3.4. Run It Anywhere

On your colleague's machine (or a cluster, or the cloud):

```bash
nextflow run main.nf -profile docker
```

```console title="Output"
Pulling biocontainers/fastqc:0.12.1 ... done
Pulling biocontainers/fastp:0.23.4 ... done
Pulling biocontainers/spades:3.15.5 ... done
Pulling biocontainers/quast:5.2.0 ... done

[... pipeline runs identically ...]
```

**No manual installation. No version conflicts. Perfect reproducibility.**

### 3.5. Mix Incompatible Tools

Ever tried to install tools that need different Python versions?
With containers, each process is isolated:

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

They run in separate containers. No conflicts.

### 3.6. Effortless Collaboration

With containers, sharing your analysis is trivial:

```bash
# You send your colleague:
git clone https://github.com/your-lab/bacterial-assembly
cd bacterial-assembly
nextflow run main.nf -profile docker --samples their_data.csv
```

**Same workflow. Same containers. Same results.**
No "can you send me your installation notes?" emails.

This enables:

- **Reproducible science** - Reviewers can verify your analysis
- **Collaborative research** - Teams get identical results
- **Shared pipelines** - Communities like [nf-core](https://nf-co.re/) publish production-ready workflows

### Takeaway

Bash scripts are environment-dependent and not reproducible.
Workflow managers isolate each step in containers with exact software versions.

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
    - SPAdes needs 16GB RAM, FastQC needs 2GB. How do you prevent memory exhaustion?
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
- **Resource allocation** - Don't run 50 SPAdes jobs simultaneously if you only have 64GB RAM
- **Progress tracking** - Know exactly which samples completed
- **Failure isolation** - One sample's failure doesn't affect others

### 4.3. Declarative Resources

Each process declares what it needs:

```groovy
process SPADES {
    cpus 8
    memory '16.GB'
    time '6.h'

    // ...
}
```

Nextflow uses this to schedule intelligently:

- **On your laptop (16GB RAM)**: Run 1 SPAdes at a time
- **On a cluster (256GB RAM)**: Run 16 SPAdes jobs in parallel
- **On cloud**: Provision appropriately-sized instances

**Same code. Automatic adaptation to available resources.**

### Takeaway

Bash scripts require rewriting to scale.
Workflow managers handle scaling automatically through declarative resource management.

---

## 5. Problem: No Audit Trail

### 5.1. The Traceability Problem

Your bash script ran. Now answer these questions:

- How long did SPAdes take for sample_17?
- Which samples used more than 8GB memory?
- What were the exact commands run for each step?
- When did the pipeline run, and how long did each phase take?

With bash, you'd need to manually add logging, timing, and resource monitoring.
That's a lot of code that has nothing to do with your science.

??? question "Think about it"

    Before seeing the solution, consider:

    - How would you write the "Methods" section of a paper based on your bash script?
    - If a reviewer asks "how long did assembly take for each sample?", can you answer?
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
├── trimmed/
├── assemblies/
├── quast/
├── report.html      # Execution report
├── timeline.html    # Visual timeline
└── trace.txt        # Detailed metrics
```

### 5.3. Execution Timeline

Open `results/timeline.html` in your browser:

![Timeline showing parallel execution](images/timeline_example.png)

See exactly when each task ran, which ran in parallel, and where bottlenecks are.

### 5.4. Detailed Trace

```bash
head results/trace.txt
```

```console title="Output"
task_id  hash       name              status     exit  duration  realtime  %cpu   peak_rss  peak_vmem
1        3a/4b5c6d  FASTQC (sample_01)  COMPLETED  0     1m 23s    1m 18s    98.5%  1.8 GB    2.1 GB
2        8e/9f0a1b  FASTP (sample_01)   COMPLETED  0     3m 45s    3m 40s    390%   2.4 GB    3.1 GB
3        2c/3d4e5f  SPADES (sample_01)  COMPLETED  0     28m 12s   27m 51s   785%   12.1 GB   15.8 GB
```

Every task tracked: duration, CPU usage, memory consumption, exit status.
Perfect for:

- **Optimization** - Which step is the bottleneck?
- **Resource planning** - How much memory does assembly actually need?
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

### 7.4. The Investment

| Learning Curve | Return |
|---------------|--------|
| 1-2 days to become productive | Hours saved per analysis |
| 1 week to feel comfortable | Perfect reproducibility |
| After 2-3 workflows, it's second nature | Seamless collaboration |

### 7.5. Next Steps

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
