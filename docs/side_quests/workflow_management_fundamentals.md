# Workflow Management Fundamentals

If you're coming from a background of writing shell scripts, Python scripts, or other sequential data processing scripts, Nextflow might seem like unnecessary complexity.
Why learn a whole new framework when your bash script does the job?

This side quest takes a different approach: **you'll build a bash script yourself**, experience its limitations firsthand, and then discover how Nextflow solves each problem as it emerges.
By the end, you'll understand not just *what* workflow management does, but *why* it exists.

This tutorial uses a single exemplar throughout: a quality control and assembly pipeline for bacterial genome sequencing data.

**What you'll learn:**

- **Why parallelization matters:** Experience sequential processing delays, then see automatic parallelization
- **Why resume capability matters:** Lose progress to a failure, then see instant recovery
- **Why portability matters:** See environment-dependent code, then see write-once-run-anywhere
- **Why resource management matters:** Experience resource bottlenecks, then see declarative solutions

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Be comfortable with command-line basics (file navigation, running commands, writing simple bash scripts)
- Have basic familiarity with bioinformatics concepts (FASTQ files, quality control, assembly)
- Have completed at least Part 1 of the Hello Nextflow training course

### 0.2. The Scenario

You're analyzing bacterial genome sequences. Your workflow involves:

1. Quality control of raw sequencing reads (FastQC)
2. Adapter trimming and filtering (fastp)
3. Genome assembly (SPAdes)
4. Assembly quality assessment (QUAST)

You have 3 samples to process initially, but your PI just said there are 20 more coming next week...

### 0.3. Navigate to the Project Directory

```bash
cd side-quests/workflow_management_fundamentals
```

### 0.4. Explore the Starting Point

The directory contains sample data, starter scripts, and a `solutions/` folder with completed files:

```bash
tree
```

```console title="Output"
.
├── README.md
├── data
│   ├── reads
│   │   ├── sample_01_R1.fastq.gz
│   │   ├── sample_01_R2.fastq.gz
│   │   ├── sample_02_R1.fastq.gz
│   │   ├── sample_02_R2.fastq.gz
│   │   ├── sample_03_R1.fastq.gz
│   │   └── sample_03_R2.fastq.gz
│   └── samples.csv
├── main.nf
├── modules
│   └── fastqc.nf
├── nextflow.config
├── process_one.sh
└── solutions
    ├── main.nf
    ├── modules
    │   ├── fastp.nf
    │   ├── fastqc.nf
    │   ├── quast.nf
    │   └── spades.nf
    ├── nextflow.config
    └── process_complete.sh
```

!!! tip "Solutions available"

    The `solutions/` directory contains completed versions of all files.
    Reference them if you get stuck, but try to work through the exercises yourself first!

View the sample metadata:

```bash
cat data/samples.csv
```

```csv title="data/samples.csv"
sample_id,organism,read1,read2
sample_01,E.coli,data/reads/sample_01_R1.fastq.gz,data/reads/sample_01_R2.fastq.gz
sample_02,S.aureus,data/reads/sample_02_R1.fastq.gz,data/reads/sample_02_R2.fastq.gz
sample_03,P.aeruginosa,data/reads/sample_03_R1.fastq.gz,data/reads/sample_03_R2.fastq.gz
```

### Takeaway

You now have a project directory with sample data and starter scripts ready to build upon.

### What's next?

Experience the limitations of bash scripts for bioinformatics pipelines.

---

## 1. The Bash Script Problem

### 1.1. Start Simple: Process One Sample

We've provided a simple bash script that processes a single sample.
Let's examine it and run it.

```bash
cat process_one.sh
```

```bash title="process_one.sh" linenums="1"
#!/bin/bash
set -e  # Exit on error

echo "Processing sample_01..."

# Create output directory
mkdir -p results/fastqc

# Run FastQC (mock - just echoes for learning purposes)
echo "Running FastQC on sample_01_R1.fastq.gz..."
sleep 1
echo "Running FastQC on sample_01_R2.fastq.gz..."
sleep 1

echo "Done!"
```

Make it executable and run it:

```bash
chmod +x process_one.sh
./process_one.sh
```

```console title="Output"
Processing sample_01...
Running FastQC on sample_01_R1.fastq.gz...
Running FastQC on sample_01_R2.fastq.gz...
Done!
```

This works for one sample. But you have 3 samples...

### 1.2. Scale to Multiple Samples

Modify the script to process all samples by adding a loop.
Open `process_one.sh` in your editor and replace its contents:

=== "After"

    ```bash title="process_one.sh" linenums="1" hl_lines="4-5 8-11 13-20"
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

        echo "Running FastQC on $read1..."
        sleep 1
        echo "Running FastQC on $read2..."
        sleep 1

        echo "Completed $sample_id"
    done

    echo ""
    echo "All samples processed!"
    ```

=== "Before"

    ```bash title="process_one.sh" linenums="1"
    #!/bin/bash
    set -e  # Exit on error

    echo "Processing sample_01..."

    # Create output directory
    mkdir -p results/fastqc

    # Run FastQC (mock - just echoes for learning purposes)
    echo "Running FastQC on sample_01_R1.fastq.gz..."
    sleep 1
    echo "Running FastQC on sample_01_R2.fastq.gz..."
    sleep 1

    echo "Done!"
    ```

Run the modified script and time it:

```bash
time ./process_one.sh
```

```console title="Output"
Processing all samples...
========================

Processing sample_01 (E.coli)...
Running FastQC on data/reads/sample_01_R1.fastq.gz...
Running FastQC on data/reads/sample_01_R2.fastq.gz...
Completed sample_01

Processing sample_02 (S.aureus)...
Running FastQC on data/reads/sample_02_R1.fastq.gz...
Running FastQC on data/reads/sample_02_R2.fastq.gz...
Completed sample_02

Processing sample_03 (P.aeruginosa)...
Running FastQC on data/reads/sample_03_R1.fastq.gz...
Running FastQC on data/reads/sample_03_R2.fastq.gz...
Completed sample_03

All samples processed!

real    0m6.0s
```

### 1.3. Problem #1: It's Slow (Sequential Processing)

Notice something? The samples process **one at a time**.
If each sample takes 2 seconds, the total is 6 seconds.

But these samples are completely independent! They could run **simultaneously**.

!!! warning "The scaling problem"

    With 20 samples, you're looking at 40 seconds of sequential waiting.
    On a machine with 16 CPUs, you're using only 1 CPU at a time.
    The other 15 sit idle.

You *could* parallelize this with `&` backgrounding or GNU parallel, but that adds complexity and error handling challenges.

### 1.4. Problem #2: What If It Crashes?

Imagine your script is halfway through processing 20 samples when it crashes.
You've lost all progress and must re-run everything from the beginning.

There's no built-in way to:

- Resume from where it failed
- Skip samples that already completed
- Track what has been done

### 1.5. Problem #3: Environment Dependencies

Your script assumes that tools are installed and available.
Try running this on a colleague's machine and you might see:

```console
bash: fastqc: command not found
```

Or worse - it runs with different versions and produces different results.

### Takeaway

Bash scripts are great for simple tasks, but they have fundamental limitations:

- Sequential processing wastes resources
- No resume capability after failures
- Environment-dependent (not reproducible)
- Manual parallelization is error-prone

### What's next?

See how Nextflow solves these problems automatically.

---

## 2. Enter Nextflow

### 2.1. Examine the Starter Workflow

We've provided a starter Nextflow workflow that runs FastQC on all samples.
Let's examine it:

```bash
cat main.nf
```

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

// Include process definitions
include { FASTQC } from './modules/fastqc'

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

    // Run FastQC on all samples
    FASTQC(ch_samples)
}
```

And the FASTQC process:

```bash
cat modules/fastqc.nf
```

```groovy title="modules/fastqc.nf" linenums="1"
process FASTQC {
    tag "$meta.id"
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
    echo "Running FastQC on ${meta.id}..."
    sleep 1

    for read in ${reads}; do
        BASENAME=\$(basename \$read .fastq.gz)
        echo "<html><body><h1>FastQC: \${BASENAME}</h1></body></html>" > \${BASENAME}_fastqc.html
        echo "data" > \${BASENAME}_fastqc_data.txt
        zip -q \${BASENAME}_fastqc.zip \${BASENAME}_fastqc_data.txt
        rm \${BASENAME}_fastqc_data.txt
    done
    """
}
```

### 2.2. Run the Nextflow Pipeline

```bash
nextflow run main.nf
```

```console title="Output"
N E X T F L O W  ~  version 25.04.3
Launching `main.nf` [friendly_darwin] DSL2

executor >  local (3)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3 ✔
```

**Notice the difference:**

- **Bash script**: 6 seconds (sequential: 2 sec × 3 samples)
- **Nextflow**: ~2 seconds (all 3 samples processed in parallel!)

Check the results:

```bash
ls results/fastqc/
```

```console title="Output"
sample_01_R1_fastqc.html  sample_02_R1_fastqc.html  sample_03_R1_fastqc.html
sample_01_R1_fastqc.zip   sample_02_R1_fastqc.zip   sample_03_R1_fastqc.zip
sample_01_R2_fastqc.html  sample_02_R2_fastqc.html  sample_03_R2_fastqc.html
sample_01_R2_fastqc.zip   sample_02_R2_fastqc.zip   sample_03_R2_fastqc.zip
```

### 2.3. Solution #1: Automatic Parallelization

Nextflow analyzed the data dependencies and determined:

- All 3 FASTQC tasks are independent → **Run simultaneously**

You wrote **zero parallelization code**. Nextflow figured it out from the data flow.

### 2.4. Add the FASTP Trimming Step

Now let's extend the pipeline to include adapter trimming.
First, create the FASTP module:

```bash
cat > modules/fastp.nf << 'EOF'
process FASTP {
    tag "$meta.id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_R{1,2}.fastq.gz"), emit: reads
    path "*.{json,html}", emit: reports

    script:
    def prefix = meta.id
    """
    echo "Running fastp on ${meta.id}..."
    sleep 2

    # Pass through reads (mock trimming)
    zcat ${reads[0]} | gzip > ${prefix}_trimmed_R1.fastq.gz
    zcat ${reads[1]} | gzip > ${prefix}_trimmed_R2.fastq.gz

    # Create mock reports
    echo '{"summary": {"before": {"total_reads": 3}, "after": {"total_reads": 3}}}' > ${prefix}.json
    echo "<html><body><h1>fastp: ${prefix}</h1></body></html>" > ${prefix}.html
    """
}
EOF
```

Now update `main.nf` to include FASTP:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4 24"
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

=== "Before"

    ```groovy title="main.nf" linenums="1"
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

Run the pipeline again:

```bash
nextflow run main.nf -resume
```

```console title="Output"
executor >  local (3)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3, cached: 3 ✔
[8e/9f0a1b] FASTP (sample_01)      [100%] 3 of 3 ✔
```

**Notice:** FASTQC results were **cached**! Only FASTP ran.

### 2.5. Solution #2: Built-in Resume Capability

Nextflow automatically tracks completed work using content-based hashing.
When you use `-resume`:

- Completed tasks are skipped
- Only new or changed tasks run
- You never lose progress

### 2.6. Add the SPADES Assembly Step

Create the SPADES module that depends on FASTP output:

```bash
cat > modules/spades.nf << 'EOF'
process SPADES {
    tag "$meta.id"
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
    echo "Running SPAdes assembly on ${meta.id}..."
    sleep 3

    mkdir -p ${meta.id}
    cat > ${meta.id}/contigs.fasta << FASTA
>NODE_1_length_300_cov_50.0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>NODE_2_length_250_cov_45.5
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>NODE_3_length_180_cov_40.0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
FASTA
    """
}
EOF
```

Update `main.nf` to chain SPADES after FASTP:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5 23"
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

=== "Before"

    ```groovy title="main.nf" linenums="1"
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

**Key insight:** `SPADES(FASTP.out.reads)` creates a dependency chain.
SPADES will wait for FASTP to complete before running.

Run the pipeline:

```bash
nextflow run main.nf -resume
```

```console title="Output"
executor >  local (3)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3, cached: 3 ✔
[8e/9f0a1b] FASTP (sample_01)      [100%] 3 of 3, cached: 3 ✔
[2c/3d4e5f] SPADES (sample_01)     [100%] 3 of 3 ✔
```

### 2.7. Complete the Pipeline with QUAST

Finally, add the QUAST quality assessment step:

```bash
cat > modules/quast.nf << 'EOF'
process QUAST {
    tag "$meta.id"
    publishDir "${params.outdir}/quast", mode: 'copy'

    cpus 2
    memory '4.GB'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}/*"), emit: report

    script:
    """
    echo "Running QUAST on ${meta.id}..."
    sleep 1

    mkdir -p ${meta.id}

    CONTIG_COUNT=\$(grep -c '^>' ${assembly} || echo "0")
    TOTAL_LENGTH=\$(grep -v '^>' ${assembly} | tr -d '\\n' | wc -c)

    cat > ${meta.id}/report.tsv << TSV
Assembly\t${meta.id}
# contigs\t\${CONTIG_COUNT}
Total length\t\${TOTAL_LENGTH}
N50\t250
TSV

    cat > ${meta.id}/report.html << HTML
<html><body><h1>QUAST: ${meta.id}</h1>
<p>Contigs: \${CONTIG_COUNT}</p></body></html>
HTML
    """
}
EOF
```

Update `main.nf` to complete the pipeline:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="6 25"
    #!/usr/bin/env nextflow

    include { FASTQC } from './modules/fastqc'
    include { FASTP } from './modules/fastp'
    include { SPADES } from './modules/spades'
    include { QUAST } from './modules/quast'

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
        QUAST(SPADES.out.assembly)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
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

Run the complete pipeline:

```bash
nextflow run main.nf -resume
```

```console title="Output"
executor >  local (3)
[3a/4b5c6d] FASTQC (sample_01)     [100%] 3 of 3, cached: 3 ✔
[8e/9f0a1b] FASTP (sample_01)      [100%] 3 of 3, cached: 3 ✔
[2c/3d4e5f] SPADES (sample_01)     [100%] 3 of 3, cached: 3 ✔
[7g/8h9i0j] QUAST (sample_01)      [100%] 3 of 3 ✔
```

Verify the results:

```bash
cat results/quast/sample_01/report.tsv
```

```console title="Output"
Assembly        sample_01
# contigs       3
Total length    243
N50             250
```

### Takeaway

You've built a complete pipeline that:

- Runs quality control (FASTQC) in parallel
- Trims adapters (FASTP) in parallel
- Assembles genomes (SPADES) as trimmed reads become available
- Assesses quality (QUAST) as assemblies complete

All with automatic parallelization and resume capability.

### What's next?

Learn how to configure your pipeline for different execution environments.

---

## 3. Configuration

### 3.1. Examine the Configuration File

Nextflow separates **what** to do (workflow) from **how** to do it (configuration).

```bash
cat nextflow.config
```

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
}

// Container settings - disabled for mock processes
docker.enabled = false

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

### 3.2. Add Environment Profiles

Add profiles to run the pipeline in different environments.
Update `nextflow.config` by adding profiles at the end:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="32-60"
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
    docker.enabled = false

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

    // Profiles for different execution environments
    profiles {

        standard {
            process.executor = 'local'
            docker.enabled = false
        }

        docker {
            docker.enabled = true
            process {
                withName: 'FASTQC' { container = 'biocontainers/fastqc:0.12.1' }
                withName: 'FASTP' { container = 'biocontainers/fastp:0.23.4' }
                withName: 'SPADES' { container = 'biocontainers/spades:3.15.5' }
                withName: 'QUAST' { container = 'biocontainers/quast:5.2.0' }
            }
        }

        cluster {
            process.executor = 'slurm'
            process.queue = 'general'
            singularity.enabled = true
            singularity.cacheDir = '/shared/containers'
        }

        test {
            params.samples = 'data/samples.csv'
            process.cpus = 1
            process.memory = '2.GB'
        }
    }
    ```

=== "Before"

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
    }

    // Container settings - disabled for mock processes
    docker.enabled = false

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

### 3.3. Use Different Profiles

Switch environments with one flag:

```bash
# Local execution (default)
nextflow run main.nf -profile standard

# With Docker containers (for real bioinformatics tools)
nextflow run main.nf -profile docker

# On HPC cluster
nextflow run main.nf -profile cluster

# Quick test with minimal resources
nextflow run main.nf -profile test
```

**The workflow code never changes.** Only the configuration adapts.

### 3.4. View the Execution Reports

The configuration enables automatic report generation:

```bash
ls results/*.html results/*.txt
```

```console title="Output"
results/report.html
results/timeline.html
results/trace.txt
```

Open `results/timeline.html` in a browser to see:

- Visual timeline of when each task ran
- Parallelization patterns
- Bottlenecks and idle time

Open `results/trace.txt` to see detailed metrics:

```bash
head results/trace.txt
```

```console title="Output"
task_id  hash       name              status     exit  duration  realtime
1        3a/4b5c6d  FASTQC (sample_01)  COMPLETED  0     1.1s      1s
2        8e/9f0a1b  FASTP (sample_01)   COMPLETED  0     2.1s      2s
...
```

### Takeaway

Configuration files let you:

- Set default parameters and resources
- Override settings for specific processes
- Create profiles for different environments
- Enable automatic execution reports

### What's next?

Review what you've learned and explore next steps.

---

## 4. Side-by-Side Comparison

### 4.1. The Same Task, Two Approaches

**Bash Script:**

```bash
./process_all.sh
```

- ❌ Sequential execution
- ❌ No resume capability
- ❌ Environment-dependent
- ❌ No execution tracking
- ✅ Simple to understand initially

**Nextflow:**

```bash
nextflow run main.nf -resume
```

- ✅ Automatic parallelization
- ✅ Built-in resume capability
- ✅ Container-based software isolation
- ✅ Comprehensive provenance tracking
- ⚠️ Initial learning curve

### 4.2. When to Make the Switch

**Use Nextflow when:**

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

---

## 5. Conclusion

### 5.1. What You've Learned

You started by building a bash script and experienced its limitations:

1. ⏱️ **Sequential processing** - Samples waited unnecessarily
2. 💥 **No resume capability** - Lost progress to failures
3. 📦 **Environment dependencies** - "Works on my machine" problems
4. 📊 **No provenance** - Couldn't track or optimize

Then you saw how Nextflow solves each problem:

1. ✅ **Automatic parallelization** - 3× faster with zero code changes
2. ✅ **Built-in caching** - Resume from failures instantly
3. ✅ **Container isolation** - Perfect reproducibility
4. ✅ **Comprehensive tracking** - Complete provenance built-in

### 5.2. Next Steps

**Continue learning:**

- [Hello Nextflow](../hello_nextflow/index.md) - Full introduction to Nextflow concepts
- [Groovy Essentials](./groovy_essentials.md) - Data manipulation and advanced patterns
- [nf-core](https://nf-co.re/) - Community best practices and ready-made pipelines

**Join the community:**

- [Nextflow Slack](https://www.nextflow.io/slack-invite.html) - Active community support
- [nf-core Slack](https://nf-co.re/join) - Pipeline-specific help

### 5.3. The Mindset Shift

Moving from scripts to workflows is a new way of thinking:

**Script thinking:**
> "How do I run these commands in sequence on my data?"

**Workflow thinking:**
> "What are the data dependencies? How will this scale? How can others reproduce this?"

Once you start thinking in terms of processes, channels, and data flow, you'll wonder how you ever managed with bash scripts.

Welcome to workflow management. Your future self will thank you.

---

## Quick Reference

### Nextflow vs Bash Commands

| Task | Bash | Nextflow |
|------|------|----------|
| Process samples | `for f in *.fq; do ...` | `PROCESS(channel)` (automatic) |
| Chain commands | `cmd1; cmd2; cmd3` | `CMD1(); CMD2(CMD1.out)` |
| Parallel execution | `cmd &` + management | Automatic |
| Resume after failure | Manual checkpoints | `-resume` flag |
| Specify resources | Hardcoded | `cpus = X, memory = 'Y.GB'` |

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

# Clean up work directory
nextflow clean -f

# List previous runs
nextflow log
```
