# Workflows in the AI era

If an AI agent can run your analysis on demand, or write the code for you in seconds, why do you still need to learn a workflow tool?

Bioinformatics has been through three eras of how analyses get expressed: commands typed into a terminal, scripts that captured those commands, and workflows that captured the scripts plus the engineering work needed to make them reproducible, portable, and maintainable. Each step was a better artefact for a team to live with, and each step happened because the previous form's limits started costing real time.

We are in a fourth era now, where an agent will produce any of those forms in seconds. The agent will write bash, Python, Nextflow, all of it. So the question is not "should I bother learning Nextflow when the agent could write it for me?" The question is what artefact you want to live with after the agent is done; what gets vetted at code review, tested, updated when a tool releases a new version, extended when the lab adds a new sample type, read by the colleague who inherits it next quarter.

A bash script can be that artefact, but it is the wrong shape for the job: every engineering property is buried inline, every change requires re-deriving correctness, and a maintainer six months from now has to reconstruct what the original author (human or otherwise) was thinking. A Nextflow workflow encodes the same logic in a form built for that maintenance role, with declarative containers and resources, explicit data flow, and the engineering properties supplied by the workflow boundary rather than by individual lines an author had to remember to write.

This side quest answers the question by experience. You will see the same RNA-seq analysis as a bash script (the kind an agent might hand you), find the engineering limits of that form, then see it as a Nextflow workflow (which the same agent could just as easily produce) and watch each limit disappear. By the end you should have an artefact you would be willing to put your name on and maintain.

!!! note "Why bash?"

    Part 1 uses bash to stand in for any ad-hoc form: a chat-generated script, a quick Python equivalent, an R one-liner, a sequence of commands an agent ran for you. The limitations apply to all of them. The issue is not the language; it is running tools without a workflow boundary, regardless of who or what authored the code.

### Learning goals

By the end of this side quest, you'll be able to:

- Name the production-quality properties an ad-hoc analysis (yours, your colleague's, or an agent's) does not give you for free.
- Wrap the same analysis logic in a Nextflow process so that reproducibility, software tracking, parallelisation, resume, and portability come from the workflow boundary.
- Use the workflow output system (`publish:` + `output {}` + `-output-dir`) to produce a stable layout an operator or agent can consume.
- Articulate why the case for a workflow tool gets stronger, not weaker, when AI is writing the Nextflow.

### Prerequisites

You don't need any prior Nextflow experience for this side quest; it is a motivational tour. If you want a thorough introduction afterwards, head to [Hello Nextflow](../../hello_nextflow/index.md).

You'll need:

- The Codespace (or a local environment with Docker, conda or mamba, and Nextflow >= 25.10).
- About 90 minutes.

---

## 1. What a production-grade analysis needs

These are the properties any agent producing bash has to remember to write, every time, and the properties a workflow tool supplies structurally:

- **Reproducibility**: Same inputs produce identical outputs, every time. Results can be verified and published with confidence.
- **Software management**: Each task gets its own isolated environment. No dependency conflicts between tools, and a standardised approach across your whole pipeline.
- **Scalability**: Handles 3 samples or 3,000 without code changes.
- **Efficient parallelisation**: Independent tasks run simultaneously, so analysis completes in hours, not days.
- **Resource awareness**: Respects memory and CPU limits. No crashed jobs or killed processes.
- **Failure recovery**: Can resume from where it stopped. A single failure doesn't waste hours of completed work.
- **Portability**: Runs on laptop, cluster, or cloud with the same code.

These aren't nice-to-haves; they're requirements for serious computational research.

---

## 2. The artefact those properties live in

Whatever code achieves those properties, the code itself ends up somewhere: a script in your home directory, a Nextflow workflow in a repo, a chat session nobody saved. That somewhere is the artefact, and it has to outlive the conversation that produced it. It goes into version control. Someone reads it at code review. Someone tests it. When a tool ships a new version, when the lab adds a new sample type, when the pipeline starts producing slightly different numbers and you have to figure out why, you go back to the artefact.

Two artefacts can do exactly the same analysis but be very different to live with. One can have its parallelism, software pinning, and resume logic written into individual lines that someone has to find and update by hand. Another can have those same properties supplied by the artefact's structure, visible at a glance to anyone who reads it. The second is what a workflow tool gives you. The first is what most ad-hoc analyses, whether human-written or agent-generated, end up as.

The rest of this side quest builds the same RNA-seq analysis as both forms, so you can see the difference where it counts: in the artefact, not in the runtime.

---

## 3. What you'll build

An RNA-seq analysis pipeline that runs:

1. Quality control (FastQC)
2. Adapter trimming (fastp)
3. Transcript quantification (Salmon)

How this tutorial works:

| Part | Software | What you'll do                                       | What you'll learn                                    |
| ---- | -------- | ---------------------------------------------------- | ---------------------------------------------------- |
| 1    | Bash     | Build a pipeline, try to hit those quality standards | How much infrastructure code it takes                |
| 2    | Nextflow | Rebuild with a workflow manager                      | How the same standards are achieved with less effort |

By the end, you'll understand _why_ workflow managers exist, not from explanation but from hitting the limitations yourself.

---

## 4. Building a bash pipeline

In this part, you'll build an RNA-seq analysis pipeline using bash, aiming for the production-quality properties listed above.
We'll see how far we can get and where things get difficult.

---

### 4.1. Setup

#### 4.1.1. The scenario

You're a bioinformatician analysing RNA-seq data from a yeast gene expression study.
You have 3 samples today. Your PI mentions that 50 more samples are coming next week.

#### 4.1.2. Navigate to the working directory

```bash
cd side-quests/workflows_in_the_ai_era
```

#### 4.1.3. Explore the starter files

```bash
ls bash/
```

```console title="Output"
pipeline_parallel.sh    pipeline_sequential.sh  process_sample.sh
```

Three scripts, each with TODOs you'll fill in as you progress:

- `process_sample.sh`: single sample script
- `pipeline_sequential.sh`: multi-sample loop
- `pipeline_parallel.sh`: parallel version

You'll fill in these starter files as you progress through the tutorial.

#### 4.1.4. Examine the sample data

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

### 4.2. Kick off the tool install

Part 1 needs FastQC, fastp, and Salmon installed locally so the bash scripts can call them from `PATH`. Conda solves environments slowly, so kick off the install now and let it run in the background while you keep reading. We'll check on it before you actually need to run anything.

=== "Mamba (recommended)"

    ```bash
    mamba create -n rnaseq-bash fastqc fastp salmon -c bioconda -c conda-forge -y
    ```

=== "Conda"

    ```bash
    conda create -n rnaseq-bash fastqc fastp salmon -c bioconda -c conda-forge -y
    ```

This typically takes about 3 minutes. While it runs, work through the next section. The activation and version check come later, once you actually have a script to test.

---

### 4.3. Building your first script

Open the starter file `bash/process_sample.sh`. It has the structure ready; you just need to add the tool commands.

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

#### 4.3.1. Understanding the starter script

The script accepts three arguments: a sample ID and two URLs for paired-end FASTQ files. The `set -e` line means the script will stop immediately if any command fails.

The TODO comments mark where you'll add each analysis step:

1. **Download**: Get the raw sequencing data
2. **FastQC**: Check the quality of the raw reads
3. **fastp**: Trim adapters and filter low-quality reads
4. **Salmon**: Quantify transcript expression levels

This is a typical RNA-seq preprocessing workflow. Fill in each step below.

#### 4.3.2. Add the download step

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

#### 4.3.3. Add FastQC

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

#### 4.3.4. Add fastp

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

#### 4.3.5. Add the Salmon index download

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

!!! note "Conditional logic in scripts vs workflow managers"

    This `if [ ! -d ... ]` check is a common pattern in scripts: "only do this if it hasn't been done already." While conceptually simple, it's fragile, because it detects whether the directory exists rather than whether the download actually succeeded or the files are valid. Workflow managers handle this more robustly by tracking actual task outputs rather than checking for file or directory existence.

#### 4.3.6. Add Salmon quantification

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

#### 4.3.7. Activate the env and verify it's ready

The conda install you started in section 4.2 should be done by now. Switch back to that terminal (or the same one if you ran the install in this shell) and activate the environment:

```bash
mamba activate rnaseq-bash
```

!!! warning "Remember this step"

    Every time you open a new terminal, you must activate this environment again. Forget, and your script fails with `command not found`.

Verify the three tools are on your `PATH`:

```bash
fastqc --version
fastp --version
salmon --version
```

If any of those errors out, give the install another minute and try again. Once all three respond with versions, you are ready to run the script.

#### 4.3.8. Test your script

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

This works for one sample, but you have 3 samples and 50 more coming. Running this command manually for each one isn't practical.

#### 4.3.9. Takeaway

You now have a runnable script. There is no record of which conda environment solved for which tool versions, no resume on failure, no way to scale up without rewriting the loop, no provenance trace beyond your shell history. A bash script is what an agent might hand you for a quick analysis; everything past "it ran once on my laptop" is something the agent has to remember to write separately, and probably won't get right. Asking the same agent for a workflow instead is what fixes that, as Part 2 will show.

---

### 4.4. Processing multiple samples

Now that we can process one sample, we need to handle all samples from our CSV file. Open `bash/pipeline_sequential.sh`; it already has the loop structure that reads the CSV and iterates over each sample.

The starter script handles:

- Reading the CSV file and skipping the header row
- Parsing each line into sample ID and FASTQ URLs
- Creating output directories and downloading the shared Salmon index

Your task is to add the actual processing commands inside the loop, the same steps you just wrote, but using the loop variables (`$sample_id`, `$fastq_r1`, `$fastq_r2`) instead of the script arguments.

#### 4.4.1. Add processing logic to the loop

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

#### 4.4.2. Run and time it

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

#### 4.4.3. The problem: sequential execution

```
Timeline (sequential):
WT_REP1:          [====FastQC====][====fastp====][======Salmon======]
WT_REP2:                                                              [====FastQC====]...
```

These samples are completely independent. WT_REP2's analysis doesn't depend on WT_REP1's results at all. Yet WT_REP2 sits idle while WT_REP1 runs through all its steps. With 50 samples that's 50x the runtime for no good reason; your computer has multiple cores sitting unused.

---

### 4.5. Adding parallelisation

The sequential script works, but it's slow; each sample waits for the previous one to finish completely. Since samples are independent, they could run simultaneously.

Open `bash/pipeline_parallel.sh`. This version wraps the processing logic in a function called `process_sample()`. The starter script has everything except the parallelisation itself.

Your task is to make the loop run samples in parallel by:

1. Adding `&` after each function call to run it in the background
2. Adding `wait` after the loop to pause until all background jobs complete

#### 4.5.1. Add background execution

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

- `&` at the end of the function call runs it in the background.
- `wait` pauses until all background jobs complete.

#### 4.5.2. Run and compare

```bash
chmod +x bash/pipeline_parallel.sh
time ./bash/pipeline_parallel.sh
```

Notice the interleaved output, all samples running at once:

```
[WT_REP1] Starting...
[WT_REP2] Starting...
[RAP1_IAA_30M_REP1] Starting...
[WT_REP1] Complete!
[RAP1_IAA_30M_REP1] Complete!
[WT_REP2] Complete!
```

Faster, because all samples run concurrently.

#### 4.5.3. The hidden problem

What happens with 50 samples, or 500?

```bash
# This would launch 500 Salmon jobs simultaneously
# Each Salmon uses ~4GB RAM
# That's 500 x 4GB = 2TB RAM required
```

Your machine would crash from memory exhaustion. Bash's `&` has no concept of resource limits; it just launches everything at once. An agent producing bash here would need to add a job-throttling layer: 20-30 lines tracking running jobs, waiting at the limit, and starting new jobs as others complete. The same agent producing Nextflow doesn't have to write any of it; the executor handles throttling from the resource declarations on each process.

---

### 4.6. Adding aggregation

Real pipelines need a final aggregation step that collects QC metrics from all samples into a summary report. MultiQC does this, but integrating it with parallel bash processing is tricky.

#### 4.6.1. The problem

Your parallel script launches all samples with `&` and waits for completion with `wait`. But `wait` only tells you when jobs finish, not which ones succeeded. To run MultiQC, you need to:

1. Track which samples completed successfully
2. Collect the right output files from each tool
3. Handle partial failures gracefully

#### 4.6.2. A simple (fragile) approach

Add MultiQC after the `wait`:

```bash title="bash/pipeline_parallel.sh" linenums="64"
wait
echo "All samples complete. Running MultiQC..."

# Collect outputs - assumes everything succeeded
multiqc results/fastqc results/fastp results/salmon -o results/multiqc
```

This works if everything succeeds, but breaks silently on partial failures: MultiQC runs on whatever files exist, potentially giving you an incomplete report without warning.

#### 4.6.3. A robust approach (more code)

Proper aggregation requires tracking job status:

```bash
declare -A job_status
for sample_id in "${samples[@]}"; do
    process_sample "$sample_id" &
    pids[$sample_id]=$!
done

# Wait and track success/failure
for sample_id in "${!pids[@]}"; do
    if wait ${pids[$sample_id]}; then
        job_status[$sample_id]="success"
    else
        job_status[$sample_id]="failed"
        echo "WARNING: $sample_id failed"
    fi
done

# Only aggregate successful samples
successful_dirs=""
for sample_id in "${!job_status[@]}"; do
    if [[ ${job_status[$sample_id]} == "success" ]]; then
        successful_dirs+=" results/salmon/$sample_id"
    fi
done

multiqc results/fastqc results/fastp $successful_dirs -o results/multiqc
```

That's another 20+ lines of infrastructure code for something that should be simple: "run this after everything else finishes."

---

### 4.7. The pain points

We won't implement these, but consider what you'd need for a production-ready pipeline:

#### 4.7.1. Failure recovery

If Salmon fails on sample 47 of 50, what happens?

- With `set -e`, everything stops immediately.
- You have no record of which 46 samples completed successfully.
- To resume, you'd need to manually figure out what finished and what didn't.
- Implementing proper checkpoint logic would add 40+ lines of file-based state tracking.

#### 4.7.2. Reproducibility

When your colleague tries to run your script:

```
salmon: command not found
```

They need to replicate your exact conda environment, same tool versions, same dependencies.

Or worse, they have salmon 1.4.0 but you used 1.10.3. The script runs fine, but results differ silently. Months later, you can't reproduce your own analysis because conda updated something.

#### 4.7.3. Scaling to cluster

Your laptop worked fine for 3 samples. For 500 samples, you need the SLURM cluster:

```bash
#SBATCH --job-name=salmon
#SBATCH --cpus-per-task=4
```

Now you're rewriting the entire script with job arrays, dependency tracking, and cluster-specific syntax. Your collaborator's PBS cluster needs yet another rewrite. Each environment requires maintaining a separate codebase.

---

### 4.8. Takeaway

Checking progress against the production-quality properties:

- :warning: **Reproducibility**: Partial. Conda environment helps, but requires discipline to maintain.
- :warning: **Software management**: Partial. All tools share one environment, so conflicts are possible. Setup is manual and must be replicated by every user.
- :x: **Scalability**: Limited. Works on one machine until you hit resource limits. Scaling to cluster or cloud requires a complete rewrite.
- :warning: **Efficient parallelisation**: Partial. Works, but no resource limits.
- :x: **Resource awareness**: Missing. Would need 20-30 lines of job-limiting code.
- :x: **Failure recovery**: Missing. Would need 40+ lines of checkpoint logic.
- :x: **Portability**: Missing. Complete rewrite for each cluster type.

We achieved the basics, but the production-quality properties require significant infrastructure code that has nothing to do with our science.

The fundamental problem: scripts mix _what_ to compute with _how_ to compute it. Every quality requirement adds more infrastructure code, and this would be true whether the script came from bash, Python, or a chat session with an agent. This is where workflow managers come in.

---

## 5. Building a Nextflow pipeline

Now you'll rebuild the same pipeline using Nextflow.
At each step, we'll contrast with the scripting approach from Part 1.

!!! note "This isn't a Nextflow tutorial"

    We're skipping over a lot of Nextflow detail here to focus on illustrating how workflow managers solve the problems from Part 1. If you're sold on the idea and want to learn Nextflow properly, head to [Hello Nextflow](../../hello_nextflow/index.md) for a thorough introduction.

---

### 5.1. What is Nextflow?

Nextflow is a workflow manager. Instead of writing imperative scripts that say "do this, then do that," you _declare_ what processes exist and how data flows between them. Nextflow handles everything else: the parallelisation, scheduling, error handling, and resource management you'd otherwise write yourself.

#### 5.1.1. Key concepts

| Concept    | What it is                                         | Bash equivalent                      |
| ---------- | -------------------------------------------------- | ------------------------------------ |
| `process`  | A unit of work (like running FastQC on one sample) | A function or script                 |
| `channel`  | A queue of data flowing between processes          | Variables passed between commands    |
| `workflow` | How processes connect together                     | The order of commands in your script |

The key difference: in bash, you explicitly manage data flow with variables and file paths. In Nextflow, you declare what each process needs, and Nextflow figures out the execution order automatically.

#### 5.1.2. Software management

In Part 1, you installed FastQC, fastp, and Salmon with conda, hoping dependencies wouldn't conflict and documenting versions manually.

With Nextflow, each process declares its own software requirements. The tools are configured automatically at runtime with support for whichever software packaging tool you prefer. Your colleague runs the same pipeline and gets the exact same software environment based on the process definition, not their system.

---

### 5.2. Your first process: FastQC

In bash, you wrote FastQC commands inside loops, managed file paths with variables, and handled parallelisation with `&`. In Nextflow, each tool runs inside a `process`, a self-contained unit that declares its inputs, outputs, and command. You don't write loops; Nextflow runs the process once per input item automatically.

Copy this complete process into `nextflow/modules/fastqc.nf`:

```groovy title="nextflow/modules/fastqc.nf"
process FASTQC {
    tag "$id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*.html"), emit: html
    tuple val(id), path("*.zip"), emit: zip

    script:
    """
    fastqc --quiet --threads 2 ${reads}
    """
}
```

#### 5.2.1. Understanding the process

Each part has a purpose:

- `tag`: Labels each job with the sample ID (visible in logs).
- `container`: The Docker image containing FastQC (replaces your conda setup).
- `input`: Declares expected data, here a tuple with sample ID (`val(id)`) and file paths (`path(reads)`).
- `output`: Declares produced files with named channels (`emit: html`) for downstream processes.
- `script`: The actual command, nearly identical to your bash version.

You're not writing any loop logic, file existence checks, or error handling.

!!! tip "Contrast with the agent's script"

    The one line `container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'` handles all software dependencies. The version is locked forever. Your colleague, your cluster, and your cloud all get the exact same FastQC, whether the agent generated this process for you or you wrote it by hand.

#### 5.2.2. Call FASTQC in main.nf

Open `nextflow/main.nf` and add the FASTQC call:

=== "After"

    ```groovy title="nextflow/main.nf" linenums="27" hl_lines="1"
    FASTQC(ch_samples)
    ```

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="27"
    // TODO: Call FASTQC with ch_samples
    ```

#### 5.2.3. Test it

```bash
cd nextflow
nextflow run main.nf -output-dir results
```

```console title="Output"
N E X T F L O W  ~  version 25.10.4

executor >  local (3)
[a1/b2c3d4] FASTQC (WT_REP1)           [100%] 3 of 3 ✔
```

All 3 samples ran in parallel automatically. In Part 1, you wrote `&` and `wait` and worried about resource limits. Nextflow figures out optimal parallelisation from your process definition alone, with no extra code from you or the agent.

---

### 5.3. Adding fastp

Now add fastp following the same pattern. The key difference: fastp's output (trimmed reads) needs to flow to Salmon. In bash, you managed this with file paths like `results/fastp/${sample_id}_trimmed_R1.fastq.gz`. In Nextflow, you declare the output and let channels handle the data flow.

Open `nextflow/modules/fastp.nf` and fill in the input, output, and script blocks.

#### 5.3.1. Complete fastp.nf

Pay attention to the `emit: reads` on the trimmed reads output. This names the output channel so Salmon can consume it.

=== "After"

    ```groovy title="nextflow/modules/fastp.nf" hl_lines="6 9-11 15-21"
    process FASTP {
        tag "$id"
        container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

        input:
        tuple val(id), path(reads)

        output:
        tuple val(id), path("*_trimmed*.fastq.gz"), emit: reads
        tuple val(id), path("*.json"), emit: json
        tuple val(id), path("*.html"), emit: html

        script:
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${id}_trimmed_R1.fastq.gz \\
            --out2 ${id}_trimmed_R2.fastq.gz \\
            --json ${id}.fastp.json \\
            --html ${id}.fastp.html \\
            --thread 4
        """
    }
    ```

=== "Before"

    ```groovy title="nextflow/modules/fastp.nf"
    process FASTP {
        tag "$id"
        container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

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

#### 5.3.2. Call FASTP in main.nf

Now connect FASTP to the sample channel in the workflow:

=== "After"

    ```groovy title="nextflow/main.nf" linenums="29" hl_lines="1"
    FASTP(ch_samples)
    ```

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="29"
    // TODO: Call FASTP with ch_samples
    ```

!!! tip "Contrast with the agent's script"

    In Part 1, all tools shared one conda environment and conflicts were a constant risk. Here, fastp uses a completely different container than FastQC. They could require incompatible Python versions or conflicting libraries; doesn't matter. Each process gets its own isolated environment, regardless of what the agent picked when it set up your conda env.

---

### 5.4. Connecting to Salmon

This is where Nextflow's data flow model pays off. In bash, you manually coordinated file paths; Salmon needed to know where fastp wrote its output. In Nextflow, you declare that Salmon needs fastp's output, and the framework handles the rest.

Salmon needs two inputs:

1. **Trimmed reads**: Different for each sample (from FASTP).
2. **Reference index**: Same for all samples (from UNTAR).

This is a common pattern: per-sample data combined with a shared reference. In bash, you'd check if the index exists with `if [ ! -d ... ]`. In Nextflow, the channel system handles this automatically.

#### 5.4.1. Complete salmon.nf

The UNTAR process (which extracts the pre-built Salmon index) is already complete. Your task is to fill in SALMON_QUANT, which has two input declarations:

=== "After"

    ```groovy title="nextflow/modules/salmon.nf" linenums="17" hl_lines="9-10 13 17-23"
    process SALMON_QUANT {
        tag "$id"
        container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

        cpus 4
        memory '8.GB'

        input:
        tuple val(id), path(reads)
        path index

        output:
        tuple val(id), path("${id}"), emit: results

        script:
        """
        salmon quant \\
            --index ${index} \\
            --libType A \\
            --mates1 ${reads[0]} \\
            --mates2 ${reads[1]} \\
            --output ${id} \\
            --threads ${task.cpus}
        """
    }
    ```

=== "Before"

    ```groovy title="nextflow/modules/salmon.nf" linenums="17"
    process SALMON_QUANT {
        tag "$id"
        container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

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

- `${task.cpus}`: Nextflow makes declared resources available as variables. In bash, you hardcoded `--threads 2`. Here, you declare resources once (`cpus 4`) and reference them in the script. Change the declaration, and the script adapts automatically.
- Two inputs: the process receives the trimmed reads tuple and the index path separately. In bash, you coordinated these with file paths. Here, the data flow is explicit.

#### 5.4.2. Wire it up in main.nf

Now connect the processes. SALMON_QUANT needs two arguments. Then add `publish:` and `output {}` blocks so each tool's results land in a predictable layout.

=== "After"

    ```groovy title="nextflow/main.nf" linenums="30" hl_lines="1 3-9 12-18"
    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    fastp_reads = FASTP.out.reads
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
    salmon_quant = SALMON_QUANT.out.results
    }

    output {
        fastqc_html { path 'fastqc' }
        fastqc_zip { path 'fastqc' }
        fastp_reads { path 'fastp' }
        fastp_json { path 'fastp' }
        fastp_html { path 'fastp' }
        salmon_quant { path 'salmon' }
    }
    ```

=== "Before"

    ```groovy title="nextflow/main.nf" linenums="30"
    // TODO: Call SALMON_QUANT with FASTP output and the index
    // Hint: Use FASTP.out.reads and UNTAR.out.index.first()

    publish:
    // TODO: Publish the outputs you want kept after the run.
    }

    output {
        // TODO: Configure output paths for each published channel.
    }
    ```

Two things to notice:

- `FASTP.out.reads`: the trimmed reads channel (one item per sample).
- `UNTAR.out.index.first()`: `.first()` converts the channel to a single value that gets reused for every sample. Without it, Nextflow would try to pair each sample with a separate index item.

The `publish:` block names the channels you want kept; the top-level `output {}` block tells Nextflow what subdirectory to put each one in. The actual root directory comes from `-output-dir` on the command line.

#### 5.4.3. Run and watch

Run the complete pipeline:

```bash
nextflow run main.nf -output-dir results
```

Watch what happens. Nextflow automatically determines the execution order from the data flow you defined:

- FastQC and FASTP can run in parallel; both only need the raw reads, no dependency between them.
- SALMON_QUANT must wait for FASTP; it needs the trimmed reads output.
- Each sample runs independently; sample 2's Salmon can start as soon as sample 2's fastp finishes, even if sample 1 is still running.

!!! tip "Contrast with the agent's script"

    In Part 1, the bash version needed `&` and `wait` and worried about memory limits at 500 samples. The same agent producing this Nextflow process gets parallelisation and resource-aware throttling for free; producing the bash version, it has to remember to write the throttling itself, and probably won't.

---

### 5.5. Aggregation with MultiQC

Real pipelines don't just run processes per sample; they often need to aggregate results across all samples. MultiQC collects QC metrics from FastQC, fastp, Salmon, and other tools into a single summary report.

In bash, aggregation is awkward. You'd need to:

1. Wait for all parallel jobs to complete
2. Track which output files exist
3. Handle the case where some samples failed
4. Construct file lists for the aggregation tool

With channels, this is straightforward. Add a MULTIQC process:

```groovy title="nextflow/modules/multiqc.nf"
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data

    script:
    """
    multiqc . --filename multiqc_report
    """
}
```

Then wire it up to collect outputs from all processes:

```groovy title="nextflow/main.nf"
// Collect all QC outputs for MultiQC
ch_multiqc = FASTQC.out.zip
    .map { id, files -> files }
    .mix(FASTP.out.json.map { id, files -> files })
    .mix(SALMON_QUANT.out.results.map { id, files -> files })
    .collect()

MULTIQC(ch_multiqc)
```

Add `multiqc_report = MULTIQC.out.report` to the workflow's `publish:` block and `multiqc_report { path 'multiqc' }` to the `output {}` block.

The `.collect()` operator waits for all upstream processes to complete, then passes everything to MultiQC as a single batch. Nextflow tracks which files to collect automatically.

!!! tip "Contrast with the agent's script"

    In Part 1, you'd need to track job completion status, build file lists manually, and handle partial failures. Here, the channel operators handle all of that. The `.collect()` naturally waits for all inputs, and `.mix()` combines outputs from different processes; an agent generating bash from scratch would have to reinvent each of these every time.

---

### 5.6. Resume and caching

In Part 1, failure recovery was a pain point. If sample 47 failed, you had no record of what completed, and implementing checkpoint logic would take 40+ lines of state tracking.

Run the Nextflow pipeline, then modify something and run with `-resume`:

```bash
nextflow run main.nf -output-dir results -resume
```

```console title="Output"
[a1/b2c3d4] UNTAR              [100%] 1 of 1, cached: 1 ✔
[e5/f6g7h8] FASTQC (WT_REP1)   [100%] 3 of 3, cached: 3 ✔
[i9/j0k1l2] FASTP (WT_REP1)    [100%] 3 of 3, cached: 3 ✔
[m3/n4o5p6] SALMON_QUANT       [100%] 3 of 3 ✔  <- Only this re-ran
```

Nextflow automatically tracks what completed successfully. Failed tasks can be fixed and re-run without repeating successful work. One flag (`-resume`) replaces 40+ lines of bash checkpoint logic you'd have to write and debug, whether you wrote it or asked an agent to.

---

### 5.7. Configuration profiles

In Part 1, scaling to clusters required rewriting the entire script with SLURM job arrays and cluster-specific syntax. Your collaborator on PBS would need yet another version.

With Nextflow, the `nextflow.config` file lets you run the same workflow anywhere:

```bash
# On your laptop
nextflow run main.nf -output-dir results -profile docker

# On SLURM cluster
nextflow run main.nf -output-dir results -profile slurm

# On AWS Batch
nextflow run main.nf -output-dir results -profile awsbatch
```

Same workflow code. Same scientific results. Different infrastructure. You write your analysis once, and configuration profiles adapt it to any environment.

---

### 5.8. Takeaway

Checking progress against the same properties:

- :white_check_mark: **Reproducibility**: Containers lock exact tool versions. Same results everywhere.
- :white_check_mark: **Software management**: Each process gets its own isolated container. No conflicts, no manual setup.
- :white_check_mark: **Scalability**: Same code for 3 or 3,000 samples. Manages resources on one machine, or distributes to cluster or cloud with a config change.
- :white_check_mark: **Efficient parallelisation**: Automatic from data flow declarations.
- :white_check_mark: **Resource awareness**: Declarative `cpus` and `memory` per process.
- :white_check_mark: **Failure recovery**: `-resume`. One flag.
- :white_check_mark: **Portability**: Same code, different `-profile`.

Every property we struggled with in Part 1 is built into Nextflow. You didn't write any infrastructure code; you just declared what each process needs and produces.

The fundamental difference: scripts mix _what_ to compute with _how_ to compute it. Workflow managers separate these concerns. You (or your agent) declare the science, the framework handles the infrastructure.

---

## 6. Why workflows still matter when AI does the writing

Look back at the four eras the introduction named. Commands gave way to scripts, and scripts gave way to workflows, because each step was a better artefact for a team to live with: easier to vet, easier to validate, easier to maintain when tools changed and labs grew. The workflow form did not appear because writing scripts was hard. It appeared because keeping scripts around was.

We are in the fourth era now. An agent producing bash and the same agent producing Nextflow take roughly the same number of seconds. The bottleneck has moved. It is no longer how long it takes to write code; it is how long it takes a team to verify, ship, and maintain what was written. The artefact's role has not changed.

The case for the workflow tool gets stronger, not weaker, the more code AI is producing. Volumes of generated code that nobody can read or vet are a liability. A workflow puts the engineering structure where humans can see it: container pins as one-line declarations, parallelism as data flow, resume as a property of the cache rather than a property the author remembered to write. That is what makes review tractable, tests writable, and maintenance possible. AI accelerates the authoring; the workflow tool makes the operate-debug-maintain loop tractable.

You still read what the AI produced. You still verify container pins, validate cache behaviour, lock the contract with tests; these are the durable skills. The [`nf-test` side quest](../nf_test/index.md) is the right next step for locking what got authored.

The answer to "why workflows when AI can write the code for me" is that the artefact has to outlive the conversation that made it.

---

## Summary

You started with a question: _why learn a whole new framework when an AI agent can do my analysis or build my pipeline for me?_

You then built the same RNA-seq analysis twice. Once as a bash script of the kind an agent might produce on demand, where every production-quality property cost extra infrastructure code that someone, human or otherwise, had to remember to write. Once as a Nextflow workflow (which the same agent could just as easily produce), where reproducibility, software tracking, scalability, parallelisation, resource awareness, failure recovery, and portability came from the workflow boundary itself.

The pipeline is the durable artefact. The agent will write either form for you on demand. The form that makes the result trustworthy across time and across team members is the one shaped for the maintenance role: vetted at code review, validated by tests, updated when tools change, read by the next person who has to live with it.

The artefact has to outlive the conversation that made it.

### Key patterns

Pin every tool with a container, not a manual install:

```groovy
container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
```

Capture outputs with the workflow output system, so the run produces a stable layout:

```groovy
publish:
fastqc_html = FASTQC.out.html

}

output {
    fastqc_html { path 'fastqc' }
}
```

Pick the run's output root at the command line:

```bash
nextflow run main.nf -output-dir results
```

Resume after a failure with one flag:

```bash
nextflow run main.nf -output-dir results -resume
```

Move the same workflow between laptop, cluster, and cloud with a profile:

```bash
nextflow run main.nf -output-dir results -profile slurm
```

### Additional resources

- [Hello Nextflow](../../hello_nextflow/index.md): the thorough introduction to Nextflow if you want to keep going.
- [Testing with nf-test](../nf_test/index.md): lock the contract on the artefact, whoever or whatever wrote it.
- [Nextflow documentation](https://www.nextflow.io/docs/latest/): reference for the workflow output system, profiles, and `nextflow log`.
- [nf-core](https://nf-co.re/): community-maintained, production-grade pipelines that already use the patterns you saw here.
- [Seqera AI](https://seqera.io/ask-ai/): an AI assistant trained on Nextflow and nf-core resources for help drafting and debugging workflows.

### What's next?

Return to the [menu of Side Quests](../index.md) or click the button in the bottom right of the page to move on to the next topic.
