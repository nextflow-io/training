# Part 2: Single-sample implementation

In this part of the course, we're going to write the simplest possible workflow that wraps all the commands we ran in Part 1 to automate running them, and we'll just aim to process one sample at a time.

We'll do this in three stages:

1. Write a single-stage workflow that runs the initial QC step
2. Add adapter trimming and post-trimming QC
3. Add alignment to the reference genome

!!! warning "Prerequisite"

    You must work through Part 1 of the course before starting this lesson.
    Specifically, working through section 1.2.3 creates the genome index file (`data/genome_index.tar.gz`) required for the alignment step in this lesson.

As a starting point, we provide you with a workflow file, `rnaseq.nf`, that outlines the main parts of the workflow, as well as four module files in the `modules/` directory — `fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf`, and `multiqc.nf` — that outline the structure of each process.
These files are not functional; their purpose is just to serve as scaffolds for you to fill in with the interesting parts of the code.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Module INCLUDE statements

/*
 * Pipeline parameters
 */

// Primary input

workflow {

    // Create input channel

    // Call processes

}
```

---

## 1. Write a single-stage workflow that runs the initial QC

This first step focuses on the basics: loading a FASTQ file and running quality control on it.

Recall the `fastqc` command from [Part 1](01_method.md):

```bash
fastqc <reads>
```

The command takes a FASTQ file as input and produces a quality control report as a `.zip` archive and an `.html` summary.
The container URI was `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

We're going to take this information and wrap it in Nextflow in two stages:

1. Set up the input
2. Write the QC process and call it in the workflow

### 1.1. Set up the input

We need to declare an input parameter with a default value and create an input channel.

#### 1.1.1. Add an input parameter declaration

In `rnaseq.nf`, under the `Pipeline parameters` section, declare a parameter called `reads` with a default path pointing to one of the test FASTQ files.

=== "After"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

#### 1.1.2. Set up the input channel

In the workflow block, create an input channel from the parameter value using the `.fromPath` channel factory (as used in [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "After"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="3-4"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Call processes

    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        // Create input channel

        // Call processes

    }
    ```

Now we need to create the process to run QC on this input.

### 1.2. Write the QC process and call it in the workflow

We need to fill in the process definition in the module file, import it into the workflow using an include statement, and call it on the input.

#### 1.2.1. Fill in the module for the QC process

Open `modules/fastqc.nf` and examine the outline of the process definition.
You should recognize the main structural elements; if not, consider reading through [Hello Nextflow](../../hello_nextflow/01_hello_world.md) for a refresher.

Go ahead and fill in the process definition by yourself using the information provided above, then check your work against the solution in the "After" tab below.

=== "Before"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container
        publishDir

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 9 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
        publishDir "results/fastqc", mode: 'symlink'

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

The `publishDir` directive tells Nextflow to publish the process outputs to a results directory.
We use `mode: 'symlink'` to create symbolic links instead of copying files.
The `simpleName` accessor strips all extensions from the filename, so `ENCSR000COQ1_1.fastq.gz` becomes `ENCSR000COQ1_1`.
We use the `emit:` syntax to assign names to each output channel, which will be useful later for wiring outputs between processes.

!!! note

    Even though the data files we're using here are very small, in genomics they can get very large.
    For the purposes of demonstration in the teaching environment, we're using the 'symlink' publishing mode to avoid unnecessary file copies.
    You shouldn't do this in your final workflows, since you'll lose results when you clean up your `work` directory.

Once you've completed this, the process is complete.
To use it in the workflow, you'll need to import the module and add a process call.

#### 1.2.2. Include the module

In `rnaseq.nf`, add an `include` statement to make the process available to the workflow:

=== "After"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    ```

The process is now available in the workflow scope.

#### 1.2.3. Call the QC process on the input

Add a call to `FASTQC` in the workflow block, passing the input channel as an argument.

=== "After"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="6-7"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Initial quality control
        FASTQC(read_ch)

    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Call processes

    }
    ```

### 1.3. Run the workflow

At this point, we have a one-step QC workflow that should be fully functional.

We could use the `--reads` parameter to specify an input from the command line, but during development we can be lazy and just use the default we set up in the parameter declaration.

```bash
nextflow run rnaseq.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

This should run very quickly if you worked through Part 1 and have already pulled the container.
If you skipped it, Nextflow will pull the container for you; you don't have to do anything for it to happen, but you may need to wait up to a minute.

You can check the outputs under `results/fastqc` as specified by the `publishDir` directive.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

### Takeaway

You know how to create a module containing a process, import it into a workflow, call it with an input channel, and publish the results.

### What's next?

Add adapter trimming with post-trimming QC as a second step in the workflow.

---

## 2. Add adapter trimming and post-trimming QC

Now that we have the initial QC in place, we can add the adapter trimming step with its built-in post-trimming QC.

Recall the `trim_galore` command from [Part 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

The command trims adapters from a FASTQ file and runs FastQC on the trimmed output.
It produces trimmed reads, a trimming report, and FastQC reports for the trimmed reads.
The container URI was `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

We just need to write the process definition, import it, and call it in the workflow.

### 2.1. Write the trimming process and call it in the workflow

#### 2.1.1. Fill in the module for the trimming process

Open `modules/trim_galore.nf` and examine the outline of the process definition.

Go ahead and fill in the process definition by yourself using the information provided above, then check your work against the solution in the "After" tab below.

=== "Before"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container
        publishDir

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 9 12 15 16 17 21"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
        publishDir "results/trimming", mode: 'symlink'

        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    }
    ```

This process has three named outputs: the trimmed reads that feed into the alignment step, the trimming report, and the post-trimming FastQC reports.
The `--fastqc` flag tells Trim Galore to automatically run FastQC on the trimmed output.

#### 2.1.2. Include the module

Update `rnaseq.nf` to import the new module:

=== "After"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

#### 2.1.3. Call the trimming process on the input

Add the process call in the workflow block:

=== "After"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="9-10"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Initial quality control
        FASTQC(read_ch)

        // Adapter trimming and post-trimming QC
        TRIM_GALORE(read_ch)

    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Initial quality control
        FASTQC(read_ch)

    }
    ```

### 2.2. Run the workflow

```bash
nextflow run rnaseq.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

This should run very quickly too, since we're running on such a small input file.

You can find the outputs under `results/trimming` as specified by the `publishDir` directive.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

### Takeaway

You know how to add a second processing step that runs independently on the same input, producing multiple named outputs.

### What's next?

Add the alignment step that chains off the trimmed reads output.

---

## 3. Add alignment to the reference genome

Finally we can add the genome alignment step using HISAT2.

Recall the alignment command from [Part 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

The command aligns reads to a reference genome and converts the output to BAM format.
It requires a pre-built genome index archive and produces a BAM file and an alignment summary log.
The container URI was `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

This process requires an additional input (the genome index archive), so we need to set that up first, then write and wire the process.

### 3.1. Set up the inputs

#### 3.1.1. Add a parameter for the genome index

Add a parameter declaration for the genome index archive in `rnaseq.nf`:

=== "After"

    ```groovy title="rnaseq.nf" linenums="9" hl_lines="5-6"
    params {
        // Primary input
        reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

        // Reference genome archive
        hisat2_index_zip: Path = "data/genome_index.tar.gz"
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="9"
    params {
        // Primary input
        reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
    }
    ```

### 3.2. Write the alignment process and call it in the workflow

#### 3.2.1. Fill in the module for the alignment process

Open `modules/hisat2_align.nf` and examine the outline of the process definition.

Go ahead and fill in the process definition by yourself using the information provided above, then check your work against the solution in the "After" tab below.

=== "Before"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container
        publishDir

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 9 12 13 16 17 21 22 23 24"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
        publishDir "results/align", mode: 'symlink'

        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    }
    ```

This process takes two inputs: the reads and the genome index archive.
The script block first extracts the index from the archive, then runs the HISAT2 alignment piped into `samtools view` to convert the output to BAM format.
The `simpleName` accessor on `index_zip` extracts the base name of the archive (`genome_index`) to use as the index prefix.

#### 3.2.2. Include the module

Update `rnaseq.nf` to import the new module:

=== "After"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

#### 3.2.3. Call the alignment process

The trimmed reads are in the `TRIM_GALORE.out.trimmed_reads` channel output by the previous step.
We use `#!groovy file(params.hisat2_index_zip)` to provide the genome index archive.

=== "After"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="12-13"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Initial quality control
        FASTQC(read_ch)

        // Adapter trimming and post-trimming QC
        TRIM_GALORE(read_ch)

        // Alignment to a reference genome
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        // Create input channel from a file path
        read_ch = channel.fromPath(params.reads)

        // Initial quality control
        FASTQC(read_ch)

        // Adapter trimming and post-trimming QC
        TRIM_GALORE(read_ch)

    }
    ```

### 3.3. Run the workflow

```bash
nextflow run rnaseq.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

You can find the outputs under `results/align` as specified by the `publishDir` directive.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

This completes the basic processing we need to apply to each sample.

_We'll add MultiQC report aggregation in Part 3, after we've made the workflow accept multiple samples at a time._

---

### Takeaway

You know how to wrap all the core steps to process single-end RNAseq samples individually.

### What's next?

Learn how to modify the workflow to process multiple samples in parallel, aggregate QC reports across all steps for all samples, and enable running the workflow on paired-end RNAseq data.
