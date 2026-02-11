# Part 2: Single-sample implementation

In this part of the course, we're going to write the simplest possible workflow that wraps all the commands we ran in Part 1 to automate running them, and we'll just aim to process one sample at a time.

!!! warning "Prerequisite"

    You must work through [Part 1: Method overview](./01_method.md) before starting this lesson.
    Specifically, working through section 1.2.3 creates the genome index file (`data/genome_index.tar.gz`) required for the alignment step in this lesson.

## Assignment

In this part of the course, we're going to develop a workflow that does the following:

1. Run quality control (FastQC) on input reads
2. Trim adapters and run post-trimming QC (Trim Galore)
3. Align trimmed reads to a reference genome (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

This automates the steps from the first section of [Part 1: Method overview](./01_method.md#1-single-sample-processing), where you ran these commands manually in their containers.

As a starting point, we provide you with a workflow file, `rnaseq.nf`, that outlines the main parts of the workflow, as well as four module files in the `modules/` directory (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf`, and `multiqc.nf`) that outline the structure of each process.

??? full-code "Scaffold files"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Module INCLUDE statements

    /*
     * Pipeline parameters
     */

    // Primary input

    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }

    output {
        // Configure publish targets
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

These files are not functional; their purpose is just to serve as scaffolds for you to fill in with the interesting parts of the code.

## Lesson plan

In order to make the development process more educational, we've broken this down into three stages:

1. **Write a single-stage workflow that runs the initial QC step.**
   This covers setting up a CLI parameter, creating an input channel, writing a process module, and configuring output publishing.
2. **Add adapter trimming and post-trimming QC.**
   This introduces chaining processes by connecting one process's output to another's input.
3. **Add alignment to the reference genome.**
   This covers handling additional reference inputs and working with compressed archives.

Each step focuses on a specific aspect of workflow development.

!!! tip

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Write a single-stage workflow that runs the initial QC

This first step focuses on the basics: loading a FASTQ file and running quality control on it.

Recall the `fastqc` command from [Part 1](01_method.md):

```bash
fastqc <reads>
```

The command takes a FASTQ file as input and produces a quality control report as a `.zip` archive and an `.html` summary.
The container URI was `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

We're going to take this information and wrap it in Nextflow in three stages:

1. Set up the input
2. Write the QC process and call it in the workflow
3. Configure the output handling

### 1.1. Set up the input

We need to declare an input parameter, create a test profile to provide a convenient default value, and create an input channel.

#### 1.1.1. Add an input parameter declaration

In `rnaseq.nf`, under the `Pipeline parameters` section, declare a parameter called `reads` with the type `Path`.

=== "After"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        input: Path
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

That sets up the CLI parameter, but we don't want to type out the file path every time we run the workflow during development.
There are multiple options for providing a default value; here we use a test profile.

#### 1.1.2. Create a test profile with a default value in `nextflow.config`

A test profile provides convenient default values for trying out a workflow without specifying inputs on the command line.
This is a common convention in the Nextflow ecosystem (see [Hello Config](../../hello_nextflow/06_hello_config.md) for more detail).

Add a `profiles` block to `nextflow.config` with a `test` profile that sets the `reads` parameter to one of the test FASTQ files.

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Here, we're using `#!groovy ${projectDir}`, a built-in Nextflow variable that points to the directory where the workflow script is located.
This makes it easy to reference data files and other resources without hardcoding absolute paths.

The parameter now has a convenient default. Next, we need to create a channel from it.

#### 1.1.3. Set up the input channel

In the workflow block, create an input channel from the parameter value using the `.fromPath` channel factory (as used in [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "After"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
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

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

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

The `simpleName` accessor strips all extensions from the filename, so `ENCSR000COQ1_1.fastq.gz` becomes `ENCSR000COQ1_1`.
We use the `emit:` syntax to assign names to each output channel, which will be useful for wiring outputs into the publish block.

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

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Initial quality control
        FASTQC(read_ch)

        publish:
        // Declare outputs to publish
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

The workflow now loads the input and runs the QC process on it.
Next, we need to configure how the output is published.

### 1.3. Configure the output handling

We need to declare which process outputs to publish and specify where they should go.

#### 1.3.1. Declare outputs in the `publish:` section

The `publish:` section inside the workflow block declares which process outputs should be published.
Assign the outputs of `FASTQC` to named targets.

=== "After"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Declare outputs to publish
    }
    ```

Now we need to tell Nextflow where to put the published outputs.

#### 1.3.2. Configure the output targets in the `output {}` block

The `output {}` block sits outside the workflow and specifies where each named target is published.
Configure both targets to publish into a `fastqc/` subdirectory.

=== "After"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configure publish targets
    }
    ```

!!! note

    By default, Nextflow publishes output files as symbolic links, which avoids unnecessary duplication.
    Even though the data files we're using here are very small, in genomics they can get very large.
    Symlinks will break when you clean up your `work` directory, so for production workflows you may want to override the default publish mode to `'copy'`.

### 1.4. Run the workflow

At this point, we have a one-step QC workflow that should be fully functional.

We run with `-profile test` to use the default value set up in the test profile, avoiding the need to write the path on the command line.

```bash
nextflow run rnaseq.nf -profile test
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

You can check the outputs in the results directory.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

The QC reports for the sample are now published in the `fastqc/` subdirectory.

### Takeaway

You know how to create a module containing a process, import it into a workflow, call it with an input channel, and publish the results using the workflow-level output block.

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

We just need to write the process definition, import it, call it in the workflow, and update the output handling.

### 2.1. Write the trimming process and call it in the workflow

As before, we need to fill in the process definition, import the module, and add the process call.

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

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

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

Now add the process call to the workflow.

#### 2.1.3. Call the trimming process on the input

Add the process call in the workflow block:

=== "After"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Initial quality control
        FASTQC(read_ch)

        // Adapter trimming and post-trimming QC
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Initial quality control
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

The trimming process is now wired into the workflow.

### 2.2. Update the output handling

We need to add the trimming outputs to the publish declaration and configure where they go.

#### 2.2.1. Add publish targets for the trimming outputs

Add the trimming outputs to the `publish:` section:

=== "After"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Now we need to tell Nextflow where to put these outputs.

#### 2.2.2. Configure the new output targets

Add entries for the trimming targets in the `output {}` block, publishing them into a `trimming/` subdirectory:

=== "After"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

The output configuration is complete.

### 2.3. Run the workflow

The workflow now includes both initial QC and adapter trimming.

```bash
nextflow run rnaseq.nf -profile test
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

You can find the trimming outputs in the results directory.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

The trimming outputs and post-trimming QC reports are now in the `trimming/` subdirectory.

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

We need to declare a parameter for the genome index archive.

#### 3.1.1. Add a parameter for the genome index

Add a parameter declaration for the genome index archive in `rnaseq.nf`:

=== "After"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Primary input
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primary input
        input: Path
    }
    ```

#### 3.1.2. Add the genome index default to the test profile

Just as we did for `reads` in section 1.1.2, add a default value for the genome index to the test profile in `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

The parameter is ready; now we can create the alignment process.

### 3.2. Write the alignment process and call it in the workflow

As before, we need to fill in the process definition, import the module, and add the process call.

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

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

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

Now add the process call to the workflow.

#### 3.2.3. Call the alignment process

The trimmed reads are in the `TRIM_GALORE.out.trimmed_reads` channel output by the previous step.
We use `#!groovy file(params.hisat2_index_zip)` to provide the genome index archive.

=== "After"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Initial quality control
        FASTQC(read_ch)

        // Adapter trimming and post-trimming QC
        TRIM_GALORE(read_ch)

        // Alignment to a reference genome
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)

        // Initial quality control
        FASTQC(read_ch)

        // Adapter trimming and post-trimming QC
        TRIM_GALORE(read_ch)
    ```

The alignment process is now wired into the workflow.

### 3.3. Update the output handling

We need to add the alignment outputs to the publish declaration and configure where they go.

#### 3.3.1. Add publish targets for the alignment outputs

Add the alignment outputs to the `publish:` section:

=== "After"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

Now we need to tell Nextflow where to put these outputs.

#### 3.3.2. Configure the new output targets

Add entries for the alignment targets in the `output {}` block, publishing them into an `align/` subdirectory:

=== "After"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

The output configuration is complete.

### 3.4. Run the workflow

The workflow now includes all three processing steps: QC, trimming, and alignment.

```bash
nextflow run rnaseq.nf -profile test
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

You can find the alignment outputs in the results directory.

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
