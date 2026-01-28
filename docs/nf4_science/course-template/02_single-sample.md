# Part 2: Single-sample implementation

In this part of the course, we're going to write the simplest possible workflow that wraps all the commands we ran in Part 1 to automate running them.

We'll do this in three stages:

1. Write a single-stage workflow that runs the initial [analysis] step
2. Add [intermediate steps]
3. Add [final step]

Before we dive into development, let's review how the code is going to be structured.

---

## 0. Warmup: Review the code structure

[TODO: intro text. Include note that we're going to write processes as standalone modules and use a configuration file; refer to relevant Hello Nextflow parts for details]

[TODO: image showing elements?]

[TODO: note that we provide templates for each part]

### 0.1. The workflow template

We provide you with a workflow template file, `[course-name].nf`, that outlines the main parts of the workflow.
You'll develop your workflow by editing this file directly.

```groovy title="[course-name].nf" linenums="1"
#!/usr/bin/env nextflow

// List module import statements
include { EXAMPLE } from './modules/module-template.nf'

// Define the workflow
workflow {

    main:
    // Create input channel

    // Call processes

    publish:
    // Declare outputs that we care about
}

// Define output destinations
output {

    // List outputs declared under 'publish'

}
```

This workflow code is technically correct (you can run it) but it won't _do_ anything; its purpose is just to serve as a template for writing the actual workflow.

### 0.2. The process module template

By convention, we store process modules in a directory called `modules/`, which already contains a basic module template.

```groovy title="modules/module-template.nf" linenums="1"
#!/usr/bin/env nextflow

// Define the process
process EXAMPLE {

    // specify container directive

    input:
    // declare input(s)

    output:
    // declare output(s) (use 'emit' to tag if >1)

    script:
    """
    # command to run
    """
}
```

Just like for the workflow, you will use this as a template to write the actual processes, but you'll make a copy of the template to create each module.

### 0.3. The configuration file

Finally, we include a basic `nextflow.config` configuration file.

```groovy title="nextflow.config" linenums="1"
#!/usr/bin/env nextflow

docker.enable = true
workflow.output.mode = 'copy'

params {
  // declare default pipeline parameter values
  example = "this is just an example"
}

```

This is not a template, it's a fully functional configuration file that currently enables use of Docker containers, sets the workflow output mode and declares an example parameter.

You will add default parameter values to this file to make it easier to run commands without a lot of typing.

Now that we've reviewed these building blocks, let's start building!

---

## 1. Write a single-stage workflow that runs the [first step]

As noted above, we're going to start our pipeline project by writing a simple workflow that runs [tool] on [input] containing [type-of-data] and produces [output].

We are going to go through this in five steps:

1. Write the first process and import the module
2. Set up the input (create an input channel and configure a default value)
3. Call the process on the input
4. Set up the workflow outputs
5. Test that the workflow runs correctly

[TODO: any notes about naming conventions?]

Let's go!

### 1.1. Write the first process and import the module

[TODO: templated comment for summary]

#### 1.1.1. Create the [analysis] module

Let's create a module file called `modules/[module-name].nf` by copying the module template as follows:

```bash
cp modules/module-template.nf modules/[module-name].nf
```

Then, open the new file and edit it as follows.

```groovy title="modules/[module-name].nf" linenums="1"
#!/usr/bin/env nextflow

process [module-name] {

    container "community.wave.seqera.io/library/[container]"

    input:
    [inputs]

    output:
    [outputs (with emit if >1)]

    script:
    """
    [command]]
    """
}
```

[Add explanations as appropriate, focusing on anything specific to the analysis domain]

[TODO: templated comment for recap]

#### 1.1.2. Import the [analysis] module into the workflow

Open the workflow file and edit the example module import statement to import the [module-name] module.

[TODO: correct syntax for before/after display]

[AFTER]

```groovy title="[course-name].nf" linenums="3"
// // List module import statements
include { [module-name] } from './modules/[module-name].nf'
```

[BEFORE]

```groovy title="[course-name].nf" linenums="3"
// // List module import statements
include { EXAMPLE } from './modules/module-template.nf'
```

Now the process is available to be called in the workflow, but we need an input to run it on.

### 1.2. Set up the input

[TODO: templated comment for summary]

#### 1.2.1. Create an input channel in the workflow block

[TODO: templated comment for summary -- should explain we use the `channel.fromPath()` channel factory to create the input channel]

[Adapt as appropriate if using a different channel factory]

```groovy title="[course-name].nf" hl_lines="6" linenums="6"
// Define the workflow
workflow {

    main:
    // Create input channel
    input_ch = channel.fromPath(params.input)

    // Call processes

    publish:
    // Declare outputs that we care about
}
```

[TODO: templated comment for summary -- say we want to set a default value so we don't have to include it on the command line every time we want to test the workflow]

#### 1.2.2. Configure a default value

Open the `nextflow.config` configuration file and replace the example parameter declaration with default input value declaration.

[TODO: correct syntax for before/after display]

[AFTER]

```groovy title="nextflow.config" linenums="6"
params {
  // declare default pipeline parameter values
  input = "[default path for single input file]"
}
```

[BEFORE]

```groovy title="nextflow.config" linenums="6"
params {
  // declare default pipeline parameter values
  example = "this is just an example"
}
```

[TODO: recap what this does]

### 1.3. Call the `[process-name]` process on the input

[TODO: templated comment for summary]

```groovy title="[course-name].nf" hl_lines="9" linenums="6"
// Define the workflow
workflow {

    main:
    // Create input channel
    input_ch = channel.fromPath(params.input)

    // Call processes
    [process-name](input_ch)

    publish:
    // Declare outputs that we care about
}
```

### 1.4. Set up the workflow outputs

[TODO: templated comment for summary]

#### 1.4.1. Declare the outputs that we care about

[TODO: templated comment for summary]

```groovy title="[course-name].nf" hl_lines="13" linenums="6"
// Define the workflow
workflow {

    main:
    // Create input channel
    input_ch = channel.fromPath(params.input)

    // Call processes
    [process-name](input_ch)

    publish:
    // Declare outputs that we care about
    [output-name] = [process-name].out[.emit-tag]
}
```

That completes the definition of our single-step workflow as such.
However, we have one more thing to do for the workflow outputs to be handled appropriately, because any outputs listed under `publish:` must be listed in the `output { }` block.

#### 1.4.2. Define output destinations

By default, Nextflow will write workflow outputs under a `results/` directory.
We have to

```groovy title="[course-name].nf" hl_lines="5" linenums="18"
// Define output destinations
output {

    // List outputs declared under 'publish'
    [output-name] {
        path '.'
    }
}
```

[TODO: templated comment for recap]

### 1.5. Test that the workflow runs correctly

We could use the `--input` parameter to specify an input from command line, but during development we can be lazy and just use the test default we set up.

```bash
nextflow run [course-name].nf
```

This should run very quickly if you worked through Part 1 and have already pulled the container.
If you skipped it, Nextflow will pull the container for you; you don't have to do anything for it to happen, but you may need to wait up to a minute.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `[course-name].nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

executor >  local (1)
[d6/d94c3a] [process-name] (1) [100%] 1 of 1 ✔
```

You can find the outputs under `results/` as determined by the workflow output configuration.

```bash
ls results/
```

```console title="Output"
[console output]
```

---

## [TODO: NOT ADAPTED YET BEYOND THIS POINT]

## 2. Add adapter trimming and post-trimming quality control

We're going to use the Trim_Galore wrapper, which bundles Cutadapt for the trimming itself and FastQC for the post-trimming quality control.

### 2.1. Create a module for the trimming and QC process

Let's create a module file called `modules/trim_galore.nf` to house the `TRIM_GALORE` process:

```bash
touch modules/trim_galore.nf
```

Open the file in the code editor and copy the following code into it:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

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
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Import the module into the workflow file

Add the statement `include { TRIM_GALORE } from './modules/trim_galore.nf'` to the `rnaseq.nf` file:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Call the process on the input channel

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. Run the workflow to test that it works

```bash
nextflow run rnaseq.nf
```

This should run very quickly too, since we're runnng on such a small input file.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

executor >  local (1)
[d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
[c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
```

You can find the outputs under `results/trimming` as specified in the `TRIM_GALORE` process by the `publishDir` directive.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Align the reads to the reference genome

Finally we can run the genome alignment step using Hisat2, which will also emit FastQC-style quality control metrics.

### 3.1. Create a module for the HiSat2 process

Let's create a module file called `modules/hisat2_align.nf` to house the `HISAT2_ALIGN` process:

```bash
touch modules/hisat2_align.nf
```

Open the file in the code editor and copy the following code into it:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

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
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Import the module into the workflow file

Add the statement `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` to the `rnaseq.nf` file:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Add a parameter declaration to provide the genome index

Declare an input parameter with a default value:

```groovy title="rnaseq.nf" linenums="8"
/*
 * Pipeline parameters
 */
params.hisat2_index_zip = "data/genome_index.tar.gz"
```

### 3.4. Call the `HISAT2_ALIGN` process on the trimmed reads output by `TRIM_GALORE`

The trimmed reads are in the `TRIM_GALORE.out.trimmed_reads` channel output by the previous step.

In addition, we use `file (params.hisat2_index_zip)` to provide the Hisat2 tool with the gzipped genome index tarball.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.4. Run the workflow to test that it works

```bash
nextflow run rnaseq.nf
```

This should run very quickly too, since we're runnng on such a small input file.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

executor >  local (3)
[e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
[c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
[c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
```

You can find the outputs under `results/align` as specified in the `HISAT2_ALIGN` process by the `publishDir` directive.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

This completes the basic processing we need to apply to each sample.

_We'll add MultiQC report aggregation in Part 2, after we've made the workflow accept multiple samples at a time._

---

### Takeaway

You know how to wrap all the core steps to process single-end RNAseq samples individually.

### What's next?

Learn how to modify the workflow to process multiple samples in parallel, aggregate QC reports across all steps for all samples, and enable running the workflow on paired-end RNAseq data.
