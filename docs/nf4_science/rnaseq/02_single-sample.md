# Part 2: Single-sample implementation

In this part of the course, we're going to write the simplest possible workflow that wraps all the commands we ran in Part 1 to automate running them, and we'll just aim to process one sample at a time.

We'll do this in three stages:

1. Write a single-stage workflow that runs the initial QC step
2. Add adapter trimming and post-trimming QC
3. Add alignment to the reference genome

---

## 1. Write a single-stage workflow that runs the initial QC

Let's start by writing a simple workflow that runs the FastQC tool on a FASTQ file containing single-end RNAseq reads.

We provide you with a workflow file, `rnaseq.nf`, that outlines the main parts of the workflow.

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

Keep in mind this workflow code is correct but it's not functional; its purpose is just to serve as a skeleton that you'll use to write the actual workflow.

### 1.1. Create a directory to store modules

We'll create standalone modules for each process to make it easier to manage and reuse them, so let's create a directory to store them.

```bash
mkdir modules
```

### 1.2. Create a module for the QC metrics collection process

Let's create a module file called `modules/fastqc.nf` to house the `FASTQC` process:

```bash
touch modules/fastqc.nf
```

Open the file in the code editor and copy the following code into it:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

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
    fastqc $reads
    """
}
```

You should recognize all the pieces from what you learned in Part 1 & Part 2 of this training series; the only notable change is that this time we're using `mode: symlink` for the `publishDir` directive, and we're using a parameter to define the `publishDir`.

!!! note

    Even though the data files we're using here are very small, in genomics they can get very large. For the purposes of demonstration in the teaching environment, we're using the 'symlink' publishing mode to avoid unnecessary file copies. You shouldn't do this in your final workflows, since you'll lose results when you clean up your `work` directory.

### 1.2. Import the module into the workflow file

Add the statement `include { FASTQC } from './modules/fastqc.nf'` to the `rnaseq.nf` file:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
```

### 1.3. Add an input declaration

Declare an input parameter with a default value:

```groovy title="rnaseq.nf" linenums="10"
// Primary input
params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"
```

### 1.4. Create an input channel in the workflow block

Use a basic `.fromPath()` channel factory to create the input channel:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Call processes

}
```

### 1.5. Call the `FASTQC` process on the input channel

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
```

### 1.6. Run the workflow to test that it works

We could use the `--reads` parameter to specify an input from command line, but during development we can be lazy and just use the test default we set up.

```bash
nextflow run rnaseq.nf
```

This should run very quickly if you worked through Part 1 and have already pulled the container.
If you skipped it, Nextflow will pull the container for you; you don't have to do anything for it to happen, but you may need to wait up to a minute.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

executor >  local (1)
[d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
```

You can find the outputs under `results/fastqc` as specified in the `FASTQC` process by the `publishDir` directive.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

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
