# Part 2: Single-sample implementation of RNAseq processing

TODO: short blurb

- Write a single-stage workflow that runs the initial QC
- Add adapter trimming and post-trimming QC
- Add alignment to the reference genome
- Add comprehensive QC report generation

---

## 1. Write a single-stage workflow that runs the initial QC

Let's start by writing a simple workflow that runs the FastQC tool on a fastq file containing single-end RNAseq reads.

We provide you with a workflow file, `rnaseq-1.nf`, that outlines the main parts of the workflow.
It's not functional; its purpose is just to serve as a skeleton that you'll use to write the actual workflow.

We'll store the code for each process in a standalone module to make it easier to reuse.

### 1.1. Describe the QC metrics collection process

Let's start by writing a process, which we'll call `FASTQC`, describing the `fastqc` operation.

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'copy'

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

<!-- TODO Add a note about assuming single-end reads -->

### 1.2. Import the module into the workflow file

TODO

### 1.3. Add an input declaration

TODO

### 1.4. Create an input channel in the workflow block

TODO

### 1.5. Call the process on the input channel

TODO

### 1.6. Run the workflow to test that it works

TODO

---

## 2. Add adapter trimming and post-trimming quality control

TODO

### 2.1. Describe the trimming and QC process

Let's write a process, which we'll call `TRIM_GALORE`, that trims adapter sequences with Cutadapt and runs post-trimming quality control with FastQC (both tools are bundled into the Trim_Galore wrapper).

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'copy'

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

<!-- TODO Add a note about assuming single-end reads -->

### 2.2. Import the module into the workflow file

TODO

### 2.3. Call the process on the input channel

TODO

### 2.4. Run the workflow to test that it works

TODO

---

## 3. Align the reads to the reference genome

TODO

### 3.1. Describe the HiSat2 process

Let's write a process, which we'll call `HISAT2_ALIGN`, that aligns the trimmed reads to the reference genome.

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'copy'

    input:
    path reads
    path index

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    hisat2 -x ${index.simpleName} -U $reads --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
    samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

<!-- TODO Add a note about assuming single-end reads -->

### 3.2. Import the module into the workflow file

TODO

### 3.3. Call the process on the trimmed reads

TODO

### 3.4. Run the workflow to test that it works

TODO

---

## 4. Aggregate pre-processing QC metrics into a single MultiQC report

TODO

### 4.1. Describe the MultiQC process

Let's write a process, which we'll call `MULTIQC`, that collect QC metrics with MultiQC in a generic way.

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c"
    publishDir "results/multiqc", mode: 'copy'

    input:
    path '*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}
```

<!-- TODO Add a note about assuming single-end reads -->

### 4.2. Import the module into the workflow file

TODO

### 4.3. Call the process on the outputs of the previous QC steps

TODO

### 4.4. Run the workflow to test that it works

TODO

---

### Takeaway

TODO

### What's next?

Modify the workflow to take a CSV file of file paths as inputs to parallelize the processing for multiple samples.

---
