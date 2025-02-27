# Part 1: RNAseq processing method overview

There are multiple valid methods for processing and analyzing bulk RNAseq data.
For this course, we are following the method described [here](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) by Drs. Simon Andrews and Laura Biggins at the [Babraham Institute](https://www.babraham.ac.uk/).

Our goal is to develop a workflow that implements the following processing steps: run initial quality control on reads in a bulk RNAseq sample, trim adapter sequences from the reads, align the reads to a reference genome, and produce a comprehensive quality control (QC) report.

<figure class="excalidraw">
--8<-- "docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC_RAW:** Perform QC on the read data before trimming using FastQC
- **TRIM_GALORE:** Trim adapter sequences and perform QC after trimming using Trim Galore (bundles Cutadapt and FastQC)
- **HISAT2_ALIGN:** Align reads to the reference genome using Hisat2
- **MULTIQC:** Generate a comprehensive QC report using MultiQC

However, before we dive into writing any workflow code, we are going to try out the commands manually on some test data.
The tools we need are not installed in the GitHub Codespaces environment, so we'll use them via containers (see [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

     Make sure you're in the `nf4-science/rnaseq` directory so that the last part of the path shown when you type `pwd` is `rnaseq`.

---

## 1. Initial QC and adapter trimming

We're going to pull down a container that has both `fastqc` and `trim_galore` installed, spin it up interactively and run the trimming and QC commands on one of the example data files.

### 1.1. Pull the container

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

### 1.2. Spin up the container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

### 1.3. Run the first `fastqc` command

The [`fastqc` documentation](URL) gives us the command line to run to collect quality control metrics on the read data.

```bash
TODO
```

This should complete immediately, and you should now see output files in the working directory.

```console title="Directory contents"
data/TODO
```

### 1.4. Trim adapter sequences with `trim_galore`

The [`trim_galore`](URL) gives us the command line to use:

```bash
TODO
```

The `--fastqc` flag causes the command to automatically run a QC collection step after trimming is complete.

TODO

```console title="[OUTPUT]" linenums=""
TODO
```

TODO

### 1.6. Exit the container

```bash
exit
```

---

## 2. Align the reads to the reference genome

We're going to pull down a container that has `hisat2` installed, spin it up interactively and run the alignment command to align the RNAseq data to a reference genome.

### 2.1. Pull the `hisat2` container

```bash
docker pull <CONTAINER URI>
```

### 2.2. Spin up the `hisat2` container interactively

```bash
docker run -it -v ./data:/data <CONTAINER URI>
```

### 2.3. Run the `hisat2` command

```bash
TODO
```

TODO

```console title="[OUTPUT]" linenums=""
TODO
```

TODO

### 2.4. Exit the container

```bash
exit
```

---

## 3. Generate a comprehensive QC report

We're going to pull down a container that has `multiqc` installed, spin it up interactively and run a report generation command on the before/after FastQC report files.

### 3.1. Pull the `multiqc` container

```bash
docker pull <CONTAINER URI>
```

### 3.2. Spin up the `multiqc` container interactively

```bash
docker run -it -v ./data:/data <CONTAINER URI>
```

### 3.3. Run the `multiqc` command

```bash
TODO
```

TODO

```console title="[OUTPUT]" linenums=""
TODO
```

TODO

### 3.4. Exit the container

```bash
exit
```

---

### Takeaway

You know how to test the individual commands interactively in the relevant containers.

### What's next?

Learn how to wrap those same commands into a multi-step workflow that uses containers to execute the work.
