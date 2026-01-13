# Part 1: Method overview and manual testing

There are multiple valid methods for processing and analyzing bulk RNAseq data.
For this course, we are following the method described [here](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) by Drs. Simon Andrews and Laura Biggins at the [Babraham Institute](https://www.babraham.ac.uk/).

Our goal is to develop a workflow that implements the following processing steps: run initial quality control on reads in a bulk RNAseq sample, trim adapter sequences from the reads, align the reads to a reference genome, and produce a comprehensive quality control (QC) report.

<figure class="excalidraw">
--8<-- "docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** Perform QC on the read data before trimming using FastQC
- **TRIM_GALORE:** Trim adapter sequences and perform QC after trimming using Trim Galore (bundles Cutadapt and FastQC)
- **HISAT2_ALIGN:** Align reads to the reference genome using Hisat2
- **MULTIQC:** Generate a comprehensive QC report using MultiQC

However, before we dive into writing any workflow code, we are going to try out the commands manually on some test data.
The tools we need are not installed in the GitHub Codespaces environment, so we'll use them via containers (see [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

     Make sure you're in the `nf4-science/rnaseq` directory. The last part of the path shown when you type `pwd` should be `rnaseq`.

---

## 1. Initial QC and adapter trimming

We're going to pull a container image that has both `fastqc` and `trim_galore` installed, spin it up interactively and run the trimming and QC commands on one of the example data files.

### 1.1. Pull the container

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

This gives you the following console output as the system downloads the image:

```console title="Output"
0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
dafa2b0c44d2: Pull complete
dec6b097362e: Pull complete
f88da01cff0b: Pull complete
4f4fb700ef54: Pull complete
92dc97a3ef36: Pull complete
403f74b0f85e: Pull complete
10b8c00c10a5: Pull complete
17dc7ea432cc: Pull complete
bb36d6c3110d: Pull complete
0ea1a16bbe82: Pull complete
030a47592a0a: Pull complete
32ec762be2d0: Pull complete
d2cb90387285: Pull complete
Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

### 1.2. Spin up the container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Your prompt will change to something like `(base) root@b645838b3314:/tmp#`, which indicates that you are now inside the container.

The `-v ./data:/data` part of the command will enable us to access the contents of the `data/` directory from inside the container.

```bash
ls /data/reads
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
```

### 1.3. Run the first `fastqc` command

Let's run `fastqc` to collect quality control metrics on the read data.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

This should run very quickly:

```console title="Output"
application/gzip
Started analysis of ENCSR000COQ1_1.fastq.gz
Approx 5% complete for ENCSR000COQ1_1.fastq.gz
Approx 10% complete for ENCSR000COQ1_1.fastq.gz
Approx 15% complete for ENCSR000COQ1_1.fastq.gz
Approx 20% complete for ENCSR000COQ1_1.fastq.gz
Approx 25% complete for ENCSR000COQ1_1.fastq.gz
Approx 30% complete for ENCSR000COQ1_1.fastq.gz
Approx 35% complete for ENCSR000COQ1_1.fastq.gz
Approx 40% complete for ENCSR000COQ1_1.fastq.gz
Approx 45% complete for ENCSR000COQ1_1.fastq.gz
Approx 50% complete for ENCSR000COQ1_1.fastq.gz
Approx 55% complete for ENCSR000COQ1_1.fastq.gz
Approx 60% complete for ENCSR000COQ1_1.fastq.gz
Approx 65% complete for ENCSR000COQ1_1.fastq.gz
Approx 70% complete for ENCSR000COQ1_1.fastq.gz
Approx 75% complete for ENCSR000COQ1_1.fastq.gz
Approx 80% complete for ENCSR000COQ1_1.fastq.gz
Approx 85% complete for ENCSR000COQ1_1.fastq.gz
Approx 90% complete for ENCSR000COQ1_1.fastq.gz
Approx 95% complete for ENCSR000COQ1_1.fastq.gz
Analysis complete for ENCSR000COQ1_1.fastq.gz
```

You can find the output files in the same directory as the original data:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

```console title="Output"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Trim adapter sequences with `trim_galore`

Now let's run `trim_galore`, which bundles Cutadapt and FastQC, to trim the adapter sequences and collect post-trimming QC metrics.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

The `--fastqc` flag causes the command to automatically run a QC collection step after trimming is complete.

_The output is very verbose so what follows is abbreviated._

```console title="Output"
Multicore support not enabled. Proceeding with single-core trimming.
Path to Cutadapt set as: 'cutadapt' (default)
Cutadapt seems to be working fine (tested command 'cutadapt --version')
Cutadapt version: 4.9
single-core operation.
igzip command line interface 2.31.0
igzip detected. Using igzip for decompressing

<...>

Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
```

You can find the output files in the working directory:

```bash
ls ENCSR000COQ1_1*
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Move the output files to the filesystem outside the container

Anything that remains inside the container will be inaccessible to future work so let's move these to a new directory.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Exit the container

```bash
exit
```

---

## 2. Align the reads to the reference genome

We're going to pull a container image that has `hisat2` installed, spin it up interactively and run the alignment command to align the RNAseq data to a reference genome.

### 2.1. Pull the `hisat2` container

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

```console title="Output"
Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
5e49f68a37dc010e: Pulling from library/hisat2_samtools
dafa2b0c44d2: Already exists
dec6b097362e: Already exists
f88da01cff0b: Already exists
4f4fb700ef54: Already exists
92dc97a3ef36: Already exists
403f74b0f85e: Already exists
10b8c00c10a5: Already exists
17dc7ea432cc: Already exists
bb36d6c3110d: Already exists
0ea1a16bbe82: Already exists
030a47592a0a: Already exists
e74ed5dd390b: Pull complete
abfcf0185e51: Pull complete
Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

### 2.2. Spin up the `hisat2` container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

The command is the same as before, with the relevant container URI swapped in.

### 2.3. Create the Hisat2 genome index files

Hisat2 requires the genome reference to be provided in a very specific format, and can't just consume the `genome.fa` FASTA file that we provide, so we're going to take this opportunity to create the relevant resources.

```bash
hisat2-build /data/genome.fa genome_index
```

The output is very verbose so the following is abbreviated:

```console title="Output"
Settings:
  Output files: "genome_index.*.ht2"
<...>
Total time for call to driver() for forward index: 00:00:16
```

This creates multiple genome index files, which you can find in the working directory.

```bash
ls genome_index.*
```

```console title="Output"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

We'll use these in a moment, but first let's generate a gzipped tarball with these genome index files; we'll need them later and generating these is not typically something we want to do as part of a workflow.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

This stores a `genome_index.tar.gz` tarball containing the genome index files in the `data/` directory on our filesystem, which wil come in handy in Part 2 of this course.

### 2.4. Run the `hisat2` command

Now we can run the alignment command, which performs the alignment step with `hisat2` then pipes the output to `samtools` to write the output out as a BAM file.

The read data input is the `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` file we generated with `trim_galore` in the previous step.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

This runs almost instantly because it's a very small test file.
At real scale this could take a lot longer.

```console title="Output"
HISAT2 summary stats:
        Total reads: 27816
                Aligned 0 time: 1550 (5.57%)
                Aligned 1 time: 25410 (91.35%)
                Aligned >1 times: 856 (3.08%)
        Overall alignment rate: 94.43%
```

Once again you can find the output files in the working directory:

```bash
ls ENCSR000COQ1_1*
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Move the output files to the filesystem outside the container

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Exit the container

```bash
exit
```

---

## 3. Generate a comprehensive QC report

We're going to pull a container image that has `multiqc` installed, spin it up interactively and run a report generation command on the before/after FastQC report files.

### 3.1. Pull the `multiqc` container

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

```console title="Output"
ad8f247edb55897c: Pulling from library/pip_multiqc
dafa2b0c44d2: Already exists
dec6b097362e: Already exists
f88da01cff0b: Already exists
4f4fb700ef54: Already exists
92dc97a3ef36: Already exists
403f74b0f85e: Already exists
10b8c00c10a5: Already exists
17dc7ea432cc: Already exists
bb36d6c3110d: Already exists
0ea1a16bbe82: Already exists
030a47592a0a: Already exists
3f229294c69a: Pull complete
5a5ad47fd84c: Pull complete
Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.2. Spin up the `multiqc` container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. Run the `multiqc` command

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

MultiQC is able to search through directories for compatible QC reports and will aggregate everything it finds.

```console title="Output"

/// MultiQC 🔍 v1.27.1

       file_search | Search path: /data/reads
       file_search | Search path: /data/trimmed
       file_search | Search path: /data/aligned
         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
            hisat2 | Found 1 reports
          cutadapt | Found 1 reports
            fastqc | Found 1 reports
     write_results | Data        : ENCSR000COQ1_1_QC_data
     write_results | Report      : ENCSR000COQ1_1_QC.html
           multiqc | MultiQC complete
```

Here we see the tool found all three QC reports we generated: the initial QC we did with `fastqc`, the post-trimming report from `cutadapt` (made via `trim_galore`) and the post-alignment QC produced by `hisat2`.

The output files are once again in the working directory:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Output"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Move the output files to the filesystem outside the container

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Exit the container

```bash
exit
```

---

### Takeaway

You have tested all the individual commands interactively in the relevant containers.

### What's next?

Learn how to wrap those same commands into a multi-step workflow that uses containers to execute the work.
