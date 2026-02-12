# Part 1: Method overview

There are multiple valid methods for processing and analyzing bulk RNAseq data.
For this course, we are following the method described [here](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) by Drs. Simon Andrews and Laura Biggins at the [Babraham Institute](https://www.babraham.ac.uk/).

Our goal is to develop a workflow that implements the following processing steps: run initial quality control on reads in a bulk RNAseq sample, trim adapter sequences from the reads, align the reads to a reference genome, and produce a comprehensive quality control (QC) report.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Perform QC on the read data before trimming using FastQC
- **TRIM_GALORE:** Trim adapter sequences and perform QC after trimming using Trim Galore (bundles Cutadapt and FastQC)
- **HISAT2_ALIGN:** Align reads to the reference genome using Hisat2
- **MULTIQC:** Generate a comprehensive QC report using MultiQC

### Methods

We're going to show you how to apply these processing steps in two phases.
First we'll start with **single-sample processing** that runs the QC, trimming and alignment tools on one sample.
Then we'll extend to **multi-sample processing** that runs the same tools on multiple samples and generates an aggregated quality control report.

Before we dive into writing any workflow code for either approach, we are going to try out the commands manually on some test data.

### Dataset

We provide the following data and related resources:

- **RNAseq data** (`reads/`): FASTQ files from six samples, subset to a small region to keep file sizes down. Each sample has paired-end reads (two files per sample), though we start by working with single-end reads only.
- **A reference genome** (`genome.fa`): a small region of the human chromosome 20 (from hg19/b37).
- **CSV samplesheets** (`single-end.csv` and `paired-end.csv`): files listing the IDs and paths of the example data files.

### Software

The four main tools involved are [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control metrics collection, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for adapter trimming (bundles Cutadapt and FastQC for post-trimming QC), [HISAT2](http://daehwankimlab.github.io/hisat2/) for spliced alignment to a reference genome, and [MultiQC](https://multiqc.info/) for aggregated QC report generation.

These tools are not installed in the GitHub Codespaces environment, so we'll use them via containers retrieved via the Seqera Containers service (see [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip

     Make sure you're in the `nf4-science/rnaseq` directory. The last part of the path shown when you type `pwd` should be `rnaseq`.

---

## 1. Single-sample processing

In this section we test the commands that process a single RNAseq sample: quality control, adapter trimming, and alignment to a reference genome.
These are the commands we'll wrap into a Nextflow workflow in Part 2 of this course.

1. Run initial QC on a FASTQ file using FastQC
2. Trim adapter sequences and run post-trimming QC using Trim Galore
3. Align the trimmed reads to the reference genome using HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

We start by testing these commands on just one sample.

### 1.1. QC and adapter trimming

First, we want to run the QC and trimming commands on one of the example data files.

#### 1.1.1. Pull the container

Let's pull a container image that has both `fastqc` and `trim_galore` installed:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Command output"

    ```console
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

If you haven't downloaded this image before, it may take a minute to complete.
Once it's done, you have a local copy of the container image.

#### 1.1.2. Spin up the container interactively

To run the container interactively, use `docker run` with the `-it` flags.
The `-v ./data:/data` option mounts our local `data/` directory so we can access the input files from inside the container.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Command output"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Your prompt will change to something like `(base) root@b645838b3314:/tmp#`, which indicates that you are now inside the container.

Verify that you can see the sequence data files under `/data/reads`:

```bash
ls /data/reads
```

??? abstract "Directory contents"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

With that, you are ready to try your first command.

#### 1.1.3. Run the FastQC command

The method referenced above gives us the command line to run QC on a single file.
We only need to provide the input file; the tool will automatically generate output files in the same directory as the original data.

Run the `fastqc` command on one data file:

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Command output"

    ```console
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

This should run very quickly.
You can find the output files in the same directory as the original data:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "Directory contents"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

You should see an HTML report and a ZIP archive containing the QC metrics.
That completes the testing of the first step.

#### 1.1.4. Trim adapter sequences with Trim Galore

Now let's run `trim_galore`, which bundles Cutadapt and FastQC, to trim the adapter sequences and collect post-trimming QC metrics.
As noted above, the software is included in the same container, so no change needed there.

The command is straightforward; we simply need to add the `--fastqc` flag to automatically run a QC collection step after trimming is complete.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Command output"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

The output is very verbose, so we've highlighted the most relevant lines in the example above.
You can find the output files in the working directory:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Directory contents"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

This includes the trimmed reads, trimming report and post-trimming QC files.

#### 1.1.5. Move the output files

Anything that remains inside the container will be inaccessible to future work, so we need to move these files to a directory on the mounted filesystem.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Directory contents"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

The files are now accessible in your normal filesystem.

#### 1.1.6. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should go back to normal; that completes the testing of the first two steps.

### 1.2. Align reads to the reference genome

Next, we want to run the alignment command to align the trimmed RNAseq reads to a reference genome.

#### 1.2.1. Pull the container

Let's pull a container image that has `hisat2` and `samtools` installed:

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Command output"

    ```console
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

You'll notice that some layers show `Already exists` because they are shared with the Trim Galore container image we pulled earlier.
As a result, this pull should go faster than the first one.

#### 1.2.2. Spin up the container interactively

Spin up the container interactively, using the same approach as before with the relevant container URI swapped in.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Your prompt will change again to indicate you are inside the container.

#### 1.2.3. Create the genome index files

HISAT2 requires the genome reference to be provided in a very specific format, and can't just consume the `genome.fa` FASTA file that we provide, so we're going to take this opportunity to create the relevant resources.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Command output"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

The output is very verbose, so we've highlighted some relevant lines in the example above.

This creates multiple genome index files, which you can find in the working directory.

```bash
ls genome_index.*
```

??? abstract "Directory contents"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

We'll need these files later, and generating them is not typically something we want to do as part of a workflow, so we're going to generate a gzipped tarball containing the genome index files that we can easily pass around as needed.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Command output"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

We'll move the resulting `genome_index.tar.gz` tarball containing the genome index files to the `data/` directory on our filesystem in a few minutes.
That will come in handy in Part 2 of this course.

#### 1.2.4. Run the alignment command

Now we can run the alignment command, which performs the alignment step with `hisat2` then pipes the output to `samtools` to write the output out as a BAM file.

The read data input is the `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` file we generated with `trim_galore` in the previous step.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Command output"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

This runs almost instantly because it's a very small test file.
At full scale, this could take a lot longer.

Once again you can find the output files in the working directory:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Directory contents"

    ```console title="Output"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

The alignment produced a BAM file and a log file with alignment statistics.

#### 1.2.5. Move the output files

As before, move the output files to a directory on the mounted filesystem so they remain accessible after we exit the container.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

With that done, we have everything we need.

#### 1.2.6. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should go back to normal.
That concludes the single-sample processing testing run.

!!! example "Write it as a workflow!"

    Feel free to move on to [Part 2](./02_single-sample.md) right away if you'd like to get started implementing this analysis as a Nextflow workflow.
    You'll just need to come back to complete the second round of testing before moving on to Part 3.

---

## 2. Multi-sample QC aggregation

The commands we just tested process one sample at a time.
In practice, we typically need to process many samples and then aggregate QC results across all of them to evaluate the quality of the overall dataset.

[MultiQC](https://multiqc.info/) is a tool that searches through directories for QC reports from many common bioinformatics tools and aggregates them into a single comprehensive HTML report.
It can recognize output from FastQC, Cutadapt (via Trim Galore) and HISAT2, among many others.

Here we process two additional samples through the same per-sample tools, then use MultiQC to aggregate QC reports across all three samples.
These are the commands we'll wrap into a Nextflow workflow in Part 3 of this course.

1. Run QC and trimming on additional samples using Trim Galore
2. Run alignment on additional samples using HISAT2
3. Aggregate all QC reports into a comprehensive report using MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC and trim additional samples

The per-sample QC and trimming commands are identical to what we ran in section 1.1.
We already pulled the container image, so we can spin it up directly.

#### 2.1.1. Spin up the container

We already pulled this container image in section 1.1, so we can spin it up directly:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Your prompt changes to indicate you are inside the container.

#### 2.1.2. Run QC and trimming on additional samples

Run FastQC and Trim Galore on two more samples, one after another.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Once this completes, you should have Trim Galore output files for both samples in the working directory.

#### 2.1.3. Move the output files

Move the Trim Galore output files to the same directory we used in section 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Directory contents"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

The files are now accessible in your normal filesystem.

#### 2.1.4. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should go back to normal.

### 2.2. Align additional samples

The alignment commands are identical to what we ran in section 1.2.
We need to extract the genome index from the tarball we saved earlier, since the original index files were created inside a container that no longer exists.

#### 2.2.1. Spin up the container

We already pulled this container image in section 1.2, so we can spin it up directly:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Your prompt changes to indicate you are inside the container.

#### 2.2.2. Extract the genome index

Extract the genome index files from the tarball we saved to the mounted filesystem:

```bash
tar -xzf /data/genome_index.tar.gz
```

This restores the `genome_index.*` files in the working directory.

#### 2.2.3. Run alignment on additional samples

Run the HISAT2 alignment on the two newly trimmed samples, one after another.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Command output"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "Command output"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Once this completes, you should have BAM and log files for both samples in the working directory.

#### 2.2.4. Move the output files

Move the alignment output files to the same directory we used in section 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Directory contents"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

The files are now accessible in your normal filesystem.

#### 2.2.5. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should go back to normal.

### 2.3. Generate a comprehensive QC report

Now that we have QC, trimming and alignment output for three samples, we can use MultiQC to aggregate them into a single report.
MultiQC searches through directories for compatible QC reports and aggregates everything it finds.

#### 2.3.1. Pull the container

Let's pull a container image that has `multiqc` installed:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Command output"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
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
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

You'll notice that some layers show `Already exists` because they are shared with the container images we pulled earlier.
As a result, this pull should go faster than the previous ones.

#### 2.3.2. Spin up the container interactively

Spin up the container interactively with the data directory mounted, as before.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Your prompt will change to indicate you are inside the container.

#### 2.3.3. Run the MultiQC command

Run `multiqc`, pointing it at the directories where we stored QC-related output files for all three samples.
The `-n` flag sets the name of the output report.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Command output"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

Here we see the tool found QC reports for all three samples: the initial QC from `fastqc`, the post-trimming reports from `cutadapt` (via `trim_galore`) and the alignment summaries from `hisat2`.

The output files are in the working directory:

```bash
ls all_samples_QC*
```

??? abstract "Directory contents"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

The main output is the `all_samples_QC.html` report, accompanied by a data directory containing the underlying metrics.

#### 2.3.4. Move the output files

Move the report and its data directory to the mounted filesystem.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

The files are now accessible in your normal filesystem.

#### 2.3.5. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should go back to normal.
That concludes the testing of all the RNAseq processing commands.

---

### Takeaway

You know how to run the FastQC, Trim Galore, HISAT2 and MultiQC commands in their respective containers, including how to process multiple samples and aggregate QC reports.

### What's next?

Take a break, then head on to [Part 2](./02_single-sample.md) to learn how to wrap those same commands into workflows that use containers to execute the work.
