# Orientation

The training environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please the [Environment Setup](../../envsetup/) mini-course before going any further.

## Materials provided

Throughout this training course, we'll be working in the `nf4-science/rnaseq/` directory, which you need to move into when you open the training workspace.
This directory contains all the code files, test data and accessory files you will need.

Feel free to explore the contents of this directory; the easiest way to do so is to use the file explorer on the left-hand side of the training workspace in the VSCode interface.
Alternatively, you can use the `tree` command.
Throughout the course, we use the output of `tree` to represent directory structure and contents in a readable form, sometimes with minor modifications for clarity.

Here we generate a table of contents to the second level down:

```bash
tree . -L 3
```

If you run this inside `nf4-science/rnaseq`, you should see the following output:

```console title="Directory contents"
.
├── data
│   ├── genome.fa
│   ├── paired-end.csv
│   ├── reads
│   │   ├── ENCSR000COQ1_1.fastq.gz
│   │   ├── ENCSR000COQ1_2.fastq.gz
│   │   ├── ENCSR000COQ2_1.fastq.gz
│   │   ├── ENCSR000COQ2_2.fastq.gz
│   │   ├── ENCSR000COR1_1.fastq.gz
│   │   ├── ENCSR000COR1_2.fastq.gz
│   │   ├── ENCSR000COR2_1.fastq.gz
│   │   ├── ENCSR000COR2_2.fastq.gz
│   │   ├── ENCSR000CPO1_1.fastq.gz
│   │   ├── ENCSR000CPO1_2.fastq.gz
│   │   ├── ENCSR000CPO2_1.fastq.gz
│   │   └── ENCSR000CPO2_2.fastq.gz
│   └── single-end.csv
├── nextflow.config
├── rnaseq.nf
└── solutions
    ├── modules
    │   ├── fastqc.nf
    │   ├── fastqc_pe.nf
    │   ├── hisat2_align.nf
    │   ├── hisat2_align_pe.nf
    │   ├── multiqc.nf
    │   ├── trim_galore.nf
    │   └── trim_galore_pe.nf
    ├── rnaseq-2.1.nf
    ├── rnaseq-2.2.nf
    ├── rnaseq-2.3.nf
    ├── rnaseq-3.1.nf
    ├── rnaseq-3.2.nf
    └── rnaseq_pe-3.3.nf

```

!!!note

    Don't worry if this seems like a lot; we'll go through the relevant pieces at each step of the course.
    This is just meant to give you an overview.

**Here's a summary of what you should know to get started:**

- **The `rnaseq.nf` file** is the outline if the workflow script we will work to develop.

- **The file `nextflow.config`** is a configuration file that sets minimal environment properties. You can ignore it for now.

- **The `data` directory** contains input data and related resources:

  - _A reference genome_ called `genome.fa` consisting of a small region of the human chromosome 20 (from hg19/b37).
  - _RNAseq data_ that has been subset to a small region to keep the file sizes down, in the `reads/` directory.
  - _CSV files_ listing the IDs and paths of the example data files, for processing in batches.

- **The `solutions` directory** contains the completed workflow scripts and modules that result from each step of the course.
  They are intended to be used as a reference to check your work and troubleshoot any issues.
  The number in the filename corresponds to the step of the relevant part of the course.

!!!tip

    If for whatever reason you move out of this directory, you can always run this command to return to it:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Now, to begin the course, click on the arrow in the bottom right corner of this page.
