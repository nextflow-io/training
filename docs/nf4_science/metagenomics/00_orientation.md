# Orientation

The training environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please the [Environment Setup](../../envsetup/) mini-course before going any further.

## Materials provided

For the purpose of the course, we'll be working in the `nf4-science/metagenomics/` directory, where you will find all the code files, test data and accessory files you will need.
To move into it, run the following command:

```bash
cd nf4-science/metagenomics/
```

Before we go any further, we are going to download some files that are too large to be permanently stored within the GitHub repository.
Specifically, this is a set of files that constitute the database required by Kraken2 and Bracken.

Run the following commands in that exact order and wait until all of them are finished:

```bash
mkdir -p data/viral_db && cd "$_"
wget --no-check-certificate --no-proxy 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20241228.tar.gz'
tar -xvzf k2_viral_20241228.tar.gz
rm -r k2_viral_20241228.tar.gz
cd -
```

Briefly, this creates a directory called `viral_db` under `data/` and moves into it.
Then, it downloads an archive file with `wget`, unpacks its contents with `tar`, and deletes the original archive file.
Finally, it moves you back up to the original `nf4-science/metagenomics/` directory.

Now, let's take a look of the files contained in this directory with the `tree` command:

```bash
tree . -L 3
```

Here you should see the following directory structure:

```console title="Directory contents"
.
├── bin
│   └── report.Rmd
├── data
│   ├── samples
│   │   ├── ERR2143768
│   │   │   ├── ERR2143768_1.fastq
│   │   │   └── ERR2143768_2.fastq
│   │   ├── ERR2143769
│   │   │   ├── ERR2143769_1.fastq
│   │   │   └── ERR2143769_2.fastq
│   │   ├── ERR2143770
│   │   │   ├── ERR2143770_1.fastq
│   │   │   └── ERR2143770_2.fastq
│   │   └── ERR2143774
│   │       ├── ERR2143774_1.fastq
│   │       └── ERR2143774_2.fastq
│   ├── samplesheet.csv
│   └── yeast
│       ├── yeast.1.bt2
│       ├── yeast.2.bt2
│       ├── yeast.3.bt2
│       ├── yeast.4.bt2
│       ├── yeast.rev.1.bt2
│       └── yeast.rev.2.bt2
├── main.nf
├── modules
│   ├── bowtie2.nf
│   ├── bracken.nf
│   ├── kReport2Krona.nf
│   ├── knit_phyloseq.nf
│   ├── kraken2.nf
│   ├── kraken_biom.nf
│   └── ktImportText.nf
├── nextflow.config
└── workflow.nf
```

**This a summarized description of the files and directories found:**

- **`main.nf`** is the file we are going to invoke with the world-famous `nextflow run` command.
- **`workflow.nf`** is where all the magic happens, it stores the order of execution of tasks and how data should be handled.
- **`nextflow.config`**: you should know what this file does right? JK, with it we can manage different directives for workflow execution.
- **`modules`** is a really important folder since here we find dedicated files per each process of the pipeline.
- **`bin`** is the directory where we store customized scripts that can be run within a given process.
- **`data`** contains input data and related resources:
  - An indexed genome within the `yeast` folder representing the host genome to which we want to map the reads for contamination removal.
  - _viral_db_ is a directory that contains Kraken2 database necessary for both taxonomic annotation and species abundance re-estimation.
  - _samplesheet.csv_ lists the IDs and paths of the example data files, for processing in batches.
  - _samples_ directory is where the raw sequences are stored.
    The names correspond to accession numbers that you can search on the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)

!!!note

    Don't panic if this feels like a lot.
    This is just a glimpse of the material, and we are going to dig into each necessary file for the analysis in due time.

Now, to begin the course, click on the arrow in the bottom right corner of this page.
