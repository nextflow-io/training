# Orientation

The training environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please the [Environment Setup](../../envsetup/) mini-course before going any further.

## Materials provided

For the purpose of the course, we'll be working in the `nf4-science/metagenomics/` directory, therefore let's move into it since here you will find all the code files, test data and accessory files you will need. Run the command:

```bash
cd nf4-science/metagenomics/
```

Now, let's take a look of the files contained in the directory with the command:

```bash
tree . -L 3
```

Here you should see the following directory structure:

```console title="Directory contents"
.
├── bin
│   └── report.Rmd
├── data
│   ├── oryza
│   │   ├── oryza.1.bt2
│   │   ├── oryza.2.bt2
│   │   ├── oryza.3.bt2
│   │   ├── oryza.4.bt2
│   │   ├── oryza.rev.1.bt2
│   │   └── oryza.rev.2.bt2
│   ├── samples
│   │   ├── ERR2143758
│   │   │   ├── ERR2143758_1.fastq
│   │   │   └── ERR2143758_2.fastq
│   │   ├── ERR2143763
│   │   │   ├── ERR2143763_1.fastq
│   │   │   └── ERR2143763_2.fastq
│   │   ├── ERR2143774
│   │   │   ├── ERR2143774_1.fastq
│   │   │   └── ERR2143774_2.fastq
│   │   └── ERR2143791
│   │       ├── ERR2143791_1.fastq
│   │       └── ERR2143791_2.fastq
│   └── viral_db
│       ├── database100mers.kmer_distrib
│       ├── database150mers.kmer_distrib
│       ├── database200mers.kmer_distrib
│       ├── database250mers.kmer_distrib
│       ├── database300mers.kmer_distrib
│       ├── database50mers.kmer_distrib
│       ├── database75mers.kmer_distrib
│       ├── hash.k2d
│       ├── inspect.txt
│       ├── ktaxonomy.tsv
│       ├── library_report.tsv
│       ├── names.dmp
│       ├── nodes.dmp
│       ├── opts.k2d
│       ├── seqid2taxid.map
│       └── taxo.k2d
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
├── samplesheet.csv
└── workflow.nf

```

!!!note

    Don't panic. This is just a glipmse of the material, and we are going to dig into each necessary file for the analysis.

**This a summarized description of the files and directories found:**

- **`main.nf`** is the file we are going to invoke with the worldwide famous `nextflow run` command.
- **`workflow.nf`** is where all the magic happens, it stores the order of execution of tasks and how data should be handled.
- **`nextflow.config`**, you should know what this file does right?
- **`modules`**, this is a really important folder since here we found dedicated files per each process of the pipeline.
- **`bin`**, in this directory we store customized scripts that can be run within a given process.
- **`data`** contains input data and related resources:
  - _An indexed genome_ within the `oryza` folder representing the host genome to which we want to map the reads for contamination removal.
  - _viral_db_ is directory that contains Kraken2 database necessary for both taxonomic annotation and species abundance re-estimation.
  - _samplesheet.csv_ listing the IDs and paths of the example data files, for processing in batches.
  - _samples_ directory is where the raw sequences are stored. The names correspond to accession numbers that you can search on the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)

Now, to begin the course, click on the arrow in the bottom right corner of this page.
