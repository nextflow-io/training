# Orientation

The Gitpod environment contains some test data that will be used in this workshop.

!!! note

    Follow [this link](../envsetup/00_introduction.md) if you have not yet setup your Gitpod environment.

## Getting started

You will complete this module in the `nf-customize/` folder.

In this folder you will find three pairs of zipped fastq files (`*.fastq.gz`) in a `data/` folder and an example samplesheet (`samplesheet.csv`) in a `scripts/` folder.

```console
.
├── data
│   ├── gut_1.fastq.gz
│   ├── gut_2.fastq.gz
│   ├── liver_1.fastq.gz
│   ├── liver_2.fastq.gz
│   ├── lung_1.fastq.gz
│   └── lung_2.fastq.gz
└── scripts
    └── samplesheet.csv
```

Each file will be used in this training module.

!!! question "Exercise"

    Open the [Gitpod training environment](https://gitpod.io/#https://github.com/nextflow-io/training) and use the following command to switch to the `nf-customize` folder. View the files in this folder using the `tree` command:

    ```bash
    cd /workspace/gitpod/nf-customize
    tree .
    ```
