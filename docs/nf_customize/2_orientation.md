# Orientation

The Gitpod environment contains some test data that will be used in this workshop.

You will complete this module in the `nf-customize/` folder.

In this folder you will find three pairs of zipped fastq files (`*.fq.gz`) in a `data/` folder and an example samplesheet (`samplesheet.csv`) in a `scripts/` folder. Each will be used later in this module.

```console
.
├── data
│   ├── gut_1.fq.gz
│   ├── gut_2.fq.gz
│   ├── liver_1.fq.gz
│   ├── liver_2.fq.gz
│   ├── lung_1.fq.gz
│   └── lung_2.fq.gz
└── scripts
    └── samplesheet.csv
```

!!! question "Exercise"

    Open the [Gitpod training environment](https://gitpod.io/#https://github.com/nextflow-io/training) and use the following command to switch to the `nf-customize` folder. View the files in this folder using the `tree` command:

    ```bash
    cd /workspace/gitpod/nf-customize
    tree .
    ```
