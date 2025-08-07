# Orientation

The GitHub Codespaces environment contains some test data that will be used in this workshop.

!!! note

    Follow [this link](../envsetup/index.md) if you have not yet setup your GitHub Codespaces environment.

## Getting started

You will complete this module in the `nf-training-advanced/` folder.

In this folder you will find a series of folders that will be used during different sections of this training.

```console
nf-training-advanced
├── groovy
│   ├── main.nf
│   ├── modules
│   │   └── local
│   │       └── fastp
│   │           └── main.nf
│   └── nextflow.config
├── grouping
│   ├── data
│   │   ├── genome.fasta
│   │   ├── genome.fasta.fai
│   │   ├── intervals.bed
│   │   ├── reads
│   │   │   ├── treatmentA
│   │   │   │   └── <data files>
│   │   │   └── treatmentB
│   │   │       └── <data files>
│   │   ├── samplesheet.csv
│   │   └── samplesheet.ugly.csv
│   └── main.nf
├── metadata
│   ├── data
│   │   ├── reads
│   │   │   ├── treatmentA
│   │   │   │   └── <data files>
│   │   │   └── treatmentB
│   │   │       └── <data files>
│   │   ├── samplesheet.csv
│   │   └── samplesheet.ugly.csv
│   └── main.nf
├── operators
│   ├── data
│   │   ├── reads
│   │   │   └── <data files>
│   │   ├── samplesheet.csv
│   │   └── samplesheet.ugly.csv
│   └── main.nf
└── structure
    ├── lib
    │   └── Food.groovy
    ├── main.nf
    └── templates
        ├── adder.py
        └── demo_script.sh
```

## Selecting a Nextflow version

By default, Nextflow will pull the latest stable version into your environment.

However, Nextflow is constantly evolving as we make improvements and fix bugs.

The latest releases can be viewed on GitHub [here](https://github.com/nextflow-io/nextflow).

If you want to use a specific version of Nextflow, you can set the `NXF_VER` variable as shown below:

```bash
export NXF_VER=25.04.6
```

!!! Note

    This tutorial workshop requires `NXF_VER=23.10.0`, or later.

Run `nextflow -version` again to confirm that the change has taken effect.
