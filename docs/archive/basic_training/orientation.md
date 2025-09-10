---
title: Orientation
description: How to set up a development environment to run Nextflow
---

# Orientation

The GitHub Codespaces environment contains some test data that will be used in this workshop.

!!! note

    Follow [this link](../envsetup/index.md) if you have not yet setup your GitHub Codespaces environment.

## Getting started

You will complete this module in the `nf-training/` folder.

In this folder you will find a series of data files (`ggal`, `index`, `meta`...) and several script and configuration files.

```console
.
├── data
│   ├── ggal
│   │   └── <data files>
│   ├── index
│   │   └── <data files>
│   ├── meta
│   │   └── <data files>
│   ├── prots
│   │   └── <data files>
│   ├── reads
│   │   └── <data files>
│   └── test
│       └── <data files>
├── env.yml
├── hello.nf
├── hello_py.nf
├── modules.hello.nf
├── nextflow.config
├── script1.nf
├── script2.nf
├── script3.nf
├── script4.nf
├── script5.nf
├── script6.nf
├── script7.nf
└── snippet.nf
```

Each file will be used in this training module.

## Selecting a Nextflow version

By default, Nextflow will pull the latest stable version into your environment.

However, Nextflow is constantly evolving as improvements are made.

The latest releases can be viewed on GitHub [here](https://github.com/nextflow-io/nextflow/releases).

If you want to use a specific version of Nextflow, you can set the `NXF_VER` variable as shown below:

```bash
export NXF_VER=23.10.1
```

You can double-check `NXF_VER` by running:

```bash
nextflow -version
```

!!! question "Exercise"

    Open the GitHub Codespaces training environment and use the following command to switch to the `nf-customize` folder. View the files in this folder using the `tree` command:

    ```bash
    cd /workspaces/training/nf-training
    tree .
    ```
