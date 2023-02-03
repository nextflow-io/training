---
description: How to set up a development environment to run Nextflow locally
---

# Local installation

Nextflow can be used on any POSIX-compatible system (Linux, macOS, Windows Subsystem for Linux, etc.).

### Requirements

-   Bash
-   [Java 11 (or later, up to 18)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
-   Git
-   [Docker](https://docs.docker.com/get-docker/)

### Optional requirements for this tutorial

-   [Singularity](https://github.com/sylabs/singularity) 2.5.x (or later)
-   [Conda](https://conda.io/) 4.5 (or later)
-   [Graphviz](http://www.graphviz.org/)
-   [AWS CLI](https://aws.amazon.com/cli/)
-   A configured AWS Batch computing environment

## Download Nextflow

Enter this command in your terminal:

```bash
wget -qO- https://get.nextflow.io | bash
```

Or, if you prefer curl:

```bash
curl -s https://get.nextflow.io | bash
```

Then ensure that the downloaded binary is executable:

```bash
chmod +x nextflow
```

AND put the nextflow executable into your `$PATH` (e.g. `/usr/local/bin` or `/bin/`)

## Docker

Ensure you have Docker Desktop running on your machine. Download Docker [here](https://docs.docker.com/get-docker/).

## Training material

You can view the training material here: <https://training.nextflow.io/>

To download the material use the command:

```bash
git clone https://github.com/nextflow-io/training.git
```

Then `cd` into the `nf-training` directory.

## Checking your installation

Check the correct installation of `nextflow` by running the following command:

```bash
nextflow info
```

This should show the current version, system and runtime.
