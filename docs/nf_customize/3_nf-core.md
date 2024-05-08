# What is nf-core?

![Gitpod welcome](img/nf-core-logo.png)

nf-core is a community effort to collect a curated set of analysis pipelines built using Nextflow.

nf-core provides a standardised set of best practices, guidelines, and templates for building and sharing bioinformatics pipelines. These pipelines are designed to be modular, scalable, and portable, allowing researchers to easily adapt and execute them using their own data and compute resources.

The community is a diverse group of bioinformaticians, developers, and researchers from around the world who collaborate on developing and maintaining a growing collection of high-quality pipelines. These pipelines cover a range of applications, including transcriptomics, proteomics, and metagenomics.

One of the key benefits of nf-core is that it promotes open development, testing, and peer review, ensuring that the pipelines are robust, well-documented, and validated against real-world datasets. This helps to increase the reliability and reproducibility of bioinformatics analyses and ultimately enables researchers to accelerate their scientific discoveries.

nf-core is published in Nature Biotechnology: [Nat Biotechnol 38, 276–278 (2020). Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x)

**Key Features of nf-core pipelines**

-   **Documentation**
    -   nf-core pipelines have extensive documentation covering installation, usage, and description of output files to ensure that you won't be left in the dark.
-   **CI Testing**
    -   Every time a change is made to the pipeline code, nf-core pipelines use continuous-integration testing to ensure that nothing has broken.
-   **Stable Releases**
    -   nf-core pipelines use GitHub releases to tag stable versions of the code and software, making pipeline runs totally reproducible.
-   **Packaged software**
    -   Pipeline dependencies are automatically downloaded and handled using Docker, Singularity, Conda, or other software management tools. There is no need for any software installations.
-   **Portable and reproducible**
    -   nf-core pipelines follow best practices to ensure maximum portability and reproducibility. The large community makes the pipelines exceptionally well-tested and easy to execute.
-   **Cloud-ready**
    -   nf-core pipelines are tested on AWS after every major release. You can even browse results live on the website and use outputs for your own benchmarking.

It is important to remember all nf-core pipelines are open-source, community driven, and are under active community development.

## nf-core tools

This workshop will make use of nf-core tools, a set of helper tools for use with Nextflow pipelines. These tools have been developed to provide a range of additional functionality for using, developing, and testing pipelines.

??? tip "How to download nf-core tools"

    nf-core tools is written in Python and is available from the [Python Package Index (PyPI)](https://pypi.org/project/nf-core/):

    ```bash
    pip install nf-core
    ```

The nf-core `--version` option can be used to print your version of nf-core tools:

```bash
nf-core --version
```

!!! question "Exercise"

    Find out what version of nf-core tools you have available using the `nf-core --version` command.

nf-core tools are for everyone and has commands to help both users and developers.

For users, the tools make it easier to execute pipelines.

For developers, the tools make it easier to develop and test your pipelines using best practices.

!!! question "Exercise"

    Find out what nf-core commands and options are available using the `nf-core --help` option.

## nf-core pipelines

There are currently >100 pipelines (April 2024) available as part of nf-core.

These pipelines are at various stages of development with 60 released, 34 under development, and 11 archived.

The [nf-core website](https://nf-co.re/) has a full list of pipelines, as well as their documentation, which can be explored.

Each pipeline has a dedicated page that includes expansive documentation that is split into 6 sections:

-   **Introduction:** An introduction and overview of the pipeline
-   **Usage:** Descriptions of how to execute the pipeline
-   **Parameters:** Grouped pipeline parameters with descriptions
-   **Output:** Descriptions and examples of the expected output files
-   **Results:** Example output files generated from the full test dataset
-   **Releases & Statistics:** pipeline version history and statistics

## Pulling an nf-core pipeline

Unless you intend to develop an nf-core pipeline independently, you do not need to clone a copy of a pipeine.

Instead, you can use Nextflow’s built-in functionality to `pull` a pipeline:

```bash
nextflow pull nf-core/demo
```

!!! note "nextflow run"

    The `nextflow run` command will also automatically `pull` the pipeline if it was not already available locally.

Nextflow will `pull` the default git branch if a pipeline version is not specified. This will be the master branch for nf-core pipelines with a stable release.

nf-core pipelines use GitHub releases to tag stable versions of the code and software. You will always be able to execute a previous version of a pipeline once it is released using the `-revision` or `-r` flag.

!!! question "Exercise"

    Pull the dev of the `nf-core/demo` pipeline directly from GitHub:

    ```bash
    nextflow pull nf-core/demo -r dev
    ```

    Check that it has been pulled by listing your cached pipelines:

    ```bash
    nextflow list
    ```
