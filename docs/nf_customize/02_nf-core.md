# What is nf-core?

![nf-core logo](img/nf-core-logo.png)

nf-core is a community effort to collect a curated set of analysis pipelines built using Nextflow.

nf-core provides a standardised set of best practices, guidelines, and templates for building and sharing bioinformatics pipelines. These pipelines are designed to be modular, scalable, and portable, allowing researchers to easily adapt and execute them using their own data and compute resources.

One of the key benefits of nf-core is that it promotes open development, testing, and peer review, ensuring that the pipelines are robust, well-documented, and validated against real-world datasets. This helps to increase the reliability and reproducibility of bioinformatics analyses and ultimately enables researchers to accelerate their scientific discoveries.

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

nf-core is published in Nature Biotechnology: [Nat Biotechnol 38, 276–278 (2020). Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x). An updated preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.05.10.592912v1).

## nf-core pipelines

There are currently >100 nf-core pipelines. These pipelines are at various stages of development with 60 released, 34 under development, and 11 archived (April 2024).

The [nf-core website](https://nf-co.re/) contains a full list of pipelines, as well as their documentation, which can be explored.

![nf-core logo](img/pipelines.png)

Each released pipeline has a dedicated page that includes 6 documentation sections:

-   **Introduction:** An introduction and overview of the pipeline
-   **Usage:** Descriptions of how to execute the pipeline
-   **Parameters:** Grouped pipeline parameters with descriptions
-   **Output:** Descriptions and examples of the expected output files
-   **Results:** Example output files generated from the full test dataset
-   **Releases & Statistics:** pipeline version history and statistics

Each section should be explored by a user to understand what the pipeline does and how it can be configured.

## Pulling an nf-core pipeline

Unless you intend to develop an nf-core pipeline independently, you do not need to clone a copy of a pipeline. Instead, you can use Nextflow’s `pull` command:

```bash
nextflow pull nf-core/demo
```

!!! note "The `nextflow run` command"

    The `nextflow run` command will also automatically `pull` the pipeline if it had not been pulled.

Nextflow will `pull` the pipelines default GitHub branch if a pipeline version is not specified. This will be the master branch for nf-core pipelines with a stable release.

nf-core pipelines use GitHub releases to tag stable versions of the code and software. You will always be able to execute different versions of a pipeline using the `-revision` or `-r` option.

Similarly, you can use the `-r` option to specify a specific GitHub branch. For example, the `dev` branch of the `nf-core/demo` pipeline could be pulled with the command:

```
nextflow pull nf-core/demo -r dev
```

If updates to a remote pipeline have been made, the pull command can be used to update or revery your local copy.

!!! question "Exercise"

    Use nextflow to pull the `nf-core/demo` pipeline:

    ```bash
    nextflow pull nf-core/demo
    ```

    Use the list command to view your cached pipelines:

    ```bash
    nextflow list
    ```

    Pulled pipelines are stored in a hidden assets folder:

    ```bash
    ls $HOME/.nextflow/assets/
    ```
