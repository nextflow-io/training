# Adding modules

The nf-core pipeline template is a working pipeline and comes pre-configured with two modules:

-   [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): A tool that performs quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses that can be used to give a quick impression of your data.
-   [MultiQC](https://multiqc.info/): A modular tool to aggregate results from bioinformatics analyses across many samples into a single report.

## Testing your pipeline

The `test` profile can be used to check if your pipeline is still working during your development cycle. It is also used as a part of GitHub Actions to test your pipeline during pull requests.

The default template `test` profile leverages small test files that are stored in the nf-core [test data GitHub repository](https://github.com/nf-core/test-datasets) as inputs for the pipeline.

Additionally, the template comes with profiles for the management of software dependencies (e.g., `docker`, `singularity`, and `conda`). For pipelines with processes are shipped with containers/images/recipes, these profiles can be used to change the way dependencies are handled when you execute your pipeline.

!!! warning

    If `-profile` for managing software dependencies is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. **This is not recommended.**

Additional test profiles can be created to test different parts of your pipeline and can also be added to GitHub actions.

!!! question "Exercise"

    Run your pipeline with the `test` and `docker` profile.

    ```bash
    cd /workspace/gitpod/nf-develop
    nextflow run nf-core-myfirstpipeline -profile test,docker --outdir results
    ```

## Adding a new tool to your pipeline

Here, as an example, you want to trim low-quality bases from both ends of your fastq files using the [`seqtk trim`](https://github.com/lh3/seqtk) command.

The `seqtk trim` module will take fastq files from the sample sheet as inputs and will produce trimmed fastq files that can be used as an input for other tools and version information about the `seqtk` tools that will be mixed into the inputs for the `MultiQC` process.

<figure class="excalidraw">
--8<-- "docs/nf_template/img/pipeline.excalidraw.svg"
</figure>

While you could develop a module for this tool independently, you can save a lot of time and effort by leveraging nf-core modules and subworkflows.

nf-core modules and subworkflows are written and maintained by the nf-core community. They are designed to be flexible but may require additional configuration to suit different use cases. Currently, there are more than [1200 nf-core modules](https://nf-co.re/modules) and [60 nf-core subworkflows](https://nf-co.re/subworkflows) (April 2024) available.

Modules and subworkflows can be listed, installed, updated, removed, and patched using nf-core tooling.

### Working with branches

GitHub branches are used to isolate development work without affecting other branches in a repository. Each repository has one default branch, and can have multiple other branches.

You can merge updates from one branch into another branch using a pull request.

The `nf-core create` command initiated three branches: `main`, `dev`, and `TEMPLATE`.

In nf-core, the `main` branch is for stable releases and the `dev` branch is for merging feature branches together. This enables the `main` branch to remain fully functional while new features are developed in feature branches, collected in the `dev` branch, and then merged into `main` once they are ready.

<figure class="excalidraw">
--8<-- "docs/nf_template/img/branches.excalidraw.svg"
</figure>

Feature branches should be checked out from the `dev` branch.

!!! question "Exercise"

    Checkout a new feature branch named `myFeature` from the dev branch

    ```
    git checkout -b myFeature dev
    ```

You can find out more about working collaboratively with branches on the [GitHub documentation](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests).

!!! note "Executing revisions"

    Remote GitHub branches can be executed with Nextflow using the revision flag (e.g., `-r dev`).

### The `TEMPLATE` branch

The `TEMPLATE` branch is used by the `nf-core sync` command to integrate template changes to your pipeline. You should **never** modify the `TEMPLATE` branch as any changes will likely disrupt the syncing functionality.

### Installing the `fastp` module

The `nf-core modules list` command can be used to show the modules in your local pipeline or the nf-core remote repository.

```
nf-core modules list remote
```

The `nf-core modules install` command can be used to install the `fastp` module directly from the nf-core repository:

```
cd nf-core-myfirstpipeline
nf-core modules install
```

You can follow the prompts to find and install the module you are interested in:

```console
? Tool name: seqtk_trim
```

Once selected, the tooling will install the module in the `modules/nf-core/` folder and suggest code that you can add to your main workflow file (`workflows/mypipeline.nf`).

```console
INFO     Installing 'seqtk_trim'
INFO     Use the following statement to include this module:

include { SEQTK_TRIM } from '../modules/nf-core/seqtk/trim/main'

```

!!! question "Exercise"

    Run the `nf-core modules install` command to add the `seqtk_trim` module to your pipeline.

To enable reporting and reproducibility, modules and subworkflows from the nf-core repository are tracked using hashes in the `modules.json` file. When modules are installed or removed using the nf-core tooling the `modules.json` file will be automatically updated.

!!! question "Exercise"

    Open your `modules.json` file and see if the `seqtk_trim` module is being tracked.

### Adding a module to your pipeline

Although the module has been installed in your local pipeline repository, it is not yet added to your pipeline.

The suggested `include` statement needs to be added to your `workflows/mypipeline.nf` file and the process call (with inputs) needs to be added to the workflow block.

```groovy title="workflows/mypipeline.nf" linenums="7"
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { SEQTK_TRIM             } from '../modules/nf-core/seqtk/trim/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
```

!!! question "Exercise"

    Add the suggested `include` statement to your `mypipeline.nf` file.

    ```groovy title="workflows/mypipeline.nf" linenums="8"
    include { SEQTK_TRIM             } from '../modules/nf-core/seqtk/trim/main'
    ```

To add the `SEQTK_TRIM` module to your workflow you will need to check what inputs are required.

You can view the input channels for the module by opening the `./modules/nf-core/seqtk/trim/main.nf` file.

```groovy title="/modules/nf-core/seqtk/trim/main.nf" linenums="10"
input:
tuple val(meta), path(reads)
```

Each nf-core module also has a `meta.yml` file which describes the inputs and outputs. This meta file is rendered on the [nf-core website](https://nf-co.re/modules/seqtk_trim), or can be viewed using the `nf-core modules info` command.

!!! question "Exercise"

    Use the `nf-core modules info` command to view information for the `seqtk_trim` module

    ```
    nf-core modules info seqtk_trim
    ```

Using this module information you can work out what inputs are required for the `SEQTK_TRIM` process:

1.  `tuple val(meta), path(reads)`

    -   A tuple with a meta _map_ and a list of fastq _files_
    -   The channel `ch_samplesheet` used by the `FASTQC` process can be used as the reads input.

As only one input channel required and it already exists it can be added to your `mypipeline.nf` file without any additional channel creation or modifications.

!!! question "Exercise"

    Add the `SEQTK_TRIM` process to your `mypipeline.nf` file.

    ```groovy title="workflows/mypipeline.nf" linenums="40"
    //
    // MODULE: Run SEQTK_TRIM
    //
    SEQTK_TRIM (
        ch_samplesheet
    )
    ```

As with the inputs, you can view the outputs for the module by opening the `/modules/nf-core/seqtk/trim/main.nf` file and viewing the module metadata.

```groovy title="/modules/nf-core/seqtk/trim/main.nf" linenums="13"
output:
tuple val(meta), path("*.fastq.gz"), emit: reads
path "versions.yml"                , emit: versions
```

To help with organization and readability it is beneficial to create named output channels.

For `SEQTK_TRIM`, the `reads` output could be put into a channel named `ch_trimmed`.

```groovy title="workflows/mypipeline.nf" linenums="46"
ch_trimmed  = SEQTK_TRIM.out.reads
```

Similarly, it is beneficial immediately mix the versions of tools into the `ch_versions` channel so they can be used as an input for the `MULTIQC` process.

```groovy title="workflows/mypipeline.nf" linenums="47"
ch_versions = ch_versions.mix(SEQTK_TRIM.out.versions.first())
```

!!! question "Exercise"

    Create a channel named `ch_trimmed` from the `SEQTK_TRIM.out.reads` output mix the `SEQTK_TRIM.out.versions` output with the `ch_versions` channel.

    ```groovy title="workflows/mypipeline.nf" linenums="46"
    ch_trimmed  = SEQTK_TRIM.out.reads
    ch_versions = ch_versions.mix(SEQTK_TRIM.out.versions.first())
    ```

!!! note

    The `first` operator is used to emit the first item from `SEQTK_TRIM.out.versions` to avoid duplication.

### Additional configuration options

To prevent changing the nf-core modules, additional configuration options can be applied to a module using scopes within configuration files.

The configuration of modules is commonly added to the `modules.conf` file in the `conf` folder. Process selectors (e.g., `withName`) are used to apply configuration to modules selectively. Process selectors must be used within the `process` scope.

Extra configuration may also be applied as directives by using `args`. You can find many examples of how arguments are added to modules in nf-core pipelines, for example, the nf-core/rnaseq [modules.config](https://github.com/nf-core/rnaseq/blob/master/conf/modules.config) file.

!!! question "Exercise"

    Add this snippet to your `conf/modules.config` file to save the trimmed fastq files reports in folders named using `meta.id`.

    ```console title="conf/modules.config" linenums="31"
    withName: 'SEQTK_TRIM' {
        publishDir = [
            path: { "${params.outdir}/fq/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{fastq.gz}"
        ]
    }
    ```

!!! note "Closures"

    Closures can be used in configuration files to inject code evaluated at runtime.
