# Part 8: Hello nf-core

nf-core is a community effort to develop and maintain a curated set of analysis pipelines built using Nextflow.

![nf-core logo](img/nf-core-logo.png)

nf-core provides a standardized set of best practices, guidelines, and templates for building and sharing scientific pipelines.
These pipelines are designed to be modular, scalable, and portable, allowing researchers to easily adapt and execute them using their own data and compute resources.

One of the key benefits of nf-core is that it promotes open development, testing, and peer review, ensuring that the pipelines are robust, well-documented, and validated against real-world datasets.
This helps to increase the reliability and reproducibility of scientific analyses and ultimately enables researchers to accelerate their scientific discoveries.

nf-core is published in Nature Biotechnology: [Nat Biotechnol 38, 276–278 (2020). Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x). An updated preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.05.10.592912v1).

## nf-core pipelines and other components

The nf-core collection currently offers [over 100 pipelines](https://nf-co.re/pipelines/) in various stages of development, [72 subworkflows](https://nf-co.re/subworkflows/) and [over 1300 modules](https://nf-co.re/modules/) that you can use to build your own pipelines.

Each released pipeline has a dedicated page that includes 6 documentation sections:

-   **Introduction:** An introduction and overview of the pipeline
-   **Usage:** Descriptions of how to execute the pipeline
-   **Parameters:** Grouped pipeline parameters with descriptions
-   **Output:** Descriptions and examples of the expected output files
-   **Results:** Example output files generated from the full test dataset
-   **Releases & Statistics:** Pipeline version history and statistics

You should read the pipeline documentation carefully to understand what a given pipeline does and how it can be configured before attempting to run it.

### Pulling an nf-core pipeline

One really cool aspect of how Nextflow manages pipelines is that you can pull a pipeline from a GitHub repository without cloning the repository.
This is really convenient if you just want to run a pipeline without modifying the code.

So if you want to try out an nf-core pipeline with minimal effort, you can start by pulling it using the `nextflow pull` command:

```bash
nextflow pull nf-core/demo
```

Nextflow will `pull` the pipeline's default GitHub branch
For nf-core pipelines with a stable release, that will be the master branch.
You select a specific branch with `-r`; we'll cover that later.

```console title="Output"
Checking nf-core/demo ...
 downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
```

To be clear, you can do this with any Nextflow pipeline that is appropriately set up in GitHub, not just nf-core pipelines.
However nf-core is the largest open collection of Nextflow pipelines.

!!!tip

    One detail that sometimes trips people up is that the pipelines you pull this way are stored in a hidden assets folder:

    ```bash
    ls $HOME/.nextflow/assets/
    ```

    So you don't actually see them listed in your working directory.
    However, you can view a list of your cached pipelines using the `nextflow list` command:

    ```bash
    nextflow list
    ```

    ```console title="Output"
    nf-core/demo
    ```

Now that we've got the pipeline pulled, we can try running it!

### Trying out an nf-core pipeline with the test profile

Conveniently, every nf-core pipeline comes with a `test` profile.
This is a minimal set of configuration settings for the pipeline to run using a small test dataset that is hosted on the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository. It's a great way to try out a pipeline at small scale.

The `test` profile for `nf-core/demo` is shown below:

```groovy title="conf/test.config" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Limit resources so that this can run on GitHub Actions
    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

This tells us that the `nf-core/demo` `test` profile already specifies the input parameter, so you don't have to provide any input yourself.
However, the `outdir` parameter is not included in the `test` profile, so you have to add it to the execution command using the `--outdir` flag.

Here, we're also going to specify `-profile docker`, which by nf-core convention enables the use of Docker.

Try running it yourself:

```bash
nextflow run nf-core/demo -profile docker,test --outdir results
```

Annoyingly, we get the following output:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `https://github.com/nf-core/demo` [serene_stallman] DSL2 - revision: 04060b4644 [master]

Downloading plugin nf-schema@2.1.1
Nextflow version 24.02.0-edge does not match workflow required version: >=24.04.2
```

Apparently we are using a slightly older version of Nextflow than what the workflow requires.

This is annoying because we might not want to update Nextflow without checking that the update doesn't affect our current work. However, Nextflow makes it up to us by letting us request a specific version on a one-time basis!

You just have to add `NXF_VER=<version>` to the start of your command, like this:

```bash
NXF_VER=24.09.2-edge nextflow run nf-core/demo -profile docker,test --outdir results
```

And boom, that gets us past the version mismatch without committing us to any big changes.

Here's our pipeline output, by the way:

```console title="Output"
 N E X T F L O W   ~  version 24.09.2-edge

Launching `https://github.com/nf-core/demo` [naughty_bell] DSL2 - revision: 04060b4644 [master]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/demo 1.0.1
------------------------------------------------------
Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
  outdir                    : results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Core Nextflow options
  revision                  : master
  runName                   : naughty_bell
  containerEngine           : docker
  launchDir                 : /workspace/gitpod/hello-nextflow
  workDir                   : /workspace/gitpod/hello-nextflow/work
  projectDir                : /home/gitpod/.nextflow/assets/nf-core/demo
  userName                  : gitpod
  profile                   : docker,test
  configFiles               :

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------* The pipeline
  https://doi.org/10.5281/zenodo.12192442

* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/demo/blob/master/CITATIONS.md

executor >  local (7)
[0a/e694d8] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     [100%] 3 of 3 ✔
[85/4198c1] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) [100%] 3 of 3 ✔
[d8/fe153e] NFCORE_DEMO:DEMO:MULTIQC                 [100%] 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
Completed at: 28-Oct-2024 03:24:58
Duration    : 1m 13s
CPU hours   : (a few seconds)
Succeeded   : 7
```

Isn't that neat?

And that's all you need to know for now. Congratulations!
You have now run your first nf-core pipeline.

### Takeaway

You have a general idea of what nf-core offers and you know how to run an nf-core pipeline using its built-in test profile.

### What's next?

Celebrate and take another break! Next, we'll show you how to take advantage of Seqera Platform to launch and monitor your workflows more conveniently and efficiently on any compute infrastructure.
