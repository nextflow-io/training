# BioContainers

Another useful resource linking together Bioconda and containers is the [BioContainers](https://biocontainers.pro) project. BioContainers is a community initiative that provides a registry of container images for every Bioconda recipe.

So far, we’ve seen how to install packages with conda and micromamba, both locally and within containers. With BioContainers, you don’t need to create your own container image for the tools you want, and you don’t need to use conda or micromamba to install the packages. It already provides you with a Docker image containing the programs you want installed. For example, you can get the container image of fastqc using BioContainers with:

```bash
docker pull biocontainers/fastqc:v0.11.5
```

You can check the registry for the packages you want in [BioContainers official website](https://biocontainers.pro/registry).

Contrary to other registries that will pull the latest image when no tag (version) is provided, you must specify a tag when pulling BioContainers (after a colon `:`, e.g. fastqc:v0.11.5). Check the tags within the registry and pick the one that better suits your needs.

!!! tip

    You can have more complex definitions within your process block by letting the appropriate container image or conda package be used depending on if the user selected singularity, Docker or conda to be used. You can click [here](https://nf-co.re/docs/contributing/modules#software-requirements) for more information and [here](https://github.com/nf-core/modules/blob/61f68913fefc20241ceccb671b104230b2d775d7/modules/bowtie2/align/main.nf#L6-L9) for an example.

## :material-progress-question: Exercises

!!! exercise

    During the earlier RNA-Seq tutorial (script2.nf), we created an index with the salmon tool. Given we do not have salmon installed locally in the machine provided by Gitpod, we had to either run it with `-with-conda` or `-with-docker`. Your task now is to run it again `-with-docker`, but without having to create your own Docker container image. Instead, use the BioContainers image for salmon 1.7.0.


    ??? result

        ```bash
        nextflow run script2.nf -with-docker quay.io/biocontainers/salmon:1.7.0--h84f40af_0
        ```

!!! exercise "Bonus Exercise"

    Change the process directives in `script5.nf` or the `nextflow.config` file to make the pipeline automatically use BioContainers when using salmon, or fastqc.

    !!! tip "Hint"

        Temporarily comment out the line `#!groovy process.container = 'nextflow/rnaseq-nf'` in the `nextflow.config` file to make sure the processes are using the BioContainers that you set, and not the container image we have been using in this training.

    ??? result

        With these changes, you should be able to run the pipeline with BioContainers by running the following in the command line:

        ```bash
        nextflow run script5.nf
        ```

        with the following container directives for each process:

        ```groovy
        process FASTQC {
            container 'biocontainers/fastqc:v0.11.5'
            tag "FASTQC on $sample_id"
        ...
        ```

        and

        ```groovy
        process QUANTIFICATION {
            tag "Salmon on $sample_id"
            container 'quay.io/biocontainers/salmon:1.7.0--h84f40af_0'
            publishDir params.outdir, mode:'copy'
        ...
        ```

        Check the `.command.run` file in the work directory and ensure that the run line contains the correct Biocontainers.
