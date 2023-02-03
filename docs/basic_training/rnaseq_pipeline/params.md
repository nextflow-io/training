---
title: Pipeline parameters
---

# Define the pipeline parameters

Parameters are inputs and options that can be changed when the pipeline is run.

The script `script1.nf` defines the pipeline input parameters.

```groovy
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"

println "reads: $params.reads"
```

Run it by using the following command:

```bash
nextflow run script1.nf
```

Try to specify a different input parameter in your execution command, for example:

```bash
nextflow run script1.nf --reads '/workspace/training/nf-training/data/ggal/lung_{1,2}.fq'
```

## :material-progress-question: Exercises

!!! exercise

    Modify the `script1.nf` by adding a fourth parameter named `outdir` and set it to a default path that will be used as the pipeline output directory.

    ??? result

        ```groovy
        params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.multiqc = "$projectDir/multiqc"
        params.outdir = "results"
        ```

!!! exercise

    Modify `script1.nf` to print all of the pipeline parameters by using a single `log.info` command as a [multiline string](https://www.nextflow.io/docs/latest/script.html#multi-line-strings) statement.

    !!! tip ""

        :material-lightbulb: See an example [here](https://github.com/nextflow-io/rnaseq-nf/blob/3b5b49f/main.nf#L41-L48).


    ??? result

        Add the following to your script file:

        ```groovy
        log.info """\
                    R N A S E Q - N F   P I P E L I N E
                    ===================================
                    transcriptome: ${params.transcriptome_file}
                    reads        : ${params.reads}
                    outdir       : ${params.outdir}
                    """
                    .stripIndent()
        ```

## :material-check-all: Summary

In this step you have learned:

1. How to define parameters in your pipeline script
2. How to pass parameters by using the command line
3. The use of `$var` and `${var}` variable placeholders
4. How to use multiline strings
5. How to use `log.info` to print information and save it in the log execution file
