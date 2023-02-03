# Expression quantification

`script4.nf` adds a gene expression `QUANTIFICATION` process and call within the workflow scope. Quantification requires the index transcriptome and RNA-Seq read pair fastq files.

In the workflow scope, note how the `index_ch` channel is assigned as output in the `INDEX` process.

Next, note that the first input channel for the `QUANTIFICATION` process is the previously declared `index_ch`, which contains the `path` to the `salmon_index`.

Also, note that the second input channel for the `QUANTIFICATION` process, is the `read_pair_ch` we just created. This being a `tuple` composed of two elements (a value: `sample_id` and a list of paths to the fastq reads: `reads`) in order to match the structure of the items emitted by the `fromFilePairs` channel factory.

Execute it by using the following command:

```bash
nextflow run script4.nf -resume
```

You will see the execution of the `QUANTIFICATION` process.

When using the `-resume` option, any step that has already been processed is skipped.

Try to execute the same script again with more read files, as shown below:

```bash
nextflow run script4.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

You will notice that the `QUANTIFICATION` process is executed multiple times.

Nextflow parallelizes the execution of your pipeline simply by providing multiple sets of input data to your script.

!!! tip

    It may be useful to apply optional settings to a specific process using [directives](https://www.nextflow.io/docs/latest/process.html#directives) by specifying them in the process body.

## :material-progress-question: Exercises

!!! exercise

    Add a [tag](https://www.nextflow.io/docs/latest/process.html#tag) directive to the `QUANTIFICATION` process to provide a more readable execution log.

    ??? result

        Add the following before the input declaration:

        ```groovy
        tag "Salmon on $sample_id"
        ```

!!! exercise

    Add a [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directive to the `QUANTIFICATION` process to store the process results in a directory of your choice.

    ??? result

        Add the following before the `input` declaration in the `QUANTIFICATION` process:

        ```groovy
        publishDir params.outdir, mode:'copy'
        ```

## :material-check-all: Summary

In this step you have learned:

1. How to connect two processes together by using the channel declarations
2. How to `resume` the script execution and skip cached steps
3. How to use the `tag` directive to provide a more readable execution output
4. How to use the `publishDir` directive to store a process results in a path of your choice
