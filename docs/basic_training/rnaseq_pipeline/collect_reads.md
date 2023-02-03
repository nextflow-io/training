---
title: Collect reads
---

# Collect read files by pairs

This step shows how to match **read** files into pairs, so they can be mapped by **Salmon**.

Edit the script `script3.nf` by adding the following statement as the last line of the file:

```groovy
read_pairs_ch.view()
```

Save it and execute it with the following command:

```bash
nextflow run script3.nf
```

It will print something similar to this:

```bash
[gut, [/.../data/ggal/gut_1.fq, /.../data/ggal/gut_2.fq]]
```

The above example shows how the `read_pairs_ch` channel emits tuples composed of two elements, where the first is the read pair prefix and the second is a list representing the actual files.

Try it again specifying different read files by using a glob pattern:

```bash
nextflow run script3.nf --reads 'data/ggal/*_{1,2}.fq'
```

!!! warning

    File paths that include one or more wildcards ie. `*`, `?`, etc., MUST be wrapped in single-quoted characters to avoid Bash expanding the glob.

## :material-progress-question: Exercises

!!! exercise

    Use the [set](https://www.nextflow.io/docs/latest/operator.html#set) operator in place of `=` assignment to define the `read_pairs_ch` channel.

    ??? result

        ```groovy
        Channel
            .fromFilePairs( params.reads )
            .set { read_pairs_ch }
        ```

!!! exercise

    Use the `checkIfExists` option for the [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) method to check if the specified path contains file pairs.

    ??? result

        ```groovy
        Channel
            .fromFilePairs( params.reads, checkIfExists: true )
            .set { read_pairs_ch }
        ```

## :material-check-all: Summary

In this step you have learned:

1. How to use `fromFilePairs` to handle read pair files
2. How to use the `checkIfExists` option to check for the existence of input files
3. How to use the `set` operator to define a new channel variable

!!! info

    The declaration of a channel can be before the workflow scope of within it. As long as it is upstream of the process that requires the specific channel.
