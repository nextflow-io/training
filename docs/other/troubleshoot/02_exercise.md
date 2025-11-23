# Exercise 2

Move into the exercise 2 directory and execute the `hello-gatk.nf` script.

```bash
cd /workspaces/training/troubleshoot/exercise2
nextflow run hello-gatk.nf
```

!!! warning "Error message"

    Although there is no error message, the `GATK_HAPLOTYPECALLER` process is only being run once not three times as expected.

    ```console
    executor >  local (5)
    [d2/839095] process > SAMTOOLS_INDEX (1)       [100%] 3 of 3 ✔
    [1a/2dafee] process > GATK_HAPLOTYPECALLER (1) [100%] 1 of 1 ✔
    [59/278fc4] process > GATK_JOINTGENOTYPING (1) [100%] 1 of 1 ✔
    ```

??? Solution

    A common mistake by developers is having the wrong channel type.

    Nextflow distinguishes two different kinds of channels: queue channels and value channels.

    - A **queue** channel is an asynchronous unidirectional FIFO queue that connects two processes or operators. Channel elements are consumed till the channel is empty.
    - A **value** channel (a.k.a. a singleton channel) is bound to a single value and it can be read unlimited times without consuming its contents.

    It is likely that one or more channel type has been converted to a *queue* channel and is only being read once.

    In this example, the `GATK_HAPLOTYPECALLER` process is only executing one task.

    You can see that the channel factory `channel.of` has been used.

    While parameters are treated as value channels, the `channel.of` channel factory has converted these inputs to `queue` channels and are only being read once.

    By removing these the inputs will be read multiple times.

    Currently:
    ```console title="hello-gatk.nf" linenums="111"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        channel.of(params.genome_reference),
        channel.of(params.genome_reference_index),
        channel.of(params.genome_reference_dict),
        channel.of(params.calling_intervals)
    )
    ```

    After:
    ```console
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        params.genome_reference,
        params.genome_reference_index,
        params.genome_reference_dict,
        params.calling_intervals
    )
    ```
