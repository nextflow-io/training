# Splitting and Grouping: Challenges

The bonus exercise at the end of the [Splitting and Grouping](index.md) side quest can be solved with what the side quest covers.
The challenges on this page go further on purpose: each one asks you to change the side quest's workflow with a Nextflow feature that the side quest does not cover.

Each challenge starts from the workflow as it stands at the end of the side quest.
From the `side-quests/splitting_and_grouping` directory, copy the solution file to `splitting-challenge.nf`, as shown in each challenge, and make your changes in that file.
Your own `main.nf` stays as it is.

The hint for each challenge introduces the new feature and links to the Nextflow documentation.
Each solution ends with a pointer to the training material that covers the feature in depth.

| Challenge                                                                                | New feature             | Covered in depth in                                                                             |
| ---------------------------------------------------------------------------------------- | ----------------------- | ----------------------------------------------------------------------------------------------- |
| [1. Split the samples in a single pass](#challenge-1-split-the-samples-in-a-single-pass) | The `branch()` operator | [Essential Scripting Patterns](../essential_scripting_patterns/index.md#51-routing-with-branch) |

The solution files are in `solutions/splitting_and_grouping/`, one directory per challenge.
To run a solution directly, launch it from the `side-quests/splitting_and_grouping` directory:

```bash
nextflow run ../solutions/splitting_and_grouping/challenge-1/main.nf
```

---

## Challenge 1: Split the samples in a single pass

Start from the workflow as it stands at the end of the side quest:

```bash
cp ../solutions/splitting_and_grouping/main.nf splitting-challenge.nf
```

The workflow goes through the samples twice to separate them: one `filter` keeps the normal samples, another keeps the tumor samples, and each is followed by the same `map`.
Change `splitting-challenge.nf` so that a single operator splits the samples into normal and tumor channels in one pass, and reshapes them as it does.
The rest of the workflow, and its output, should stay the same.

??? tip "Hint"

    The [`branch()`](https://www.nextflow.io/docs/latest/reference/operator.html#branch) operator sends each item of a channel to one of several named outputs, based on conditions you define.
    You can then use each output as a channel with a dot, as in `ch_split_samples.normal`.
    A condition can also be followed by a `return` statement that sets the value sent to that output, so you can reshape the items as you route them.

??? solution

    Replace the named closure and the two `filter` and `map` chains with a single `branch`:

    === "After"

        ```groovy title="splitting-challenge.nf" linenums="8" hl_lines="2-8"
            ch_split_samples = ch_samples
                .branch { meta, file ->
                    normal: meta.type == 'normal'
                        return [meta.subMap(['id', 'repeat']), file]
                    tumor: meta.type == 'tumor'
                        return [meta.subMap(['id', 'repeat']), file]
                }
            ch_joined_samples = ch_split_samples.normal.join(ch_split_samples.tumor)
        ```

    === "Before"

        ```groovy title="splitting-challenge.nf" linenums="8"
            getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }

            ch_normal_samples = ch_samples
                .filter { meta, file -> meta.type == 'normal' }
                .map(getSampleIdAndReplicate)
            ch_tumor_samples = ch_samples
                .filter { meta, file -> meta.type == 'tumor' }
                .map(getSampleIdAndReplicate)
            ch_joined_samples = ch_normal_samples.join(ch_tumor_samples)
        ```

    The `normal` and `tumor` conditions route each sample to an output, and the `return` after each condition reshapes it, so the separate `map` steps and the named closure are no longer needed.
    `ch_split_samples.normal` and `ch_split_samples.tumor` take the place of the two filtered channels.

    Run the workflow:

    ```bash
    nextflow run splitting-challenge.nf
    ```

    The output is the same as at the end of the side quest:

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `splitting-challenge.nf` [condescending_liskov] revision: 1424555bd0

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

    To learn more about `branch()`, including how to route samples to different processes, see [Routing with `.branch()`](../essential_scripting_patterns/index.md#51-routing-with-branch) in the Essential Scripting Patterns side quest.
