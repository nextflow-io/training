# Exercise 4

Move into the exercise 4 directory and execute the `hello-gatk.nf` script.

```bash
cd /workspaces/training/troubleshoot/exercise4
nextflow run hello-gatk.nf
```

!!! warning "Error message"

    ```console
    ERROR ~ Script compilation error
    - file : /workspaces/training/troubleshoot/exercise4/hello-gatk.nf
    - cause: Unexpected input: '{' @ line 100, column 10.
    workflow {
                ^

    1 error

    NOTE: If this is the beginning of a process or workflow, there may be a syntax error in the body, such as a missing or extra comma, for which a more specific error message could not be produced.

    -- Check '.nextflow.log' file for details
    ```

??? Solution

    There is no easy way to identify a syntax error in the body.

    As the error message suggests, this could be as simple as a missing comma.

    While there is no quick way to find the syntax error it is advisable to test regularly when you are developing your pipeline to find errors when you only have smaller sections of changes to check.

    After reviewing the workflow block, you will see that a comma is missing from the end of line 115.

    ```console title="hello-gatk.nf" linenums="111"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        params.genome_reference,
        params.genome_reference_index,
        params.genome_reference_dict
        params.calling_intervals
    )
    ```
