# Part 6: Troubleshoot input validation errors

In [Part 2](./02_configure_execution.md#14-input-validation), you saw that nf-core pipelines validate their input samplesheet at launch and stop if it contains errors.
In [Part 5](./05_samplesheets.md), you wrote a samplesheet that passed validation on the first try.

Real samplesheets rarely do.
They are often assembled by hand or exported from a spreadsheet, and they collect typos, missing values and inconsistencies along the way.

In this part, you diagnose and fix two broken samplesheets for `nf-core/rnaseq`.
The first contains errors that the input schema catches, and the second contains an error that only the pipeline itself can detect.

---

## 1. Fix errors caught by the input schema

The input schema you examined in [Part 5](./05_samplesheets.md#12-examine-the-input-schema) defines rules for each column.
At launch, the pipeline checks every row of the samplesheet against those rules and reports every violation it finds.

### 1.1. Run the pipeline with a broken samplesheet

Your working directory contains a samplesheet named `broken_samplesheet_1.csv`.

```csv title="broken_samplesheet_1.csv" linenums="1"
--8<-- "nfcore-use/broken_samplesheet_1.csv"
```

Run the pipeline with it, using `-preview` as in Part 5 so that no tasks run:

```bash
nextflow run nf-core/rnaseq -profile test,docker --input broken_samplesheet_1.csv --outdir results -preview
```

```console title="Command output (excerpt)"
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (broken_samplesheet_1.csv): Validation of file failed:
    -> Entry 2: Error for field 'fastq_1' (data/reads/WT_1_S1_L02_R1_001.fastq.gz): the file or directory 'data/reads/WT_1_S1_L02_R1_001.fastq.gz' does not exist (FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq', '.fastq', '.fq.gz' or '.fastq.gz')
    -> Entry 3: Error for field 'sample' (WT 2): "WT 2" does not match regular expression [^\S+$] (Sample name must be provided and cannot contain spaces)
    -> Entry 4: Error for field 'strandedness' (reversed): Expected any of [[forward, reverse, unstranded, auto]] (Strandedness must be provided and be one of 'auto', 'forward', 'reverse' or 'unstranded')
    -> Entry 5: Missing required field(s): strandedness
```

The pipeline stops before running anything and lists all four problems at once.

### 1.2. Interpret the error messages

Each line identifies the problem by **entry**, which is the row number counting from the first row after the header.
So entry 1 is line 2 of the file.

Each message has three parts: the field that failed, the value it contained, and the explanation.
The explanation in parentheses is the `errorMessage` defined in the input schema, while the text before it describes which schema rule was violated:

| Entry | Field          | Rule violated                               | Schema rule    |
| ----- | -------------- | ------------------------------------------- | -------------- |
| 2     | `fastq_1`      | The file does not exist                     | `exists: true` |
| 3     | `sample`       | The value contains a space                  | `pattern`      |
| 4     | `strandedness` | The value is not one of the accepted values | `enum`         |
| 5     | `strandedness` | The value is empty for a required column    | `required`     |

The message for entry 2 includes the full `errorMessage` for `fastq_1`, which mentions spaces and file extensions.
The rule that actually failed is stated before it: the file does not exist.
Compare the path with the files in `data/reads/` to find the cause.

### 1.3. Fix the samplesheet

??? exercise

    Fix the four errors in `broken_samplesheet_1.csv`, then re-run the command from 1.1 to confirm that the samplesheet passes validation.

??? solution

    The fixed rows are highlighted:

    ```csv title="broken_samplesheet_1.csv" linenums="1" hl_lines="3-6"
    sample,fastq_1,fastq_2,strandedness
    WT_1,data/reads/WT_1_S1_L001_R1_001.fastq.gz,data/reads/WT_1_S1_L001_R2_001.fastq.gz,auto
    WT_1,data/reads/WT_1_S1_L002_R1_001.fastq.gz,data/reads/WT_1_S1_L002_R2_001.fastq.gz,auto
    WT_2,data/reads/WT_2_S2_L001_R1_001.fastq.gz,data/reads/WT_2_S2_L001_R2_001.fastq.gz,auto
    KO_1,data/reads/KO_1_S3_L001_R1_001.fastq.gz,data/reads/KO_1_S3_L001_R2_001.fastq.gz,reverse
    KO_2,data/reads/KO_2_S4_L001_R1_001.fastq.gz,data/reads/KO_2_S4_L001_R2_001.fastq.gz,auto
    ```

    - Line 3: the lane in the `fastq_1` path was `L02` instead of `L002`.
    - Line 4: the space in `WT 2` is replaced by an underscore.
    - Line 5: `reversed` is replaced by `reverse`.
    - Line 6: the empty strandedness value is set to `auto`.

    Re-running the command now reports `Pipeline completed successfully`.

### Takeaway

Schema validation reports every error in the samplesheet at once, identifying each by row and column.
Reading the rule that failed, rather than only the custom error message, tells you what to fix.

### What's next?

Some errors involve more than one row and are not visible to the schema.

---

## 2. Fix errors caught by the pipeline

The input schema checks each row on its own.
It cannot express rules that involve several rows, such as "all runs of the same sample must have the same strandedness".
Pipelines that need such rules implement them in their own code, which runs after schema validation.

### 2.1. Validate a samplesheet that hides a cross-row error

Your working directory also contains `broken_samplesheet_2.csv`.

```csv title="broken_samplesheet_2.csv" linenums="1"
--8<-- "nfcore-use/broken_samplesheet_2.csv"
```

Check it with `-preview`:

```bash
nextflow run nf-core/rnaseq -profile test,docker --input broken_samplesheet_2.csv --outdir results -preview
```

```console title="Command output (excerpt)"
-[nf-core/rnaseq] Pipeline completed successfully-
```

Every row is individually valid, so the samplesheet passes schema validation.

### 2.2. Run the pipeline's own checks

The pipeline's own checks run as part of the workflow, which `-preview` skips.
To reach them without spending time on real computation, use the `-stub` option.
It runs the pipeline with placeholder commands that produce empty output files instead of running the actual tools.

```bash
nextflow run nf-core/rnaseq -profile test,docker --input broken_samplesheet_2.csv --outdir results -stub
```

```console title="Command output (excerpt)"
ERROR ~ Please check input samplesheet -> Multiple runs of a sample must have the same strandedness!: WT_1
```

The pipeline stops shortly after launch, before any task starts.

!!! note

    `-stub` only works for pipelines whose processes define a stub, which is required for nf-core modules.
    It is a convenient way to test the logic of a pipeline run, but the outputs it produces are empty placeholders.

### 2.3. Fix the samplesheet

??? exercise

    Fix `broken_samplesheet_2.csv` and verify the fix with `-stub`.

??? solution

    The two runs of `WT_1` must declare the same strandedness:

    ```csv title="broken_samplesheet_2.csv" linenums="1" hl_lines="2"
    sample,fastq_1,fastq_2,strandedness
    WT_1,data/reads/WT_1_S1_L001_R1_001.fastq.gz,data/reads/WT_1_S1_L001_R2_001.fastq.gz,auto
    WT_1,data/reads/WT_1_S1_L002_R1_001.fastq.gz,data/reads/WT_1_S1_L002_R2_001.fastq.gz,auto
    WT_2,data/reads/WT_2_S2_L001_R1_001.fastq.gz,data/reads/WT_2_S2_L001_R2_001.fastq.gz,auto
    KO_1,data/reads/KO_1_S3_L001_R1_001.fastq.gz,data/reads/KO_1_S3_L001_R2_001.fastq.gz,auto
    KO_2,data/reads/KO_2_S4_L001_R1_001.fastq.gz,data/reads/KO_2_S4_L001_R2_001.fastq.gz,auto
    ```

    With the fix in place, the stub run proceeds past the samplesheet checks and starts submitting tasks.
    You can stop it with ++ctrl+c++ once tasks appear.

`nf-core/rnaseq` applies several checks of this kind to samples with multiple rows.
Runs of the same sample must also agree on whether they are single-end or paired-end, and on whether their files are compressed.
The pipeline's usage documentation describes the rules that apply to its samplesheet.

??? exercise "Apply it to the pipeline you chose in Part 4"

    Look at the input schema of the pipeline you chose in [Part 4](./04_choose_pipeline.md) and list the rules it enforces: required columns, patterns, accepted values.
    Which mistakes in your draft samplesheet from Part 5 would it catch?

    Then check the pipeline's usage documentation for rules that span multiple rows, such as consistency between samples that belong to the same individual or family.
    These are the checks that `-preview` will not reach.

### Takeaway

Schema validation catches errors within a row; pipeline-specific checks catch inconsistencies across rows.
`-preview` triggers only the first kind, while `-stub` reaches both without running real computation.

### What's next?

Use the reports generated by the pipeline to tune its compute resource allocations.

---

## Summary

In this part you learned to:

- Interpret schema validation errors by row, field and rule
- Fix common samplesheet errors: missing files, invalid names, invalid or missing values
- Use `-stub` to trigger pipeline-specific checks that span multiple rows
