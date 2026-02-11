# Part 3: Multi-sample paired-end implementation

Previously, you built a per-sample variant calling pipeline that processed each sample's data independently.
In this part of the course, we're going to take our simple workflow to the next level by turning it into a powerful batch automation tool to handle arbitrary numbers of samples.
And while we're at it, we're also going to update it to expect paired-end data, which is more common in newer studies.

??? info "How to begin from this section"

    This section of the course assumes you have completed [Part 1: Method Overview](./01_method.md), [Part 2: Single-sample implementation](./02_single-sample.md) and have a working `rnaseq.nf` pipeline with filled-in module files.

    If you did not complete Part 2 or want to start fresh for this part, you can use the Part 2 solution as your starting point.
    Run these commands from inside the `nf4-science/rnaseq/` directory:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    This gives you a complete single-sample processing workflow.
    You can test that it runs successfully:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Assignment

In this part of the course, we're going to extend the workflow to do the following:

1. Read sample information from a CSV samplesheet
2. Run per-sample QC, trimming, and alignment on all samples in parallel
3. Aggregate all QC reports into a comprehensive MultiQC report

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

This automates the steps from the second section of [Part 1: Method overview](./01_method.md#2-multi-sample-qc-aggregation), where you ran these commands manually in their containers.

## Lesson plan

We've broken this down into three stages:

1. **Make the workflow accept multiple input samples.**
   This covers switching from a single file path to a CSV samplesheet, parsing it with `splitCsv()`, and running all existing processes on multiple samples.
2. **Add comprehensive QC report generation.**
   This introduces the `collect()` operator to aggregate outputs across samples, and adds a MultiQC process to produce a combined report.
3. **Switch to paired-end RNAseq data.**
   This covers adapting processes for paired-end inputs (using tuples), creating paired-end modules, and setting up a separate test profile.

This implements the method described in [Part 1: Method Overview](./01_method.md) (second section covering the multi-sample use case) and builds directly on the workflow produced by Part 2.

!!! tip

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Make the workflow accept multiple input samples

To run on multiple samples, we need to change how we manage the input: instead of providing a single file path, we'll read sample information from a CSV file.

We provide a CSV file containing sample IDs and FASTQ file paths in the `data/` directory.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

This CSV file includes a header line that names the columns.

Note that this is still single-end read data.

!!! warning

    The file paths in the CSV are absolute paths that must match your environment.
    If you are not running this in the training environment we provide, you will need to update the paths to match your system.

### 1.1. Change the primary input to a CSV of file paths in the test profile

First, we need to update the test profile in `nextflow.config` to provide the CSV file path instead of the single FASTQ path.

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Now we need to update the channel creation to read from this CSV.

### 1.2. Update the channel factory to parse CSV input

We need to load the contents of the file into the channel instead of just the file path itself.

We can do this using the same pattern we used in [Part 2 of Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): applying the [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) operator to parse the file, then a `map` operation to extract the FASTQ file path from each row.

=== "After"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Create input channel from the contents of a CSV file
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Create input channel from a file path
        read_ch = channel.fromPath(params.input)
    ```

One thing that is new compared to what you encountered in the Hello Nextflow course is that this CSV has a header line, so we add `#!groovy header: true` to the `splitCsv()` call.
That allows us to reference columns by name in the `map` operation: `#!groovy row.fastq_path` extracts the file path from the `fastq_path` column of each row.

The input handling is updated and the workflow is ready to test.

### 1.3. Run the workflow

The workflow now reads sample information from a CSV file and processes all samples in parallel.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

This time each step gets run 6 times, once for each sample in the CSV file.

That's all it took to get the workflow to run on multiple files.
Nextflow handles all the parallelism for us.

### Takeaway

You know how to switch from a single-file input to CSV-based multi-sample input that Nextflow processes in parallel.

### What's next?

Add a QC report aggregation step that combines metrics from all samples.

---

## 2. Aggregate pre-processing QC metrics into a single MultiQC report

All this produces a lot of QC reports, and we don't want to have to dig through individual reports.
This is the perfect point to put in a MultiQC report aggregation step.

Recall the `multiqc` command from [Part 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

The command scans the current directory for recognized QC output files and aggregates them into a single HTML report.
The container URI was `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

We need to set up an additional parameter, prepare the inputs, write the process, wire it in, and update the output handling.

### 2.1. Set up the inputs

The MultiQC process needs a report name parameter and the collected QC outputs from all previous steps bundled together.

#### 2.1.1. Add a `report_id` parameter

Add a parameter to name the output report.

=== "After"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Primary input
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path

        // Report ID
        report_id: String
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primary input
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

Add the report ID default to the test profile:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Now we need to prepare the inputs for the MultiQC process.

#### 2.1.2. Collect and combine QC outputs from previous steps

We need to give the `MULTIQC` process all the QC-related outputs from previous steps bundled together.

For that, we use the `.mix()` operator, which aggregates multiple channels into a single one.
We start from `channel.empty()` and mix in all the output channels we want to combine.
This is cleaner than chaining `.mix()` onto one of the output channels directly, because it treats all inputs symmetrically.

In our workflow, the QC-related outputs to aggregate are:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

We mix them into a single channel, then use `.collect()` to aggregate the reports across all samples into a single list.

Add these lines to the workflow body after the `HISAT2_ALIGN` call:

=== "After"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Alignment to a reference genome
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Comprehensive QC report generation
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="38"
        // Alignment to a reference genome
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

Using intermediate variables makes each step clear: `multiqc_files_ch` contains all individual QC files mixed into one channel, and `multiqc_files_list` is the collected bundle ready to pass to MultiQC.

### 2.2. Write the QC aggregation process and call it in the workflow

As before, we need to fill in the process definition, import the module, and add the process call.

#### 2.2.1. Fill in the module for the QC aggregation process

Open `modules/multiqc.nf` and examine the outline of the process definition.

Go ahead and fill in the process definition by yourself using the information provided above, then check your work against the solution in the "After" tab below.

=== "Before"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

This process uses `#!groovy path '*'` as the input qualifier for the QC files.
The `'*'` wildcard tells Nextflow to stage all the collected files into the working directory without requiring specific names.
The `val output_name` input is a string that controls the report filename.

The `multiqc .` command scans the current directory (where all the staged QC files are) and generates the report.

Once you've completed this, the process is ready to use.

#### 2.2.2. Include the module

Add the import statement to `rnaseq.nf`:

=== "After"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Now add the process call to the workflow.

#### 2.2.3. Add the process call

Pass the collected QC files and the report ID to the `MULTIQC` process:

=== "After"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

The MultiQC process is now wired into the workflow.

### 2.3. Update the output handling

We need to add the MultiQC outputs to the publish declaration and configure where they go.

#### 2.3.1. Add publish targets for the MultiQC outputs

Add the MultiQC outputs to the `publish:` section:

=== "After"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

Now we need to tell Nextflow where to put these outputs.

#### 2.3.2. Configure the new output targets

Add entries for the MultiQC targets in the `output {}` block, publishing them into a `multiqc/` subdirectory:

=== "After"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Before"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

The output configuration is complete.

### 2.4. Run the workflow

We use `-resume` so that the previous processing steps are cached and only the new MultiQC step runs.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

A single call to MULTIQC has been added after the cached process calls.

You can find the MultiQC outputs in the results directory.

```bash
tree -L 2 results/multiqc
```

```console title="Output"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

That last `all_single-end.html` file is the full aggregated report, conveniently packaged into one easy to browse HTML file.

### Takeaway

You know how to collect outputs from multiple channels, bundle them with `.mix()` and `.collect()`, and pass them to an aggregation process.

### What's next?

Adapt the workflow to handle paired-end RNAseq data.

---

## 3. Enable processing paired-end RNAseq data

Right now our workflow can only handle single-end RNAseq data.
It's increasingly common to see paired-end RNAseq data, so we want to be able to handle that.

Making the workflow completely agnostic of the data type would require using slightly more advanced Nextflow language features, so we're not going to do that here, but we can make a paired-end processing version to demonstrate what needs to be adapted.

### 3.1. Copy the workflow and update the inputs

We start by copying the single-end workflow file and updating it for paired-end data.

#### 3.1.1. Copy the workflow file

Create a copy of the workflow file to use as a starting point for the paired-end version.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Now update the parameters and input handling in the new file.

#### 3.1.2. Add a paired-end test profile

We provide a second CSV file containing sample IDs and paired FASTQ file paths in the `data/` directory.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Add a `test_pe` profile to `nextflow.config` that points to this file and uses a paired-end report ID.

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

The test profile for paired-end data is ready.

#### 3.1.3. Update the channel factory

The `.map()` operator needs to grab both FASTQ file paths and return them as a list.

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Create input channel from the contents of a CSV file
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Create input channel from the contents of a CSV file
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

The input handling is configured for paired-end data.

### 3.2. Adapt the FASTQC module for paired-end data

Copy the module to create a paired-end version:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

The FASTQC process input doesn't need to change — when Nextflow receives a list of two files, it stages both and `reads` expands to both filenames.
The only change needed is in the output block: since we now get two FastQC reports per sample, we switch from `simpleName`-based patterns to wildcards.

=== "After"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Before"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

This generalizes the process in a way that makes it able to handle either single-end or paired-end data.

Update the import in `rnaseq_pe.nf` to use the paired-end version:

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

The FASTQC module and its import are updated for paired-end data.

### 3.3. Adapt the TRIM_GALORE module for paired-end data

Copy the module to create a paired-end version:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

This module needs more substantial changes:

- The input changes from a single path to a tuple of two paths
- The command adds the `--paired` flag and takes both read files
- The output changes to reflect Trim Galore's paired-end naming conventions, producing separate FastQC reports for each read file

=== "After"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Before"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

Update the import in `rnaseq_pe.nf`:

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

The TRIM_GALORE module and its import are updated for paired-end data.

### 3.4. Adapt the HISAT2_ALIGN module for paired-end data

Copy the module to create a paired-end version:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

This module needs similar changes:

- The input changes from a single path to a tuple of two paths
- The HISAT2 command changes from `-U` (unpaired) to `-1` and `-2` (paired) read arguments
- All uses of `reads.simpleName` change to `read1.simpleName` since we now reference a specific member of the pair

=== "After"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Before"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

Update the import in `rnaseq_pe.nf`:

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

The HISAT2_ALIGN module and its import are updated for paired-end data.

### 3.5. Update the MultiQC aggregation for paired-end outputs

The paired-end `TRIM_GALORE` process now produces two separate FastQC report channels (`fastqc_reports_1` and `fastqc_reports_2`) instead of one.
Update the `.mix()` block in `rnaseq_pe.nf` to include both:

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

The MultiQC aggregation now includes both sets of paired-end FastQC reports.

### 3.6. Update the output handling for paired-end outputs

The `publish:` section and `output {}` block also need to reflect the two separate FastQC report channels from the paired-end `TRIM_GALORE` process.

Update the `publish:` section in `rnaseq_pe.nf`:

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Update the corresponding entries in the `output {}` block:

=== "After"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Before"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

The paired-end workflow is now fully updated and ready to run.

### 3.7. Run the workflow

We don't use `-resume` since this wouldn't cache, and there's twice as much data to process than before, but it should still complete in under a minute.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Now we have two slightly divergent versions of our workflow, one for single-end read data and one for paired-end data.
The next logical step would be to make the workflow accept either data type on the fly, which is out of scope for this course, but we may tackle that in a follow-up.

---

### Takeaway

You know how to adapt a single-sample workflow to parallelize processing of multiple samples, generate a comprehensive QC report, and adapt the workflow to use paired-end read data.

### What's next?

Give yourself a big pat on the back! You have completed the Nextflow for RNAseq course.

Head on to the final [course summary](./next_steps.md) to review what you learned and find out what comes next.
