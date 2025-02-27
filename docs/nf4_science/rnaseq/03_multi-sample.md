# Part 3: Multi-sample paired-end implementation of RNAseq processing

In this final part of the course, we're going to take our simple workflow to the next level by turning it into a powerful batch automation tool to handle arbitrary numbers of samples.
And while we're at it, we're also going to switch it to expect paired-end data, which is more common in newer studies.

We'll do this in three stages:

1. Make the workflow accept multiple input samples and parallelize execution
2. Add comprehensive QC report generation
3. Switch to paired-end RNAseq data

---

## 1. Make the workflow accept multiple input samples and parallelize execution

We're going to need to change how we manage the input.

## 1.1. Change the primary input to be a CSV of file paths instead of a single file

We provide a CSV file containing sample IDs and FASTQ file paths in the `data/` directory.
This CSV file includes a header line.
Note that the FASTQ file paths are absolute paths.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Let's rename the primary input parameter to `input_csv` and change the default to be the path to the `single-end.csv` file.

```groovy title="rnaseq.nf" linenums="13"
// Primary input
params.input_csv = "data/single-end.csv"
```

## 1.2. Update the input channel factory to handle a CSV as input

We're going to want to load the contents of the file into the channel instead of just the file path itself, so we use the `.splitCsv()` operator to parse the CSV format, then the `.map()` operator to grab the specific piece of information we want (the FASTQ file path).

```groovy title="rnaseq.nf" linenums="16"
    // Create input channel from the contents of a CSV file
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

## 1.3. Run the workflow to test that it works

```bash
nextflow run rnaseq.nf
```

This time we see each step gets run 6 times, on each of the 6 data files we provided.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

executor >  local (18)
[07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
[cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
[68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
```

That's all it took to get the workflow to run on multiple files!
Nextflow handles all the parallelism for us.

---

## 2. Aggregate pre-processing QC metrics into a single MultiQC report

All this produces a lot of QC reports, and we don't want to have to dig through individual reports.
This is the perfect point to put in a MultiQC report aggregation step!

### 2.1. Create a module for the QC aggregation process

Let's create a module file called `modules/multiqc.nf` to house the `MULTIQC` process:

```bash
touch modules/multiqc.nf
```

Open the file in the code editor and copy the following code into it:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c"
    publishDir "results/multiqc", mode: 'copy'

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

### 2.2. Import the module into the workflow file

Add the statement `include { MULTIQC } from './modules/multiqc.nf'` to the `rnaseq.nf` file:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Add a `report_id` parameter and give it a sensible default

```groovy title="rnaseq.nf" linenums="9"
/*
 * Pipeline parameters
 */
params.hisat2_index_zip = "data/genome_index.tar.gz"
params.report_id = "all_single-end"
```

### 2.4. Call the process on the outputs of the previous steps

We need to give the `MULTIQC` process all the QC-related outputs from previous steps.

For that, we're going to use the `.mix()` operator, which aggregates multiple channels into a single one.

If we had four processes called A, B, C and D with a simple `.out` channel each, the syntax would look like this: `A.out.mix( B.out, C.out, D.out )`. As you can see, you apply it to the first of the channels you want to combine (doesn't matter which) and just add all the others, separated by commas, in the parenthesis that follows.

In the case of our workflow, we have the following outputs to aggregate:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

So the syntax example becomes:

```groovy title="Applying .mix() in the MULTIQC call"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

To this, we need to append the `collect()` operator in order to pull the reports for all the samples into a single call to `MULTIQC`, and also give it the `report_id` parameter:

```groovy title="The completed MULTIQC call" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

In the context of the full workflow block, it ends up looking like this:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Create input channel from the contents of a CSV file
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Run the workflow to test that it works

```bash
nextflow run rnaseq.nf -resume
```

This time we see a single call to MULTIQC added after the cached process calls:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

executor >  local (1)
[07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
[2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
[a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
[56/e1f102] MULTIQC          [100%] 1 of 1 ✔
```

You can find the outputs under `results/trimming` as specified in the `TRIM_GALORE` process by the `publishDir` directive.

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

---

<!-- Bonus if there's time -->

## 3. Switch to paired-end RNAseq data

TODO

### 3.1 Modify the default `input_csv` to use the paired-end data

We provide a second CSV file containing sample IDs and paired FASTQ file paths in the `data/` directory

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Let's change the `input_csv` default to be the path to the `paired-end.csv` file.

```groovy title="rnaseq.nf" linenums="13"
// Primary input
params.input_csv = "data/paired-end.csv"
```

## 3.1. Modify channel factory to grab the second file of reads

TODO

## 3.2. Generalize the FASTQC process

TODO

## 3.3. Adapt the TRIM_GALORE process to expect paired-end reads

TODO

## 3.4. Update the call to the MULTIQC process to expect two reports from TRIM_GALORE

TODO

## 3.5. Adapt the HISAT2 process to expect paired-end reads

TODO

## 3.6. Run the workflow to test that it works

TODO

---

### Takeaway

You know how to adapt a single-sample workflow to parallelize processing of multiple samples, generate a comprehensive QC report and switch to using paired-end read data if needed.

_Making the workflow accept either data type on the fly is out of scope for this training._

### What's next?

Congratulations, you've completed the Nextflow For RNAseq mini-course! Celebrate your success and take a well deserved break!

Next, we ask you to complete a very short survey about your experience with this training course, then we'll take you to a page with links to further training resources and helpful links.
