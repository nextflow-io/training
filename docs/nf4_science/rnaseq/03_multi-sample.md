# Part 3: Multi-sample paired-end implementation

In this final part of the course, we're going to take our simple workflow to the next level by turning it into a powerful batch automation tool to handle arbitrary numbers of samples.
And while we're at it, we're also going to switch it to expect paired-end data, which is more common in newer studies.

We'll do this in three stages:

1. Make the workflow accept multiple input samples and parallelize execution
2. Add comprehensive QC report generation
3. Switch to paired-end RNAseq data

---

## 1. Make the workflow accept multiple input samples and parallelize execution

We're going to need to change how we manage the input.

### 1.1. Change the primary input to be a CSV of file paths instead of a single file

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

### 1.2. Update the input channel factory to handle a CSV as input

We're going to want to load the contents of the file into the channel instead of just the file path itself, so we use the `.splitCsv()` operator to parse the CSV format, then the `.map()` operator to grab the specific piece of information we want (the FASTQ file path).

```groovy title="rnaseq.nf" linenums="16"
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Run the workflow to test that it works

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

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

That will collect QC reports per sample.
But since we want to aggregate them across all samples, we need to add the `collect()` operator in order to pull the reports for all the samples into a single call to `MULTIQC`.
And we also need to give it the `report_id` parameter.

This gives us the following:

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
    read_ch = channel.fromPath(params.input_csv)
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

## 3. Enable processing paired-end RNAseq data

Right now our workflow can only handle single-end RNAseq data.
It's increasingly common to see paired-end RNAseq data, so we want to be able to handle that.

Making the workflow completely agnostic of the data type would require using slightly more advanced Nextflow language features, so we're not going to do that here, but we can make a paired-end processing version to demonstrate what needs to be adapted.

### 3.1. Make a copy of the workflow called `rnaseq_pe.nf`

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Modify the default `input_csv` to point to the paired-end data

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

```groovy title="rnaseq_pe.nf" linenums="15"
// Primary input
params.input_csv = "data/paired-end.csv"
```

### 3.3. Update the channel factory

We need to tell the `.map()` operator to grab both FASTQ file paths now.

So `row -> file(row.fastq_path)` becomes `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Make a paired-end version of the FASTQC process

Let's make a copy of the module so we can have both version on hand.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Open up the new `fastqc_pe.nf` module file in the code editor and make the following code changes:

- Change `fastqc $reads` to `fastqc ${reads}` in the `script` block (line 17) so that the `reads` input will be unpacked, since it's now a tuple of two paths instead of a single path.
- Replace `${reads.simpleName}` with a wildcard (`*`) to avoid having to handle the output files individually.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Technically this generalizes the `FASTQC` process in a way that makes it able to handle either single-end or paired-end RNAseq data.

Finally, update the module import statement to use the paired-end version of the module.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Make a paired-end version of the TRIM_GALORE process

Make a copy of the module so we can have both version on hand.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Open up the new `trim_galore_pe.nf` module file in the code editor and make the following code changes:

- Change the input declaration from `path reads` to `tuple path(read1), path(read2)`
- Update the command in the `script` block, replacing `$reads` with `--paired ${read1} ${read2}`
- Update the output declarations to reflect the added files and different naming conventions, using wildcards to avoid having to list everything.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Finally, update the module import statement to use the paired-end version of the module.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Update the call to the MULTIQC process to expect two reports from TRIM_GALORE

The `TRIM_GALORE` process now produces an additional output channel, so we need to feed that to MultiQC.

Replace `TRIM_GALORE.out.fastqc_reports,` with `TRIM_GALORE.out.fastqc_reports_1,` plus `TRIM_GALORE.out.fastqc_reports_2,`:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

While we're on MultiQC, let's also update the `report_id` parameter default from `"all_single-end"` to `"all_paired-end"`.

```groovy title="rnaseq_pe.nf" linenums="9"
/*
 * Pipeline parameters
 */
params.hisat2_index_zip = "data/genome_index.tar.gz"
params.report_id = "all_paired-end"
```

### 3.7. Make a paired-end version of the HISAT2_ALIGN process

Make a copy of the module so we can have both version on hand.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Open up the new `hisat2_align_pe.nf` module file in the code editor and make the following code changes:

- Change the input declaration from `path reads` to `tuple path(read1), path(read2)`
- Update the command in the `script` block, replacing `-U $reads` with `-1 ${read1} -2 ${read2}`
- Replace all instances of `${reads.simpleName}` with `${read1.simpleName}` in the command in the `script` block as well as in the output declarations.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Finally, update the module import statement to use the paired-end version of the module.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Run the workflow to test that it works

```bash
nextflow run rnaseq_pe.nf
```

We don't use `-resume` since this wouldn't cache, and there's twice as much data to process than before, but it should still complete in under a minute.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

executor >  local (19)
[c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
[e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
[3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
[e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
```

And that's it! Now we have two slightly divergent versions of our workflow, one for single-end read data and one for paired-end data.
The next logical step would be to make the workflow accept either data type on the fly, which is out of scope for this course, but we may tackle that in a follow-up.

---

### Takeaway

You know how to adapt a single-sample workflow to parallelize processing of multiple samples, generate a comprehensive QC report and adapt the workflow to use paired-end read data if needed.

### What's next?

Congratulations, you've completed the Nextflow For RNAseq mini-course! Celebrate your success and take a well deserved break!

Next, we ask you to complete a very short survey about your experience with this training course, then we'll take you to a page with links to further training resources and helpful links.
