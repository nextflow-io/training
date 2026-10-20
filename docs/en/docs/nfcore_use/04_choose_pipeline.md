# Part 4: Choose a pipeline for a scientific question

In Parts 1 to 3, you ran pipelines that were chosen for you.
In practice, the first challenge is deciding which pipeline to use, or whether any existing pipeline fits your question at all.

The nf-core catalog contains over a hundred pipelines, several of which can look relevant to the same question.
A pipeline that processes the right data type does not necessarily produce the result you need.

In this part, you take a research scenario, find candidate pipelines, and evaluate them using their documentation and the results of a completed run.
Because you examine runs that have already completed, you can look at real outputs and execution reports without waiting for a pipeline to run.

This part works for individual learners and for groups in a workshop.
In a group, discuss each step together and present your decision to the room at the end.
Working alone, write down your answers as you go and compare them with the worked answer for your scenario.

---

## 1. Pick a scenario and define the goal

### 1.1. Pick a scenario

Choose the scenario closest to your own research interests.

#### 1.1.1. Structural variants in an unsolved rare disease trio

A clinical research lab is studying a child with a suspected genetic disorder.
Exome sequencing of the child and both unaffected parents did not find a causal variant, so the family was re-sequenced with long-read whole-genome sequencing, which detects structural variants that short reads miss.
The lab wants a short list of candidate structural variants, prioritized by their relevance to the child's clinical features, which have been recorded as HPO terms.

- **Data**: Oxford Nanopore whole-genome FASTQ files for three individuals, and HPO terms for the child
- **Reference**: GRCh38

#### 1.1.2. Differential expression in bulk RNA-seq

A neurogenetics group has bulk RNA-seq from iPSC-derived neurons carrying a risk variant and isogenic controls, with four replicates each.
They want to know which genes are differentially expressed between the two groups.

- **Data**: paired-end Illumina FASTQ files for eight samples
- **Reference**: GRCh38

#### 1.1.3. Single-cell RNA-seq of immune cells

An immunology group has 10x Genomics single-cell RNA-seq of peripheral blood mononuclear cells from patients and healthy controls.
They want to identify cell types and compare their proportions between patients and controls.

- **Data**: 10x Chromium 3' FASTQ files for six samples
- **Reference**: GRCh38

#### 1.1.4. DNA methylation at imprinted loci

An epigenetics group has whole-genome bisulfite sequencing from patients with a suspected imprinting disorder and matched controls.
They want methylation levels at known imprinted loci for each sample.

- **Data**: paired-end Illumina bisulfite FASTQ files for six samples
- **Reference**: GRCh38

### 1.2. Define the goal

Before looking at any pipeline, write down in one sentence the result the analysis must produce.
For example: "a list of candidate causal variants, ranked by likely pathogenicity".

This sentence is what you evaluate each pipeline against.
Without it, it is easy to settle on a pipeline that accepts your data but stops short of the result you need.

### Takeaway

You have a scenario and a one-sentence definition of the result you need.

### What's next?

Find pipelines that could produce that result.

---

## 2. Find candidate pipelines

### 2.1. Search the catalog

Open the [nf-core pipelines catalog](https://nf-co.re/pipelines) and search using keywords from your scenario: the data type (`long-read`, `bisulfite`, `single-cell`), the assay or platform (`nanopore`, `10x`), or the analysis (`structural variants`, `differential expression`).
Try several keywords, since a pipeline's name does not always reflect everything it does.

Shortlist up to two candidates.
If you haven't found a candidate after about five minutes, open the hint for your scenario.

??? tip "Hint: Structural variants in an unsolved rare disease trio"

    Compare [nf-core/raredisease](https://nf-co.re/raredisease) and [nf-core/longraredisease](https://nf-co.re/longraredisease).
    Both address rare disease; which one fits your data?

??? tip "Hint: Differential expression in bulk RNA-seq"

    Look at [nf-core/rnaseq](https://nf-co.re/rnaseq) and [nf-core/differentialabundance](https://nf-co.re/differentialabundance).

??? tip "Hint: Single-cell RNA-seq of immune cells"

    Look at [nf-core/scrnaseq](https://nf-co.re/scrnaseq).

??? tip "Hint: DNA methylation at imprinted loci"

    Look at [nf-core/methylseq](https://nf-co.re/methylseq).

### 2.2. Check the pipeline status

For each candidate, open its page and note:

- **What it does**: the introduction summarizes the pipeline's purpose and its main steps, often with a diagram.
- **Release status**: a stable release means the pipeline is considered ready for production use.
  Pipelines still under development are listed in the catalog without a stable release; they may change significantly from one version to the next, so results may not be reproducible across updates.
- **Last update**: a pipeline that is actively maintained is more likely to support current tool versions and reference data.

### Takeaway

You have a shortlist of up to two candidate pipelines.

### What's next?

Evaluate whether each candidate produces the result you defined in 1.2.

---

## 3. Evaluate the candidates

For each candidate, work through the questions below, using the pipeline's **Usage** and **Output** documentation pages and the results of a completed run.

### 3.1. Open the results of a completed run

For every release, the nf-core project runs each pipeline on a full-size dataset in the cloud and publishes the results.
You can browse them in the **Results** tab of the pipeline's page on the nf-core website.

The results are organized in the same way as when you run the pipeline yourself:

- The output directories contain the files described in the output documentation, including the MultiQC report in `multiqc/`.
- The `pipeline_info/` directory contains the execution report and trace, as in [Part 3](./03_run_production_pipeline.md).

Some pipelines run several test configurations and organize their results in one subdirectory per configuration, for example per aligner in `nf-core/scrnaseq`.
For some releases, only the `pipeline_info/` directory is published; in that case, select an earlier release from the version menu, or rely on the output documentation.

!!! tip "Viewing runs on Seqera Platform"

    The Community Showcase workspace on [Seqera Platform](https://cloud.seqera.io) contains completed runs of several nf-core pipelines, including `nf-core/rnaseq` and `nf-core/scrnaseq`.
    A run page shows the parameters, the status and resource usage of every task, and the reports the pipeline produced.
    Viewing it requires a free Seqera account; see [Scale with Seqera](../seqera_scale/index.md).

### 3.2. Inputs

Find the samplesheet format in the usage documentation, as you did for `nf-core/demo` in [Part 2](./02_configure_execution.md#14-input-validation).

- What columns does the samplesheet require?
- Does your scenario provide that information, or would you need to collect more?

### 3.3. Reference data

- Does the pipeline support GRCh38?
- Does it need additional reference files beyond the genome, such as annotation databases?

### 3.4. Outputs

Read the output documentation, then look at the output files and the MultiQC report in the published results.

- Does the pipeline produce the result you defined in 1.2, or an intermediate step toward it?
- If it stops short, what would you need to do next, and is there another pipeline for that step?

### 3.5. Compute

Open the execution report in the `pipeline_info/` directory of the published results.

- Which processes took the longest, and which used the most memory?
- How might that scale to the size of your dataset?

### Takeaway

You have evaluated each candidate on its inputs, reference requirements, outputs and compute needs.

### What's next?

Make a decision.

---

## 4. Decide

### 4.1. Fill in the evaluation card

Decide which pipeline to use, or conclude that none of the candidates is sufficient on its own, then fill in the evaluation card:

```markdown title="Evaluation card"
Scenario:
Pipeline chosen:
Why it fits:
What it does not cover (and what you would do about it):
Inputs you would need to prepare:
Compute concern:
One open question:
```

If you're working in a group, choose one person to present the card in one minute.

### 4.2. Compare with the worked answer

Compare your conclusions with the worked answer for your scenario.
There is room for judgment, so a different answer with good reasons is not necessarily wrong.

??? solution "Structural variants in an unsolved rare disease trio"

    **Pipeline chosen**: `nf-core/longraredisease`.

    **Why not `nf-core/raredisease`**: it is an obvious match when you search the catalog for "rare disease", and it addresses the same clinical question, but it only accepts Illumina short-read data.
    Checking the data type first rules it out.

    **Why it fits**: `nf-core/longraredisease` is designed for long-read data from Oxford Nanopore or PacBio, and focuses on structural variants.
    It calls them with three callers (Sniffles, CuteSV and SVIM) and merges the results into a consensus set, supports trio analysis including detection of de novo variants, annotates the variants with AnnotSV, and ranks them against the patient's HPO terms with SvAnna.
    It can also call SNVs, copy number variants, short tandem repeats and methylation, which long reads provide from the same data.

    **Inputs to prepare**: the samplesheet requires `sample` and `file_path`, and uses `hpo_terms`, `sex`, `phenotype`, `family_id`, `maternal_id` and `paternal_id` to describe the trio and the child's clinical features.

    **What to keep in mind**: check the release history; `nf-core/longraredisease` is much younger than `nf-core/raredisease`, so it has less of a track record.
    The published test results cover one chromosome of a single sample, and do not include the phenotype-based prioritization step; use the output documentation to see what those results look like.

    **Compute concern**: the published execution report comes from that small test, so it greatly understates the resources a trio of whole genomes requires.
    Plan for read alignment and variant calling on whole genomes to dominate the run.

??? solution "Differential expression in bulk RNA-seq"

    **Pipelines chosen**: `nf-core/rnaseq` followed by `nf-core/differentialabundance`.

    **Why**: neither pipeline answers the question alone.
    `nf-core/rnaseq` performs QC, alignment and quantification, producing gene count matrices, but its documentation states that it does not compare samples statistically.
    The DESeq2 step it runs is a QC step (sample PCA and clustering), not a differential expression test.
    `nf-core/differentialabundance` takes the count matrix from `nf-core/rnaseq`, a sample sheet describing the groups, and a contrasts file defining the comparison, and runs DESeq2 (or limma) to produce differential expression tables, plots and a report.

    **Inputs to prepare**: a samplesheet for `nf-core/rnaseq` (as in Part 5), then a sample sheet with a `condition` column and a contrasts file for `nf-core/differentialabundance`.

    **Reference data**: the `nf-core/rnaseq` documentation recommends providing a recent GENCODE or Ensembl FASTA and GTF rather than the built-in `--genome` references, whose annotations are out of date.

    **Compute concern**: alignment and quantification dominate the run time and memory; differential analysis on a count matrix is light.

??? solution "Single-cell RNA-seq of immune cells"

    **Pipeline chosen**: `nf-core/scrnaseq`, with a caveat.

    **Why it fits partially**: it takes 10x Chromium data to count matrices, with a choice of tools (Alevin-fry by default, STARsolo, kallisto | bustools or Cell Ranger), and outputs them in formats ready for analysis in Python (`.h5ad`) or R (`.rds` for Seurat).

    **What it does not cover**: clustering, cell type annotation and comparing proportions between groups.
    The nf-core community is developing `nf-core/scdownstream` for these steps, but at the time of writing it has no stable release.
    Until then, the group would do the downstream analysis interactively, for example with Scanpy or Seurat.

    **Inputs to prepare**: a samplesheet with `sample`, `fastq_1` and `fastq_2`, and optionally `expected_cells`; you also need to know the 10x chemistry version for the `--protocol` parameter.

    **Compute concern**: index building and alignment are memory-intensive for the human genome.

??? solution "DNA methylation at imprinted loci"

    **Pipeline chosen**: `nf-core/methylseq`.

    **Why it fits**: it aligns bisulfite sequencing reads (with Bismark by default, or bwa-meth) and produces per-CpG methylation calls for each sample.
    It can also restrict the methylation output to regions of interest: `--run_targeted_sequencing` with `--target_regions_file` filters the results to the regions in a BED file, which could contain the imprinted loci.

    **What it does not cover**: statistical comparison between patients and controls.
    The scenario only asks for methylation levels per sample, so this is sufficient here, but a differential methylation analysis would need a separate tool such as methylKit or DSS.

    **Inputs to prepare**: a simple samplesheet with `sample`, `fastq_1` and `fastq_2`, and a BED file of the imprinted loci.

    **Compute concern**: aligning whole-genome bisulfite data is slow and memory-intensive; the alignment step dominates the run.

### Takeaway

Choosing a pipeline means checking that it produces the result you need from the inputs you have, not only that it accepts your data type.
Completed runs let you check outputs and compute needs before running anything yourself.

### What's next?

The next parts cover practical challenges of running a pipeline on your own data: preparing the samplesheet, troubleshooting validation errors, and tuning compute resources.
They use `nf-core/rnaseq` as the example; at the end of each part, you can apply the same steps to the pipeline you chose here.

---

## Summary

In this part you learned to:

- Define the result an analysis must produce before choosing a pipeline
- Find candidate pipelines in the nf-core catalog and check their status
- Evaluate a pipeline's inputs, reference requirements, outputs and compute needs using its documentation and the results of a completed run
