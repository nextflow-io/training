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

#### 1.1.1. Germline variants in a rare disease trio

A clinical research lab has whole-exome sequencing data for a family trio: an affected child and both unaffected parents.
They suspect a recessive disorder and want a short list of candidate causal variants for review.

- **Data**: paired-end Illumina exome FASTQ files for three individuals
- **Reference**: GRCh38
- **Starting points**: [nf-core/raredisease](https://nf-co.re/raredisease), [nf-core/sarek](https://nf-co.re/sarek)

#### 1.1.2. Differential expression in bulk RNA-seq

A neurogenetics group has bulk RNA-seq from iPSC-derived neurons carrying a risk variant and isogenic controls, with four replicates each.
They want to know which genes are differentially expressed between the two groups.

- **Data**: paired-end Illumina FASTQ files for eight samples
- **Reference**: GRCh38
- **Starting points**: [nf-core/rnaseq](https://nf-co.re/rnaseq), [nf-core/differentialabundance](https://nf-co.re/differentialabundance)

#### 1.1.3. Single-cell RNA-seq of immune cells

An immunology group has 10x Genomics single-cell RNA-seq of peripheral blood mononuclear cells from patients and healthy controls.
They want to identify cell types and compare their proportions between patients and controls.

- **Data**: 10x Chromium 3' FASTQ files for six samples
- **Reference**: GRCh38
- **Starting point**: [nf-core/scrnaseq](https://nf-co.re/scrnaseq)

#### 1.1.4. DNA methylation at imprinted loci

An epigenetics group has whole-genome bisulfite sequencing from patients with a suspected imprinting disorder and matched controls.
They want methylation levels at known imprinted loci for each sample.

- **Data**: paired-end Illumina bisulfite FASTQ files for six samples
- **Reference**: GRCh38
- **Starting point**: [nf-core/methylseq](https://nf-co.re/methylseq)

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

Open the [nf-core pipelines catalog](https://nf-co.re/pipelines) and search using keywords from your scenario, such as the data type (`exome`, `bisulfite`, `single-cell`) or the analysis (`variant calling`, `differential expression`).

Your scenario lists one or two starting points, but don't stop there.
Look for other pipelines that could be relevant, and shortlist up to two candidates.

### 2.2. Check the pipeline status

For each candidate, open its page and note:

- **What it does**: the introduction summarizes the pipeline's purpose and its main steps, often with a diagram.
- **Release status**: a pipeline with a stable release has been reviewed by the nf-core community.
  Pipelines still under development are listed in the catalog but have no release yet.
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

??? solution "Germline variants in a rare disease trio"

    **Pipeline chosen**: `nf-core/raredisease`.

    **Why it fits**: it is built for this question.
    It calls SNVs, indels, structural variants and other variant types from exome or genome data, annotates them with VEP, and ranks them with GENMOD, taking the family structure and inheritance models into account (including compound heterozygotes).
    Its outputs include ranked VCFs per family, which is the short list the lab needs.

    **Why not `nf-core/sarek`**: Sarek calls and annotates germline variants very well, and can genotype samples jointly, but it does not use pedigree information or rank variants.
    You would have to do the prioritization yourself.

    **Inputs to prepare**: the samplesheet requires `sex`, `phenotype` and `case_id` columns, and uses `paternal_id` and `maternal_id` to describe the trio.
    For exome data, you set `--analysis_type wes` and provide a BED file of the capture targets.

    **Compute concern**: beyond the genome, the pipeline needs several reference resources (VEP cache, variant catalogs, scoring configurations) that you need to download or build before the first run.

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
