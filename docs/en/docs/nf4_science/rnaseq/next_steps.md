# Course summary

Congratulations on completing the Nextflow for RNAseq training course!

## Your journey

You started by running RNAseq processing tools manually in the terminal to understand the methodology.
Then you built a single-sample Nextflow pipeline to automate the process, scaled it to handle multiple samples in parallel, and extended it to handle paired-end data and aggregate QC reports across samples.

### What you built

- An RNAseq processing pipeline that takes FASTQ files as input and produces trimmed reads, alignments and aggregated QC reports as output.
- Processes for trimming (Trim Galore), alignment (HISAT2), quality control (FastQC) and report aggregation (MultiQC) stored in separate module files.
- The pipeline automatically parallelizes the processing of input samples using Nextflow's dataflow paradigm.
- The final pipeline handles paired-end sequencing data.

### Skills acquired

Through this hands-on course, you've learned how to:

- Write a linear workflow to apply basic RNAseq processing and QC methods
- Handle domain-specific files such as FASTQ and reference genome resources appropriately
- Handle single-end and paired-end sequencing data
- Leverage Nextflow's dataflow paradigm to parallelize per-sample RNAseq processing
- Aggregate QC reports across multiple steps and samples using relevant channel operators

You're now equipped to start applying Nextflow to RNAseq analysis workflows in your own work.

## Next steps to build your skills

Here are our top suggestions for what to do next:

- Apply Nextflow to other scientific analysis use cases with [Nextflow for Science](../index.md)
- Get started with nf-core with [Hello nf-core](../../hello_nf-core/index.md)
- Explore more advanced Nextflow features with the [Side Quests](../../side_quests/index.md)

Finally, we recommend you have a look at [**Seqera Platform**](https://seqera.io/), a cloud-based platform developed by the creators of Nextflow that makes it even easier to launch and manage your workflows, as well as manage your data and run analyses interactively in any environment.

## Getting help

For help resources and community support, see the [Help page](../../help.md).

## Feedback survey

Before you move on, please take a minute to complete the course survey! Your feedback helps us improve our training materials for everyone.

[Take the survey :material-arrow-right:](survey.md){ .md-button .md-button--primary }
