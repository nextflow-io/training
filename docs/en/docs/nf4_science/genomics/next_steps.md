# Course summary

Congratulations on completing the Nextflow for Genomics training course! 🎉

## Your journey

You started by running variant calling tools manually in the terminal to understand the methodology.
Then you built a single-sample Nextflow pipeline to automate the process, scaled it to handle multiple samples in parallel, and added multi-sample joint genotyping using channel operators.

### What you built

- A variant calling pipeline that takes BAM files as input and produces joint-called VCFs as output.
- Three processes (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER`, and `GATK_JOINTGENOTYPING`) stored in separate module files.
- The pipeline scales automatically to any number of input samples using Nextflow's dataflow paradigm.
- The results are published to a directory called `results/`.

### Skills acquired

Through this hands-on course, you've learned how to:

- Write a linear workflow to apply variant calling to a single sample
- Handle accessory files such as index files and reference genome resources appropriately
- Leverage Nextflow's dataflow paradigm to parallelize per-sample variant calling
- Implement multi-sample joint calling using relevant channel operators
  You're now equipped to start applying Nextflow to genomics analysis workflows in your own work.

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
