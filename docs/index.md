---
title: Nextflow Training
description: Welcome to the Nextflow community training portal!
hide:
  - toc
  - footer
---

# Nextflow Training

Welcome to the Nextflow community training portal!

We have several distinct training courses available on this website. Scroll down to find the one that's right for you!

The training courses listed below are designed to be usable as a self-service resource; you can work through them on your own at any time (see Environment Setup for practical details). However, you may get even more out of them by joining a group training event.

- Free online events are run regularly by the nf-core community, see the [nf-core events page](https://nf-co.re/events) for more.
- Seqera (the company that develops Nextflow) runs a variety of training events, see the [Seqera Events](https://seqera.io/events/) page and look for 'Seqera Sessions' and 'Nextflow Summit'.
- Our Community team also regularly teaches trainings hosted by third party organizations; announcements and signups for those are typically managed by the third-party hosts.

When you're ready to get down to work, click on the 'Open in GitHub Codespaces' button, either on this page or on the index page of the course you chose, to open a web-based training environment (requires a free GitHub account).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Training Environment Setup

!!! exercise "Environment Setup"

    !!! tip inline end ""

        :material-laptop:{.nextflow-primary} Set up your environment for the first time.

    Instructions for setting up your environment to work through training materials (all courses). Provides an orientation to GitHub Codespaces as well as alternate installation instructions for working on your own local machine.

    [Start the Environment Setup training :material-arrow-right:](envsetup/index.md){ .md-button .md-button--primary }

## Nextflow for Newcomers

These are foundational, domain-agnostic courses intended for those who are completely new to Nextflow. Each course consists of a series of training modules that are designed to help learners build up their skills progressively.

!!! exercise "Hello Nextflow"

    !!! tip inline end ""

        :material-run-fast:{.nextflow-primary} Learn to develop pipelines in Nextflow.

        :fontawesome-brands-youtube:{.youtube} Video material available.

    This is a course for newcomers who wish to learn how to develop their own pipelines. The course covers the core components of the Nextflow language in enough detail to enable developing simple but fully functional pipelines. It also covers key elements of pipeline design, development and configuration practices.

    The course is calibrated to take a full day to cover in group trainings.

    [Start the Hello Nextflow training :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--primary }

!!! exercise "Nextflow Run"

    !!! tip inline end ""

        :material-run-fast:{.nextflow-primary} Learn to run pipelines in Nextflow.

    This is an abridged version of the Hello Nextflow course intended for newcomers who wish to learn how to run their own pipelines but do not necessarily plan to develop pipelines themselves. The course covers the core components of the Nextflow language in enough detail to understand the basic structure of Nextflow pipelines and relate that to the experience of running and monitoring pipeline executions on the command-line. It also covers the basics of tool management and pipeline configuration practices.

    The course is calibrated to take a half day to cover in group trainings.

    [Start the Nextflow Run training :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--primary }

!!! exercise "Hello nf-core"

    !!! tip inline end ""

        :material-run-fast:{.nextflow-primary} Learn to develop nf-core compliant pipelines.

    This is a course for newcomers who wish to learn run and develop [nf-core](https://nf-co.re/) compliant pipelines. The course covers the structure of nf-core pipelines in enough detail to enable developing simple but fully functional pipelines that follow the nf-core template and development best practices.

    The course is calibrated to take a half day to cover in group trainings.

    [Start the Hello nf-core training :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--primary }

## Nextflow for Science

These are courses that demonstrate how to apply the concepts and components presented in 'Hello Nextflow' (see above) to specific scientific use cases. Each course consists of a series of training modules that are designed to help learners build up their skills progressively.

!!! exercise "Nextflow for Genomics"

    !!! tip inline end ""

        :material-dna:{.nextflow-primary} Learn to develop a pipeline for genomics in Nextflow.

    This is a course for researchers who wish to learn how to develop their own genomics pipelines. The course uses a variant calling use case to demonstrate how to develop a simple but functional genomics pipeline.

    [Start the Nextflow for Genomics training :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--primary }

!!! exercise "Nextflow for RNAseq"

    !!! tip inline end ""

        :material-flask:{.nextflow-primary} Learn to develop a pipeline for RNAseq data processing in Nextflow.

    This is a course for researchers who wish to learn how to develop their own RNAseq pipelines. The course uses a bulk RNAseq processing use case to demonstrate how to develop a simple but functional RNAseq pipeline.

    [Start the Nextflow for RNAseq training :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--primary }

!!! exercise "Nextflow for Imaging"

    !!! tip inline end ""

        :material-microscope:{.nextflow-primary} Learn to run pipelines for spatial transcriptomics imaging data in Nextflow.

    This is a course for researchers in imaging and spatial omics who wish to learn how to run and customize analysis pipelines. The course uses the nf-core/molkart pipeline to provide a biologically-relevant pipeline demonstrate how to run, configure, and manage inputs for Nextflow pipelines workflows.

    [Start the Nextflow for Imaging training :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--primary }

Let us know what other domains and use cases you'd like to see covered here by posting in the [Training section](https://community.seqera.io/c/training/) of the community forum.

## Advanced Nextflow Training

These are materials that cover advanced concepts and mechanisms for developing and deploying Nextflow pipelines to address real-world use cases. These materials are organized into Side Quests that cover individual topics, and Training Collections that combine multiple Side Quests in order to provide a comprehensive learning experience around a particular around a common theme or use case.

!!! exercise "Side Quests"

    !!! tip inline end ""

        :material-compass:{.nextflow-primary} Training modules for a variety of topics of interest.

    Side Quests are individual training modules intended for Nextflow developers who wish to widen their range and/or deepen their skills on particular topics. Although the modules are presented linearly, learners are welcome to pick and choose topics in any order. Any dependencies on components/skills that go beyond the scope of the 'Hello Nextflow' course are indicated in the corresponding module overview. For structured learning paths combining multiple Side Quests, see Training Collections below.

    [Start the Side Quests training :material-arrow-right:](side_quests/){ .md-button .md-button--primary }

!!! exercise "Training Collections"

    !!! tip inline end ""

        :material-compass:{.nextflow-primary} Collections of Side Questions grouped around a particular theme or use case.

    Training Collections combine multiple Side Quests in order to provide a comprehensive learning experience around a particular theme or use case.

    [Browse the Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--primary }

## Archived materials

These are the original Nextflow training materials that were developed at the start of the project. We are in the process of deprecating them in favor of the newer materials listed above. However, some topics covered in the original trainings are not yet represented in the newer material, so we are keeping these around for reference, with the caveat that they are no longer maintained and some exercises may no longer work.

!!! exercise "Fundamentals Training"

    !!! quote inline end ""

        :octicons-mortar-board-16:{.nextflow-primary} Comprehensive training material for exploring the full scope of Nextflow's capabilities.

    The fundamentals training material covers all things Nextflow. Intended as a reference material for anyone looking to build complex workflows with Nextflow.

    [Start the Fundamentals Training :material-arrow-right:](archive/basic_training){ .md-button .md-button--primary }

!!! exercise "Advanced Training"

    !!! quote inline end ""

        :fontawesome-solid-hat-wizard:{.nextflow-primary} Advanced training material for mastering Nextflow.

    Advanced material exploring the more advanced features of the Nextflow language and runtime, and how to use them to write efficient and scalable data-intensive workflows.

    [Start the Advanced Training :material-arrow-right:](archive/advanced){ .md-button .md-button--primary }

## Other/Experimental

These are other training courses that are not being actively taught/maintained and that we may repurpose elsewhere or delete in the near future.
The corresponding materials are not available within the training environment.
You can still find the materials in the GitHub repository and download them for local use.

- **nf-customize** — Configuring nf-core pipelines ([docs](other/nf_customize) / [code](https://github.com/nextflow-io/training/tree/master/other/nf-customize))

- **troubleshoot** — Troubleshooting exercises ([docs](other/troubleshoot) / [code](https://github.com/nextflow-io/training/tree/master/other/troubleshoot))

- **hands-on (rnaseq)** — Developing a pipeline for bulk RNAseq (deprecated) ([docs](other/hands_on) / [code](https://github.com/nextflow-io/training/tree/master/other/hands-on))

## Resources

Quick reference to some handy links:

| Reference                                                   |  Community                                                   |
| ----------------------------------------------------------- | ------------------------------------------------------------ |
| [Nextflow Docs](https://nextflow.io/docs/latest/index.html) | [Nextflow Slack](https://www.nextflow.io/slack-invite.html)  |
| [Nextflow Homepage](https://nextflow.io/)                   | [nf-core](https://nf-co.re/)                                 |
| [Seqera](https://seqera.io/)                                | [Seqera Community Forum](https://community.seqera.io)        |

Not sure where to go? Check out the [Getting help](help.md) page.

## Credits and contributions

[![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0](assets/img/cc_by-nc-nd.svg){ align=right }](https://creativecommons.org/licenses/by-nc-nd/4.0/)

This training material is developed and maintained by [Seqera](https://seqera.io) and released under an open-source license ([CC BY-NC-ND](https://creativecommons.org/licenses/by-nc-nd/4.0/)) for the benefit of the community. You are welcome to reuse these materials according to the terms of the license. If you are an instructor running your own trainings, we'd love to hear about how it goes and what we could do to make it easier. If you wish to use this material in a way that falls outside the scope of the license, please contact us at community@seqera.io to discuss your request.

We welcome fixes and improvements from the community. Every page has a :material-file-edit-outline: icon in the top right of the page, which will take you to GitHub where you can propose changes to the training source material via a pull request.

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
