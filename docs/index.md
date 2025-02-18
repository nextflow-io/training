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

        :material-lightbulb: Set up your environment for the first time.

    Instructions for setting up your environment to work through training materials (all courses). Provides an orientation to GitHub Codespaces as well as alternate installation instructions for working on your own local machine.

    [Launch the Environment Setup training :material-arrow-right:](envsetup/index.md){ .md-button .md-button--primary }

## Nextflow for Newcomers

These are foundational, domain-agnostic courses intended for those who are completely new to Nextflow. Each course consists of a series of training modules that are designed to help learners build up their skills progressively.

!!! exercise "Hello Nextflow"

    !!! tip inline end ""

        :material-run-fast: Learn to develop pipelines in Nextflow.

    This is a course for newcomers who wish to learn how to develop their own pipelines. The course covers the core components of the Nextflow language in enough detail to enable developing simple but fully functional pipelines. It also covers key elements of pipeline design, development and configuration practices.

    [Launch the Hello Nextflow training :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--primary }

!!! info ""

    **Coming soon:** "Nextflow Run" — Learn to run Nextflow pipelines (run only, no code development)

<!-- COMMENTED OUT UNTIL THIS IS READY
!!! exercise "Nextflow Run"

    !!! tip inline end ""

        :material-run-fast: Learn to run Nextflow pipelines.

    This is a course for newcomers who wish to learn how to run existing pipelines. The course covers the bare essentials of the Nextflow language in order to enable interpretation of existing pipelines, as well as the mechanics for configuring and running Nextflow pipelines from a command-line environment. It also covers important components of the Nextflow ecosystem, including the nf-core project, which offers a large number of community-curated pipelines, and the Seqera platform for managing pipeline execution at scale (operated by the creators of Nextflow).

    [Launch the Nextflow Run training :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--primary }
-->

## Nextflow for Science

These are courses that demonstrate how to apply the concepts and components presented in 'Hello Nextflow' (see above) to specific scientific use cases. Each course consists of a series of training modules that are designed to help learners build up their skills progressively.

!!! exercise "Nextflow for Genomics"

    !!! tip inline end ""

        :material-run-fast: Learn to develop a pipeline for genomics in Nextflow.

    This is a course for researchers who wish to learn how to develop their own genomics pipelines. The course uses a variant calling use case to demonstrate how to develop a simple but functional genomics pipeline.

    [Launch the Nextflow for Genomics training :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--primary }

!!! info ""

    **Coming soon:** "Nextflow for RNAseq" — Learn to develop a pipeline for bulk RNAseq analysis in Nextflow

<!-- COMMENTED OUT UNTIL THIS IS READY
!!! exercise "Nextflow for RNAseq"

    !!! tip inline end ""

        :material-run-fast: Learn to develop a pipeline for bulk RNAseq analysis in Nextflow.

    This is a course for researchers who wish to learn how to develop their own RNAseq pipelines. The course uses a bulk RNAseq processing use case to demonstrate how to develop a simple but functional RNAseq pipeline.

    [Launch the Nextflow for RNAseq training :material-arrow-right:]((nf4_science/rnaseq/){ .md-button .md-button--primary }
-->

## In-depth Nextflow Training

These are courses that demonstrate how to use Nextflow features in more detail or at a more advanced level. Each course consists of one or more training modules that are designed to help learners hone their skills on the corresponding topics.

<!-- COMMENTED OUT UNTIL THE FIRST ONE IS READY
!!! exercise "Side Quests"

    !!! tip inline end ""

        :material-run-fast: Training modules for a variety of topics of interest.

    This is a course for Nextflow developers who wish to widen their range and/or deepen their skills. Although the modules are presented linearly, learners are welcome to pick and choose topics in any order. Any dependencies on components/skills that go beyond the scope of the 'Hello Nextflow' course are indicated in the corresponding module overview.

    [Launch the Side Quests training :material-arrow-right:](side_quests/index.md){ .md-button .md-button--primary }
-->

!!! exercise "Fundamentals Training"

    !!! quote inline end ""

        :material-lightbulb: Comprehensive training material for exploring the full scope of Nextflow's capabilities.

    The fundamentals training material covers all things Nextflow. Intended as a reference material for anyone looking to build complex workflows with Nextflow.

    [Launch the Fundamentals Training :material-arrow-right:](basic_training/index.md){ .md-button .md-button--primary }

!!! exercise "Advanced Training"

    !!! quote inline end ""

        :material-lightbulb: Advanced training material for mastering Nextflow.

    Advanced material exploring the more advanced features of the Nextflow language and runtime, and how to use them to write efficient and scalable data-intensive workflows.

    [Launch the Advanced Training :material-arrow-right:](advanced/index.md){ .md-button .md-button--primary }

## Other/Experimental

These are training courses that are not being actively taught/maintained and that we may repurpose elsewhere or delete in the near future.
The corresponding materials are not available within the training environment.
You can still find the materials in the GitHub repository and download them for local use.

- **nf-customize** — Configuring nf-core pipelines ([docs](other/nf_customize) / [code](../other/nf-customize))

- **nf-develop** — Developing a pipeline with the nf-core template ([docs](other/nf_develop) / [code](../other/nf-develop))

- **troubleshoot** — Troubleshooting exercises ([docs](other/troubleshoot) / [code](../other/troubleshoot))

- **hands-on (rnaseq)** — Developing a pipeline for bulk RNAseq (deprecated) ([docs](other/hands_on) / [code](../other/hands-on))

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

This training material is developed and maintained by [Seqera](https://seqera.io) and released under an open-source license ([CC BY-NC-ND](https://creativecommons.org/licenses/by-nc-nd/4.0/)) for the benefit of the community. You are welcome to reuse these materials according to the terms of the license. If you are an instructor running your own trainings, we'd love to hear about how it goes and what we could do to make it easier.

We welcome fixes and improvements from the community. Every page has a :material-file-edit-outline: icon in the top right of the page, which will take you to GitHub where you can propose changes to the training source material via a pull request.

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
