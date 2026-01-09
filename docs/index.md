---
title: Nextflow Training
description: Welcome to the Nextflow community training portal!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Self-service courses__

    ---

    **Welcome to the Nextflow community training portal!**

    The training courses listed below are designed to be usable as a self-service resource.
    You can work through them on your own at any time either in the web-based environment we provide via Github Codespaces or in your own environment.

    [Explore the courses :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Additional information__

    ---

    ??? warning "Version compatibility"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **As of January 2026, all of our Nextflow training courses require Nextflow version 25.10.2 or later, with strict syntax activated, unless otherwise noted.**

        For more information about version requirements and strict syntax, please see the [Nextflow docs migration guide](https://nextflow.io/docs/latest/strict-syntax.html).

        Older versions of the training material corresponding to prior syntax are available via the version selector in the menu bar of this webpage.

    ??? terminal "Environment options"

        We provide a web-based training environment where everything you need to take the training is preinstalled, available through Github Codespaces (requires a free GitHub account).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        If this does not suit your needs, please see the other [Environment options](./envsetup/index.md).

    ??? learning "Training events"

        If you'd prefer to take Nextflow training as part of a structured event, there are many opportunities to do so. We recommend checking out the following options:

        - **[Training Weeks]()** organized quarterly by the Community team
        - **[Seqera Events](https://seqera.io/events/)** include in-person training events organized by Seqera (search for 'Seqera Sessions' and 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organize events for their local community
        - **[nf-core events](https://nf-co.re/events)** include community hackathons

    ??? people "Information for trainers"

        If you are an instructor running your own trainings, you are welcome to use our materials directly from the training portal as long as you attribute proper credit. See 'Credits and contributions' below for details.

        In addition, we'd love to hear from you on how we could better support your training efforts! Please contact us at [community@seqera.io](mailto:community@seqera.io) or on the community forum (see [Help](help.md) page).

    ??? licensing "Open-source license and contribution policy"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0](assets/img/cc_by-nc-nd.svg){ align=right }](https://creativecommons.org/licenses/by-nc-nd/4.0/)

        This training material is developed and maintained by [Seqera](https://seqera.io) and released under an open-source license ([CC BY-NC-ND](https://creativecommons.org/licenses/by-nc-nd/4.0/)) for the benefit of the community. If you wish to use this material in a way that falls outside the scope of the license (note the limitations on commercial use and redistribution), please contact us at [community@seqera.io](mailto:community@seqera.io) to discuss your request.

        We welcome improvements, fixes and bug reports from the community. Every page has a :material-file-edit-outline: icon in the top right of the page linking to the code repository, where you can report issues or propose changes to the training source material via a pull request. See the `README.md` in the repository for more details.

</div>

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Introductory track__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow for Newcomers {.mt-1}

    Domain-agnostic courses intended for those who are completely new to Nextflow. Each course consists of a series of training modules that are designed to help learners build up their skills progressively.

    ??? courses "**Hello Nextflow:** Learn to develop your own pipelines"

        This course covers the core components of the Nextflow language in enough detail to enable developing simple but fully functional pipelines, plus key elements of pipeline design, development and configuration practices.

        [Start the Hello Nextflow training :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Learn to run existing pipelines"

        This is an abridged version of the Hello Nextflow course, covering the core components of the Nextflow language in enough detail to understand the basic structure of Nextflow pipelines and relate that to the experience of running and monitoring pipeline executions on the command-line, plus the basics of tool management and pipeline configuration practices.

        [Start the Nextflow Run training :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow for Science {.mt-1}

    Learn to apply the concepts and components presented in 'Hello Nextflow' to specific scientific use cases.

    ??? courses "**Nextflow for Genomics** (variant calling)"

        For researchers who wish to learn how to develop their own genomics pipelines. The course uses a variant calling use case to demonstrate how to develop a simple but functional genomics pipeline.

        [Start the Nextflow for Genomics training :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        For researchers who wish to learn how to develop their own RNAseq pipelines. The course uses a bulk RNAseq processing use case to demonstrate how to develop a simple but functional RNAseq pipeline.

        [Start the Nextflow for RNAseq training :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (spatial omics)"

        For researchers in imaging and spatial omics who wish to learn how to run and customize analysis pipelines. The course uses the nf-core/molkart pipeline to provide a biologically-relevant pipeline demonstrate how to run, configure, and manage inputs for Nextflow pipelines workflows.

        [Start the Nextflow for Imaging training :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Advanced track__

    ---

    ### :material-bridge:{.nextflow-primary} From Nextflow to nf-core {.mt-1}

    Learn to utilize code and best practices from the [nf-core](https://nf-co.re/) community project.

    These courses help you go from Nextflow fundamentals to nf-core best practices.
    Understand how and why the nf-core community builds pipelines, and how you can contribute and reuse these techniques.

    ??? courses "**Hello nf-core:** Get started with nf-core"

        For developers who wish to learn run and develop [nf-core](https://nf-co.re/) compliant pipelines. The course covers the structure of nf-core pipelines in enough detail to enable developing simple but fully functional pipelines that follow the nf-core template and development best practices, as well as use existing nf-core modules.

        [Start the Hello nf-core training :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Advanced Nextflow Training {.mt-1}

    Learn advanced concepts and mechanisms for developing and deploying Nextflow pipelines to address real-world use cases.

    ??? courses "**Side Quests:** Deep dives into standalone topics"

        Standalone mini-courses intended for Nextflow developers who wish to widen their range and/or deepen their skills on particular topics. They are presented linearly but can be taken in any order (see dependencies in each mini-course overview).

        [Browse the Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Recommended learning paths through the Side Quests"

        Training Collections combine multiple Side Quests in order to provide a comprehensive learning experience around a particular theme or use case.

        [Browse the Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

## Archived materials

### Original Nextflow training materials

These were developed at the start of the project. We are in the process of deprecating them in favor of the newer materials listed above. However, some topics covered in the original trainings are not yet represented in the newer material, so we are keeping these around for reference, with the caveat that they are no longer maintained and some exercises may no longer work.

- **Fundamentals Training** — Intended as a reference material for anyone looking to build complex workflows with Nextflow ([docs](archive/basic_training/) / [code](https://github.com/nextflow-io/training/tree/master/nf-training))

- **Advanced Training** — More advanced features of the Nextflow language and runtime, and how to use them to write efficient and scalable data-intensive workflows ([docs](archive/advanced/) / [code](https://github.com/nextflow-io/training/tree/master/nf-training-advanced))

### Other/Experimental

These are other training courses that are not being actively taught/maintained and that we may repurpose elsewhere or delete in the near future.
The corresponding materials are not available within the training environment.
You can still find the materials in the GitHub repository and download them for local use.

- **nf-customize** — Configuring nf-core pipelines ([docs](other/nf_customize) / [code](https://github.com/nextflow-io/training/tree/master/other/nf-customize))

- **troubleshoot** — Troubleshooting exercises ([docs](other/troubleshoot) / [code](https://github.com/nextflow-io/training/tree/master/other/troubleshoot))

- **hands-on (rnaseq)** — Developing a pipeline for bulk RNAseq (deprecated) ([docs](other/hands_on) / [code](https://github.com/nextflow-io/training/tree/master/other/hands-on))

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
