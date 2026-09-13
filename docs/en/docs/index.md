---
title: Home
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

    Work through the courses below at your own pace, in our web-based environment or your own.
    Each course is hands-on, with goal-oriented exercises you can complete independently.

    [Explore the courses :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Training Events__

    ---

    **Looking for something beyond self-service?**

    Find structured training events, guidance for running your own trainings, and our open-source license and contribution policy.

    [See training events :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __For users__

    ---

    ### :material-play-circle:{.nextflow-primary} Run pipelines {.mt-1}

    Learn to run existing pipelines without writing any code.

    ??? courses "**Nextflow Run:** Run pipelines with Nextflow"

        A fast-track introduction to running Nextflow pipelines that does not require understanding code. Covers launching pipelines, retrieving outputs, using containers, and configuring execution at a basic level.

        [View the training :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Find and run community-curated pipelines"

        A fast-track introduction to finding, running, and configuring pipelines from the nf-core community project, starting with a minimal demo pipeline then scaling up to a production-scale analysis pipeline.

        [View the training :material-arrow-right:](nfcore_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Launch and monitor pipelines at scale"

        A hands-on introduction to launching and monitoring Nextflow pipelines with Seqera Platform, from both the web interface and the command line.

        [View the training :material-arrow-right:](seqera_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Manage execution {.mt-1}

    Learn to manage pipeline execution effectively.

    ??? courses "**Execution Config:** Configure pipelines like a pro"

        A hands-on introduction to configuring Nextflow pipeline execution: adapting to different compute environments, controlling resource allocations and retries, and switching between preset configuration profiles.

        [View the training :material-arrow-right:](nextflow_config/index.md){ .md-button .md-button--secondary }

    !!! info compact "More topics coming"

        Performance tuning, HPC/cloud execution, and more are planned for this section.
        Vote on what to cover next in our short interest poll.
        <!-- TODO: link to interest poll once available -->

-   :material-code-tags:{ .lg .middle } __For developers__

    ---

    ### :material-wrench:{.nextflow-primary} Write pipelines {.mt-1}

    Learn to develop your own Nextflow pipelines.

    ??? courses "**Hello Nextflow:** Develop your own pipelines from scratch"

        This course covers the core components of the Nextflow language in enough detail to enable developing simple but fully functional pipelines, plus key elements of pipeline design, development and configuration practices.

        [View the training :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** Use the nf-core tools and rules"

        For Nextflow developers who wish to learn to develop [nf-core](https://nf-co.re/) compliant pipelines.
        The course covers the structure of nf-core pipelines in enough detail to enable developing simple but fully functional pipelines that leverage the nf-core template and development best practices, as well as use existing nf-core modules.

        [View the training :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Dive into advanced Nextflow topics"

        A collection of standalone mini-courses intended for Nextflow developers who wish to widen their range and/or deepen their skills on particular topics.
        They are presented linearly but can be taken in any order (see dependencies in each mini-course overview).

        [Browse the Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow for Science {.mt-1}

    Learn to develop Nextflow pipelines for specific scientific applications.

    ??? courses "**Genomics:** Develop a variant calling pipeline"

        A course for researchers who wish to learn how to develop their own genomics pipelines, using a variant calling use case to demonstrate essential Nextflow development patterns.

        [View the training :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Develop a bulk RNAseq processing pipeline"

        A course for researchers who wish to learn how to develop their own RNAseq pipelines, using a bulk RNAseq processing use case to demonstrate essential Nextflow development patterns.

        [View the training :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Run and configure imaging pipelines"

        A course for researchers who wish to learn how to run and configure bioimaging pipelines, using nf-core/molkart to demonstrate essential Nextflow usage patterns.

        [View the training :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Setup & Help

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Training Environment__

    ---

    Options for setting up your environment for the Nextflow trainings.

    [View the training environments :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Nextflow versions__

    ---

    Understanding and managing the evolution of Nextflow's syntax versions.

    [Check version requirements :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __The Hello pipeline__

    ---

    Recap of what the Hello pipeline does and how it is structured.

    [Read the recap :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Getting help__

    ---

    Helpful resources when you have a problem with Nextflow training.

    [Find help :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
