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

    The training courses listed below are designed to be usable as a self-service resource.
    You can work through them on your own at any time, either in the web-based environment we provide via Github Codespaces or in your own environment.

    [Explore the courses :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Additional information__

    ---

    ??? warning "Version compatibility"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **As of January 2026, all of our Nextflow training courses require Nextflow version 25.10.2 or later, with strict syntax activated, unless otherwise noted.**

        For more information about version requirements and strict syntax, please see the [Nextflow versions](info/nxf_versions.md) guide.

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

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        This training material is developed and maintained by [Seqera](https://seqera.io) and released under an open-source license ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) for the benefit of the community. If you wish to use this material in a way that falls outside the scope of the license (note the limitations on commercial use and redistribution), please contact us at [community@seqera.io](mailto:community@seqera.io) to discuss your request.

        We welcome improvements, fixes and bug reports from the community. Every page has a :material-file-edit-outline: icon in the top right of the page linking to the code repository, where you can report issues or propose changes to the training source material via a pull request. See the `README.md` in the repository for more details.

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

    ### :material-tune:{.nextflow-primary} Configure execution {.mt-1}

    Learn to configure your pipelines effectively.

    ??? courses "**Nextflow Config:** Configure pipelines like a pro"

        A hands-on introduction to configuring Nextflow pipeline execution: managing inputs and outputs, adapting to different compute environments, and switching between preset configuration profiles.

        [View the training :material-arrow-right:](nextflow_config/index.md){ .md-button .md-button--secondary }

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

    ??? courses "**Nextflow for Genomics:** Develop a variant calling pipeline"

        A course for researchers who wish to learn how to develop their own genomics pipelines, using a variant calling use case to demonstrate essential Nextflow development patterns.

        [View the training :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq:** Develop an RNAseq processing pipeline"

        A course for researchers who wish to learn how to develop their own RNAseq pipelines, using a bulk RNAseq processing use case to demonstrate essential Nextflow development patterns.

        [View the training :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Bioimaging:** Run and configure imaging pipelines"

        A course for researchers who wish to learn how to run and configure bioimaging pipelines, using nf-core/molkart to demonstrate essential Nextflow usage patterns.

        [View the training :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

!!! info "Looking for archived training materials?"

    Older training materials (Fundamentals Training, Advanced Training, and other experimental courses) have been removed from the training portal as they are incompatible with Nextflow 3.0 strict syntax.
    If you need access to these materials, they are available in the [git history](https://github.com/nextflow-io/training) prior to January 2026.

    The previous version of the Nextflow Run course, superseded by the current one, is still browsable at [archive/nextflow_run](archive/nextflow_run/index.md).

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
