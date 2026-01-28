---
title: Environment options
description: Options for setting up your environment for the Nextflow trainings
hide:
  - toc
  - footer
---

# Environment options

We aim to provide a consistent and thoroughly tested environment that allows learners to focus on learning Nextflow without having to spend time and effort on managing software.
To that end, we have developed a containerized environment that contains all the necessary software, code files and example data to work through all of our courses.

This containerized environment can be run out of the box on Github Codespaces or locally in VS Code with the Devcontainers extension.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces is a web-based service that allows us to provide a pre-built environment for training, with all tools and data included, backed by virtual machines in the cloud. It is accessible for free to anyone with a Github account.

    [Use Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Local Devcontainers__

    ---

    VS Code with Devcontainers provides a locally-run containerized development environment with all training tools pre-configured. It offers the same pre-built environment as Codespaces but running entirely on your local hardware.

    [Use Devcontainers locally :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instructions for manual installation

If neither of the options above suit your needs, you can replicate this environment on your own local system by installing the software dependencies manually and cloning the training repository.

[Manual installation :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Deprecation of Gitpod"

    Nextflow Training used to use [Gitpod](https://gitpod.io) until February 2025.
    However, the makers of Gitpod decided to retire the free functionality in favor of the [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) system.
    For that reason, we switched to using GitHub Codespaces, which also offer a one-click developer environment with no prior setup.

    Depending on when you signed up to Gitpod and when exactly they retire the service, you may still be able to launch the training in their old cloud IDE, though we cannot guarantee reliable access going forward:
    [Open in Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
