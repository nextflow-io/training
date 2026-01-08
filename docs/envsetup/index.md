# Environment options

Our goal is to provide a consistent and thoroughly tested environment that allows learners to focus on learning Nextflow without having to spend time and effort on managing software.
To that end, we have developed a containerized environment that contains all the necessary software, code files and example data to work through all of our courses.

This containerized environment can be run out of the box on [Github Codespaces](01_setup.md) or in any other custom environment that supports [devcontainers](./03_devcontainer.md).
Alternatively, you can replicate this environment on your own local system by installing the software dependencies and cloning the training repository, as detailed [here](02_local.md).

For more details about the three options, please see the relevant articles.

- [Using GitHub Codespaces for Nextflow training](01_setup.md)
- [Setting up VS Code with Devcontainers](03_devcontainer.md)
- [Installing everything on your local system](02_local.md)

!!! warning

    Many of our courses require a working internet connection by default, but in most cases it is possible to pre-download or cache everything locally if necessary.

---

!!! info "Deprecation of Gitpod"

    Nextflow Training used to use [Gitpod](https://gitpod.io) until February 2025.
    However, the makers of Gitpod decided to retire the free functionality in favor of the [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) system.
    For that reason, we switched to using GitHub Codespaces, which also offer a one-click developer environment with no prior setup.

    Depending on when you signed up to Gitpod and when exactly they retire the service, you may still be able to launch the training in their old cloud IDE, though we cannot guarantee reliable access going forward:
    [Open in Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
