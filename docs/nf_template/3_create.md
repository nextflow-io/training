# Getting started

The nf-core pipeline template is a standardized framework designed to streamline the development of Nextflow-based bioinformatics pipelines.

Creating a pipeline using the nf-core template is greatly simplified by the nf-core tooling, which will help you create a pipeline using the set framework, which you can then modify to suit your own purposes.

Here, you will use the nf-core template to kickstart your pipeline development.

## The `nf-core create` command

The `nf-core create` command makes a new pipeline using the nf-core base template with pipeline name, description, and author. It is the first and most important step for creating a pipeline that will seamlessly integrate to the rest of the wider Nextflow ecosystem.

```bash
nf-core create
```

Although you can provide options on the command line, it’s easiest to use the interactive prompts.

```bash
? Workflow name mypipeline
? Description My first pipeline
? Author Chris
? Do you want to customize which parts of the template are used? (y/N) <n>
```

!!! question "Exercise"

    Create a pipeline named `mypipeline` using the `nf-core create` command.

!!! note

    If the command has run successfully, you will see a new folder in your current directory that has been given the name `nf-core-mypipeline`.

### Customizing the template

The nf-core pipeline comes packed with features. However, you may not want to include all of these in your pipeline.

Instead of manually removing these features once the template has been created, you can customize certain areas of the template when it is being created.

The following template areas can be customized:

-   **GitHub hosting:** Files required for GitHub hosting of the pipeline, e.g., `.github/` and `.gitignore`.
-   **GitHub CI:** Files required for GitHub continuous integration tests, e.g., `.github/workflows/`.
-   **GitHub badges:** GitHub badges in the `README.md` file.
-   **iGenomes config:** Pipeline options related to iGenomes, e.g., `conf/igenomes.config`.
-   **nf-core/configs:** Repository options that integrate nf-core config profiles.

If you choose to customize the template, the nf-core tooling will provide a series of interactive prompts help guide your choices.

```bash
? Workflow name mypipeline
? Description My first pipeline
? Author Chris
? Do you want to customize which parts of the template are used? (y/N) y
? Pipeline prefix myorg
Skip template areas?
   ○ GitHub hosting
   ○ GitHub CI
   ● GitHub badges
   ○ iGenomes config
   ○ nf-core/configs
```

!!! note

    GitHub badges were skipped in the example above.

## Working with GitHub

GitHub is a web-based interface allowing real-time collaboration. It encourages teams to work together in developing code. With GitHub, you can easily track changes, navigate revisions, and automate some of the more mundane tasks, such as testing.

It is recommended that you work with GitHub (or another code repository) when you are developing your pipelines.

### Getting started with GitHub

To get started with GitHub, you'll need to create a free personal account on [github.com](https://github.com/) and verify your email address.

Every person who uses github.com signs in to a personal account. Your personal account is your identity on github.com and has a username and profile.

Signing up for a new personal account:

1. Navigate to [https://github.com/](https://github.com/).
2. Click **Sign up**.
3. Follow the prompts to create your personal account.

!!! warning "Verified accounts"

    You won't be able to complete some basic GitHub tasks, such as creating a repository, without a verified email address.

### Submit your code to GitHub

The `nf-core create` command suggests commands for submitting this to GitHub.

Before you can use these commands, you will need to create an empty repository for your template to be pushed to.

When you are logged into GitHub, you can use the `New` repository button or navigate to [https://github.com/new](https://github.com/new) and follow the prompts to make a new repository.

Once you have created the repository you can push your template to GitHub.

```bash
cd /workspace/gitpod/nf-template/nf-core-mypipeline
git remote add origin git@github.com:<USERNAME>/<REPO>.git
git push --all origin
```

!!! question "Exercise"

    Push your new nf-core template pipeline to GitHub using the commands above.

### Working with branches

GitHub branches are used to isolate development work without affecting other branches in a repository. Each repository has one default branch, and can have multiple other branches. You can merge a branch into another branch using a pull request.

The `nf-core create` command will create three branches that will all be pushed to GitHub using the commands above.

```console
remote: Resolving deltas: 100% (10/10), done.
To https://github.com/<USERNAME>/<REPO>.git
 * [new branch]      TEMPLATE -> TEMPLATE
 * [new branch]      dev -> dev
 * [new branch]      main -> main
```

In nf-core the `main` branch for stable releases and the `dev` branch is for merging feature branches.

<figure class="excalidraw">
--8<-- "docs/nf_template/img/branches.excalidraw.svg"
</figure>

Feature branches should be checked out from the `dev` branch.

!!! question "Exercise"

    Checkout a new feature branch from the dev branch

    ```
    git checkout -b myFeature dev
    ```

You can find out more about working collaboratively with branches on the [GitHub documentation](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests).

!!! note

    You can execute different remote GitHub branches using the revision flag (`-r`).

### The `TEMPLATE` branch

The `TEMPLATE` branch is used by the `nf-core sync` command to integrate template changes to your pipeline. You should **never** modify the `TEMPLATE` branch as any changes will likely disrupt the syncing functionality.

You will learn more about the `TEMPLATE` branch in a later section.

## Template tour

The nf-core pipeline template comes packed with a lot of files and folders.

While the template may feel overwhelming a complete understanding isn't required start developing your pipeline.

### Workflows, subworkflows, and modules

The nf-core pipeline template has a `main.nf` file that calls a `mypipeline.nf` file from the `workflows` folder. The `mypipeline.nf` file is the main pipeline file and is used to bring everything else together. Instead of having one large monolithic script, it's broken up into smaller components, namely, modules and subworkflows:

-   **Modules** are wrappers around a single process.
-   **Subworkflows** are two or more modules that are packaged together.

<figure class="excalidraw">
--8<-- "docs/nf_template/img/nested.excalidraw.svg"
</figure>

Modules and subworkflows are stored within `local` and `nf-core` folders. nf-core folders are for components that have come from nf-core while local folders are for components that have been developed independently.

```
modules/
├── local
│   └── <toolname>
│       └── main.nf
└── nf-core
    ├── <toolname>
    │   ├── environment.yml
    │   ├── main.nf
    │   ├── meta.yml
    │   └── tests
    │       ├── main.nf.test
    │       ├── main.nf.test.snap
    │       └── tags.yml
```

### Configuration files

The nf-core pipeline template utilizes Nextflows flexible customization options and has a series of configuration files throughout the template.

In the template, the `nextflow.config` file is a central configuration file and is used to set default values for parameters and other configuration options.

There are several additional configuration files that are stored in the `conf` folder and are either added to the configuration scope by default or optionally as profiles:

-   `base.config`: sensible defaults for pipeline resource requests.
-   `igenomes.config`: configuration settings required to access the igenomes registry.
-   `modules.config`: additional module directives amd arguments.
-   `test.config`: a profile to run the pipeline with minimal test data.
-   `test_full.config`: a profile to run the pipeline with a full-sized test dataset.

### `.nf-core.yml`

The `.nf-core.yml` file is used to specify the repository type and manage linting tests.

```yml title=".nf-core.yml"
repository_type: pipeline
```

By default, the `.nf-core.yml` file will only show the repository is a pipeline. However, if the template is customized and parts of the template are removed, this file needs to be modified for linting tests to pass.

### `nextflow_schema.json`

The `nextflow_schema.json` is a file used to store parameter related information including type, description and help text in a machine readable format. This file is used for various purposes, including automated parameter validation, help text generation, and interactive parameter form rendering in UI interfaces.

### GitHub actions workflows

Automated workflows are an important part of the nf-core pipeline template.

By default, the template comes with several automated tests that utilize GitHub Actions, each of which are configured in the `.github/workflows` folder. Brief descriptions about each test and file are outlined below:

-   `branch.yml` (_Configured for nf-core repo only_)
    -   Sets the branch protection for the nf-core repository
-   `ci.yml` (_Configured for all pipelines_)
    -   Run small pipeline tests with the small test datasets
-   `clean-up.yml` (_Configured for nf-core repo only_)
    -   Automated testing for stale and closed GitHub issues and PRs in nf-core repo
-   `download_pipeline.yml` (_Configured for nf-core repo only_)
    -   Test a pipeline download with `nf-core download`.
-   `fix-linting.yml` (_Configured for nf-core repo only_)
    -   Fix linting by adding a comment to a PR
-   `linting_comment.yml` (_Configured for nf-core repo only_)
    -   Triggered after the linting action and posts an automated comment to the PR, even if the PR is coming from a fork
-   `linting.yml` (_Configured for all pipelines_)
    -   Triggered on pushes and PRs to the repository and runs `nf-core lint` and markdown lint tests to ensure that the code meets the nf-core guidelines
-   `release-announcements.yml` (_Configured for nf-core repo only_)
    -   Automatic release toot and tweet announcements for nf-core pipeline releases

Notable, many of these tests are only configured for the nf-core repo. However, they can be modified to better suit your needs or ignored if they are superfluous to your requirements.

You can read more about creating and modifying workflows on the [GitHub Actions documentation webpage](https://docs.github.com/en/actions).

!!! note

    To enable these workflows you need to click `Enable Actions on this Repository` under the `Actions` tab in your GitHub repository.

!!! warning "Deleting Workflows"

    Even though many of these action workflows are not relevant for private repositories, it is recommended to keep them in place to prevent `nf-core lint` from throwing errors.

## Changes to the template structure

Occasionally, the structure of the nf-core pipeline template is updated during a new release of the nf-core tooling. Most of the time these changes are minor. However, sometimes, larger structural changes are adopted to align with changes in the wider ecosystem.

For example, as of nf-core tools 2.13, the groovy code that once lived in the `lib` folder has been moved to `subworkflows/`. Moving this code has made it easier to find, modify, and test the code. Importantly, it's modular and is paving the way for a more flexible template in the future.

!!! note

    The `TEMPLATE` branch is essential for adopting these changes.
