# Getting started

The nf-core pipeline template is a standardized framework designed to streamline the development of Nextflow-based bioinformatics pipelines.

Creating a pipeline using the nf-core template is greatly simplified by the nf-core tooling. It will help you create a pipeline using the set framework that can be modified to suit your own purposes.

Here, you will use the nf-core template to kickstart your pipeline development using the latest version of Nextflow and the nf-core tooling.

## The `nf-core create` command

nf-core tooling has commands for pipeline users and developers.

```
nf-core --help
```

Here we will focus on the tooling to assist pipeline developers.

The `nf-core create` command makes a new pipeline using the nf-core base template with a pipeline name, description, and author. It is the first and most important step for creating a pipeline that will integrate with the wider Nextflow ecosystem.

```bash
nf-core create
```

When executing this command you will be given a series of prompts to setup your pipeline.

```bash
? Workflow name <name>
? Description <description>
? Author <name>
? Do you want to customize which parts of the template are used? (y/N) <n>
```

If the command has run successfully, you will see a new folder in your current directory.

!!! question "Exercise"

    Create a pipeline named `mypipeline` using the `nf-core create` command.

!!! note

    If run successfully you will see a new folder in your current directory that has been given the name `nf-core-mypipeline`.

### Customizing the template

The nf-core pipeline comes packed with features. However, you may not want to include all of these in your pipeline.

Instead of manually removing these features once the template has been created, you can customize certain areas of the template when it is being created.

The following template areas can be customized:

-   **GitHub hosting:** Files required for GitHub hosting of the pipeline, e.g., `.github/` and `.gitignore`.
-   **GitHub CI:** Files required for GitHub continuous integration tests, e.g., `.github/workflows/`.
-   **GitHub badges:** GitHub badges in the `README.md` file.
-   **iGenomes config:** Pipeline options related to iGenomes, e.g., `conf/igenomes.config`.
-   **nf-core/configs:** Repository options that integrate nf-core config profiles.

If you choose to customize the template, the nf-core tooling will provide a series of interactive prompts to help guide your choices.

```bash
? Workflow name <name>
? Description <description>
? Author <name>
? Do you want to customize which parts of the template are used? (y/N) <y>
? Pipeline prefix <organisation>
Skip template areas?
   ○ GitHub hosting
   ○ GitHub CI
   ● GitHub badges
   ○ iGenomes config
   ○ nf-core/configs
```

!!! note

    GitHub badges were skipped in the example above.

### Submit your code to GitHub

The `nf-core create` command suggests commands for submitting this to GitHub.

You will first need to create an empty repository on GitHub for your template to be pushed to.

When you are logged into GitHub, you can use the `New` repository button or navigate to [https://github.com/new](https://github.com/new) and follow the prompts to make a new repository.

![GitHub new repo](img/github.new.png)

Once you have created the repository you can push your template to GitHub.

```bash
cd /workspace/gitpod/nf-develop/nf-core-mypipeline
git remote add origin https://github.com/<USERNAME>/<REPO>.git
git push --all origin
```

!!! warning "https vs ssh"

    When using Gitpod you must use **https** not **ssh** when adding your remote:

    For example, you should use:

    ```
    git remote add origin https://github.com/<USERNAME>/<REPO>.git
    ```

    Not:

    ```
    git remote add origin git@github.com:<USERNAME>/<REPO>.git
    ```

By default, three branches will all be pushed to GitHub using the commands above.

```console
remote: Resolving deltas: 100% (10/10), done.
To https://github.com/<USERNAME>/myfirstpipeline.git
 * [new branch]      TEMPLATE -> TEMPLATE
 * [new branch]      dev -> dev
 * [new branch]      main -> main
```

The role each of these branches have in pipeline development will be explained in subsequent sections.

!!! question "Exercise"

    Create a new GitHub repository named `myfirstpipeline` and push your new pipeline using the commands above. You will need to replace `<USERNAME>` and `<REPO>` with your GitHub username and `myfirstpipeline`, respectively.

## Template tour

The nf-core pipeline template comes packed with a lot of files and folders.

While the template may feel overwhelming, a complete understanding isn't required to start developing your pipeline.

### Workflows, subworkflows, and modules

The nf-core pipeline template has a `main.nf` file that calls `mypipeline.nf` from the `workflows` folder. The `mypipeline.nf` file is the central pipeline file that is used to bring everything else together. Instead of having one large monolithic pipeline script, it's broken up into smaller script components, namely, modules and subworkflows:

-   **Modules** are wrappers around a single process.
-   **Subworkflows** are two or more modules that are packaged together.

<figure class="excalidraw">
--8<-- "docs/nf_template/img/nested.excalidraw.svg"
</figure>

Within your pipeline repository, `modules` and `subworkflows` are stored within `local` and `nf-core` folders. The `nf-core` folder is for components that have come from the online nf-core repository while the `local` folder is for components that have been developed independently.

```
modules/
├── local
│   └── <toolname>
│   │   └── main.nf
│   .
│
└── nf-core
    ├── <toolname>
    │   ├── environment.yml
    │   ├── main.nf
    │   ├── meta.yml
    │   └── tests
    │       ├── main.nf.test
    │       ├── main.nf.test.snap
    │       └── tags.yml
    .
```

Modules from nf-core follow the same structure and contain a small number of additional files that are used for testing and documentation.

### Configuration files

The nf-core pipeline template utilizes Nextflows flexible customization options and has a series of configuration files throughout the template.

In the template, the `nextflow.config` file is a central configuration file and is used to set default values for parameters and other configuration options. The majority of these configuration options are applied by default while others (e.g., software dependency profiles) are included as optional profiles.

There are several configuration files that are stored in the `conf` folder and are added to the configuration by default or optionally as profiles:

-   `base.config`: sensible defaults for pipeline resource requests.
-   `igenomes.config`: configuration settings required to access the igenomes registry.
-   `modules.config`: additional module directives and arguments.
-   `test.config`: a profile to run the pipeline with minimal test data.
-   `test_full.config`: a profile to run the pipeline with a full-sized test dataset.

### `.nf-core.yml`

The `.nf-core.yml` file is used to specify the repository type and manage linting tests.

```yml title=".nf-core.yml"
repository_type: pipeline
```

By default, the `.nf-core.yml` file will only show the repository is a pipeline. However, if the template is customized and parts of the template are removed, this file needs to be modified for linting tests to pass.

### `nextflow_schema.json`

The `nextflow_schema.json` is a file used to store parameter related information including type, description and help text in a machine readable format. The schema is used for various purposes, including automated parameter validation, help text generation, and interactive parameter form rendering in UI interfaces.

### GitHub actions workflows

Automated workflows are an important part of the nf-core pipeline template.

By default, the template comes with several automated tests that utilize GitHub Actions, each of which are configured in the `.github/workflows` folder:

-   `branch.yml`: Sets the branch protection for the nf-core repository
-   `ci.yml`: Run small pipeline tests with the small test datasets
-   `clean-up.yml`: Automated testing for stale and closed GitHub issues and PRs in the nf-core repo
-   `download_pipeline.yml`: Test a pipeline download with `nf-core download`.
-   `fix-linting.yml`: Fix linting by adding a comment to a PR
-   `linting_comment.yml`: Triggered after the linting action and posts an automated comment to the PR, even if the PR is coming from a fork
-   `linting.yml`: Triggered on pushes and PRs to the repository and runs `nf-core lint` and markdown lint tests to ensure that the code meets the nf-core guidelines
-   `release-announcements.yml`: Automatic release toot and tweet announcements for nf-core pipeline releases

Notably, many of these tests are only configured for the nf-core repo. However, they can be modified for your repository or ignored if they are superfluous to your requirements.

You can read more about creating and modifying workflows on the [GitHub Actions documentation webpage](https://docs.github.com/en/actions).

Even though many of these action workflows are not relevant for private repositories, it is recommended to keep them in place to prevent `nf-core lint` from throwing errors.

!!! note

    To enable these workflows you need to click `Enable Actions on this Repository` under the `Actions` tab in your GitHub repository.

    ![GitHub actions](img/github.actions.png)
