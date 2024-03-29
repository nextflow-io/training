# Creating a pipeline using the nf-core template

The nf-core pipeline template is a standardized framework designed to streamline the development of Nextflow-based bioinformatics pipelines.

Creating a pipeline using the nf-core template is greatly simplified by the nf-core tooling, which will help you create a pipeline using the set framework, which you can then modify to suit your own purposes.

## The `nf-core create` command

The `nf-core create` command makes a new pipeline using the nf-core base template with pipeline name, description and author. It is the first and most important step for creating a pipeline that will follow nf-core best practices and seamlessly plugin to the rest of the nf-core ecosystem.

```bash
nf-core create
```

Although you can provide options on the command line, it’s easiest to use the interactive prompts.

```bash
? Workflow name <pipeline name>
? Description <pipeline description>
? Author <your name>
? Do you want to customize which parts of the template are used? (y/N) n
```

If the command has run successfully, you will see a new folder in your current directory that has been given the name `nf-core-<pipeline name>`.

!!! note

    If you do not customize the template, all features of the template will be included by default.

### Customizing the template

The nf-core pipeline comes packed with features. However, you may not want to include all of these in your pipeline.

Instead of manually removing these features once the template has been creates, you can customize certain areas of the template when it is being created.

The following template areas can be customized:

-   **GitHub hosting:** Files required for GitHub hosting of the pipeline. E.g., `.github/` and `.gitignore`.
-   **GitHub CI:** Files required for GitHub continuous integration tests. E.g., `.github/workflows/`.
-   **GitHub badges:** GitHub badges in the `README.md` file.
-   **iGenomes config:** Pipeline options related to iGenomes. E.g., `conf/igenomes.config`.
-   **nf-core/configs:** Repository options that integrate nf-core config profiles.

If you choose to customize the template, the nf-core tooling will provide a series of interactive prompts help guide your choices.

```bash
? Workflow name mypipeline
? Description My pipeline
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

!!! question "Exercise"

    Create a pipeline named `mypipeline` using the `nf-core create` command.

    Do **not** modify the template.

## Working with GitHub

GitHub is a web-based interface allowing real-time collaboration. It encourages teams to work together in developing code. With GitHub, you can easily track changes, navigate revisions, and automate some of the more mundane tasks, such as testing.

It is recommended that you work with GitHub (or another code repository) when you are developing your pipelines.

### Getting started with GitHub

To get started with GitHub, you'll need to create a free personal account on [github.com](https://github.com/) and verify your email address.

Every person who uses github.com signs in to a personal account.
Your personal account is your identity on github.com and has a username and profile.

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
cd /workspace/gitpod/docs/nf_templates/nf-core-mypipeline
git remote add origin git@github.com:<USERNAME>/<REPO>.git
git push --all origin
```

!!! note

    `<USERNAME>` is your GitHub handle and `<REPO>` is the repository name.

!!! question "Exercise"

    Push your new nf-core template pipeline to GitHub using the commands above.

### Working with branches

GitHub branches are used to isolate development work without affecting other branches in the repository. Each repository has one default branch, and can have multiple other branches. You can merge a branch into another branch using a pull request.

The `nf-core create` command will create three branches that will all be pushed to GitHub.

```console
remote: Resolving deltas: 100% (10/10), done.
To https://github.com/<USERNAME>/<REPO>.git
 * [new branch]      TEMPLATE -> TEMPLATE
 * [new branch]      dev -> dev
 * [new branch]      main -> main
```

It is recommended that you use the `main` branch for stable releases and the `dev` branch for merging `feature` branches.

<figure class="excalidraw">
--8<-- "docs/nf_templates/img/3_branches.excalidraw.svg"
</figure>

!!! note "Feature branches"

    Feature branches should be checked out from the `dev` branch.

You can find out more about working collaboratively with branches on the [GitHub documentation](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests) website.

!!! note "What is the `TEMPLATE` branch?"

    The `TEMPLATE` branch is used by the `nf-core sync` command to integrate template changes to your pipeline. You should **never** modify the `TEMPLATE` branch as any changes will likely disrupt the sync functionality.

    You will learn more about the TEMPLATE branch in a later section.
