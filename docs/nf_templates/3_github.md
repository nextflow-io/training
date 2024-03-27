# Working with GitHub

GitHub is a web-based interface allowing real-time collaboration. It encourages teams to work together in developing code. With GitHub, you can easily track changes, navigate revisions, and automate some of the more mundane tasks, such as testing.

It is recommended that you work with GitHub (or another code repository) when you are developing your pipelines.

## Getting started with GitHub

To get started with GitHub, you'll need to create a free personal account on [github.com](https://github.com/) and verify your email address.

Every person who uses github.com signs in to a personal account.
Your personal account is your identity on github.com and has a username and profile.

Signing up for a new personal account:

1. Navigate to https://github.com/.
2. Click Sign up.
3. Follow the prompts to create your personal account.

!!! warning "Verified accounts"

    You won't be able to complete some basic GitHub tasks, such as creating a repository, without a verified email address.

## Submit your code to GitHub

The `nf-core create` command suggests commands for submitting this to GitHub.

Before you can use these commands, you will need to create an empty repository for your template to be pushed to. When you are logged into GitHub, you can use the `New` repository button or navigate to [https://github.com/new](https://github.com/new) and follow the prompts to make a new repository.

Once you have created the repository you can push your template to GitHub.

```bash
cd /Users/workspace/nf-core-mypipeline
git remote add origin git@github.com:<USERNAME>/<REPO>.git
git push --all origin
```

!!! note

    `<USERNAME>` is your GitHub handle and `<REPO>` is the repository name.

!!! question "Exercise"

    Push your new nf-core template pipeline to GitHub using the commands above.


## Working with branches

GitHub branches are used to isolate development work without affecting other branches in the repository. Each repository has one default branch, and can have multiple other branches. You can merge a branch into another branch using a pull request.

The nf-core template will create three branches that will all be pushed to GitHub.

```console
remote: Resolving deltas: 100% (10/10), done.
To https://github.com/<USERNAME>/<REPO>.git
 * [new branch]      TEMPLATE -> TEMPLATE
 * [new branch]      dev -> dev
 * [new branch]      main -> main
```

It is recommended that you use the `main` branch for stable releases and the `dev` branch for merging `feature` branches.

Feature branches should be checked out from the `dev` branch.

<figure class="excalidraw">
--8<-- "docs/nf_templates/img/3_branches.excalidraw.svg"
</figure>

!!! note "What is the `TEMPLATE` branch?"

    The `TEMPLATE` branch is used by the `nf-core sync` command to integrate template changes to your pipeline. You should **never** modify the `TEMPLATE` branch as any changes will likely disrupt the sync functionality.
