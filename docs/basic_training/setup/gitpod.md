---
description: How to set up a development environment to run Nextflow using Gitpod
---

# Gitpod

A preconfigured Nextflow development environment is available using Gitpod.

### Requirements

-   A GitHub account
-   Web browser (Google Chrome, Firefox)
-   Internet connection

## Gitpod quick start

To run Gitpod:

-   Click the following URL: <https://gitpod.io/#https://github.com/nextflow-io/training>
    -   This is our GitHub repository URL, prefixed with `https://gitpod.io/#`
-   Log in to your GitHub account (and allow authorization).

Once you have signed in, Gitpod should load (skip prebuild if asked).

## Explore your Gitpod IDE

You should now see something similar to the following:

![Gitpod welcome](../img/gitpod.welcome.png)

-   **The sidebar** allows you to customize your Gitpod environment and perform basic tasks (copy, paste, open files, search, git, etc.). Click the Explorer button to see which files are in this repository.
-   **The terminal** allows you to run all the programs in the repository. For example, both `nextflow` and `docker` are installed and can be executed.
-   **The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. You should also see the nf-training material browser (<https://training.nextflow.io/>).

To test that the environment is working correctly, type the following into the terminal:

```bash
nextflow info
```

This should come up with the Nextflow version and runtime information:

```
Version: 22.09.7-edge build 5806
Created: 28-09-2022 09:36 UTC
System: Linux 5.15.0-47-generic
Runtime: Groovy 3.0.13 on OpenJDK 64-Bit Server VM 11.0.1+13-LTS
Encoding: UTF-8 (UTF-8)
```

## Gitpod resources

-   Gitpod gives you 500 free credits per month, which is equivalent to 50 hours of free environment runtime using the standard workspace (up to 4 cores, 8 GB RAM and 30GB storage).
-   There is also a large workspace option that gives you up to 8 cores, 16GB RAM, and 50GB storage for each workspace (up to 4 parallel workspaces in the free version). However, the large workspace uses your free credits quicker so your account will have fewer hours of access to this space.
-   Gitpod will time out after 30 minutes. However, changes are saved for up to 2 weeks (see the next section for reopening a timed out session).

See [gitpod.io](https://www.gitpod.io) for more details.

## Reopening a Gitpod session

You can reopen an environment from <https://gitpod.io/workspaces>. Find your previous environment in the list, then select the ellipsis (three dots icon) and select Open.

If you have saved the URL for your previous Gitpod environment, you can simply open it your browser to open the previous environment.

Alternatively, you can start a new workspace by following the Gitpod URL: <https://gitpod.io/#https://github.com/nextflow-io/training>

If you have lost your environment, you can find the main scripts used in this tutorial in the `nf-training` directory to resume with a new environment.

## Saving files from Gitpod to your local machine

To save any file from the explorer panel, right-click the file and select `Download`.

## Training material

The training course can be accessed in your browser from <https://training.nextflow.io/>
