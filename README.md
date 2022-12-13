# Nextflow Training Guide.

Welcome to the Nextflow training repo. 

We are excited to have you on the path to writing reproducible and scalable scientific workflows using Nextflow. To access the material for training, click [here](https://training.seqera.io).

Our training course complements the full Nextflow documentation - if you ever have any doubts, head over to the docs located [here](https://www.nextflow.io/docs/latest/).

There are two main ways to get started with Seqera's Nextflow training course.

1. Install Locally - best if you are already confident with Git and Docker, or working offline. Follow the instructions [here](https://training.seqera.io/#_local_installation), section 1.1.

2. Gitpod - (recommended), is a containerized environment with all the programs and data pre-installed. Simply click the link and login via a GitHub account to start the tutorial. The full instructions are below.

## Gitpod requirements:

- A GitHub account
- Web browser (Google Chrome, Firefox)
- Internet connection

## Gitpod quick start

To run Gitpod:

- Click the following URL:

https://gitpod.io/#https://github.com/seqeralabs/nf-training-public

(which is our Github repository URL, prefixed with https://gitpod.io/#).

- Log in to your Github account (and allow authorization).

Once you have signed in, Gitpod should load (skip prebuild if asked).

## Explore your Gitpod IDE

You should now see something similar to the following:

![PNG](/asciidocs/img/gitpod.welcome.png)

**The sidebar** allows you to customize your Gitpod environment and perform basic tasks (copy, paste, open files, search, git, etc.). Click the Explorer button to see which files are in this repository.

**The terminal** allows you to run all the programs in the repository. For example, both `nextflow` and `docker` are installed and can be executed.

**The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. You should also see the nf-training material browser (https://training.seqera.io/).

To test that the environment is working correctly, type the following into the terminal:

	nextflow info

This should come up with the Nextflow version and runtime information:

    Version: 22.04.2 build 5701
    Created: 16-05-2022 17:52 UTC
    System: Linux 5.16.20-051620-generic
    Runtime: Groovy 3.0.10 on OpenJDK 64-Bit Server VM 11.0.13+8-LTS
    Encoding: UTF-8 (UTF-8)

### Gitpod resources

- Gitpod gives you up to 50 hours per month to run the environment for free.
- It includes up to 16 cpus and 30GB of workspace.
- Gitpod will timeout after 30 minutes. However any changes are saved for up to two week (see next section for reopening a timed out session).

See www.gitpod.io for more details.

### Reopening a Gitpod session

You can reopen an environment by going to https://gitpod.io/workspaces and finding your previous environment, then clicking the button with three dots and selecting Open.

If you save the URL from your previous Gitpod environment, you can just paste this into your browser to open the previous environment.

Alternatively, you can start a new workspace by following the Gitpod URL:
https://gitpod.io/#https://github.com/seqeralabs/nf-training-public

This tutorial provides all the scripts, so don't worry if you have lost your environment. In the `nf-training` directory, you can find the main scripts used in the tutorial.

### Saving files from Gitpod to your local machine.

To save your files, select your file of interest from the explorer panel, then right click the file to click `Download`


### Copyright
<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a>

Copyright 2020-2022, Seqera Labs. All examples and descriptions are licensed under the <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
