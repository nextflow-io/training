# Orientation

This workshop should be completed using the [Nextflow Training Gitpod Environment](https://gitpod.io/#https://github.com/nextflow-io/training).

Gitpod offers a virtual machine with everything already set up for you, accessible from your web browser or built into your code editor (e.g., VSCode).

To start, click on the button below.

[![Open in Gitpod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/nextflow-io/training)

## Running Gitpod

To run Gitpod:

-   Click the following URL: <https://gitpod.io/#https://github.com/nextflow-io/training>
    -   This is our GitHub repository URL, prefixed with `https://gitpod.io/#`
-   Log in with your GitHub account.

### Getting Started with GitHub

    [GitHub](https://github.com/) is a web-based interface allowing real-time collaboration. It encourages teams to work together in developing code. With GitHub, you can easily track changes, navigate revisions, and automate some of the more mundane tasks, such as testing.

    To get started with GitHub, you will need to create a free personal account on [github.com](https://github.com/) and verify your email address.

    Every person who uses github.com signs in to a personal account. Your account is your identity and has a username and profile.

    Signing up for a new personal account:

    1. Navigate to [https://github.com/](https://github.com/).
    2. Click `Sign up` at the top right corner.
    3. Follow the prompts to create your account.

!!! warning "Verified accounts"

    You won't be able to complete some basic GitHub tasks, such as creating a repository, without a verified email address.

When you first sign up to Gitpod you may need to authenticate your account.

### Explore your Gitpod IDE

You should now see something similar to the following:

![Gitpod welcome](img/gitpod.welcome.png)

-   **The sidebar** allows you to customize your Gitpod environment and perform basic tasks (copy, paste, open files, search, git, etc.). You can click the explorer button to see which files are in this repository.
-   **The terminal** allows you to run all the programs in the repository. For example, both `nextflow` and `docker` are installed and can be executed.
-   **The file explorer** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window.
-   **The Simple Browser** lets you view the nf-training material browser (<https://training.nextflow.io/>). If you close it by accident, you can load the simple browser again by typing the following in the terminal: `gp preview https://training.nextflow.io`.

### Gitpod resources

-   Gitpod gives you 500 free credits per month, which is equivalent to 50 hours of free environment runtime using the standard workspace (up to 4 cores, 8 GB RAM, and 30 GB storage).
-   There is also a large workspace option that gives you up to 8 cores, 16GB RAM, and 50GB storage. However, the large workspace will use your free credits quicker and you will have fewer hours of access to this space.
-   Gitpod will time out after 30 minutes of inactivity and will save your changes for up to 2 weeks.

More information about Gitpod is available at [gitpod.io](https://www.gitpod.io).

### Reopening a Gitpod session

You can reopen an environment from <https://gitpod.io/workspaces>. Previous environments will be listed. You can select the ellipsis (three dots icon) and then select `Open` to reopen a previous environment.

If you have saved the URL for your previous Gitpod environment, you can simply open it in your browser.

Alternatively, you can open a new training workspace by following the Gitpod URL: <https://gitpod.io/#https://github.com/nextflow-io/training>

### Saving files from Gitpod to your local machine

To save any file from the explorer panel, right-click the file and select `Download`.

## Getting started

You will need to move into the `troubleshoot` folder, copy, and unzip the data required for this training.

!!! question "Exercise"

    Use the following command to switch to the empty `troubleshoot` folder:

    ```bash
    cd /workspace/gitpod/troubleshoot
    cp -r /workspace/gitpod/hello-nextflow/data/ .
    tar -zxvf data/ref.tar.gz -C data/
    ```

To check everything is working as expected you can run the `hello-gatk.nf` script located in the `troubleshoot` folder.

If all of the data has been copied and unzipped correctly you should see the pipeline execute three processes:

- `SAMTOOLS_INDEX`
- `GATK_HAPLOTYPECALLER`
- `GATK_JOINTGENOTYPING`
