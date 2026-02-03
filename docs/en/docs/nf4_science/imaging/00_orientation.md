# Orientation

This orientation assumes you have already opened the training environment by clicking on the "Open in GitHub Codespaces" button.
If not, please do so now, ideally in a second browser window or tab so you can refer back to these instructions.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Machine size requirement"

    Make sure to select an **8-core machine** when creating your Codespace for this training course. The bioimaging workflows require additional compute resources.

## GitHub Codespaces

The GitHub Codespaces environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) GitHub account to log in, and if you are unfamiliar with the interface you should take a few minutes to familiarize yourself with it by completing the [GitHub Codespaces Orientation](../../envsetup/index.md) mini-course.

## Pre-download Docker images

Once you've opened your Codespace, let's pre-download all the Docker images we'll need for this training course.
This will save time later and ensure smooth execution of the workflows.

Open a new terminal tab and run the following command:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

This command will download all necessary Docker images in the background.
You can continue with the rest of the orientation while this runs.

!!!tip

    The `-stub` flag allows the pipeline to run quickly without processing real data, which is perfect for downloading images. You can monitor the progress in the terminal tab.

## Working directory

Throughout this training course, we'll be working in the `nf4-science/imaging/` directory.

Change directory now by running this command in the terminal:

```bash
cd nf4-science/imaging/
```

!!!tip

    If for whatever reason you move out of this directory, you can always use the full path to return to it, assuming you're running this within the GitHub Codespaces training environment:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Now, to begin the course, click on the arrow in the bottom right corner of this page.**
