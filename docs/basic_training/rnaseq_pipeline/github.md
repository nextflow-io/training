# Run a project from GitHub

Nextflow allows the execution of a pipeline project directly from a GitHub repository (or similar services, e.g., BitBucket and GitLab).

This simplifies the sharing and deployment of complex projects and tracking changes in a consistent manner.

The following GitHub repository hosts a complete version of the workflow introduced in this tutorial: <https://github.com/nextflow-io/rnaseq-nf>

You can run it by specifying the project name and launching each task of the execution as a Docker container run command:

```bash
nextflow run nextflow-io/rnaseq-nf -with-docker
```

It automatically downloads the container and stores it in the `$HOME/.nextflow` folder.

Use the command `info` to show the project information:

```bash
nextflow info nextflow-io/rnaseq-nf
```

Nextflow allows the execution of a specific revision of your project by using the `-r` command line option. For example:

```bash
nextflow run nextflow-io/rnaseq-nf -r v2.1 -with-docker
```

Revision are defined by using Git tags or branches defined in the project repository.

Tags enable precise control of the changes in your project files and dependencies over time.
