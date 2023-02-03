---
title: Conda
---

# Conda / Bioconda packages

Conda is a popular package and environment manager. The built-in support for Conda allows Nextflow pipelines to automatically create and activate the Conda environment(s), given the dependencies specified by each process.

In this Gitpod environment, conda is already installed.

## Using conda

A Conda environment is defined using a YAML file, which lists the required software packages. The first thing you need to do is to initiate conda for shell interaction, and then open a new terminal by running bash.

```bash
conda init
bash
```

Then write your YAML file (to `env.yml`). For example:

```yaml
name: nf-tutorial
channels:
    - conda-forge
    - defaults
    - bioconda
dependencies:
    - bioconda::salmon=1.5.1
    - bioconda::fastqc=0.11.9
    - bioconda::multiqc=1.12
    - conda-forge::tbb=2020.2
```

Given the recipe file, the environment is created using the command shown below:

```bash
conda env create --file env.yml
```

You can check the environment was created successfully with the command shown below:

```bash
conda env list
```

This should look something like this:

```bash
 conda environments:

base                  *  /opt/conda
nf-tutorial              /opt/conda/envs/nf-tutorial
```

To enable the environment, you can use the `activate` command:

```bash
conda activate nf-tutorial
```

Nextflow is able to manage the activation of a Conda environment when its directory is specified using the `-with-conda` option (using the same path shown in the `list` function. For example:

```bash
nextflow run script7.nf -with-conda /opt/conda/envs/nf-tutorial
```

!!! info

    When creating a Conda environment with a YAML recipe file, Nextflow automatically downloads the required dependencies, builds the environment and activates it.

This makes easier to manage different environments for the processes in the workflow script.

See the [docs](https://www.nextflow.io/docs/latest/conda.html) for details.

## Create and use conda-like environments using micromamba

Another way to build conda-like environments is through a `Dockerfile` and [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

`micromamba` is a fast and robust package for building small conda-based environments.

This saves having to build a conda environment each time you want to use it (as outlined in previous sections).

To do this, you simply require a `Dockerfile` and you use micromamba to install the packages. However, a good practice is to have a YAML recipe file like in the previous section, so we’ll do it here too.

```yaml
name: nf-tutorial
channels:
    - conda-forge
    - defaults
    - bioconda
dependencies:
    - bioconda::salmon=1.5.1
    - bioconda::fastqc=0.11.9
    - bioconda::multiqc=1.12
    - conda-forge::tbb=2020.2
```

Then, we can write our Dockerfile with micromamba installing the packages from this recipe file.

```dockerfile
FROM mambaorg/micromamba:0.25.1

LABEL image.author.name "Your Name Here"
LABEL image.author.email "your@email.here"

COPY --chown=$MAMBA_USER:$MAMBA_USER micromamba.yml /tmp/env.yml

RUN micromamba create -n nf-tutorial

RUN micromamba install -y -n nf-tutorial -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/nf-tutorial/bin:$PATH
```

The above `Dockerfile` takes the parent image _mambaorg/micromamba_, then installs a `conda` environment using `micromamba`, and installs `salmon`, `fastqc` and `multiqc`.

Try executing the RNA-Seq pipeline from earlier (script7.nf). Start by building your own micromamba `Dockerfile` (from above), save it to your docker hub repo, and direct Nextflow to run from this container (changing your `nextflow.config`).

!!! warning

    Building a Docker container and pushing to your personal repo can take &gt;10 minutes.

??? example "For an overview of steps to take, click here:"

    1. Make a file called `Dockerfile` in the current directory (with the code above).

    2. Build the image: `docker build -t my-image .` (don’t forget the _._).

    3. Publish the docker image to your online docker account.

        Something similar to the following, with `<myrepo>` replaced with your own Docker ID, without _&lt;_ and _&gt;_ characters!

        `my-image` could be replaced with any name you choose. As good practice, choose something memorable and ensure the name matches the name you used in the previous command.

        ```bash
        docker login
        docker tag my-image <myrepo>/my-image
        docker push <myrepo>/my-image
        ```

    4. Add the image file name to the `nextflow.config` file.

        e.g. remove the following from the `nextflow.config`:

        ```groovy
        process.container = 'nextflow/rnaseq-nf'
        ```

        and replace with:

        ```groovy
        process.container = '<myrepo>/my-image'
        ```

    5. Trying running Nextflow, e.g.:

        ```bash
        nextflow run script7.nf -with-docker
        ```

    Nextflow should now be able to find `salmon` to run the process.
