---
description: Basic Nextflow Training Workshop
---

# Manage dependencies and containers

Computational workflows are rarely composed of a single script or tool. More often, they depend on dozens of software components or libraries.

Installing and maintaining such dependencies is a challenging task and the most common source of irreproducibility in scientific applications.

To overcome these issues, we use containers that enable software dependencies, i.e. tools and libraries required by a data analysis application, to be encapsulated in one or more self-contained, ready-to-run, immutable Linux container images, that can be easily deployed in any platform that supports the container runtime.

Containers can be executed in an isolated manner from the hosting system. Having its own copy of the file system, processing space, memory management, etc.

!!! info

    Containers were first introduced with kernel 2.6 as a Linux feature known as _Control Groups_ or [Cgroups](https://en.wikipedia.org/wiki/Cgroups).

## Docker

Docker is a handy management tool to build, run and share container images.

These images can be uploaded and published in a centralized repository known as [Docker Hub](https://hub.docker.com), or hosted by other parties like [Quay](https://quay.io).

### Run a container

A container can be run using the following command:

```bash
docker run <container-name>
```

Try for example the following publicly available container (if you have docker installed):

```bash
docker run hello-world
```

### Pull a container

The pull command allows you to download a Docker image without running it. For example:

```bash
docker pull debian:stretch-slim
```

The above command downloads a Debian Linux image. You can check it exists by using:

```bash
docker images
```

### Run a container in interactive mode

Launching a BASH shell in the container allows you to operate in an interactive mode in the containerized operating system. For example:

```
docker run -it debian:stretch-slim bash
```

Once launched, you will notice that it is running as root (!). Use the usual commands to navigate the file system. This is useful to check if the expected programs are present within a container.

To exit from the container, stop the BASH session with the `exit` command.

### Your first Dockerfile

Docker images are created by using a so-called `Dockerfile`, which is a simple text file containing a list of commands to assemble and configure the image with the software packages required.

Here, you will create a Docker image containing cowsay and the Salmon tool.

!!! warning

    The Docker build process automatically copies all files that are located in the current directory to the Docker daemon in order to create the image. This can take a lot of time when big/many files exist. For this reason, it’s important to _always_ work in a directory containing only the files you really need to include in your Docker image. Alternatively, you can use the `.dockerignore` file to select paths to exclude from the build.

Use your favorite editor (e.g., `vim` or `nano`) to create a file named `Dockerfile` and copy the following content:

```dockerfile
FROM debian:stretch-slim

LABEL image.author.name "Your Name Here"
LABEL image.author.email "your@email.here"

RUN apt-get update && apt-get install -y curl cowsay

ENV PATH=$PATH:/usr/games/
```

### Build the image

Build the Dockerfile image by using the following command:

```bash
docker build -t my-image .
```

Where "my-image" is the user-specified name tag for the Dockerfile, present in the current directory.

!!! tip

    Don’t miss the dot in the above command.

When it completes, verify that the image has been created by listing all available images:

```bash
docker images
```

You can try your new container by running this command:

```bash
docker run my-image cowsay Hello Docker!
```

### Add a software package to the image

Add the Salmon package to the Docker image by adding the following snippet to the `Dockerfile`:

```dockerfile
RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz | tar xz \
&& mv /salmon-*/bin/* /usr/bin/ \
&& mv /salmon-*/lib/* /usr/lib/
```

Save the file and build the image again with the same command as before:

```bash
docker build -t my-image .
```

You will notice that it creates a new Docker image with the same name **but** with a different image ID.

### Run Salmon in the container

Check that Salmon is running correctly in the container as shown below:

```bash
docker run my-image salmon --version
```

You can even launch a container in an interactive mode by using the following command:

```bash
docker run -it my-image bash
```

Use the `exit` command to terminate the interactive session.

### File system mounts

Create a genome index file by running Salmon in the container.

Try to run Salmon in the container with the following command:

```bash
docker run my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

The above command fails because Salmon cannot access the input file.

This happens because the container runs in a completely separate file system and it cannot access the hosting file system by default.

You will need to use the `--volume` command-line option to mount the input file(s) e.g.

```bash
docker run --volume $PWD/data/ggal/transcriptome.fa:/transcriptome.fa my-image \
    salmon index -t /transcriptome.fa -i transcript-index
```

!!! warning

    The generated `transcript-index` directory is still not accessible in the host file system.

An easier way is to mount a parent directory to an identical one in the container, this allows you to use the same path when running it in the container e.g.

```bash
docker run --volume $PWD:$PWD --workdir $PWD my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Or set a folder you want to mount as an environmental variable, called `DATA`:

```bash
DATA=/workspace/gitpod/nf-training/data
docker run --volume $DATA:$DATA --workdir $PWD my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Now check the content of the `transcript-index` folder by entering the command:

```bash
ls -la transcript-index
```

!!! note

    Note that the permissions for files created by the Docker execution is `root`.

### Upload the container in the Docker Hub (bonus)

Publish your container in the Docker Hub to share it with other people.

Create an account on the <https://hub.docker.com> website. Then from your shell terminal run the following command, entering the user name and password you specified when registering in the Hub:

```bash
docker login
```

Tag the image with your Docker user name account:

```bash
docker tag my-image <user-name>/my-image
```

Finally push it to the Docker Hub:

```bash
docker push <user-name>/my-image
```

After that anyone will be able to download it by using the command:

```bash
docker pull <user-name>/my-image
```

Note how after a pull and push operation, Docker prints the container digest number e.g.

```console
Digest: sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266
Status: Downloaded newer image for nextflow/rnaseq-nf:latest
```

This is a unique and immutable identifier that can be used to reference a container image in a univocally manner. For example:

```bash
docker pull nextflow/rnaseq-nf@sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266
```

### Run a Nextflow script using a Docker container

The simplest way to run a Nextflow script with a Docker image is using the `-with-docker` command-line option:

```bash
nextflow run script2.nf -with-docker my-image
```

As seen in the last section, you can also configure the Nextflow config file (`nextflow.config`) to select which container to use instead of having to specify it as a command-line argument every time.

## Singularity

[Singularity](http://singularity.lbl.gov) is a container runtime designed to work in high-performance computing data centers, where the usage of Docker is generally not allowed due to security constraints.

Singularity implements a container execution model similar to Docker. However, it uses a completely different implementation design.

A Singularity container image is archived as a plain file that can be stored in a shared file system and accessed by many computing nodes managed using a batch scheduler.

!!! warning

    Singularity will not work with Gitpod. If you wish to try this section, please do it locally, or on an HPC.

### Create a Singularity images

Singularity images are created using a `Singularity` file in a similar manner to Docker but using a different syntax.

```singularity
Bootstrap: docker
From: debian:stretch-slim

%environment
export PATH=$PATH:/usr/games/

%labels
AUTHOR <your name>

%post

apt-get update && apt-get install -y locales-all curl cowsay
curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.0.0/salmon-1.0.0_linux_x86_64.tar.gz | tar xz \
&& mv /salmon-*/bin/* /usr/bin/ \
&& mv /salmon-*/lib/* /usr/lib/
```

Once you have saved the `Singularity` file, you can create the image with these commands:

```bash
sudo singularity build my-image.sif Singularity
```

Note: the `build` command requires `sudo` permissions. A common workaround consists of building the image on a local workstation and then deploying it in the cluster by copying the image file.

### Running a container

Once done, you can run your container with the following command

```bash
singularity exec my-image.sif cowsay 'Hello Singularity'
```

By using the `shell` command you can enter in the container in interactive mode. For example:

```bash
singularity shell my-image.sif
```

Once in the container instance run the following commands:

```bash
touch hello.txt
ls -la
```

!!! info

    Note how the files on the host environment are shown. Singularity automatically mounts the host `$HOME` directory and uses the current work directory.

### Import a Docker image

An easier way to create a Singularity container without requiring `sudo` permission and boosting the containers interoperability is to import a Docker container image by pulling it directly from a Docker registry. For example:

```bash
singularity pull docker://debian:stretch-slim
```

The above command automatically downloads the Debian Docker image and converts it to a Singularity image in the current directory with the name `debian-jessie.simg`.

### Run a Nextflow script using a Singularity container

Nextflow allows the transparent usage of Singularity containers as easy as with Docker.

Simply enable the use of the Singularity engine in place of Docker in the Nextflow configuration file by using the `-with-singularity` command-line option:

```bash
nextflow run script7.nf -with-singularity nextflow/rnaseq-nf
```

As before, the Singularity container can also be provided in the Nextflow config file. We’ll see how to do this later.

### The Singularity Container Library

The authors of Singularity, [SyLabs](https://www.sylabs.io/) have their own repository of Singularity containers.

In the same way that we can push Docker images to Docker Hub, we can upload Singularity images to the Singularity Library.

## Conda/Bioconda packages

Conda is a popular package and environment manager. The built-in support for Conda allows Nextflow pipelines to automatically create and activate the Conda environment(s), given the dependencies specified by each process.

In this Gitpod environment, conda is already installed.

### Using conda

A Conda environment is defined using a YAML file, which lists the required software packages. The first thing you need to do is to initiate conda for shell interaction, and then open a new terminal by running bash.

```bash
conda init
bash
```

Then write your YAML file (to `env.yml`). There is already a file named `env.yml` in the `nf-training` folder as an example. Its content is shown below.

```yaml
--8<-- "nf-training/env.yml"
```

Given the recipe file, the environment is created using the command shown below. The `conda env create` command may take several minutes, as conda tries to resolve dependencies of the desired packages at runtime, and then downloads everything that is required.

```bash
conda env create --file env.yml
```

You can check the environment was created successfully with the command shown below:

```bash
conda env list
```

This should look something like this:

```bash
# conda environments:
#
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

### Create and use conda-like environments using micromamba

Another way to build conda-like environments is through a `Dockerfile` and [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

`micromamba` is a fast and robust package for building small conda-based environments.

This saves having to build a conda environment each time you want to use it (as outlined in previous sections).

To do this, you simply require a `Dockerfile` and you use micromamba to install the packages. However, a good practice is to have a YAML recipe file like in the previous section, so we’ll do it here too, using the same `env.yml` as before.

```yaml
--8<-- "nf-training/env.yml"
```

Then, we can write our Dockerfile with micromamba installing the packages from this recipe file.

```dockerfile
FROM mambaorg/micromamba:0.25.1

LABEL image.author.name "Your Name Here"
LABEL image.author.email "your@email.here"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

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

    4. Add the container image name to the `nextflow.config` file.

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

## BioContainers

Another useful resource linking together Bioconda and containers is the [BioContainers](https://biocontainers.pro) project. BioContainers is a community initiative that provides a registry of container images for every Bioconda recipe.

So far, we’ve seen how to install packages with conda and micromamba, both locally and within containers. With BioContainers, you don’t need to create your own container image for the tools you want, and you don’t need to use conda or micromamba to install the packages. It already provides you with a Docker image containing the programs you want installed. For example, you can get the container image of fastqc using BioContainers with:

```bash
docker pull biocontainers/fastqc:v0.11.5
```

You can check the registry for the packages you want in [BioContainers official website](https://biocontainers.pro/registry).

Contrary to other registries that will pull the latest image when no tag (version) is provided, you must specify a tag when pulling BioContainers (after a colon `:`, e.g. fastqc:v0.11.5). Check the tags within the registry and pick the one that better suits your needs.

!!! tip

    You can have more complex definitions within your process block by letting the appropriate container image or conda package be used depending on if the user selected singularity, Docker or conda to be used. You can click [here](https://nf-co.re/docs/contributing/modules#software-requirements) for more information and [here](https://github.com/nf-core/modules/blob/61f68913fefc20241ceccb671b104230b2d775d7/modules/bowtie2/align/main.nf#L6-L9) for an example.

### :material-progress-question: Exercises

!!! exercise

    During the earlier RNA-Seq tutorial (script2.nf), we created an index with the salmon tool. Given we do not have salmon installed locally in the machine provided by Gitpod, we had to either run it with `-with-conda` or `-with-docker`. Your task now is to run it again `-with-docker`, but without having to create your own Docker container image. Instead, use the BioContainers image for salmon 1.7.0.


    ??? result

        ```bash
        nextflow run script2.nf -with-docker quay.io/biocontainers/salmon:1.7.0--h84f40af_0
        ```

!!! exercise "Bonus Exercise"

    Change the process directives in `script5.nf` or the `nextflow.config` file to make the pipeline automatically use BioContainers when using salmon, or fastqc.

    !!! tip "Hint"

        Temporarily comment out the line `#!groovy process.container = 'nextflow/rnaseq-nf'` in the `nextflow.config` file to make sure the processes are using the BioContainers that you set, and not the container image we have been using in this training.

    ??? result

        With these changes, you should be able to run the pipeline with BioContainers by running the following in the command line:

        ```bash
        nextflow run script5.nf
        ```

        with the following container directives for each process:

        ```groovy
        process FASTQC {
            container 'biocontainers/fastqc:v0.11.5'
            tag "FASTQC on $sample_id"
        ...
        ```

        and

        ```groovy
        process QUANTIFICATION {
            tag "Salmon on $sample_id"
            container 'quay.io/biocontainers/salmon:1.7.0--h84f40af_0'
            publishDir params.outdir, mode:'copy'
        ...
        ```

        Check the `.command.run` file in the work directory and ensure that the run line contains the correct Biocontainers.
