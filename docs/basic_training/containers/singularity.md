# Singularity

[Singularity](http://singularity.lbl.gov) is a container runtime designed to work in high-performance computing data centers, where the usage of Docker is generally not allowed due to security constraints.

Singularity implements a container execution model similar to Docker. However, it uses a completely different implementation design.

A Singularity container image is archived as a plain file that can be stored in a shared file system and accessed by many computing nodes managed using a batch scheduler.

!!! warning

    Singularity will not work with Gitpod. If you wish to try this section, please do it locally, or on an HPC.

## Create a Singularity images

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

Once you have saved the `Singularity` file. You can create the image with these commands:

```bash
sudo singularity build my-image.sif Singularity
```

Note: the `build` command requires `sudo` permissions. A common workaround consists of building the image on a local workstation and then deploying it in the cluster by copying the image file.

## Running a container

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

## Import a Docker image

An easier way to create a Singularity container without requiring `sudo` permission and boosting the containers interoperability is to import a Docker container image by pulling it directly from a Docker registry. For example:

```bash
singularity pull docker://debian:stretch-slim
```

The above command automatically downloads the Debian Docker image and converts it to a Singularity image in the current directory with the name `debian-jessie.simg`.

## Run a Nextflow script using a Singularity container

Nextflow allows the transparent usage of Singularity containers as easy as with Docker.

Simply enable the use of the Singularity engine in place of Docker in the Nextflow configuration file by using the `-with-singularity` command-line option:

```bash
nextflow run script7.nf -with-singularity nextflow/rnaseq-nf
```

As before, the Singularity container can also be provided in the Nextflow config file. We’ll see how to do this later.

## The Singularity Container Library

The authors of Singularity, [SyLabs](https://www.sylabs.io/) have their own repository of Singularity containers.

In the same way that we can push Docker images to Docker Hub, we can upload Singularity images to the Singularity Library.
