# Docker

Docker is a handy management tool to build, run and share container images.

These images can be uploaded and published in a centralized repository known as [Docker Hub](https://hub.docker.com), or hosted by other parties like [Quay](https://quay.io).

## Run a container

A container can be run using the following command:

```bash
docker run <container-name>
```

Try for example the following publicly available container (if you have docker installed):

```bash
docker run hello-world
```

## Pull a container

The pull command allows you to download a Docker image without running it. For example:

```bash
docker pull debian:stretch-slim
```

The above command downloads a Debian Linux image. You can check it exists by using:

```bash
docker images
```

## Run a container in interactive mode

Launching a BASH shell in the container allows you to operate in an interactive mode in the containerized operating system. For example:

```
docker run -it debian:stretch-slim bash
```

Once launched, you will notice that it is running as root (!). Use the usual commands to navigate the file system. This is useful to check if the expected programs are present within a container.

To exit from the container, stop the BASH session with the `exit` command.

## Your first Dockerfile

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

## Build the image

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

## Add a software package to the image

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

## Run Salmon in the container

Check that Salmon is running correctly in the container as shown below:

```bash
docker run my-image salmon --version
```

You can even launch a container in an interactive mode by using the following command:

```bash
docker run -it my-image bash
```

Use the `exit` command to terminate the interactive session.

## File system mounts

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
DATA=/workspace/training/nf-training/data
docker run --volume $DATA:$DATA --workdir $PWD my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Now check the content of the `transcript-index` folder by entering the command:

```bash
ls -la transcript-index
```

!!! note

    Note that the permissions for files created by the Docker execution is `root`.

## Upload the container in the Docker Hub (bonus)

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

## Run a Nextflow script using a Docker container

The simplest way to run a Nextflow script with a Docker image is using the `-with-docker` command-line option:

```bash
nextflow run script2.nf -with-docker my-image
```

As seen in the last section, you can also configure the Nextflow config file (`nextflow.config`) to select which container to use instead of having to specify it as a command-line argument every time.
