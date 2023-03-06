---
description: Material de treinamento básico do Nextflow
---

# Gerencie dependencias e containers

Fluxos de trabalhos computacionais são raramente compostos de um só script ou ferramenta.  Muitas vezes, eles dependem de duzias de componentes de softwares ou bibliotecas.

Instalar e manter tais dependências é uma tarefa desafiadora e uma fonte comum de ireprodutibilidade em aplicações científicas.

Para superar esses problemas, nós utilizamos containers que habilitam essas dependências de softwares, isto é ferramentas e bibliotecas necessárias para uma análise de dados, para estar encapsuladas em um ou mais independente, pronto para executar, imutável imagens de containers Linux, que facilmente podem ser implementado em qualquer plataforma que suporta o motor de conteinerização.

Containers podem ser executados de uma forma isolada pelo sistema do hospedeiro. Tendo sua própria cópia do sistema de arquivos, espaço de processamento, gerencimento de memória, etc.

!!! info

    Containers forma introduzidos com o kernel 2.6 como um recurso do Linux conhecido como _Control Groups_ or [Cgroups](https://en.wikipedia.org/wiki/Cgroups).

## Docker

Docker é uma ferramenta útil para o gerenciamento, execução e compartilhamento de imagens de containers.

Essa imagens can podem ser carregadas e publicadas em um repositório centralizado conhecido como [Docker Hub](https://hub.docker.com), ou hospedada por outros grupos como [Quay](https://quay.io).

### Exectute um container

Um container pode ser executado com o seguinte comando:

```bash
docker run <container-name>
```

Tente, por exemplo, o seguinte container público (se você tem docker instalado):

```bash
docker run hello-world
```

### Baixar um container

O comando pull possibilita que você baixe uma imagem Docker sem que execute-a. Por exemplo:

```bash
docker pull debian:stretch-slim
```

O comando acima baixa uma imagem Debian Linux. Você pode checar se ela existe usando:

```bash
docker images
```

### Executar um container em mode interativo

Inciar um BASH shell em um container permite que você opere em modo interativo no sistema operacional containerizado. Por exemplo:

```
docker run -it debian:stretch-slim bash
```

Uma vez iniciado, você vai notar que está como root (!). Use os comandos usuais para navegar pelo sistema de arquivos. Isso é útil para checar se os programas necessários estão presentes no container.

Para sair do container, para a seção BASH com o comando `exit`.

### Seu primeiro Dockerfile

Imagens docker são criadas utilizando um arquivo chamado `Dockerfile`, que é um simples arquivo de texto contendo uma lista de comandos para montar e configurar uma imagem com os pacotes de programas necessários.

Aqui, você criará uma imagem Docker contendo cowsay e Salmon tool

!!! warning

    O processo de montagem do Docker automaticamente copia todos arquivos que estão no diretório atual para o Docker daemon para que possa criar a imagem. Isso pode custar muito tempo quando existem vários ou grandes arquivos. Por essa razão, é importante que _sempre_ trabalhe em um diretório contendo apenas os arquviso que você realmente precisa incluir em sua imagem Docker. Alternativamente, você pode usar o arquivo `.dockerignore` para selecionar os aquivos que serão excluidos da montagem.

Use seu editor favorito (ex.: `vim` ou `nano`) para criar um arquivo chamado `Dockerfile` e copiar o seguinte conteúdo:

```dockerfile
FROM debian:stretch-slim

LABEL image.author.name "Your Name Here"
LABEL image.author.email "your@email.here"

RUN apt-get update && apt-get install -y curl cowsay

ENV PATH=$PATH:/usr/games/
```

### Monte a imagem

Monte a imagem Dockerfile utilizando o seguinte comando:

```bash
docker build -t my-image .
```

Onde "my-image" é o nome que o usuário especificou para Dockerfile, presente no diretório atual.

!!! tip

    Não esqueça do ponto no comando acima.

Quando completo, verifique se a imagem foi criada listando todas imagens disponíveis:

```bash
docker images
```

Você pode testar seu novo container executando esse comando:

```bash
docker run my-image cowsay Hello Docker!
```

### Adicione um pacote de programa a imagem

Adicione o pacote Salmon para a imagem Docker adicionando o seguinte trecho para o `Dockerfile`:

```dockerfile
RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz | tar xz \
&& mv /salmon-*/bin/* /usr/bin/ \
&& mv /salmon-*/lib/* /usr/lib/
```

Salve o arquivo e monte a imagem novamente utilizando o mesmo comando anterior:

```bash
docker build -t my-image .
```

Você perceberá que isso cria uma nova imagem Docker **mas** com um ID de imagem diferente.

### Execute Salmon no container

Cheque se Salmon está executando corretamente no container como mostrado abaixo:

```bash
docker run my-image salmon --version
```

Você pode até iniciar o container no modo interativo utilizando o seguinte comando:

```bash
docker run -it my-image bash
```

Use o comando `exit` para finalizar a seção interativa.

### Montagem do sitema de arquivos

Crie um arquivo index de genoma utilizando Salmon no container.

Tente executar Salmon no container com o seguinte comando:

```bash
docker run my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

O comando acima falha porque Salmon não pode acessar o arquivo de entrada.

Isso acontece porque o container executa em um sistema de arquivos totalmente diferente e não pode acessar o arquivo no sistema de arquivo do host por default.

Você precisará usar a opção de linha de comando `--volume` para montar o(s) arquivo(s) de entrada por exemplo.

```bash
docker run --volume $PWD/data/ggal/transcriptome.fa:/transcriptome.fa my-image \
    salmon index -t /transcriptome.fa -i transcript-index
```

!!! warning

    O `transcript-index` diretório gerado ainda está inacessível no sistema de arquivo do host.

Um jeito mais fácil é montar o diretório original em um indêntico no container, isso permite que você utilize o mesmo caminho durante a execução dentro do container por exemplo.

```bash
docker run --volume $PWD:$PWD --workdir $PWD my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Ou definir uma pasta que você queira montar como uma variável de ambiente, chamada `DATA`:

```bash
DATA=/workspace/gitpod/nf-training/data
docker run --volume $DATA:$DATA --workdir $PWD my-image \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Agora cheque o conteúdo da pasta `transcript-index` utilizando o comando:

```bash
ls -la transcript-index
```

!!! note

    Note que as permissões para criação dos arquivos utilizado pelo Docker são `root`.

### Disponibilize o container no Docker Hub (bonus)

Publique seu container no Docker Hub para compartilha-lo com outras pessoas.

Crie uma conta no site <https://hub.docker.com>. Entao no seu terminal shell execute o seguinte comando, utilizando seu usuário e senha que criou quando se registrou no Hub:

```bash
docker login
```

Marque a imagem com seu nome de usuário Docker:

```bash
docker tag my-image <user-name>/my-image
```

Finalmente mande para o Docker Hub:

```bash
docker push <user-name>/my-image
```

Depois qualquer um conseguira baixar a imagem utilizando o comando:

```bash
docker pull <user-name>/my-image
```

Note que depois de uma operação push e pull, o Docker printa o número de registro do container, por exemplo:

```console
Digest: sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266
Status: Downloaded newer image for nextflow/rnaseq-nf:latest
```

Isso é um identificador imutável e único que pode ser usado para referênciar a imagem de container de uma forma única. Por exemplo:

```bash
docker pull nextflow/rnaseq-nf@sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266
```

### Execute um script do Nextflow utilizando um container Docker

A maneira mais simples de rodar um script Nextflow é usando a opção de linha de comando `-with-docker`:

```bash
nextflow run script2.nf -with-docker my-image
```

Como visto na última seção, você também pode configurar o arquivo config (`nextflow.config`) para selecionar qual container utilizar invés de ter que especificar como um argumento de linha de comando toda vez.

## Singularity

[Singularity](http://singularity.lbl.gov) é um motor de conteinerização desenvolvido para trabalhar com computação de alta performace em centro de dados, onde geralmente o Docker não é permitido por motivos de restrições de segurança.

Singularity implementa um modelo de execução de container similar ao Docker. Entretanto, ele usa um design de implementação completamente diferente.

Uma imagem container do Singularity é arquivada como um arquivo plain file que pode ser armazenado em um sistema de arquivo compartilhado e acessado por muitos nós computacionais gerenciados usando um escalonador de lote.

!!! warning

    Singularity não irá funcionar com Gitpod. Se você quer testar essa seção, por favor faça localmente, ou em um HPC.

### Crie imagens do Singularity

Imagens Singularity são criadas utilizando um arquivo `Singularity` de uma forma similar ao Docker mas utilizando uma sintase diferente.

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

Uma vez que você salvou o arquivo `Singularity`. Você pode criar uma imagem utilizando esses comandos:

```bash
sudo singularity build my-image.sif Singularity
```

Note: o comando `build` requer permissões `sudo`. Uma forma de contorna isso consiste em criar a imagem em uma estação de trabalho local e então implementar no cluster copiando o arquivo imagem.

### Executando um container

Quando terminar, você pode executar o container com o seguinte comando:

```bash
singularity exec my-image.sif cowsay 'Hello Singularity'
```

Utilizando o comando `shell` você pode entrar no container utilizando o modo interativo. Por exemplo:

```bash
singularity shell my-image.sif
```

Uma vez dentro da instância do container execute o comando:

```bash
touch hello.txt
ls -la
```

!!! info

    Note como os arquivos do hospedeiro são mostrados. Singularity automaticamente monta o diretório do hospedeiro `$HOME` e usa como diretório de trabalho.

### Importe uma imagem do Docker

Uma forma mais fácil de criar um container Singularity não necessitando da permissão `sudo` permission e melhorando a interoperabilidade dos container é não importando um imagem de container Docker puxando diretamente do registro Docker. Por exemplo:

```bash
singularity pull docker://debian:stretch-slim
```

O comando acima automaticamente baixa uma imagem Debian Docker e converte para uma imagem Singularity no diretório atual com o nome`debian-jessie.simg`.

### Execute um script do Nextflow utilizando um container Singularity

Nextflow permite um uso transparente de containers Singularity facilmente como com o Docker.

Simplesmente ative o uso do motor Singularity no lugar do Docker no arquivo de configurações do Nextflow utilizando a opção de linha de comando `-with-singularity`:

```bash
nextflow run script7.nf -with-singularity nextflow/rnaseq-nf
```

Com anteriomente, o container Singularity também pode ser disponibilizado no arquivo de confirguração do Nextflow. Nós iremos ver como funciona isso mais tarde.

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
