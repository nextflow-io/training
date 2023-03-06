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

Crie um arquivo índice de genoma utilizando Salmon no container.

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

Anteriormente, o container Singularity também pode ser disponibilizado no arquivo de confirguração do Nextflow. Nós iremos ver como funciona isso mais tarde.

### A Biblioteca de Containers Singularity

Os autores do Singularity, [SyLabs](https://www.sylabs.io/) tem o seu próprio repositório de containers Singularity.

Da mesma forma que disponibilizamos as imagens Docker no Docker Hub, nós podemos disponibilizar as imagens Singularity na Singularity Library.

## Pacotes Conda/Bioconda

Conda é um popular gerenciador de pacotes e ambientes. O suporte a Conda permite pipelines Nextflow automaticamente criarem e ativarem ambientes Conda and activate the Conda, dadas as dependências especificadas de cada processo.

Neste ambiente Gitpod, conda já está instalado.

### Usando conda

Um ambiente Conda é definido utilizando um arquivo YAML, que lista todos os pacotes de programas. A primeira coisa que você precisa fazer é inciar o conda para uma interação shell, e daí abrir um novo terminal utilizando bash.

```bash
conda init
bash
```

Com isso escrever seu arquivo YAML (para `env.yml`). Já existe um arquivo `env.yml` na pasta `nf-training` como um exemplo. O seu conteúdo é mostrado abaixo.

```yaml
--8<-- "nf-training/env.yml"
```

Dado a receita do arquivo, o ambiente é criado utilizando o comando abaixo. O comando `conda env create` deve demorar vários minutos, pois o conda tenta resolver todas dependências dos desejados pacotes durante a execução, e então baixa tudo que é requerido.

```bash
conda env create --file env.yml
```

Você pode checar se o ambiente foi criado com êxito com o comando abaixo:

```bash
conda env list
```

Isso deve parecer com algo assim:

```bash
# conda environments:
#
base                  *  /opt/conda
nf-tutorial              /opt/conda/envs/nf-tutorial
```

Para habilitar o ambiente, você pode usar o comando `activate`:

```bash
conda activate nf-tutorial
```

Nextflow consegue gerenciar a ativação de um ambiente Conda quando seu diretório é especificado utilizando a opção `-with-conda` (usando o mesmo caminho mostrado na função `list`). Por exemplo:

```bash
nextflow run script7.nf -with-conda /opt/conda/envs/nf-tutorial
```

!!! info

    Quando criar um ambiente Conda com o arquivo receita YAML, Nextflow automaticamente baixará todas dependências necessárias, montará o ambiente e ativará.

Isso deixa fácil gerenciar diferentes ambientes pra os processos no fluxo de trabalho do script.

Veja a [documentação](https://www.nextflow.io/docs/latest/conda.html) para detalhes.

### Crie e utilize ambientes conda-like utilizando micromamba

Outra forma de construir um ambiente conda-like é pelo `Dockerfile` e [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

`micromamba` é um pacote rápido e robusto para montar pequenos ambientes conda-like.

Isso previne montar um ambiente conda para cada vez que queira utiliza-lo (como delineado nas seções anteriores).

Para fazer isso, você simplesmente precisa de um `Dockerfile` e utilizar micromamba para instalar os pacotes. Porém, uma boa prática é ter o arquivo receita YAML como nas seções anteriores, então nós iremos fazer isso aqui também, utilizando o mesmo `env.yml` como antes.

```yaml
--8<-- "nf-training/env.yml"
```

Então, nós podemos escrever nosso Dockerfile com micromamba instalando os pacotes por esse arquivo receita.

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

O `Dockerfile` acima pega a imagem pai _mambaorg/micromamba_, e instala um ambiente `conda` utilizando `micromamba`, e instala `salmon`, `fastqc` e `multiqc`.

Tente executar o pipeline RNA-seq de antes (script7.nf). Comece montando seu próprio `Dockerfile`  micromamba (como mostrado acima), salve no seu repositório docker hub, e oriente o Nextflow a rodar por esse container (mudando seu `nextflow.config`).

!!! warning

    Montar um container Docker e disponibilizar no seu repósitorio pessoal pode levar &gt;10 minutos.

??? example "Para um resumo dos passos a tomar, clique aqui:"

    1. crie um arquivo chamado `Dockerfile` no diretório atual (com os código acima).

    2. Monte a imagem: `docker build -t my-image .` (não esqueça o _._).

    3. Publique a imagem docker na sua conta docker.

        Algo parecido como o seguinte, com `<myrepo>` substituido para seu próprio Docker ID, sem _&lt;_ e _&gt;_ caracteres!

        `my-image` pode ser substituido por qualquer nome que escolher. como boa prática, escolha algo memorável e certifique que o nome combine com o nome usado no comando anterior.

        ```bash
        docker login
        docker tag my-image <myrepo>/my-image
        docker push <myrepo>/my-image
        ```

    4. Adicione o arquivo imagem no arquivo `nextflow.config`.

        ex. remova o seguinte de `nextflow.config`:

        ```groovy
        process.container = 'nextflow/rnaseq-nf'
        ```

        e mude por:

        ```groovy
        process.container = '<myrepo>/my-image'
        ```

    5. Tente executar o Nextflow, ex.:

        ```bash
        nextflow run script7.nf -with-docker
        ```

    Nextflow deve conseguir achar `salmon` pra rodar o processo.

## BioContainers

Outro útilo recurso para conectar Bioconda e containers é o projeto [BioContainers](https://biocontainers.pro). BioContainers é uma inciativa da comunidade que prover um registro de imagens de container para toda receita Bioconda.

Até agora, nós vimos como instalar pacotes com conda and micromamba, ambos localmente e com container. Com BioContainers, você não precisa criar sua própria imagem container pra as ferramentas que você quer, e não precisa utilizar conda or micromamba para instalar pacotes. Ele já disponiliza uma imagem Docker contendo os programas que você quer instalado. Por exemplo, você pode adquirir a imagem container do fastqc utilizando BioContainers:

```bash
docker pull biocontainers/fastqc:v0.11.5
```

Você pode checar o registro dos pacotes que quer no [site oficial do BioContainers](https://biocontainers.pro/registry).

Contrariamente a outros registros que irão puxar a ultima imagem quando nenhuma tag (version) é disponibilizada, você precisa especificar uma tag quando baixando do BioContainers (depois de dois pontos `:`, ex. fastqc:v0.11.5). Cheque as tags com o registro e escolha a que melhor se adequa a suas necessidades.

!!! tip

    Você pode ter definições mais complexas dentro de seu bloco de processo deixando a imagem de container apropriada ou o pacote conda seria usada dependendo se o usuário seleciona singularity, Docker ou conda. Você pode clicar [aqui](https://nf-co.re/docs/contributing/modules#software-requirements) para mais informações e [aqui](https://github.com/nf-core/modules/blob/61f68913fefc20241ceccb671b104230b2d775d7/modules/bowtie2/align/main.nf#L6-L9) para um exemplo.

### :material-progress-question: Exercises

!!! exercise

    Durante o anterior tutorial RNA-Seq (script2.nf), nós criamos um índice. Dado que nós não temos salmom instalado localmente na maquina provida pelo Gitpod, nós temos que ou executar com `-with-conda` ou `-with-docker`. Sua tarefa agora é executar novamente com `-with-docker`, mas sem ter que criar sua própria imagem de container Docker. Invés disso, use a imagem BioContainers para salmon 1.7.0.


    ??? result

        ```bash
        nextflow run script2.nf -with-docker quay.io/biocontainers/salmon:1.7.0--h84f40af_0
        ```

!!! exercise "Bonus Exercise"

    Mude as diretivas do processo no `script5.nf` ou no arquivo `nextflow.config` para fazer o pipeline utilizar automaticamente BioContainers quando usando salmon, ou fastqc.

    !!! tip "Dica"

        temporariamente comente a linha `#!groovy process.container = 'nextflow/rnaseq-nf'` no arquivo `nextflow.config` para ter certeza que o processo está utilizando BioContainers que você configurou, e não a imagem de container que estavamos usando durante esse treinamento.

    ??? result

        Com essas mudanças, você deve ser capaz de executar o pipeline com BioContainers executando a seguinte linha de comando:

        ```bash
        nextflow run script5.nf
        ```

        com as sequintes diretivas de container para cada processo:

        ```groovy
        process FASTQC {
            container 'biocontainers/fastqc:v0.11.5'
            tag "FASTQC on $sample_id"
        ...
        ```

        e

        ```groovy
        process QUANTIFICATION {
            tag "Salmon on $sample_id"
            container 'quay.io/biocontainers/salmon:1.7.0--h84f40af_0'
            publishDir params.outdir, mode:'copy'
        ...
        ```

        Cheque o arquivo `.command.run` no diretório de trabalho e certifique-se que a linha de execução contém o Biocontainers correto.
