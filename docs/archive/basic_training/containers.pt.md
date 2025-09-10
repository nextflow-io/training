---
description: Material de treinamento básico do Nextflow
---

!!! warning

    Some of the translations on the training portal are out of date.
    The translated material may be incomplete or incorrect.
    We plan to update the translations later this year.
    In the meantime, please try to work through the English-language material if you can.

# Gerencie dependências e contêineres

Fluxos de trabalhos computacionais são raramente compostos de um só script ou ferramenta. Muitas vezes, eles dependem de dúzias de componentes de softwares ou bibliotecas.

Instalar e manter tais dependências é uma tarefa desafiadora e uma fonte comum de irreprodutibilidade em aplicações científicas.

Para superar esses problemas, nós utilizamos contêineres que habilitam essas dependências de softwares, isto é ferramentas e bibliotecas necessárias para uma análise de dados, para estar encapsuladas em uma ou mais imagens de contêiner Linux independentes, prontas para serem executadas e imutáveis. Elas podem facilmente ser implementadas em qualquer plataforma que suporta o motor de conteinerização.

Contêineres podem ser executados de uma forma isolada pelo sistema do hospedeiro. Elas tem sua própria cópia do sistema de arquivos, espaço de processamento, gerenciamento de memória, etc.

!!! info

    Contêineres foram introduzidos no kernel 2.6 como um recurso do Linux conhecido como _Control Groups_ ou [Cgroups](https://en.wikipedia.org/wiki/Cgroups).

## Docker

Docker é uma ferramenta útil para o gerenciamento, execução e compartilhamento de imagens de contêineres.

Essas imagens can podem ser carregadas e publicadas em um repositório centralizado conhecido como [Docker Hub](https://hub.docker.com), ou hospedadas por terceiros como o [Quay](https://quay.io).

### Execute um contêiner

Um contêiner pode ser executado com o seguinte comando:

```bash
docker run <nome-do-contêiner>
```

Tente executar o seguinte contêiner público (se você tiver o Docker instalado), por exemplo:

```bash
docker run hello-world
```

### Baixe um contêiner

O comando pull possibilita que você baixe uma imagem Docker sem que a execute. Por exemplo:

```bash
docker pull debian:bullseye-slim
```

O comando acima baixa uma imagem Debian Linux. Você pode checar se ela existe usando:

```bash
docker images
```

### Execute um contêiner em mode interativo

Iniciar uma shell BASH em um contêiner permite que você opere em modo interativo no sistema operacional conteinerizado. Por exemplo:

```
docker run -it debian:bullseye-slim bash
```

Uma vez iniciado, você vai notar que está como root (!). Use os comandos usuais para navegar pelo sistema de arquivos. Isso é útil para checar se os programas necessários estão presentes no contêiner.

Para sair do contêiner, pare a sessão BASH com o comando `exit`.

### Seu primeiro Dockerfile

Imagens Docker são criadas utilizando um arquivo chamado `Dockerfile`, que é um simples arquivo de texto contendo uma lista de comandos para montar e configurar uma imagem com os pacotes de programas necessários.

Aqui, você criará uma imagem Docker contendo o cowsay e a ferramenta Salmon

!!! warning

    O processo de montagem do Docker automaticamente copia todos os arquivos que estão no diretório atual para o Docker daemon para que ele possa criar a imagem. Isso pode custar muito tempo quando existem vários ou grandes arquivos. Por essa razão, é importante que _sempre_ se trabalhe em um diretório contendo apenas os arquivos que você realmente precisa incluir em sua imagem Docker. Alternativamente, você pode usar o arquivo `.dockerignore` para selecionar os aquivos que serão excluídos da montagem.

Use seu editor favorito (ex.: `vim` ou `nano`) para criar um arquivo chamado `Dockerfile` e copiar o seguinte conteúdo:

```dockerfile
FROM debian:bullseye-slim

LABEL image.author.name "Seu nome aqui"
LABEL image.author.email "seu@email.aqui"

RUN apt-get update && apt-get install -y curl cowsay

ENV PATH=$PATH:/usr/games/
```

### Monte a imagem

Monte a imagem do Docker com base no Dockerfile utilizando o seguinte comando:

```bash
docker build -t minha-imagem .
```

Onde "minha-imagem" é o nome que o usuário especificou para a imagem que será criada.

!!! tip

    Não esqueça do ponto no comando acima.

Quando completo, verifique se a imagem foi criada listando todas imagens disponíveis:

```bash
docker images
```

Você pode testar seu novo contêiner executando esse comando:

```bash
docker run minha-imagem cowsay Olá Docker!
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
docker build -t minha-imagem .
```

Você perceberá que isso cria uma nova imagem Docker **mas** com um ID de imagem diferente.

### Execute Salmon no contêiner

Cheque se o Salmon está executando corretamente no contêiner como mostrado abaixo:

```bash
docker run minha-imagem salmon --version
```

Você pode até iniciar o contêiner no modo interativo utilizando o seguinte comando:

```bash
docker run -it minha-imagem bash
```

Use o comando `exit` para finalizar a sessão interativa.

### Montagem do sistema de arquivos

Crie um arquivo índice de genoma utilizando o Salmon no contêiner.

Tente executar o Salmon no contêiner com o seguinte comando:

```bash
docker run minha-imagem \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

O comando acima falha porque o Salmon não pode acessar o arquivo de entrada.

Isso acontece porque o contêiner executa em um sistema de arquivos totalmente diferente e não pode acessar o arquivo no sistema de arquivo do hospedeiro por padrão.

Você precisará usar a opção de linha de comando `--volume` para montar o(s) arquivo(s) de entrada, por exemplo

```bash
docker run --volume $PWD/data/ggal/transcriptome.fa:/transcriptome.fa minha-imagem \
    salmon index -t /transcriptome.fa -i transcript-index
```

!!! warning

    O diretório `transcript-index` gerado ainda está inacessível no sistema de arquivo do sistema operacional hospedeiro.

Um jeito mais fácil é montar o diretório original em um idêntico no contêiner. Isso permite que você utilize o mesmo caminho durante a execução dentro do contêiner, por exemplo

```bash
docker run --volume $PWD:$PWD --workdir $PWD minha-imagem \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Ou definir uma pasta que você queira montar como uma variável de ambiente, chamada `DATA`:

```bash
DATA=/workspaces/training/nf-training/data
docker run --volume $DATA:$DATA --workdir $PWD minha-imagem \
    salmon index -t $PWD/data/ggal/transcriptome.fa -i transcript-index
```

Agora cheque o conteúdo da pasta `transcript-index` utilizando o comando:

```bash
ls -la transcript-index
```

!!! note

    Note que as permissões para criação dos arquivos utilizado pelo Docker são `root`.

### Disponibilize o contêiner no Docker Hub (bônus)

Publique seu contêiner no Docker Hub para compartilhá-lo com outras pessoas.

Crie uma conta no site <https://hub.docker.com>. Então no seu terminal shell execute o seguinte comando, utilizando seu usuário e senha que criou quando se registrou no Hub:

```bash
docker login
```

Renomeie a imagem para incluir seu nome de usuário Docker:

```bash
docker tag minha-imagem <nome-de-usuario>/minha-imagem
```

Finalmente mande para o Docker Hub:

```bash
docker push <nome-de-usuario>/minha-imagem
```

Depois disso, qualquer um conseguirá baixar a imagem utilizando o comando:

```bash
docker pull <nome-de-usuario>/minha-imagem
```

Note que depois de uma operação push e pull, o Docker imprime na tela o número de registro do contêiner, por exemplo:

```console
Digest: sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266
Status: Downloaded newer image for nextflow/rnaseq-nf:latest
```

Isso é um identificador imutável e único que pode ser usado para referenciar a imagem de contêiner de uma forma única. Por exemplo:

```bash
docker pull nextflow/rnaseq-nf@sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266
```

### Execute um script do Nextflow utilizando um contêiner Docker

A maneira mais simples de rodar um script Nextflow é usando a opção de linha de comando `-with-docker`:

```bash
nextflow run script2.nf -with-docker minha-imagem
```

Como visto na última seção, você também pode configurar o arquivo config (`nextflow.config`) para selecionar qual contêiner utilizar invés de ter que especificar como um argumento de linha de comando toda vez.

## Singularity

[Singularity](http://singularity.lbl.gov) é um motor de conteinerização desenvolvido para trabalhar com computação de alta performance em centro de dados, onde geralmente o Docker não é permitido por motivos de restrições de segurança.

O Singularity implementa um modelo de execução de contêiner similar ao Docker. Entretanto, ele usa um design de implementação completamente diferente.

Uma imagem de contêiner do Singularity é arquivada como um arquivo de texto simples que pode ser armazenado em um sistema de arquivo compartilhado e acessado por muitos nós computacionais gerenciados usando um escalonador de lote.

!!! warning

    O Singularity não irá funcionar com Gitpod. Se você quer testar essa seção, por favor faça localmente, ou em um cluster.

### Crie imagens do Singularity

Imagens do Singularity são criadas utilizando um arquivo `Singularity` de uma forma similar ao Docker mas utilizando uma sintaxe diferente.

```singularity
Bootstrap: docker
From: debian:bullseye-slim

%environment
export PATH=$PATH:/usr/games/

%labels
AUTHOR <seu nome>

%post

apt-get update && apt-get install -y locales-all curl cowsay
curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.0.0/salmon-1.0.0_linux_x86_64.tar.gz | tar xz \
&& mv /salmon-*/bin/* /usr/bin/ \
&& mv /salmon-*/lib/* /usr/lib/
```

Uma vez que você salvou o arquivo `Singularity`, você pode criar uma imagem utilizando esses comandos:

```bash
sudo singularity build minha-imagem.sif Singularity
```

Nota: O comando `build` requer permissões `sudo`. Uma forma de contornar isso consiste em criar a imagem em uma estação de trabalho local e então implementar no cluster copiando o arquivo imagem.

### Executando um contêiner

Quando terminar, você pode executar o contêiner com o seguinte comando:

```bash
singularity exec minha-imagem.sif cowsay 'Olá Singularity'
```

Utilizando o comando `shell` você pode entrar no contêiner utilizando o modo interativo. Por exemplo:

```bash
singularity shell minha-imagem.sif
```

Uma vez dentro da instância do contêiner execute o comando:

```bash
touch ola.txt
ls -la
```

!!! info

    Note como os arquivos do sistema de arquivos hospedeiro são mostrados. O Singularity automaticamente monta o diretório do hospedeiro `$HOME` e usa como diretório de trabalho.

### Importe uma imagem do Docker

Uma forma mais fácil de criar um contêiner com o Singularity não necessitando da permissão `sudo` e melhorando a interoperabilidade dos contêineres é importando uma imagem de contêiner do Docker puxando diretamente do repositório de imagens do Docker. Por exemplo:

```bash
singularity pull docker://debian:bullseye-slim
```

O comando acima automaticamente baixa uma imagem Docker do Debian e converte para uma imagem Singularity no diretório atual com o nome `debian-jessie.simg`.

### Execute um script do Nextflow utilizando um contêiner Singularity

O Nextflow permite o uso de contêineres Singularity de forma tão fácil e transparente quanto com o Docker.

Simplesmente ative o uso do motor Singularity no lugar do Docker na linha de comando do Nextflow utilizando a opção de linha de comando `-with-singularity`:

```bash
nextflow run script7.nf -with-singularity nextflow/rnaseq-nf
```

Como antes, o contêiner Singularity também pode ser disponibilizado no arquivo de configuração do Nextflow. Nós iremos ver como funciona isso mais tarde.

### A Biblioteca de Contêineres Singularity

Os autores do Singularity, [SyLabs](https://www.sylabs.io/) tem o seu próprio repositório de contêineres Singularity.

Da mesma forma que disponibilizamos as imagens Docker no Docker Hub, nós podemos disponibilizar as imagens Singularity na Singularity Library.

## Pacotes Conda/Bioconda

O Conda é um popular gerenciador de pacotes e ambientes. O suporte a Conda permite que fluxos de trabalho Nextflow automaticamente criem e ativem ambientes Conda, dadas as dependências especificadas de cada processo.

Neste ambiente Gitpod, o conda já está instalado.

### Usando conda

Um ambiente Conda é definido utilizando um arquivo YAML, que lista todos os pacotes de programas. A primeira coisa que você precisa fazer é iniciar o conda para uma interação shell, e daí abrir um novo terminal utilizando bash.

```bash
conda init
bash
```

Com isso, escreva seu arquivo YAML (`env.yml`). Já existe um arquivo `env.yml` na pasta `nf-training` como um exemplo. O seu conteúdo é mostrado abaixo.

```yaml
--8<-- "nf-training/env.yml"
```

Dado o arquivo de receita, o ambiente é criado utilizando o comando abaixo. O comando `conda env create` deve demorar vários minutos, pois o conda tenta resolver todas as dependências dos pacotes desejados durante a execução, e então baixa tudo que é requerido.

```bash
conda env create --file env.yml
```

Você pode checar se o ambiente foi criado com êxito com o comando abaixo:

```bash
conda env list
```

Você deve ver algo similar ao mostrado abaixo:

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

O Nextflow consegue gerenciar a ativação de um ambiente Conda quando seu diretório é especificado utilizando a opção `-with-conda` (usando o mesmo caminho mostrado com a função `list`). Por exemplo:

```bash
nextflow run script7.nf -with-conda /opt/conda/envs/nf-tutorial
```

!!! info

    Quando criar um ambiente Conda com o arquivo de receita YAML, o Nextflow automaticamente baixará todas dependências necessárias, montará o ambiente e o ativará.

Isso torna fácil gerenciar diferentes ambientes para os processos no fluxo de trabalho do script.

Veja a [documentação](https://www.nextflow.io/docs/latest/conda.html) para mais detalhes.

### Crie e utilize ambientes no estilo do conda com o micromamba

Outra forma de construir um ambiente no estilo do conda é com o `Dockerfile` e o [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

O `micromamba` é um pacote rápido e robusto para montar pequenos ambientes baseados no conda.

Isso economiza tempo ao montar um ambiente conda toda vez que quiser utilizá-lo (como delineado nas seções anteriores).

Para fazer isso, você simplesmente precisa de um `Dockerfile` e utilizar o micromamba para instalar os pacotes. Porém, uma boa prática é ter o arquivo de receita YAML como nas seções anteriores, então nós iremos fazer isso aqui também, utilizando o mesmo `env.yml` utilizado anteriormente.

```yaml
--8<-- "nf-training/env.yml"
```

Então, nós podemos escrever nosso Dockerfile com o micromamba instalando os pacotes por esse arquivo de receita.

```dockerfile
FROM mambaorg/micromamba:0.25.1

LABEL image.author.name "Seu Nome Aqui"
LABEL image.author.email "seu@email.aqui"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba create -n nf-tutorial

RUN micromamba install -y -n nf-tutorial -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/nf-tutorial/bin:$PATH
```

O `Dockerfile` acima pega a imagem pai _mambaorg/micromamba_, e instala um ambiente `conda` utilizando `micromamba`, e então instala o `salmon`, o `fastqc` e o `multiqc`.

Tente executar o fluxo de trabalho RNA-seq visto anteriormente (script7.nf). Comece montando seu próprio `Dockerfile` micromamba (como mostrado acima), salve no seu repositório no Docker Hub, e oriente o Nextflow a rodar por esse contêiner (mudando seu `nextflow.config`).

!!! warning

    Montar um contêiner Docker e disponibilizar no seu repositório pessoal pode levar &gt;10 minutos.

??? example "Para um resumo dos passos a tomar, clique aqui:"

    1. Crie um arquivo chamado `Dockerfile` no diretório atual (com os código acima).

    2. Monte a imagem: `docker build -t minha-imagem .` (não esqueça o _._).

    3. Publique a imagem Docker na sua conta do Docker Hub.

        Algo parecido como o seguinte, com `<meurepo>` substituído para seu próprio ID do Docker, sem _&lt;_ e _&gt;_ caracteres!

        `minha-imagem` pode ser substituído por qualquer nome que você escolher. Como boa prática, escolha algo memorável e certifique-se de que o nome combine com o nome usado no comando anterior.

        ```bash
        docker login
        docker tag minha-imagem <meurepo>/minha-imagem
        docker push <meurepo>/minha-imagem
        ```

    4. Adicione a imagem do contêiner no arquivo `nextflow.config`.

        ex. remova o seguinte do arquivo `nextflow.config`:

        ```groovy
        process.container = 'nextflow/rnaseq-nf'
        ```

        e adicione:

        ```groovy
        process.container = '<meurepo>/minha-imagem'
        ```

    5. Tente executar o Nextflow, por exemplo:

        ```bash
        nextflow run script7.nf -with-docker
        ```

    Agora, o Nextflow deve conseguir achar `salmon` pra rodar o processo.

## BioContainers

Outro recurso útil para conectar Bioconda e contêineres é o projeto [BioContainers](https://biocontainers.pro). BioContainers é uma iniciativa da comunidade para prover um repositório de imagens de contêiner para cada receita do Bioconda.

Até agora, nós vimos como instalar pacotes com conda e micromamba, ambos localmente e com contêiner. Com o BioContainers, você não precisa criar sua própria imagem de contêiner para as ferramentas que você quiser, e não precisa utilizar conda ou micromamba para instalar pacotes. O BioContainers já disponibiliza uma imagem Docker contendo os programas que você quer instalado. Por exemplo, você pode adquirir a imagem de contêiner do fastqc utilizando BioContainers:

```bash
docker pull biocontainers/fastqc:v0.11.5
```

Você pode checar o repositório dos pacotes que quer no [site oficial do BioContainers](https://biocontainers.pro/registry). Para encontrar imagens de container com várias ferramentas, confira a página [Multi-package images](https://biocontainers.pro/multipackage).

Diferente de outros repositórios que irão puxar a imagem mais recente quando nenhum rótulo (versão) é especificado, você precisa especificar um rótulo quando for baixar do BioContainers (depois de dois pontos `:`, por exemplo fastqc:v0.11.5). Cheque os rótulos com o registro e escolha o que melhor se adéqua a suas necessidades.

Você também pode instalar o pacote `galaxy-util-tools` e procurar por imagens de container _mulled_ através da linha de comando. Veja as instruções abaixo, usando o `conda` para instalar o pacote.

```bash
conda activate um-ambiente-conda-que-voce-ja-criou
conda install galaxy-tool-util
mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep mulled
```

!!! tip

    Você pode ter definições mais complexas dentro de seu bloco de processo deixando a imagem de contêiner apropriada ou o pacote conda para serem usadas dependendo se o usuário selecionou singularity, Docker ou conda. Você pode clicar [aqui](https://nf-co.re/docs/contributing/modules#software-requirements) para mais informações e [aqui](https://github.com/nf-core/modules/blob/61f68913fefc20241ceccb671b104230b2d775d7/modules/bowtie2/align/main.nf#L6-L9) para um exemplo.

### :material-progress-question: Exercises

!!! exercise

    Durante a seção onde construímos o fluxo de trabalho de RNA-Seq, nós criamos um índice (script2.nf). Dado que nós não temos Salmon instalado localmente na máquina provida pelo Gitpod, nós temos que ou executar com `-with-conda` ou `-with-docker`. Sua tarefa agora é executar novamente com `-with-docker`, mas sem ter que criar sua própria imagem de contêiner Docker. Invés disso, use a imagem do BioContainers para Salmon 1.7.0.


    ??? Solution

        ```bash
        nextflow run script2.nf -with-docker quay.io/biocontainers/salmon:1.7.0--h84f40af_0
        ```

!!! exercise "Bonus Exercise"

    Mude as diretivas do processo no `script5.nf` ou no arquivo `nextflow.config` para fazer o fluxo de trabalho utilizar automaticamente BioContainers quando usando salmon, ou fastqc.

    !!! tip "Dica"

        temporariamente comente a linha `#!groovy process.container = 'nextflow/rnaseq-nf'` no arquivo `nextflow.config` para ter certeza que o processo está utilizando o contêiner do BioContainers que você configurou, e não a imagem de contêiner que estávamos usando durante esse treinamento.

    ??? Solution

        Com essas mudanças, você deve ser capaz de executar o fluxo de trabalho com BioContainers executando a seguinte linha de comando:

        ```bash
        nextflow run script5.nf
        ```

        com as seguintes diretivas de contêiner para cada processo:

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
            publishDir params.outdir, mode: 'copy'
        ...
        ```

        Cheque o arquivo `.command.run` no diretório de trabalho e certifique-se que a linha de execução contém o Biocontainers correto.
