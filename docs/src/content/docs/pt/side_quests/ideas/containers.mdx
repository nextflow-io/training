# Parte 1: Mais Contêineres

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Como encontrar ou criar imagens de contêiner

Alguns desenvolvedores de software fornecem imagens de contêiner para seus softwares que estão disponíveis em registros de contêineres como o Docker Hub, mas muitos não o fazem.
Nesta seção opcional, mostraremos duas maneiras de obter uma imagem de contêiner para ferramentas que você deseja usar em seus pipelines Nextflow: usando o Seqera Containers e construindo a imagem de contêiner você mesmo.

Você obterá/construirá uma imagem de contêiner para o pacote pip `quote`, que será usado no exercício no final desta seção.

### 1.1. Obter uma imagem de contêiner do Seqera Containers

Seqera Containers é um serviço gratuito que constrói imagens de contêiner para ferramentas instaláveis via pip e conda (incluindo bioconda).
Navegue até [Seqera Containers](https://www.seqera.io/containers/) e procure pelo pacote pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Clique em "+Add" e depois em "Get Container" para solicitar uma imagem de contêiner para o pacote pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Se esta for a primeira vez que um contêiner da comunidade é construído para esta versão do pacote, pode levar alguns minutos para ser concluído.
Clique para copiar o URI (por exemplo, `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) da imagem de contêiner que foi criada para você.

Agora você pode usar a imagem de contêiner para executar o comando `quote` e obter uma frase aleatória de Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Saída:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Construir a imagem de contêiner você mesmo

Vamos usar alguns detalhes de construção do site Seqera Containers para construir a imagem de contêiner para o pacote pip `quote` nós mesmos.
Retorne ao site Seqera Containers e clique no botão "Build Details".

O primeiro item que veremos é o `Dockerfile`, um tipo de arquivo de script que contém todos os comandos necessários para construir a imagem de contêiner.
Adicionamos alguns comentários explicativos ao Dockerfile abaixo para ajudá-lo a entender o que cada parte faz.

```Dockerfile title="Dockerfile"
# Inicia a partir da imagem base docker do micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Copia o arquivo conda.yml para dentro do contêiner
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Instala vários utilitários para o Nextflow usar e os pacotes no arquivo conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Executa o contêiner como usuário root
USER root
# Define a variável de ambiente PATH para incluir o diretório de instalação do micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

O segundo item que veremos é o arquivo `conda.yml`, que contém a lista de pacotes que precisam ser instalados na imagem de contêiner.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Copie o conteúdo desses arquivos para os esboços localizados no diretório `containers/build`, depois execute o seguinte comando para construir a imagem de contêiner você mesmo.

!!! Note "Nota"

    Usamos a flag `-t quote:latest` para marcar a imagem de contêiner com o nome `quote` e a tag `latest`.
    Poderemos usar esta tag para nos referir à imagem de contêiner ao executá-la neste sistema.

```bash
docker build -t quote:latest containers/build
```

Depois que terminar de construir, você pode executar a imagem de contêiner que acabou de construir.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Conclusão

Você aprendeu duas maneiras diferentes de obter uma imagem de contêiner para uma ferramenta que deseja usar em seus pipelines Nextflow: usando o Seqera Containers e construindo a imagem de contêiner você mesmo.

### O que vem a seguir?

Você tem tudo o que precisa para continuar para o [próximo capítulo](./04_hello_genomics.md) desta série de treinamento.
Você também pode continuar com um exercício opcional para buscar citações de pioneiros da computação/biologia usando o contêiner `quote` e exibi-las usando o contêiner `cowsay`.

---

## 2. Faça a vaca citar cientistas famosos

Esta seção contém alguns exercícios de desafio, para praticar o que você aprendeu até agora.
Fazer esses exercícios _não é obrigatório_ para entender as partes posteriores do treinamento, mas fornecem uma maneira divertida de reforçar seus aprendizados descobrindo como fazer a vaca citar cientistas famosos.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. Modificar o script `hello-containers.nf` para usar um processo getQuote

Temos uma lista de pioneiros da computação e biologia no arquivo `containers/data/pioneers.csv`.
Em alto nível, para completar este exercício você precisará:

- Modificar o `params.input_file` padrão para apontar para o arquivo `pioneers.csv`.
- Criar um processo `getQuote` que usa o contêiner `quote` para buscar uma citação para cada entrada.
- Conectar a saída do processo `getQuote` ao processo `cowsay` para exibir a citação.

Para a imagem de contêiner `quote`, você pode usar tanto a que você construiu no exercício de desafio anterior quanto a que você obteve do Seqera Containers.

!!! Hint "Dica"

    Uma boa escolha para o bloco `script` do seu processo getQuote seria:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Você pode encontrar uma solução para este exercício em `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modificar seu pipeline Nextflow para permitir que ele execute nos modos `quote` e `sayHello`.

Adicione alguma lógica de ramificação ao seu pipeline para permitir que ele aceite entradas destinadas tanto para `quote` quanto para `sayHello`.
Aqui está um exemplo de como usar uma instrução `if` em um fluxo de trabalho Nextflow:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "Dica"

    Você pode usar `new_ch = processName.out` para atribuir um nome ao canal de saída de um processo.

Você pode encontrar uma solução para este exercício em `containers/solutions/hello-containers-4.2.nf`.

### Conclusão

Você sabe como usar contêineres no Nextflow para executar processos, e como construir alguma lógica de ramificação em seus pipelines!

### O que vem a seguir?

Comemore, faça uma pausa para se alongar e beba um pouco de água!

Quando estiver pronto, passe para a Parte 3 desta série de treinamento para aprender como aplicar o que você aprendeu até agora a um caso de uso de análise de dados mais realista.
