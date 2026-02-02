# Começando

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja [a playlist completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) no canal do Nextflow no YouTube.

:green_book: A transcrição do vídeo está disponível [aqui](./transcripts/00_orientation.md).
///

!!! tip "Dica"

    Os vídeos do YouTube têm alguns superpoderes!

    - :fontawesome-solid-closed-captioning: Legendas de alta qualidade (manualmente curadas). Ative-as com o ícone :material-subtitles:
    - :material-bookmark: Capítulos de vídeo na linha do tempo que correspondem aos títulos das páginas.

-->

## Inicie um ambiente de treinamento

Para usar o ambiente pré-construído que fornecemos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, consulte [Opções de ambiente](../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use clique com botão direito, ctrl-clique ou cmd-clique dependendo do seu equipamento) para que você possa continuar lendo enquanto o ambiente carrega.
Você precisará manter essas instruções abertas em paralelo para acompanhar o curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Noções básicas do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para trabalhar no curso de treinamento, então você não precisa instalar nada por conta própria.

O codespace é configurado com uma interface VSCode, que inclui um explorador de sistema de arquivos, um editor de código e um terminal shell.
Todas as instruções dadas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') se referem a essas três partes da interface do VSCode, salvo indicação em contrário.

Se você está trabalhando neste curso por conta própria, por favor, familiarize-se com as [noções básicas do ambiente](../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este treinamento é projetado para Nextflow 25.10.2 ou posterior **com o analisador de sintaxe v2 ATIVADO**.
Se você está usando um ambiente local ou personalizado, certifique-se de estar usando as configurações corretas conforme documentado [aqui](../info/nxf_versions.md).

## Prepare-se para trabalhar

Assim que seu codespace estiver em execução, há duas coisas que você precisa fazer antes de mergulhar no treinamento: definir seu diretório de trabalho para este curso específico e dar uma olhada nos materiais fornecidos.

### Defina o diretório de trabalho

Por padrão, o codespace abre com o diretório de trabalho definido na raiz de todos os cursos de treinamento, mas para este curso, trabalharemos no diretório `hello-nextflow/`.

Mude de diretório agora executando este comando no terminal:

```bash
cd hello-nextflow/
```

Você pode configurar o VSCode para focar neste diretório, de modo que apenas os arquivos relevantes apareçam na barra lateral do explorador de arquivos:

```bash
code .
```

!!! tip "Dica"

    Se por algum motivo você sair deste diretório (por exemplo, seu codespace adormecer), você pode sempre usar o caminho completo para retornar a ele, assumindo que você está executando isso dentro do ambiente de treinamento do GitHub Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Agora vamos dar uma olhada no conteúdo.

### Explore os materiais fornecidos

Você pode explorar o conteúdo deste diretório usando o explorador de arquivos no lado esquerdo do espaço de trabalho de treinamento.
Alternativamente, você pode usar o comando `tree`.

Ao longo do curso, usamos a saída do `tree` para representar a estrutura e o conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice até o segundo nível:

```bash
tree . -L 2
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Clique na caixa colorida para expandir a seção e visualizar seu conteúdo.
Usamos seções recolhíveis como esta para incluir a saída esperada do comando de forma concisa.

- **Os arquivos `.nf`** são scripts de fluxo de trabalho que são nomeados com base em qual parte do curso eles são usados.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente.
  Você pode ignorá-lo por enquanto.

- **O arquivo `greetings.csv`** em `data/` contém dados de entrada que usaremos na maior parte do curso. Ele é descrito na Parte 2 (Canais), quando o introduzimos pela primeira vez.

- **Os arquivos `test-params.*`** são arquivos de configuração que usaremos na Parte 6 (Configuração). Você pode ignorá-los por enquanto.

- **O diretório `solutions`** contém os scripts de fluxo de trabalho completos que resultam de cada etapa do curso.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.

## Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está funcionando
- [ ] Eu defini meu diretório de trabalho apropriadamente

Se você pode marcar todas as caixas, está pronto para começar.

**Para continuar para a [Parte 1: Hello World](./01_hello_world.md), clique na seta no canto inferior direito desta página.**
