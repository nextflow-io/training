# Primeiros passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Inicie um ambiente de treinamento

Para usar o ambiente pré-construído que fornecemos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, veja [Opções de ambiente](../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use clique direito, ctrl-clique ou cmd-clique dependendo do seu equipamento) para que você possa continuar lendo enquanto o ambiente carrega.
Você precisará manter estas instruções abertas em paralelo para trabalhar no curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Noções básicas do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para trabalhar no curso de treinamento, então você não precisa instalar nada.

O codespace é configurado com uma interface VSCode, que inclui um explorador de sistema de arquivos, um editor de código e um terminal shell.
Todas as instruções dadas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') referem-se a essas três partes da interface VSCode, a menos que especificado de outra forma.

Se você está fazendo este curso sozinho, por favor familiarize-se com as [noções básicas do ambiente](../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este curso requer Nextflow 25.10.2 ou posterior, com o parser de sintaxe v2 habilitado (o padrão na versão 25.10+).
Se você está usando um ambiente local ou personalizado, certifique-se de estar usando as configurações corretas conforme documentado [aqui](../info/nxf_versions.md).

## Prepare-se para trabalhar

Uma vez que seu codespace esteja rodando, há duas coisas a fazer antes de mergulhar no treinamento: definir seu diretório de trabalho e dar uma olhada nos materiais fornecidos.

### Defina o diretório de trabalho

Por padrão, o codespace abre na raiz de todos os cursos de treinamento.
Para este curso, mude para o diretório `nextflow-run/`:

```bash
cd nextflow-run/
```

Em seguida, configure o VSCode para focar neste diretório, para que apenas os arquivos relevantes apareçam na barra lateral do explorador de arquivos:

```bash
code .
```

!!! tip "Dica"

    Se por qualquer razão você sair deste diretório (por exemplo, seu codespace adormecer), você sempre pode usar o caminho completo para retornar a ele, assumindo que está executando dentro do ambiente de treinamento do GitHub Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

### Explore os materiais fornecidos

Você pode explorar os materiais do curso usando o explorador de arquivos à esquerda, ou com o comando `tree`.
Execute o seguinte no terminal para ver a estrutura completa:

```bash
tree . -L 2
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── 1-hello.nf
    ├── 2-inputs.nf
    ├── data
    │   ├── greetings-extended.csv
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Os **arquivos `.nf`** são scripts de fluxo de trabalho de complexidade crescente, usados nessa ordem ao longo do curso.

O **diretório `data/`** contém os arquivos CSV de entrada que usaremos a partir da seção 2.

O **diretório `modules/`** contém as definições de processos usadas pelo `main.nf`.

O **arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente. Você pode ignorá-lo por enquanto; vamos abordá-lo na seção 4.

## Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está funcionando
- [ ] Eu defini meu diretório de trabalho apropriadamente

Se você pode marcar todas as caixas, está pronto para começar.

**Para continuar para [Parte 1: Executar Nextflow](./01_run_nextflow.md), clique na seta no canto inferior direito desta página.**
