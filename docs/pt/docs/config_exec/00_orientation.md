# Primeiros passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


## Inicie um ambiente de treinamento

Para usar o ambiente pré-configurado que disponibilizamos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, consulte [Opções de ambiente](../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use o botão direito do mouse, ctrl+clique ou cmd+clique dependendo do seu equipamento) para que você possa continuar lendo enquanto o ambiente carrega.
Você precisará manter estas instruções abertas em paralelo para acompanhar o curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Noções básicas do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para acompanhar o curso, portanto você não precisa instalar nada por conta própria.

O codespace é configurado com uma interface VSCode, que inclui um explorador de arquivos, um editor de código e um terminal.
Todas as instruções fornecidas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') se referem a essas três partes da interface VSCode, salvo indicação em contrário.

Se você está acompanhando este curso por conta própria, familiarize-se com as [noções básicas do ambiente](../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este curso requer Nextflow 25.10.2 ou posterior, com o parser de sintaxe v2 habilitado (padrão no 25.10+).
Se você estiver usando um ambiente local ou personalizado, certifique-se de estar usando as configurações corretas conforme documentado [aqui](../info/nxf_versions.md).

## Prepare-se para trabalhar

Assim que seu codespace estiver em execução, há duas coisas a fazer antes de começar: definir seu diretório de trabalho e dar uma olhada nos materiais fornecidos.

### Defina o diretório de trabalho

Por padrão, o codespace abre na raiz de todos os cursos de treinamento.
Para este curso, mude para o diretório `config-exec/`:

```bash
cd config-exec/
```

Em seguida, configure o VSCode para focar neste diretório, de modo que apenas os arquivos relevantes apareçam na barra lateral do explorador de arquivos:

```bash
code .
```

!!! tip "Dica"

    Se por algum motivo você sair deste diretório (por exemplo, seu codespace entrar em modo de suspensão), você sempre pode usar o caminho completo para retornar a ele, assumindo que você está executando dentro do ambiente de treinamento do Github Codespaces:

    ```bash
    cd /workspaces/training/config-exec
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
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Os arquivos **`main.nf`** e **`modules/`** são o mesmo pipeline de múltiplas etapas do [Nextflow Run](../nextflow_run/index.md), e o arquivo **`nextflow.config`** é a mesma configuração que você já viu lá.
Você vai expandir ambos ao longo destes exercícios.

O diretório **`data/`** contém o arquivo de entrada CSV que o pipeline lê.

## Lista de verificação de prontidão

Acha que está pronto para começar?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está funcionando
- [ ] Defini meu diretório de trabalho adequadamente

Se você conseguir marcar todas as caixas, está pronto para continuar.

**Para continuar para a [Parte 1: Adapte ao seu ambiente de computação](./01_packaging_and_execution.md), clique na seta no canto inferior direito desta página.**
