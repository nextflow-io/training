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

Este treinamento é projetado para Nextflow 25.10.2 ou posterior **com o parser de sintaxe v2 HABILITADO**.
Se você está usando um ambiente local ou personalizado, certifique-se de estar usando as configurações corretas conforme documentado [aqui](../info/nxf_versions.md).

## Prepare-se para trabalhar

Uma vez que seu codespace esteja rodando, há duas coisas que você precisa fazer antes de mergulhar no treinamento: definir seu diretório de trabalho para este curso específico e dar uma olhada nos materiais fornecidos.

### Defina o diretório de trabalho

Por padrão, o codespace abre com o diretório de trabalho definido na raiz de todos os cursos de treinamento, mas para este curso, trabalharemos no diretório `nextflow-run/`.

Mude de diretório agora executando este comando no terminal:

```bash
cd nextflow-run/
```

Você pode configurar o VSCode para focar neste diretório, para que apenas os arquivos relevantes apareçam na barra lateral do explorador de arquivos:

```bash
code .
```

!!! tip "Dica"

    Se por qualquer razão você sair deste diretório (por exemplo, seu codespace adormecer), você sempre pode usar o caminho completo para retornar a ele, assumindo que está executando dentro do ambiente de treinamento do GitHub Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Agora vamos dar uma olhada no conteúdo.

### Explore os materiais fornecidos

Você pode explorar o conteúdo deste diretório usando o explorador de arquivos no lado esquerdo do espaço de trabalho de treinamento.
Alternativamente, você pode usar o comando `tree`.

Ao longo do curso, usamos a saída do `tree` para representar a estrutura e conteúdo de diretórios de forma legível, às vezes com pequenas modificações para clareza.

Aqui geramos uma tabela de conteúdo até o segundo nível:

```bash
tree . -L 2
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Clique na caixa colorida para expandir a seção e visualizar seu conteúdo.
Usamos seções recolhíveis como esta para exibir a saída esperada de comandos, bem como conteúdo de diretórios e arquivos de forma concisa.

- **Os arquivos `.nf`** são scripts de fluxo de trabalho numerados com base na parte do curso em que são usados.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente.
  Você pode ignorá-lo por enquanto.

- **O arquivo `greetings.csv`** em `data/` contém dados de entrada que usaremos na maior parte do curso. Ele é descrito na Parte 2 (Executar pipelines), quando o introduzimos pela primeira vez.

- **Os arquivos `test-params.*`** são arquivos de configuração que usaremos na Parte 3 (Configuração). Você pode ignorá-los por enquanto.

- **O diretório `solutions`** contém o estado final do fluxo de trabalho e seus arquivos acessórios (config e módulos) que resultam da conclusão do curso.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.

## Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está funcionando
- [ ] Eu defini meu diretório de trabalho apropriadamente

Se você pode marcar todas as caixas, está pronto para começar.

**Para continuar para [Parte 1: Operações Básicas de Execução](./01_basics.md), clique na seta no canto inferior direito desta página.**
