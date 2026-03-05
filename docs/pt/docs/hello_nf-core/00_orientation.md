# Primeiros Passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar um ambiente de treinamento

Para usar o ambiente pré-construído que fornecemos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, consulte [Opções de ambiente](../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use clique-direito, ctrl+clique ou cmd+clique dependendo do seu equipamento) para que você possa continuar lendo enquanto o ambiente carrega.
Você precisará manter estas instruções abertas em paralelo para trabalhar no curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Básico do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para trabalhar no curso de treinamento, então você não precisa instalar nada por conta própria.

O codespace é configurado com uma interface VSCode, que inclui um explorador de sistema de arquivos, um editor de código e um terminal shell.
Todas as instruções dadas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') referem-se a essas três partes da interface VSCode, salvo especificação em contrário.

Se você está trabalhando neste curso por conta própria, por favor familiarize-se com o [básico do ambiente](../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este treinamento é projetado para **Nextflow 25.10.2** ou posterior **com o analisador de sintaxe v2 DESABILITADO**.

#### Se você está usando nosso ambiente de treinamento:

Você DEVE executar o seguinte comando antes de prosseguir:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Se você está usando um ambiente local ou personalizado:

Por favor certifique-se de que está usando as configurações corretas conforme documentado [aqui](../info/nxf_versions.md).

O treinamento adicionalmente requer **nf-core tools 3.4.1**.
Se você usar uma versão diferente das ferramentas nf-core, você pode ter dificuldades em acompanhar.

Você pode verificar qual versão está instalada em seu ambiente usando o comando `nf-core --version`.

## Prepare-se para trabalhar

Uma vez que seu codespace esteja executando, há duas coisas que você precisa fazer antes de mergulhar no treinamento: configurar seu diretório de trabalho para este curso específico e dar uma olhada nos materiais fornecidos.

### Configurar o diretório de trabalho

Por padrão, o codespace abre com o diretório de trabalho configurado na raiz de todos os cursos de treinamento, mas para este curso, trabalharemos no diretório `hello-nf-core/`.

Mude de diretório agora executando este comando no terminal:

```bash
cd hello-nf-core/
```

!!! tip "Dica"

    Se por algum motivo você sair deste diretório (por exemplo, seu codespace entrar em modo de espera), você sempre pode usar o caminho completo para retornar a ele, assumindo que você está executando isto dentro do ambiente de treinamento GitHub Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Agora vamos dar uma olhada no conteúdo deste diretório.

### Explorar os materiais fornecidos

Você pode explorar o conteúdo deste diretório usando o explorador de arquivos no lado esquerdo do espaço de trabalho de treinamento.
Alternativamente, você pode usar o comando `tree`.

Durante o curso, usamos a saída de `tree` para representar a estrutura e conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice até o segundo nível:

```bash
tree . -L 2
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Clique na caixa colorida para expandir a seção e visualizar seu conteúdo.
Usamos seções recolhíveis como esta para incluir a saída esperada de comandos de forma concisa.

- **O arquivo `greetings.csv`** é um CSV contendo alguns dados colunares mínimos que usamos para fins de teste.

- **O diretório `original-hello`** contém uma cópia do código fonte produzido ao trabalhar na série completa de treinamento Hello Nextflow (com Docker habilitado).

- **O diretório `solutions`** contém os scripts de fluxo de trabalho completos que resultam de cada etapa do curso.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.

## Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está ativo e funcionando
- [ ] Eu me certifiquei de que o analisador de sintaxe está configurado para **v1**
- [ ] Eu configurei meu diretório de trabalho apropriadamente

Se você pode marcar todas as caixas, você está pronto para começar.

**Para continuar para a Parte 1, clique na seta no canto inferior direito desta página.**
