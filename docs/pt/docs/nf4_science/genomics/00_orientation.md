# Primeiros Passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Inicie um ambiente de treinamento

Para usar o ambiente pré-configurado que fornecemos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, consulte [Opções de ambiente](../../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use clique direito, ctrl-clique ou cmd-clique dependendo do seu equipamento) para que você possa continuar lendo enquanto o ambiente carrega.
Você precisará manter estas instruções abertas em paralelo para trabalhar no curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Básico do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para trabalhar no curso de treinamento, então você não precisa instalar nada por conta própria.

O codespace é configurado com uma interface VSCode, que inclui um explorador de arquivos, um editor de código e um terminal shell.
Todas as instruções dadas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') referem-se a essas três partes da interface VSCode, a menos que especificado de outra forma.

Se você está trabalhando neste curso sozinho, por favor familiarize-se com os [fundamentos do ambiente](../../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este treinamento foi projetado para Nextflow 25.10.2 ou posterior **com o parser de sintaxe v2 HABILITADO**.
Se você está usando um ambiente local ou personalizado, certifique-se de estar usando as configurações corretas conforme documentado [aqui](../../info/nxf_versions.md).

## Prepare-se para trabalhar

Uma vez que seu codespace esteja em execução, há duas coisas que você precisa fazer antes de mergulhar no treinamento: definir seu diretório de trabalho para este curso específico e dar uma olhada nos materiais fornecidos.

### Defina o diretório de trabalho

Por padrão, o codespace abre com o diretório de trabalho definido na raiz de todos os cursos de treinamento, mas para este curso, trabalharemos no diretório `nf4-science/genomics/`.

Mude de diretório agora executando este comando no terminal:

```bash
cd nf4-science/genomics/
```

Você pode configurar o VSCode para focar neste diretório, de modo que apenas os arquivos relevantes apareçam na barra lateral do explorador de arquivos:

```bash
code .
```

!!! tip "Dica"

    Se por qualquer motivo você sair deste diretório (por exemplo, seu codespace entrar em suspensão), você sempre pode usar o caminho completo para retornar a ele, assumindo que você está executando isso dentro do ambiente de treinamento do Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Agora vamos dar uma olhada no conteúdo.

### Explore os materiais fornecidos

Você pode explorar o conteúdo deste diretório usando o explorador de arquivos no lado esquerdo do workspace de treinamento.
Alternativamente, você pode usar o comando `tree`.

Ao longo do curso, usamos a saída do `tree` para representar a estrutura e o conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice de conteúdo até o segundo nível:

```bash
tree . -L 2
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── data
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

Clique na caixa colorida para expandir a seção e visualizar seu conteúdo.
Usamos seções recolhíveis como esta para exibir saídas de comando esperadas, bem como conteúdos de diretórios e arquivos de forma concisa.

- **O arquivo `genomics.nf`** é um script de fluxo de trabalho que você construirá ao longo do curso.

- **O diretório `modules`** contém arquivos de módulo esqueleto que você preencherá durante o curso.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente.
  Você pode ignorá-lo por enquanto.

- **O diretório `data`** contém dados de entrada e recursos relacionados, descritos posteriormente no curso.

- **O diretório `solutions`** contém arquivos de módulo completos e uma solução da Parte 2 que pode servir como ponto de partida para a Parte 3.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.

## Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está funcionando
- [ ] Defini meu diretório de trabalho apropriadamente

Se você pode marcar todas as caixas, está pronto para começar.

**Para continuar para a [Parte 1: Visão geral do método e teste manual](./01_method.md), clique na seta no canto inferior direito desta página.**
