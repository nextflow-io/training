# Primeiros Passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar um ambiente de treinamento

Para usar o ambiente pré-configurado que fornecemos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, consulte [Opções de ambiente](../../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use clique direito, ctrl-clique ou cmd-clique dependendo do seu equipamento) para que você possa ler enquanto o ambiente carrega.
Você precisará manter estas instruções abertas em paralelo para trabalhar no curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Noções básicas do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para trabalhar no curso de treinamento, então você não precisa instalar nada por conta própria.

O codespace é configurado com uma interface VSCode, que inclui um explorador de arquivos, um editor de código e um terminal shell.
Todas as instruções dadas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') referem-se a essas três partes da interface VSCode, a menos que especificado de outra forma.

Se você está trabalhando neste curso por conta própria, por favor familiarize-se com as [noções básicas do ambiente](../../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este treinamento é projetado para Nextflow 25.10.2 ou posterior **com o analisador de sintaxe v2 HABILITADO**.
Se você está usando um ambiente local ou personalizado, por favor certifique-se de estar usando as configurações corretas conforme documentado [aqui](../../info/nxf_versions.md).

## Prepare-se para trabalhar

Uma vez que seu codespace esteja em execução, há duas coisas que você precisa fazer antes de mergulhar no treinamento: definir seu diretório de trabalho para este curso específico e dar uma olhada nos materiais fornecidos.

### Definir o diretório de trabalho

Por padrão, o codespace abre com o diretório de trabalho definido na raiz de todos os cursos de treinamento, mas para este curso, trabalharemos no diretório `nf4-science/rnaseq/`.

Mude de diretório agora executando este comando no terminal:

```bash
cd nf4-science/rnaseq/
```

Você pode configurar o VSCode para focar neste diretório, de modo que apenas os arquivos relevantes apareçam na barra lateral do explorador de arquivos:

```bash
code .
```

!!! tip "Dica"

    Se por algum motivo você sair deste diretório (por exemplo, seu codespace entrar em suspensão), você sempre pode usar o caminho completo para retornar a ele, assumindo que você está executando isso dentro do ambiente de treinamento do Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Agora vamos dar uma olhada no conteúdo.

### Explorar os materiais fornecidos

Você pode explorar o conteúdo deste diretório usando o explorador de arquivos no lado esquerdo do espaço de trabalho de treinamento.
Como alternativa, você pode usar o comando `tree`.

Ao longo do curso, usamos a saída do `tree` para representar a estrutura e o conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice de conteúdo até o terceiro nível:

```bash
tree . -L 3
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

Clique na caixa colorida para expandir a seção e visualizar seu conteúdo.
Usamos seções recolhíveis como esta para exibir a saída esperada de comandos, bem como o conteúdo de diretórios e arquivos de forma concisa.

- **O arquivo `rnaseq.nf`** é um esboço de um script de fluxo de trabalho que você construirá conforme avança no curso.

- **O diretório `modules`** contém esboços de módulos de processo que você preencherá durante o curso.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente.
  Você pode ignorá-lo por enquanto.

- **O diretório `data`** contém dados de entrada e recursos relacionados, descritos mais adiante no curso.

- **O diretório `solutions`** contém scripts de fluxo de trabalho completos e módulos que resultam de cada etapa do curso.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.
  A solução da Parte 2 pode ser usada como ponto de partida para a Parte 3.

## Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está funcionando
- [ ] Defini meu diretório de trabalho apropriadamente

Se você pode marcar todas as caixas, está pronto para começar.

**Para continuar para [Parte 1: Visão geral do método](./01_method.md), clique na seta no canto inferior direito desta página.**
