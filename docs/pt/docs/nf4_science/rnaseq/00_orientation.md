# Orientação

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

O ambiente de treinamento contém todo o software, código e dados necessários para trabalhar neste curso de treinamento, então você não precisa instalar nada por conta própria.
No entanto, você precisa de uma conta (gratuita) para fazer login, e deve dedicar alguns minutos para se familiarizar com a interface.

Se você ainda não o fez, por favor complete o mini-curso de [Configuração do Ambiente](../../envsetup/) antes de prosseguir.

## Materiais fornecidos

Ao longo deste curso de treinamento, trabalharemos no diretório `nf4-science/rnaseq/`, para o qual você precisa se mover ao abrir o espaço de trabalho de treinamento.
Este diretório contém todos os arquivos de código, dados de teste e arquivos acessórios que você precisará.

Sinta-se à vontade para explorar o conteúdo deste diretório; a maneira mais fácil de fazer isso é usar o explorador de arquivos no lado esquerdo do espaço de trabalho de treinamento na interface do VSCode.
Como alternativa, você pode usar o comando `tree`.
Ao longo do curso, usamos a saída do `tree` para representar a estrutura e o conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice de conteúdo até o segundo nível:

```bash
tree . -L 3
```

??? success "Conteúdo do diretório"

    ```console
    rnaseq
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

!!!note

    Não se preocupe se isso parecer muito; vamos passar pelas partes relevantes em cada etapa do curso.
    Isso é apenas para dar uma visão geral.

**Aqui está um resumo do que você deve saber para começar:**

- **O arquivo `rnaseq.nf`** é o esboço do script do fluxo de trabalho que desenvolveremos.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente. Você pode ignorá-lo por enquanto.

- **O diretório `data`** contém dados de entrada e recursos relacionados:

  - _Um genoma de referência_ chamado `genome.fa` consistindo de uma pequena região do cromossomo humano 20 (de hg19/b37).
  - _Dados de RNAseq_ que foram reduzidos a uma pequena região para manter os tamanhos de arquivo menores, no diretório `reads/`.
  - _Arquivos CSV_ listando os IDs e caminhos dos arquivos de dados de exemplo, para processamento em lotes.

- **O diretório `solutions`** contém os scripts de fluxo de trabalho completos e módulos que resultam de cada etapa do curso.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.
  O número no nome do arquivo corresponde à etapa da parte relevante do curso.

!!!tip

    Se por algum motivo você sair deste diretório, você sempre pode executar este comando para retornar a ele:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Agora, para começar o curso, clique na seta no canto inferior direito desta página.
