# Orientação

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

O ambiente de treinamento contém todo o software, código e dados necessários para trabalhar neste curso de treinamento, então você não precisa instalar nada por conta própria.
No entanto, você precisa de uma conta (gratuita) para fazer login, e deve dedicar alguns minutos para se familiarizar com a interface.

Se você ainda não o fez, por favor siga [este link](../../../envsetup/) antes de prosseguir.

## Materiais fornecidos

Ao longo deste curso de treinamento, trabalharemos no diretório `nf4-science/genomics/`, para o qual você precisa se mover quando abrir o workspace de treinamento.
Este diretório contém todos os arquivos de código, dados de teste e arquivos acessórios que você precisará.

Sinta-se à vontade para explorar o conteúdo deste diretório; a maneira mais fácil de fazer isso é usar o explorador de arquivos no lado esquerdo do workspace de treinamento na interface do VSCode.
Alternativamente, você pode usar o comando `tree`.
Ao longo do curso, usamos a saída do `tree` para representar a estrutura e o conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice de conteúdo até o segundo nível:

```bash
tree . -L 2
```

Se você executar isso dentro de `nf4-science/genomics`, você deverá ver a seguinte saída:

```console title="Conteúdo do diretório"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "Nota"

    Não se preocupe se isso parecer muito; vamos passar pelas partes relevantes em cada etapa do curso.
    Isso é apenas para dar uma visão geral.

**Aqui está um resumo do que você deve saber para começar:**

- **Os arquivos `.nf`** são scripts de fluxo de trabalho que são nomeados com base em qual parte do curso eles são usados.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente.
  Você pode ignorá-lo por enquanto.

- **O diretório `data`** contém dados de entrada e recursos relacionados, descritos posteriormente no curso.

- **O diretório `solutions`** contém arquivos de módulo e configurações de teste que resultam das Partes 3 e 4 do curso.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.

!!!tip "Dica"

    Se por qualquer motivo você sair deste diretório, você sempre pode executar este comando para retornar a ele:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Agora, para começar o curso, clique na seta no canto inferior direito desta página.
