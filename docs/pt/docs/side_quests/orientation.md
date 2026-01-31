# Orientação

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

O ambiente GitHub Codespaces contém todo o software, código e dados necessários para trabalhar neste curso de treinamento, então você não precisa instalar nada por conta própria.
No entanto, você precisa de uma conta (gratuita) para fazer login, e deve dedicar alguns minutos para se familiarizar com a interface.

Se você ainda não fez isso, por favor siga [este link](../../envsetup/) antes de prosseguir.

## Materiais fornecidos

Ao longo deste curso de treinamento, trabalharemos no diretório `side-quests/`.
Este diretório contém todos os arquivos de código, dados de teste e arquivos acessórios que você precisará.

Sinta-se livre para explorar o conteúdo deste diretório; a maneira mais fácil de fazer isso é usar o explorador de arquivos no lado esquerdo do espaço de trabalho do GitHub Codespaces.
Como alternativa, você pode usar o comando `tree`.
Ao longo do curso, usamos a saída de `tree` para representar a estrutura e o conteúdo do diretório de forma legível, às vezes com pequenas modificações para maior clareza.

Aqui geramos um índice de conteúdo até o segundo nível:

```bash
tree . -L 2
```

Se você executar isso dentro de `side-quests`, deverá ver a seguinte saída:

```console title="Conteúdo do diretório"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Aqui está um resumo do que você deve saber para começar:**

- **Cada diretório corresponde a uma missão secundária individual.**
  Seus conteúdos são detalhados na página da missão secundária correspondente.

- **O diretório `solutions`** contém os scripts de fluxo de trabalho e/ou módulo completos que resultam da execução de várias etapas de cada missão secundária.
  Eles são destinados a serem usados como referência para verificar seu trabalho e solucionar quaisquer problemas.

!!!tip "Dica"

    Se por qualquer motivo você sair deste diretório, você sempre pode executar este comando para retornar a ele:

    ```bash
    cd /workspaces/training/side-quests
    ```

Agora, para começar o curso, clique na seta no canto inferior direito desta página.
