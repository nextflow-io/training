# Orientação

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Esta orientação pressupõe que você já abriu o ambiente de treinamento clicando no botão "Open in GitHub Codespaces".
Se não, por favor, faça isso agora, idealmente em uma segunda janela ou aba do navegador para que você possa consultar estas instruções.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Requisito de tamanho da máquina"

    Certifique-se de selecionar uma **máquina de 8 núcleos** ao criar seu Codespace para este curso de treinamento. Os fluxos de trabalho de bioimagem requerem recursos computacionais adicionais.

## GitHub Codespaces

O ambiente GitHub Codespaces contém todo o software, código e dados necessários para trabalhar neste curso de treinamento, então você não precisa instalar nada por conta própria.
No entanto, você precisa de uma conta GitHub (gratuita) para fazer login, e se você não está familiarizado com a interface, deve levar alguns minutos para se familiarizar com ela completando o mini-curso [Orientação do GitHub Codespaces](../../envsetup/index.md).

## Pré-download das imagens Docker

Depois de abrir seu Codespace, vamos fazer o pré-download de todas as imagens Docker que precisaremos para este curso de treinamento.
Isso economizará tempo mais tarde e garantirá a execução suave dos fluxos de trabalho.

Abra uma nova aba do terminal e execute o seguinte comando:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Este comando fará o download de todas as imagens Docker necessárias em segundo plano.
Você pode continuar com o restante da orientação enquanto isso é executado.

!!!tip "Dica"

    A flag `-stub` permite que o pipeline seja executado rapidamente sem processar dados reais, o que é perfeito para baixar imagens. Você pode monitorar o progresso na aba do terminal.

## Diretório de trabalho

Ao longo deste curso de treinamento, trabalharemos no diretório `nf4-science/imaging/`.

Mude de diretório agora executando este comando no terminal:

```bash
cd nf4-science/imaging/
```

!!!tip "Dica"

    Se por qualquer motivo você sair deste diretório, você sempre pode usar o caminho completo para retornar a ele, assumindo que você está executando isso dentro do ambiente de treinamento GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Agora, para começar o curso, clique na seta no canto inferior direito desta página.**
