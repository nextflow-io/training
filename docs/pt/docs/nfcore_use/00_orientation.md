# Primeiros passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Inicie um ambiente de treinamento

Para usar o ambiente pré-configurado que disponibilizamos no GitHub Codespaces, clique no botão "Open in GitHub Codespaces" abaixo. Para outras opções, consulte [Opções de ambiente](../envsetup/index.md).

Recomendamos abrir o ambiente de treinamento em uma nova aba ou janela do navegador (use o botão direito do mouse, ctrl+clique ou cmd+clique dependendo do seu equipamento) para que você possa continuar lendo enquanto o ambiente carrega.
Você precisará manter estas instruções abertas em paralelo para acompanhar o curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Noções básicas do ambiente

Este ambiente de treinamento contém todo o software, código e dados necessários para acompanhar o curso de treinamento, portanto você não precisa instalar nada por conta própria.

O codespace é configurado com uma interface VSCode, que inclui um explorador de arquivos, um editor de código e um terminal shell.
Todas as instruções fornecidas durante o curso (por exemplo, 'abra o arquivo', 'edite o código' ou 'execute este comando') se referem a essas três partes da interface do VSCode, salvo indicação em contrário.

Se você está acompanhando este curso por conta própria, familiarize-se com as [noções básicas do ambiente](../envsetup/01_setup.md) para mais detalhes.

### Requisitos de versão

Este treinamento funciona com Nextflow 25.10.2 ou posterior **com o parser de sintaxe v2**, que é o padrão a partir do Nextflow 26.04.
No nosso ambiente de treinamento você não precisa fazer nada: ele executa o Nextflow 26.04.4 com o parser v2. Se você estiver usando um ambiente local ou personalizado, consulte as [notas de versão](../info/nxf_versions.md).

!!! warning "nf-core/demo requer Nextflow 25.10.4 ou posterior"

    O pipeline `nf-core/demo` usado na Parte 1 impõe sua própria versão mínima do Nextflow (`>=25.10.4`), que é mais restritiva do que o requisito mínimo geral do treinamento de 25.10.2.
    Nosso ambiente de treinamento já satisfaz esse requisito; se você estiver usando um ambiente local ou personalizado, certifique-se de estar no Nextflow 25.10.4 ou posterior.

Este treinamento também requer **nf-core tools 4.0.2**.
Se você usar uma versão diferente das ferramentas nf-core, pode ter dificuldades para acompanhar.

Você pode verificar qual versão está instalada no seu ambiente usando o comando `nf-core --version`.

!!! warning "Compatibilidade com o parser v2"

    Muitos pipelines nf-core ainda não suportam o parser de sintaxe v2.
    Se você executar um pipeline nf-core diferente dos usados neste curso e encontrar erros, pode ser necessário alternar para o parser v1 definindo `export NXF_SYNTAX_PARSER=v1`.
    Consulte as [notas de versão](../info/nxf_versions.md) para mais detalhes.

## Prepare-se para trabalhar

Assim que seu codespace estiver em execução, há duas coisas que você precisa fazer antes de mergulhar no treinamento: definir seu diretório de trabalho para este curso específico e dar uma olhada nos materiais fornecidos.

### Defina o diretório de trabalho

Por padrão, o codespace abre com o diretório de trabalho definido na raiz de todos os cursos de treinamento, mas para este curso, trabalharemos no diretório `nfcore-use/`.

Mude para esse diretório agora executando este comando no terminal:

```bash
cd nfcore-use/
```

!!! tip "Dica"

    Se por algum motivo você sair deste diretório (por exemplo, seu codespace entrar em modo de suspensão), você sempre pode usar o caminho completo para retornar a ele, assumindo que você está executando dentro do ambiente de treinamento do Github Codespaces:

    ```bash
    cd /workspaces/training/nfcore-use
    ```

Em seguida, explore o conteúdo deste diretório.

### Explore os materiais fornecidos

Você pode explorar o conteúdo deste diretório usando o explorador de arquivos no lado esquerdo do espaço de trabalho de treinamento.
Como alternativa, você pode usar o comando `tree`.

```bash
tree .
```

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── custom.config
    ├── laptop.config
    ├── malformed_samplesheet.csv
    └── my_params.yml
    ```

- **O arquivo `laptop.config`** é um arquivo de configuração que usaremos na seção 4 para limitar o uso de recursos ao executar um pipeline em escala de produção localmente.
  Você pode ignorá-lo até lá.
- **Os arquivos `my_params.yml`, `malformed_samplesheet.csv` e `custom.config`** são usados na Parte 2, para demonstrar a definição de parâmetros a partir de um arquivo, validação de entrada e substituições de configuração no nível do processo.
  Você também pode ignorá-los até lá.

## Lista de verificação de prontidão

Acha que está pronto para começar?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu ambiente está em funcionamento
- [ ] Estou usando nf-core tools 4.0.2 (verifique com `nf-core --version`)
- [ ] Defini meu diretório de trabalho adequadamente

Se você conseguir marcar todas as caixas, está pronto para começar.

**Para continuar para a Parte 1, clique na seta no canto inferior direito desta página.**
