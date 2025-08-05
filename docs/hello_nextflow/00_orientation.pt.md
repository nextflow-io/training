# Orientação

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja a [playlist completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) no canal de YouTube do Nextflow.

:green_book: A transcrição desse vídeo está disponível [aqui](./transcripts/00_orientation.md).
///

O ambiente GitHub Codespaces contém todo o software, código e dados necessários para este curso. Você não precisa instalar nada por conta própria. No entanto, é necessária uma conta (gratuita) para logar - e recomendamos que você reserve alguns minutos para se familiarizar com a interface.

Caso ainda não tenha feito isso, siga [este link](../../envsetup/) antes de prosseguir.

## Materiais fornecidos

Ao longo deste curso, trabalharemos no diretório `hello-nextflow/`, que é carregado por padrão quando você abre o ambiente de trabalho do Gitpod. Este diretório contém todos os arquivos de código, dados de teste e arquivos auxiliares que você precisará.

Sinta-se à vontade para explorar o conteúdo deste diretório; a maneira mais fácil de fazer isso é usando o explorador de arquivos no lado esquerdo do ambiente de trabalho do Gitpod. Alternativamente, você pode usar o comando `tree`. Ao longo do curso, usamos a saída do `tree` para representar a estrutura e o conteúdo do diretório de forma legível - às vezes com pequenas modificações para maior clareza.

Aqui, geramos um índice até o segundo nível:

```bash
tree . -L 2
```

Se você executar isso dentro do diretório `hello-nextflow`, verá a seguinte saída:

```console title="Directory contents"
.
├── greetings.csv
├── hello-channels.nf
├── hello-config.nf
├── hello-containers.nf
├── hello-modules.nf
├── hello-workflow.nf
├── hello-world.nf
├── nextflow.config
├── solutions
│   ├── 1-hello-world
│   ├── 2-hello-channels
│   ├── 3-hello-workflow
│   ├── 4-hello-modules
│   ├── 5-hello-containers
│   └── 6-hello-config
└── test-params.json

7 directories, 9 files
```

!!!nota

    Não se preocupe se isso parecer muita informação até o momento. Nós passaremos pelas partes relevantes ao longo das etapas do curso. Isso é apenas para te dar uma visão geral.

**Aqui está um resumo do que você deveria saber para começar:**

- **Os arquivos `.nf`** são _scripts_ de fluxo de trabalho nomeados com base na parte do curso em que são utilizados.

- **O arquivo `nextflow.config`** é um arquivo de configuração que define propriedades mínimas do ambiente. Você pode ignorá-lo por enquanto.

- **The file `greetings.csv`** contains input data we'll use in most of the course. It is described in Part 1, when we introduce it for the first time.

- **The file `test-params.json`** is a file we'll use in Part 6. You can ignore it for now.

- **O diretório `solutions`** contém os _scripts_ de fluxo de trabalho completos resultantes de cada etapa do curso. Eles servem como referência para verificar seu trabalho e solucionar quaisquer problemas. As informações de nome e o número presentes no nome do arquivo correspondem à etapa da parte relevante do curso. Por exemplo, o arquivo `hello-world-4.nf` é o resultado esperado ao completar as etapas 1 a 4 da Parte 1: Hello World.

!!!dica

    Se, por algum motivo, você sair deste diretório, sempre poderá executar este comando para retornar:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Agora, para iniciar o curso, clique na seta no canto inferior direito desta página.
