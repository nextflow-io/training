# Parte 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=pt" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja [a playlist completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) no canal do Nextflow no YouTube.

:green_book: A transcrição do vídeo está disponível [aqui](./transcripts/01_hello_world.md).
///

Nesta primeira parte do curso de treinamento Hello Nextflow, começamos com um exemplo Hello World muito básico e independente de domínio, que vamos construir progressivamente para demonstrar o uso da lógica e componentes fundamentais do Nextflow.

??? info "O que é um exemplo Hello World?"

    Um "Hello World!" é um exemplo minimalista que tem como objetivo demonstrar a sintaxe e estrutura básicas de uma linguagem de programação ou framework de software.
    O exemplo geralmente consiste em imprimir a frase "Hello, World!" no dispositivo de saída, como o console ou terminal, ou escrevê-la em um arquivo.

---

## 0. Aquecimento: Execute um exemplo Hello World diretamente

Vamos demonstrar isso com um comando simples que executamos diretamente no terminal, para mostrar o que ele faz antes de encapsulá-lo no Nextflow.

!!! tip "Dica"

    Lembre-se de que você deve estar dentro do diretório `hello-nextflow/` conforme descrito na página [Primeiros Passos](00_orientation.md).

### 0.1. Faça o terminal dizer olá

Execute o seguinte comando no seu terminal.

```bash
echo 'Hello World!'
```

??? success "Saída do comando"

    ```console
    Hello World!
    ```

Isso exibe o texto 'Hello World' diretamente no terminal.

### 0.2. Escreva a saída em um arquivo

A execução de pipelines envolve principalmente a leitura de dados de arquivos e a escrita de resultados em outros arquivos, então vamos modificar o comando para escrever a saída de texto em um arquivo para tornar o exemplo um pouco mais relevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Saída do comando"

    ```console

    ```

Isso não exibe nada no terminal.

### 0.3. Encontre a saída

O texto 'Hello World' agora deve estar no arquivo de saída que especificamos, chamado `output.txt`.
Você pode abri-lo no explorador de arquivos ou pela linha de comando usando o utilitário `cat`, por exemplo.

??? abstract "Conteúdo do arquivo"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Isso é o que vamos tentar replicar com nosso primeiro fluxo de trabalho Nextflow.

### Conclusão

Você agora sabe como executar um comando simples no terminal que exibe algum texto e, opcionalmente, como fazer com que ele escreva a saída em um arquivo.

### O que vem a seguir?

Descubra como isso ficaria escrito como um fluxo de trabalho Nextflow.

---

## 1. Examine o script e execute-o

Fornecemos um script de fluxo de trabalho totalmente funcional, embora minimalista, chamado `hello-world.nf` que faz a mesma coisa que antes (escrever 'Hello World!'), mas com Nextflow.

Para começar, vamos abrir o script do fluxo de trabalho para que você tenha uma noção de como ele está estruturado.
Então vamos executá-lo e procurar suas saídas.

### 1.1. Examine o código

Você encontrará o script `hello-world.nf` no seu diretório atual, que deve ser `hello-nextflow`. Abra-o no painel do editor.

??? full-code "Arquivo de código completo"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usa echo para imprimir 'Hello World!' em um arquivo
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // emite uma saudação
        sayHello()
    }
    ```

Um script de fluxo de trabalho Nextflow normalmente inclui uma ou mais definições de [**process**](https://nextflow.io/docs/latest/process.html) e o [**workflow**](https://nextflow.io/docs/latest/workflow.html) em si, além de alguns blocos opcionais (não presentes aqui) que apresentaremos mais tarde.

Cada **process** descreve quais operações a etapa correspondente no pipeline deve realizar, enquanto o **workflow** descreve a lógica de fluxo de dados que conecta as várias etapas.

Vamos examinar primeiro o bloco **process** e depois veremos o bloco **workflow**.

#### 1.1.1. A definição de `process`

O primeiro bloco de código descreve um **process**.

A definição do processo começa com a palavra-chave `process`, seguida pelo nome do processo e finalmente o corpo do processo delimitado por chaves.
O corpo do processo deve conter um bloco script que especifica o comando a ser executado, que pode ser qualquer coisa que você seria capaz de executar em um terminal de linha de comando.

```groovy title="hello-world.nf" linenums="3"
/*
* Usa echo para imprimir 'Hello World!' em um arquivo
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Aqui temos um **process** chamado `sayHello` que escreve sua **saída** em um arquivo chamado `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Esta é uma definição de processo muito mínima que contém apenas uma definição de `output` e o `script` a ser executado.

A definição de `output` inclui o qualificador `path`, que diz ao Nextflow que isso deve ser tratado como um caminho (inclui tanto caminhos de diretório quanto arquivos).
Outro qualificador comum é `val`.

É importante notar que a definição de saída não _determina_ qual saída será criada.
Ela simplesmente _declara_ qual é a saída esperada, para que o Nextflow possa procurá-la assim que a execução estiver completa.
Isso é necessário para verificar se o comando foi executado com sucesso e para passar a saída para processos subsequentes, se necessário. A saída produzida que não corresponder ao que está declarado no bloco de saída não será passada para processos subsequentes.

!!! warning "Aviso"

    Este exemplo é frágil porque codificamos o nome do arquivo de saída em dois lugares separados (o script e os blocos de saída).
    Se mudarmos um mas não o outro, o script vai falhar.
    Mais tarde, você aprenderá maneiras de usar variáveis para mitigar esse problema.

Em um pipeline do mundo real, um processo geralmente contém blocos adicionais como diretivas e entradas, que apresentaremos em breve.

#### 1.1.2. A definição de `workflow`

O segundo bloco de código descreve o **workflow** em si.
A definição de fluxo de trabalho começa com a palavra-chave `workflow`, seguida por um nome opcional, depois o corpo do fluxo de trabalho delimitado por chaves.

Aqui temos um **workflow** que consiste em um bloco `main:` (que diz 'este é o corpo principal do fluxo de trabalho') contendo uma chamada ao processo `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emite uma saudação
    sayHello()
}
```

Esta é uma definição de **workflow** muito mínima.
Em um pipeline do mundo real, o fluxo de trabalho geralmente contém múltiplas chamadas a **processos** conectados por **canais**, e os processos esperam uma ou mais **entradas** variáveis.

Você aprenderá como adicionar entradas variáveis mais tarde neste módulo de treinamento; e aprenderá como adicionar mais processos e conectá-los por canais na Parte 3 deste curso.

!!! tip "Dica"

    Tecnicamente, a linha `main:` não é necessária para fluxos de trabalho simples como este, então você pode encontrar fluxos de trabalho que não a tenham.
    Mas precisaremos dela para aproveitar as saídas em nível de fluxo de trabalho, então podemos incluí-la desde o início.

### 1.2. Execute o fluxo de trabalho

Olhar para o código não é tão divertido quanto executá-lo, então vamos experimentar isso na prática.

#### 1.2.1. Lance o fluxo de trabalho e monitore a execução

No terminal, execute o seguinte comando:

```bash
nextflow run hello-world.nf
```

??? success "Saída do comando"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Se a saída do seu console se parece com isso, então parabéns, você acabou de executar seu primeiro fluxo de trabalho Nextflow!

A saída mais importante aqui é a última linha, que está destacada na saída acima:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Isso nos diz que o processo `sayHello` foi executado com sucesso uma vez (`1 of 1 ✔`).

É importante notar que esta linha também informa onde encontrar a saída da chamada do processo `sayHello`.
Vamos ver isso agora.

#### 1.2.2. Encontre a saída e os logs no diretório `work`

Quando você executa o Nextflow pela primeira vez em um determinado diretório, ele cria um diretório chamado `work` onde escreverá todos os arquivos (e quaisquer links simbólicos) gerados durante a execução.

Dentro do diretório `work`, o Nextflow organiza saídas e logs por chamada de processo.
Para cada chamada de processo, o Nextflow cria um subdiretório aninhado, nomeado com um hash para torná-lo único, onde preparará todas as entradas necessárias (usando links simbólicos por padrão), escreverá arquivos auxiliares e escreverá logs e quaisquer saídas do processo.

O caminho para esse subdiretório é mostrado de forma truncada entre colchetes na saída do console.
Olhando para o que obtivemos na execução mostrada acima, a linha de log do console para o processo sayHello começa com `[65/7be2fa]`. Isso corresponde ao seguinte caminho de diretório: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Vamos dar uma olhada no que há lá.

??? abstract "Conteúdo do diretório"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Não está vendo a mesma coisa?"

    Os nomes exatos dos subdiretórios serão diferentes no seu sistema.

    Se você navegar pelo conteúdo do subdiretório de tarefa no explorador de arquivos do VSCode, verá todos os arquivos imediatamente.
    No entanto, os arquivos de log estão configurados para serem invisíveis no terminal, então se você quiser usar `ls` ou `tree` para visualizá-los, precisará definir a opção relevante para exibir arquivos invisíveis.

    ```bash
    tree -a work
    ```

A primeira coisa que você deseja ver é a saída real do fluxo de trabalho, ou seja, o arquivo `output.txt` produzido pelo processo `sayHello`.
Abra-o e você encontrará a saudação `Hello World!`, que era o objetivo do nosso fluxo de trabalho minimalista.

??? abstract "Conteúdo do arquivo"

    ```console title="output.txt"
    Hello World!
    ```

Funcionou!

Concedido, pode parecer muito código de wrapper para um resultado tão pequeno, mas o valor de todo esse código de wrapper se tornará mais óbvio quando começarmos a ler arquivos de entrada e encadear várias etapas.

Dito isso, vamos também olhar os outros arquivos naquele diretório. Esses são arquivos auxiliares e de log produzidos pelo Nextflow como parte da execução da tarefa.

- **`.command.begin`**: Metadados relacionados ao início da execução da chamada do processo
- **`.command.err`**: Mensagens de erro (`stderr`) emitidas pela chamada do processo
- **`.command.log`**: Saída de log completa emitida pela chamada do processo
- **`.command.out`**: Saída regular (`stdout`) pela chamada do processo
- **`.command.run`**: Script completo executado pelo Nextflow para executar a chamada do processo
- **`.command.sh`**: O comando que foi realmente executado pela chamada do processo
- **`.exitcode`**: O código de saída resultante do comando

O arquivo `.command.sh` é especialmente útil porque informa o comando principal que o Nextflow executou, não incluindo toda a contabilidade e configuração de tarefa/ambiente.

??? abstract "Conteúdo do arquivo"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Isso corresponde ao que executamos anteriormente manualmente.

Neste caso é muito direto porque o comando do processo foi codificado, mas mais adiante no curso você verá comandos de processo que envolvem alguma interpolação de variáveis.
Isso torna especialmente valioso poder ver exatamente como o Nextflow interpretou o código e qual comando foi produzido quando você está solucionando problemas de uma execução com falha.

### 1.3. Execute o fluxo de trabalho novamente

Tente executar novamente o fluxo de trabalho algumas vezes e depois olhe os diretórios de tarefa em `work/`.

??? abstract "Conteúdo do diretório"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Você verá que um novo subdiretório com um conjunto completo de arquivos de saída e log foi criado para cada execução.
Isso mostra que executar o mesmo fluxo de trabalho várias vezes não sobrescreverá os resultados de execuções anteriores.

### Conclusão

Você sabe como decifrar um script Nextflow simples, executá-lo e encontrar a saída e arquivos de log relevantes no diretório work.

### O que vem a seguir?

Aprenda como publicar as saídas do fluxo de trabalho em um local mais conveniente.

---

## 2. Publique saídas

Como você acabou de aprender, a saída produzida pelo nosso pipeline está enterrada em um diretório de trabalho várias camadas abaixo.
Isso é feito propositalmente; o Nextflow está no controle deste diretório e não devemos interagir com ele.
No entanto, isso torna inconveniente recuperar saídas que nos interessam.

Felizmente, o Nextflow fornece uma maneira de publicar saídas em um diretório designado usando [definições de saída em nível de fluxo de trabalho](https://nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Uso básico

Isso vai envolver dois novos pedaços de código:

1. Um bloco `publish:` dentro do corpo do `workflow`, declarando saídas de processo.
2. Um bloco `output` no script especificando opções de saída como modo e localização.

#### 2.1.1. Declare a saída do processo `sayHello`

Precisamos adicionar um bloco `publish:` ao corpo do fluxo de trabalho (mesmo tipo de elemento de código que o bloco `main:`) e listar a saída do processo `sayHello()`.

No arquivo de script do fluxo de trabalho `hello-world.nf`, adicione as seguintes linhas de código:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emite uma saudação
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emite uma saudação
        sayHello()
    }
    ```

Você vê que podemos nos referir à saída do processo simplesmente fazendo `sayHello().out`, e atribuir a ela um nome arbitrário, `first_output`.

#### 2.1.2. Adicione um bloco `output:` ao script

Agora só precisamos adicionar o bloco `output:` onde o caminho do diretório de saída será especificado. Note que este novo bloco fica **fora** e **abaixo** do bloco `workflow` dentro do script.

No arquivo de script do fluxo de trabalho `hello-world.nf`, adicione as seguintes linhas de código:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emite uma saudação
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emite uma saudação
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Podemos usar isso para atribuir caminhos específicos a quaisquer saídas de processo declaradas no bloco `workflow`.
Mais tarde, você aprenderá sobre maneiras de gerar estruturas sofisticadas de diretório de saída, mas por enquanto, estamos apenas codificando um caminho mínimo para simplicidade.

#### 2.1.3. Execute o fluxo de trabalho

Agora execute o script de fluxo de trabalho modificado:

```bash
nextflow run hello-world.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

A saída do terminal deve parecer familiar. Externamente, nada mudou.

No entanto, verifique seu explorador de arquivos: desta vez, o Nextflow criou um novo diretório chamado `results/`.

??? abstract "Conteúdo do diretório"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

Dentro do diretório `results`, encontramos um link simbólico para o `output.txt` produzido no diretório work pelo comando que acabamos de executar.

Isso nos permite recuperar facilmente arquivos de saída sem ter que vasculhar o subdiretório work.

### 2.2. Defina um local personalizado

Ter um local padrão é ótimo, mas você pode querer personalizar onde os resultados são salvos e como eles são organizados.

Por exemplo, você pode querer organizar suas saídas em subdiretórios.
A maneira mais simples de fazer isso é atribuir um caminho de saída específico por saída.

#### 2.2.1. Modifique o caminho de saída

Mais uma vez, modificar o comportamento de publicação para uma saída específica é realmente direto.
Para definir um local personalizado, basta editar o `path` de acordo:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Como isso é definido no nível da saída individual, você pode especificar diferentes locais e subdiretórios para atender às suas necessidades.

#### 2.2.2. Execute o fluxo de trabalho novamente

Vamos experimentar.

```bash
nextflow run hello-world.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Desta vez o resultado é escrito no subdiretório especificado.

??? abstract "Conteúdo do diretório"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Você vê que o resultado da execução anterior ainda está lá.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Você pode usar quantos níveis de aninhamento desejar.
Também é possível usar o nome do processo ou outras variáveis para nomear os diretórios usados para organizar resultados, e é possível alterar o nome padrão do diretório de saída de nível superior (que é controlado pela flag CLI `-o` ou variável de configuração `outputDir`).
Cobriremos essas opções mais tarde no treinamento.

### 2.3. Defina o modo de publicação para copiar

Por padrão, as saídas são publicadas como links simbólicos do diretório `work`.
Isso significa que há apenas um único arquivo no sistema de arquivos.

Isso é ótimo quando você está lidando com arquivos muito grandes, para os quais você não quer armazenar várias cópias.
No entanto, se você excluir o diretório work em algum momento (abordaremos operações de limpeza em breve), você perderá o acesso ao arquivo.
Portanto, você precisa ter um plano para salvar cópias de quaisquer arquivos importantes em um local seguro.

Uma opção fácil é mudar o modo de publicação para copiar para as saídas que você se importa.

#### 2.3.1. Adicione a diretiva mode

Esta parte é realmente direta.
Basta adicionar `mode 'copy'` à definição de saída relevante em nível de fluxo de trabalho:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Isso define o modo de publicação para essa saída específica.

#### 2.3.2. Execute o fluxo de trabalho novamente

Vamos experimentar.

```bash
nextflow run hello-world.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Desta vez, se você olhar os resultados, o arquivo é uma cópia adequada em vez de apenas um link simbólico.

??? abstract "Conteúdo do diretório"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Como isso também é definido no nível da saída individual, permite que você defina o modo de publicação de forma granular.
Isso será especialmente útil mais tarde quando passarmos para pipelines de múltiplas etapas, onde você pode querer copiar apenas as saídas finais e deixar saídas intermediárias como links simbólicos, por exemplo.

Como observado anteriormente, existem outras opções mais sofisticadas para controlar como as saídas são publicadas.
Mostraremos como usá-las no devido tempo em sua jornada com o Nextflow.

### 2.4. Nota sobre diretivas `publishDir` em nível de processo

Até muito recentemente, a forma estabelecida de publicar saídas era fazê-lo no nível de cada processo individual usando uma diretiva `publishDir`.

Para conseguir o que acabamos de fazer para as saídas do processo `sayHello`, teríamos adicionado a seguinte linha à definição do processo:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Você ainda encontrará este padrão de código por toda parte em pipelines e módulos de processo Nextflow mais antigos, por isso é importante estar ciente dele.
No entanto, não recomendamos usá-lo em qualquer trabalho novo, pois eventualmente será proibido em versões futuras da linguagem Nextflow.

### Conclusão

Você sabe como publicar saídas de fluxo de trabalho em um local mais conveniente.

### O que vem a seguir?

Aprenda a fornecer uma entrada variável via parâmetro de linha de comando e utilizar valores padrão de forma eficaz.

---

## 3. Use uma entrada variável passada na linha de comando

Em seu estado atual, nosso fluxo de trabalho usa uma saudação codificada no comando do processo.
Queremos adicionar alguma flexibilidade usando uma variável de entrada, para que possamos mudar mais facilmente a saudação em tempo de execução.

Isso requer que façamos três conjuntos de mudanças em nosso script:

1. Alterar o processo para esperar uma entrada variável
2. Configurar um parâmetro de linha de comando para capturar a entrada do usuário
3. Passar a entrada para o processo no corpo do fluxo de trabalho

Vamos fazer essas mudanças uma de cada vez.

### 3.1. Altere o processo `sayHello` para esperar uma entrada variável

Precisamos editar a definição do processo para (1) aceitar uma variável de entrada e (2) usar essa variável na linha de comando.

#### 3.1.1. Adicione um bloco de entrada à definição do processo

Primeiro, vamos adaptar a definição do processo para aceitar uma entrada chamada `greeting`.

No bloco do processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

A variável `greeting` é prefixada por `val` para dizer ao Nextflow que é um valor (não um caminho).

#### 3.1.2. Edite o comando do processo para usar a variável de entrada

Agora trocamos o valor codificado original pelo valor da variável de entrada que esperamos receber.

No bloco do processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

O símbolo `$` e as chaves (`{ }`) dizem ao Nextflow que este é um nome de variável que precisa ser substituído pelo valor de entrada real (=interpolado).

!!! tip "Dica"

    As chaves (`{ }`) eram tecnicamente opcionais em versões anteriores do Nextflow, então você pode ver fluxos de trabalho mais antigos onde isso está escrito como `echo '$greeting' > output.txt`.

Agora que o processo `sayHello()` está pronto para aceitar uma entrada variável, precisamos de uma maneira de fornecer um valor de entrada para a chamada do processo em nível de fluxo de trabalho.

### 3.2. Configure um parâmetro de linha de comando para capturar a entrada do usuário

Poderíamos simplesmente codificar uma entrada diretamente fazendo a chamada do processo `sayHello('Hello World!')`.
No entanto, quando estivermos fazendo trabalho real com nosso fluxo de trabalho, vamos querer ser capazes de controlar suas entradas a partir da linha de comando, para que possamos fazer algo assim:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Felizmente, o Nextflow tem um sistema de parâmetros de fluxo de trabalho integrado chamado [`params`](https://nextflow.io/docs/latest/config.html#params), que facilita declarar e usar parâmetros CLI.

A sintaxe geral é declarar `params.<nome_do_parâmetro>` para dizer ao Nextflow para esperar um parâmetro `--<nome_do_parâmetro>` na linha de comando.

Aqui, queremos criar um parâmetro chamado `--input`, então precisamos declarar `params.input` em algum lugar no fluxo de trabalho.
Em princípio podemos escrevê-lo em qualquer lugar; mas como vamos querer passá-lo para a chamada do processo `sayHello()`, podemos conectá-lo lá diretamente escrevendo `sayHello(params.input)`.

No bloco do fluxo de trabalho, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emite uma saudação
    sayHello(params.input)
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emite uma saudação
    sayHello()
    ```

Isso diz ao Nextflow para executar o processo `sayHello` no valor fornecido através do parâmetro `--input`.

Na prática, realizamos as etapas (2) e (3) descritas no início da seção de uma só vez.

### 3.3. Execute o comando do fluxo de trabalho

Vamos executá-lo!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Se você fez todas essas edições corretamente, deve obter outra execução bem-sucedida.

Certifique-se de abrir o arquivo de saída para verificar se você agora tem a nova versão da saudação.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà!

Note como a nova execução sobrescreveu o arquivo de saída publicado no diretório `results`.
No entanto, os resultados das execuções anteriores ainda estão preservados nos diretórios de tarefa em `work`.

!!! tip "Dica"

    Você pode distinguir facilmente parâmetros em nível de Nextflow de parâmetros em nível de pipeline.

    - Parâmetros que se aplicam a um pipeline sempre levam um hífen duplo (`--`).
    - Parâmetros que modificam uma configuração do Nextflow, _por exemplo_ o recurso `-resume` que usamos anteriormente, levam um hífen simples (`-`).

### 3.4. Use valores padrão para parâmetros de linha de comando

Ok, isso foi conveniente, mas em muitos casos, faz sentido fornecer um valor padrão para um determinado parâmetro para que você não precise especificá-lo para cada execução.

#### 3.4.1. Defina um valor padrão para o parâmetro CLI

Vamos dar ao parâmetro `input` um valor padrão declarando-o antes da definição do fluxo de trabalho.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String = 'Holà mundo!'
}
```

Como você vê, podemos especificar o tipo de entrada que o fluxo de trabalho espera (Nextflow 25.10.2 e posterior).
A sintaxe é `nome: Tipo = valor_padrão`.
Os tipos suportados incluem `String`, `Integer`, `Float`, `Boolean` e `Path`.

!!! info

    Em fluxos de trabalho mais antigos, você pode ver que todo o bloco `params` está escrito como apenas `input = 'Holà mundo!'`.

À medida que você adiciona mais parâmetros ao seu pipeline, deve adicioná-los todos a este bloco, quer você precise ou não dar a eles um valor padrão.
Isso facilitará encontrar todos os parâmetros configuráveis de relance.

#### 3.4.2. Execute o fluxo de trabalho novamente sem especificar o parâmetro

Agora que você tem um valor padrão definido, pode executar o fluxo de trabalho novamente sem ter que especificar um valor na linha de comando.

```bash
nextflow run hello-world.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

A saída estará no mesmo lugar que anteriormente, mas o conteúdo deve ser atualizado com o novo texto.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

O Nextflow usou o valor padrão do parâmetro de saudação para criar a saída.

#### 3.4.3. Substitua o valor padrão

Se você fornecer o parâmetro na linha de comando, o valor CLI substituirá o valor padrão.

Experimente:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Mais uma vez, você deve encontrar a saída atualizada correspondente no seu diretório de resultados.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Nota"

    No Nextflow, há vários lugares onde você pode especificar valores para parâmetros.
    Se o mesmo parâmetro for definido com valores diferentes em vários lugares, o Nextflow determinará qual valor usar com base na ordem de precedência descrita [aqui](https://www.nextflow.io/docs/latest/config.html).

    Cobriremos isso com mais detalhes na Parte 6 (Configuração).

### Conclusão

Você sabe como usar uma entrada variável simples fornecida em tempo de execução via parâmetro de linha de comando, bem como configurar, usar e substituir valores padrão.

### O que vem a seguir?

Aprenda como gerenciar execuções de forma mais conveniente.

---

## 4. Gerencie execuções de fluxo de trabalho

Saber como lançar fluxos de trabalho e recuperar saídas é ótimo, mas você rapidamente descobrirá que há alguns outros aspectos do gerenciamento de fluxo de trabalho que tornarão sua vida mais fácil, especialmente se você estiver desenvolvendo seus próprios fluxos de trabalho.

Aqui mostramos como usar o recurso [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) para quando você precisar relançar o mesmo fluxo de trabalho, como inspecionar o log de execuções passadas com [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), e como excluir diretórios work mais antigos com [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

<!-- Any other cool options we should include? Added log -->

### 4.1. Relance um fluxo de trabalho com `-resume`

Às vezes, você vai querer executar novamente um pipeline que já lançou anteriormente sem refazer nenhuma etapa que já foi concluída com sucesso.

O Nextflow tem uma opção chamada [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) que permite fazer isso.
Especificamente, neste modo, quaisquer processos que já foram executados com exatamente o mesmo código, configurações e entradas serão ignorados.
Isso significa que o Nextflow só executará processos que você adicionou ou modificou desde a última execução, ou aos quais você está fornecendo novas configurações ou entradas.

Existem duas vantagens principais em fazer isso:

- Se você está no meio do desenvolvimento do seu pipeline, pode iterar mais rapidamente, pois só precisa executar o(s) processo(s) em que está trabalhando ativamente para testar suas alterações.
- Se você está executando um pipeline em produção e algo dá errado, em muitos casos você pode corrigir o problema e relançar o pipeline, e ele retomará a execução do ponto de falha, o que pode economizar muito tempo e computação.

Para usá-lo, basta adicionar `-resume` ao seu comando e executá-lo:

```bash
nextflow run hello-world.nf -resume
```

??? success "Saída do comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

A saída do console deve parecer familiar, mas há uma coisa que é um pouco diferente em comparação com antes.

Procure pela parte `cached:` que foi adicionada na linha de status do processo (linha 5), o que significa que o Nextflow reconheceu que já fez este trabalho e simplesmente reutilizou o resultado da execução anterior bem-sucedida.

Você também pode ver que o hash do subdiretório work é o mesmo da execução anterior.
O Nextflow está literalmente apontando para a execução anterior e dizendo "Eu já fiz isso lá."

!!! tip "Dica"

    Quando você executa novamente um pipeline com `resume`, o Nextflow não sobrescreve nenhum arquivo publicado fora do diretório work por quaisquer execuções que foram executadas com sucesso anteriormente.

### 4.2. Inspecione o log de execuções passadas

Quer você esteja desenvolvendo um novo pipeline ou executando pipelines em produção, em algum momento você provavelmente precisará consultar informações sobre execuções passadas.
Aqui está como fazer isso.

Sempre que você lança um fluxo de trabalho nextflow, uma linha é escrita em um arquivo de log chamado `history`, em um diretório oculto chamado `.nextflow` no diretório de trabalho atual.

??? abstract "Conteúdo do arquivo"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Este arquivo fornece o timestamp, nome da execução, status, ID de revisão, ID de sessão e linha de comando completa para cada execução do Nextflow que foi lançada a partir do diretório de trabalho atual.

Uma maneira mais conveniente de acessar essas informações é usar o comando `nextflow log`.

```bash
nextflow log
```

??? success "Saída do comando"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Isso exibirá o conteúdo do arquivo de log no terminal, aumentado com uma linha de cabeçalho.

Você notará que o ID de sessão muda sempre que você executa um novo comando `nextflow run`, EXCETO se você estiver usando a opção `-resume`.
Nesse caso, o ID de sessão permanece o mesmo.

O Nextflow usa o ID de sessão para agrupar informações de cache de execução no diretório `cache`, também localizado em `.nextflow`.

### 4.3. Exclua diretórios work mais antigos

Durante o processo de desenvolvimento, você normalmente executará seu rascunho de pipeline um grande número de vezes, o que pode levar a um acúmulo de muitos arquivos em muitos subdiretórios.

Felizmente, o Nextflow inclui um subcomando útil `clean` que pode excluir automaticamente os subdiretórios work de execuções passadas que você não se importa mais.

#### 4.3.1. Determine os critérios de exclusão

Existem várias [opções](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para determinar o que excluir.

Aqui mostramos um exemplo que exclui todos os subdiretórios de execuções antes de uma determinada execução, especificada usando seu nome de execução.

Procure a execução bem-sucedida mais recente onde você não usou `-resume`; no nosso caso, o nome da execução foi `golden_cantor`.

O nome da execução é a string de duas partes gerada por máquina mostrada entre colchetes na linha de saída do console `Launching (...)`.
Você também pode usar o log do Nextflow para procurar uma execução com base em seu timestamp e/ou linha de comando.

#### 4.3.2. Faça uma execução de teste

Primeiro usamos a flag de execução de teste `-n` para verificar o que será excluído dado o comando:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Saída do comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Sua saída terá nomes de diretório de tarefa diferentes e pode ter um número diferente de linhas, mas deve parecer semelhante ao exemplo.

Se você não ver nenhuma linha de saída, você não forneceu um nome de execução válido ou não há execuções passadas para excluir. Certifique-se de alterar `golden_cantor` no comando de exemplo para qualquer que seja o nome de execução mais recente correspondente no seu log.

#### 4.3.3. Prossiga com a exclusão

Se a saída parecer como esperado e você quiser prosseguir com a exclusão, execute novamente o comando com a flag `-f` em vez de `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Saída do comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

A saída deve ser semelhante à anterior, mas agora dizendo 'Removed' em vez de 'Would remove'.
Note que isso não remove os subdiretórios de dois caracteres (como `a3/` acima), mas esvazia seu conteúdo.

!!! Warning "Aviso"

    Excluir subdiretórios work de execuções passadas os remove do cache do Nextflow e exclui quaisquer saídas que foram armazenadas nesses diretórios.
    Isso significa que quebra a capacidade do Nextflow de retomar a execução sem executar novamente os processos correspondentes.

    Você é responsável por salvar quaisquer saídas que você se importa ou planeja confiar! Essa é a principal razão pela qual preferimos usar o modo `copy` em vez do modo `symlink` para a diretiva `publish`.

### Conclusão

Você sabe como publicar saídas em um diretório específico, relançar um pipeline sem repetir etapas que já foram executadas de forma idêntica e usar o comando `nextflow clean` para limpar diretórios work antigos.

De forma mais geral, você sabe como interpretar um fluxo de trabalho Nextflow simples, gerenciar sua execução e recuperar saídas.

### O que vem a seguir?

Faça uma pequena pausa, você mereceu!

Quando estiver pronto, passe para [**Parte 2: Hello Channels**](./02_hello_channels.md) para aprender como usar canais para alimentar entradas em seu fluxo de trabalho, o que permitirá que você aproveite o paralelismo de fluxo de dados integrado do Nextflow e outros recursos poderosos.

---

## Quiz

<quiz>
Quais são os componentes mínimos necessários de um processo Nextflow?
- [ ] Apenas blocos de entrada e saída
- [x] Blocos de saída e script
- [ ] Blocos de entrada, saída e script
- [ ] Apenas um bloco script

Saiba mais: [1.1.1. A definição de `process`](#111-a-definição-de-process)
</quiz>

<quiz>
Qual é o propósito do bloco de saída em um processo?
- [ ] Imprimir resultados no console
- [ ] Salvar arquivos no diretório work
- [x] Declarar saídas esperadas do processo
- [ ] Definir variáveis de ambiente

Saiba mais: [1.1.1. A definição de `process`](#111-a-definição-de-process)
</quiz>

<quiz>
Qual comando é usado para executar um fluxo de trabalho Nextflow?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Olhando para o diretório work de uma tarefa, qual arquivo contém o comando real que foi executado?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Saiba mais: [1.2.2. Encontre a saída e os logs no diretório `work`](#122-encontre-a-saída-e-os-logs-no-diretório-work)
</quiz>

<quiz>
O que a flag `-resume` faz?
- [ ] Reinicia o fluxo de trabalho do início
- [ ] Pausa o fluxo de trabalho
- [x] Ignora processos que já foram concluídos com sucesso
- [ ] Cria um backup do fluxo de trabalho

Saiba mais: [4.1. Relance um fluxo de trabalho com `-resume`](#41-relance-um-fluxo-de-trabalho-com--resume)
</quiz>

<quiz>
Qual é o modo padrão para publicar saídas de fluxo de trabalho?
- [ ] Copiar arquivos para o diretório de saída
- [x] Criar links simbólicos no diretório de saída
- [ ] Mover arquivos para o diretório de saída
- [ ] Comprimir arquivos no diretório de saída

Saiba mais: [2.3. Defina o modo de publicação para copiar](#23-defina-o-modo-de-publicação-para-copiar)
</quiz>

<quiz>
Como você passa um valor de parâmetro para um fluxo de trabalho Nextflow a partir da linha de comando?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Saiba mais: [3.2. Configure um parâmetro de linha de comando para capturar a entrada do usuário](#32-configure-um-parâmetro-de-linha-de-comando-para-capturar-a-entrada-do-usuário)
</quiz>

<quiz>
Como você referencia uma variável dentro de um bloco script do Nextflow?
- [ ] Use a sintaxe `%variable%`
- [x] Use a sintaxe `#!groovy ${variable}`
- [ ] Use a sintaxe `{{variable}}`
- [ ] Use a sintaxe `[variable]`
</quiz>
