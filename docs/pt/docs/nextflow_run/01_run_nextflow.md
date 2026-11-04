# Parte 1: Executar o Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta parte, apresentamos os conceitos fundamentais de execução de pipelines Nextflow.
Começamos com um simples fluxo de trabalho Hello World e avançamos até um pipeline completo de múltiplas etapas que processa várias entradas em paralelo usando contêineres.

---

## 1. Hello World

O fluxo de trabalho `1-hello.nf` recebe uma saudação via argumento de linha de comando e a escreve em um arquivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Iniciar o fluxo de trabalho

Execute o seguinte comando no seu terminal.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Saída do comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

A linha principal da saída é a linha de status do processo:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Isso nos informa que o processo `sayHello` foi executado com sucesso uma vez.
O prefixo `[6d/740edd]` é um caminho truncado para o diretório de trabalho da tarefa — mais detalhes sobre isso a seguir.
O bloco `Outputs:` que aparece em seguida lista todos os arquivos publicados pelo pipeline, rotulados de acordo com o bloco `output` abordado em [1.4](#14-optional-code-walkthrough) abaixo.

### 1.2. Encontrar a saída

Este fluxo de trabalho está configurado para publicar sua saída em um diretório `results`.
Após a execução, você deve encontrar a saída lá:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Abra o arquivo para confirmar que ele contém `Hello World!`.

### 1.3. Explorar o diretório `work/`

Nos bastidores, o Nextflow cria um diretório de tarefa único para cada chamada de processo dentro de um diretório chamado `work/`.
O hash exibido na saída do console (`[6d/740edd]`) é o caminho para esse diretório.

```bash
ls work/6d/740edd*
```

Dentro dele você encontrará o arquivo de saída junto com vários arquivos de log ocultos:

- **`.command.sh`**: o comando exato que o Nextflow executou
- **`.command.out`** / **`.command.err`**: stdout e stderr do processo
- **`.command.log`**: saída de log combinada
- **`.exitcode`**: o código de saída do processo

O arquivo `.command.sh` é especialmente útil para depuração — ele mostra exatamente o que foi executado.

### 1.4. Opcional: Exploração do código

Entender o código não é essencial se você apenas quer executar pipelines, mas se você tiver curiosidade, vale a pena dar uma olhada.

??? optional "Clique para explorar o código associado a este exercício"

    Vamos abrir `1-hello.nf` e examinar seus principais componentes.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Parâmetros do pipeline
     */
    params {
        input: String
    }

    workflow {

        main:
        // emite uma saudação
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Podemos observar o seguinte:

    - uma instrução `include` apontando para um módulo de `process`
    - um bloco `params` definindo os parâmetros do pipeline
    - um bloco `workflow` descrevendo o trabalho a ser realizado
    - um bloco `output` descrevendo o que fazer com as saídas

    Vamos examinar cada um deles.

    ### O módulo `process`

    A instrução `include` diz ao Nextflow para carregar algo chamado `sayHello` de um arquivo de código separado.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    Nesse arquivo, encontramos a definição de um processo chamado `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

    Um **process** define uma única etapa no pipeline.
    Ele declara suas entradas, saídas e o script a ser executado.
    O qualificador `val` significa que a entrada é um valor simples (string, número, etc.).
    O qualificador `path` significa que a saída é um caminho de arquivo.

    É possível escrever a definição do processo no arquivo principal do fluxo de trabalho, mas mantê-los em arquivos de módulo separados os torna reutilizáveis: o mesmo módulo pode ser importado por múltiplos scripts de fluxo de trabalho.

    ### O bloco `params`

    O bloco `params` declara os parâmetros de linha de comando que o fluxo de trabalho aceita:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Qualquer parâmetro declarado aqui fica disponível na linha de comando com dois traços (`--input`).
    Os tipos suportados incluem `String`, `Integer`, `Float`, `Boolean` e `Path`.

    !!! tip "Dica"

        Os parâmetros do fluxo de trabalho sempre usam dois traços (`--input`) para diferenciá-los das próprias flags de CLI do Nextflow, que usam um traço (ex.: `-resume`).

    ### O bloco `workflow`

    O bloco **workflow** define a lógica de fluxo de dados: quais processos executar e em que ordem.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // emite uma saudação
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Aqui há apenas um processo sendo chamado, então é bem simples; abordaremos exemplos mais realistas adiante.

    A seção `main:` chama o processo `sayHello` com o valor de `--input`.
    A seção `publish:` lista quais saídas devem ser copiadas para o diretório de resultados.

    ### O bloco `output`

    O bloco `output` no final do arquivo especifica o caminho de destino e o modo de cópia.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Cada entrada nomeada corresponde a um rótulo `publish:` no fluxo de trabalho e o mapeia para um subdiretório dentro de `results/`.

### Conclusão

Você sabe como executar um pipeline Nextflow e encontrar suas saídas, e sabe que o trabalho é executado em diretórios de tarefa dentro de `work/`.

### O que vem a seguir?

Descubra como o Nextflow lida com múltiplas entradas de forma eficiente.

---

## 2. Processar múltiplas entradas

Pipelines do mundo real geralmente processam muitos dados, não apenas um.
O fluxo de trabalho `2-inputs.nf` lê de um arquivo CSV e executa `sayHello` uma vez por linha, em paralelo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Vamos executar o fluxo de trabalho primeiro e depois examinar o mecanismo que o Nextflow usa para lidar com essas múltiplas entradas.

### 2.1. Executar o fluxo de trabalho

Execute o seguinte comando no seu terminal.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Saída do comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

O `3 of 3` nos informa que o processo `sayHello` foi chamado três vezes, uma por linha no CSV.

No diretório `results`, você deve ver agora três arquivos de saída, um por saudação:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Abra qualquer um dos arquivos de saída para confirmar que cada um contém uma saudação.

A saída condensada acima mostra uma única linha de resumo para `sayHello`, mas o Nextflow na verdade iniciou três execuções de tarefa separadas por trás dela, uma por linha no CSV, e as executou em paralelo assim que sua máquina tinha recursos disponíveis.

Assim como a tarefa única que você explorou em [1.3](#13-explore-the-work-directory), cada uma dessas três execuções tem seu próprio diretório de tarefa dentro de `work/`, completamente isolado das demais:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Cada `.command.sh` contém apenas o comando para aquela saudação específica:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Esse isolamento é o que torna a execução paralela segura: três tarefas sendo executadas ao mesmo tempo nunca compartilham um diretório de trabalho, então nada que uma tarefa escreve pode colidir com ou sobrescrever o que outra tarefa escreve, mesmo que produzam arquivos com o mesmo nome.
É também por isso que o `-resume` (abordado a seguir) pode fazer cache e reutilizar tarefas individuais de forma independente: as entradas, saídas e logs de cada tarefa vivem inteiramente dentro de seu próprio diretório, sem nada compartilhado entre tarefas que possa ficar fora de sincronia.

### 2.2. Executar o fluxo de trabalho novamente com `-ansi-log false`

Por padrão, o Nextflow condensa a saída em uma única linha de resumo por processo.
Para ver cada chamada de processo listada individualmente, adicione `-ansi-log false`:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Isso mostra todas as três chamadas de processo e o subdiretório de trabalho único criado para cada uma.

### 2.3. Usar `-resume` para pular o trabalho já concluído

Agora mude para o arquivo de entrada estendido, que adiciona mais duas saudações, e adicione `-resume` à linha de comando:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Saída do comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

O Nextflow executou apenas as duas novas entradas.
As três saudações processadas na execução anterior foram armazenadas em cache e reutilizadas automaticamente.

Isso também funciona para pular a execução de processos em etapas que já foram concluídas com sucesso em um pipeline de múltiplas etapas.
Por exemplo, se a execução de um pipeline foi interrompida por um erro de sistema, ou se você adicionou novas etapas a um pipeline em desenvolvimento.

A capacidade de `-resume` é especialmente valiosa em pipelines longos, onde recuperar-se de uma falha pode economizar tempo e recursos críticos.

### 2.4. Opcional: Exploração do código

Entender o código não é essencial se você apenas quer executar pipelines, mas se você tiver curiosidade, vale a pena dar uma olhada.

??? optional "Clique para explorar o código associado a este exercício"

    A principal mudança em `2-inputs.nf` está na seção `main:` do fluxo de trabalho:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // cria um canal para entradas a partir de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
    ```

    O que você vê aqui é chamado de **canal**: uma estrutura de fila que lida com dados de entrada de uma forma que facilita a paralelização de operações.

    - `channel.fromPath(params.input)` cria um canal a partir do caminho de arquivo fornecido com `--input`
    - `.splitCsv()` analisa o CSV em linhas
    - `#!groovy .map { line -> line[0] }` extrai a primeira coluna de cada linha

    O resultado é um canal contendo `Hello`, `Bonjour` e `Hola`.
    Quando passado para `sayHello(greeting_ch)`, o Nextflow automaticamente chama o processo uma vez por item, executando-os em paralelo quando os recursos permitem.

### Conclusão

Você sabe como processar múltiplas entradas de um arquivo CSV em paralelo e como usar `-resume` para evitar repetir o trabalho já concluído.

### O que vem a seguir?

Aprenda como um pipeline completo de múltiplas etapas encadeia processos usando canais e como usar contêineres para gerenciar ferramentas de análise e suas dependências.

---

## 3. Executar um pipeline de múltiplas etapas

Até agora você executou um único processo e depois o executou várias vezes em paralelo sobre um conjunto de entradas.
Pipelines reais geralmente vão além: eles encadeiam vários processos, alimentando a saída de um na entrada do próximo, e frequentemente dependem de mais de um software ao longo do caminho.
O fluxo de trabalho `main.nf` combina ambos em um pipeline completo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Cada saudação de entrada passa por todas as quatro etapas: `sayHello` a escreve em um arquivo, `convertToUpper` converte o texto para maiúsculas, `collectGreetings` mescla todos os resultados em um único arquivo, e `cowpy` gera arte ASCII a partir da saída mesclada usando uma ferramenta em contêiner.
O Nextflow conecta essas etapas com canais: a saída de um processo se torna a entrada do próximo, então toda a cadeia é executada automaticamente conforme os dados ficam disponíveis, sem que você precise orquestrar cada etapa manualmente.

Observe que este fluxo de trabalho usa módulos: cada processo é definido em seu próprio arquivo dentro de `modules/`, e `main.nf` os importa com instruções `include` em vez de defini-los diretamente.
Isso torna cada processo reutilizável em múltiplos fluxos de trabalho sem duplicar código. Para saber mais, consulte a seção de exploração de código mais abaixo.

### 3.1. Executar o fluxo de trabalho

Execute o seguinte comando no seu terminal.

```bash
nextflow run main.nf --input data/greetings.csv
```

O parâmetro `character` tem como padrão `turkey` em `nextflow.config`, então a arte ASCII usa um peru a menos que você o substitua (tente adicionar `--character tux`).

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Quatro processos foram executados, mas não o mesmo número de vezes.
`sayHello` e `convertToUpper` foram executados uma vez por entrada (3 of 3): cada saudação precisa ser escrita e convertida para maiúsculas individualmente.
`collectGreetings` e `cowpy` foram executados apenas uma vez (1 of 1): mesclar as saudações e gerar a arte ASCII só faz sentido depois que todos os resultados individuais estão prontos.
Esse formato de expansão-e-contração (fan-out-then-fan-in), em que várias tarefas paralelas alimentam um número menor de tarefas posteriores, é comum em pipelines reais.

O Nextflow não espera que uma etapa inteira termine antes de iniciar a próxima.
Assim que uma saída de `sayHello` está pronta, a tarefa correspondente de `convertToUpper` pode começar, então tarefas de diferentes processos são executadas de forma concorrente em vez de em lotes estritamente sequenciais.
`collectGreetings` e `cowpy` precisam esperar, pois cada um deles depende de todos os resultados anteriores estarem disponíveis primeiro.

O diretório `results` reflete essa contração, além do que o autor do pipeline escolheu publicar e onde: lembre-se do bloco `output` da exploração de código em 1.4, que é o que define essa estrutura.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

O diretório de nível superior recebe o nome do parâmetro `batch`, que tem como padrão `batch`; você verá ele mudar em exercícios posteriores.

Verifique `cowpy-COLLECTED-batch-output.txt` para ver o arquivo de arte ASCII.

??? abstract "Conteúdo do arquivo"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
     ---------
      \                                  ,+*^^*+___+++_
       \                           ,*^^^^              )
        \                       _+*                     ^**+_
         \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                 {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
               {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
               U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
             (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
               (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                 ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                         ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Assim como em [2.1](#21-run-the-workflow), cada uma dessas 8 execuções de tarefa, em todos os quatro processos, tem seu próprio diretório dentro de `work/`, completamente isolado das demais.
`collectGreetings` é um bom exemplo de por que isso importa: ele depende das saídas de todas as três tarefas de `convertToUpper`, que vivem em três diretórios de tarefa diferentes, então o Nextflow cria links simbólicos para esses arquivos dentro do próprio diretório de `collectGreetings` em vez de fazê-lo ler diretamente dos diretórios das tarefas anteriores:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Cada tarefa vê apenas os arquivos específicos de que precisa, independentemente de onde vieram, e nunca o conteúdo interno do diretório de outra tarefa.
Em todo um pipeline, esse mesmo isolamento que você viu com um único processo em [2.1](#21-run-the-workflow) é o que permite ao Nextflow executar cada tarefa de cada processo de forma concorrente e segura.

!!! note "Nota"

    A etapa `cowpy` é executada dentro de um contêiner Docker em vez de depender de software instalado localmente.
    Um contêiner empacota uma aplicação junto com tudo o que ela precisa para ser executada, então você não precisa instalar e gerenciar dependências por conta própria, e o pipeline se comporta da mesma forma em qualquer máquina que possa executar o contêiner.
    O Nextflow também suporta Conda como alternativa aos contêineres; consulte a [Parte 2](./02_configure_pipeline.md) para saber como alternar entre eles.

### 3.2. Opcional: Exploração do código

Entender o código não é essencial se você apenas quer executar pipelines, mas se você tiver curiosidade, vale a pena dar uma olhada.

??? optional "Clique para explorar o código associado a este exercício"

    ### Como os dados fluem de uma etapa para a próxima

    Cada processo passa seu canal de saída para o próximo:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // cria um canal para entradas a partir de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    O padrão `processName.out` se refere ao canal de saída de um processo.

    O operador `.collect()` reúne todas as saídas individuais de `convertToUpper` em um único item de canal antes de passá-las para `collectGreetings`.

    ### Usando módulos de processo

    `main.nf` não define nenhum código de processo diretamente.
    Em vez disso, ele importa cada processo de seu próprio arquivo dentro de `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Cada arquivo de módulo contém uma única definição de processo, estruturada da mesma forma que o módulo `sayHello` em [1.4](#14-optional-code-walkthrough).
    Manter os processos em arquivos separados os torna reutilizáveis em múltiplos fluxos de trabalho sem duplicar código.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Usando software em contêiner

    O processo `cowpy` é executado dentro de um contêiner Docker especificado em seu arquivo de módulo:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

    O Nextflow automaticamente baixa a imagem, executa o script dentro do contêiner e faz a limpeza depois.
    O Docker está habilitado para este projeto em `nextflow.config`:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Essa única linha habilita o Docker para qualquer processo no pipeline que tenha um contêiner especificado.

### Conclusão

Você executou um pipeline completo de múltiplas etapas que processa múltiplas entradas em paralelo usando uma ferramenta em contêiner.

### O que vem a seguir?

Siga para a [Parte 2](./02_configure_pipeline.md), onde você aprenderá como configurar o comportamento do pipeline usando `nextflow.config`.

---

## Resumo

Nesta parte você aprendeu a:

- Executar um fluxo de trabalho Nextflow e encontrar suas saídas
- Explorar o diretório `work/` e seus arquivos de log
- Processar múltiplas entradas de um arquivo CSV em paralelo
- Usar `-resume` para pular o trabalho já concluído ao adicionar novas entradas
- Executar um pipeline de múltiplas etapas que usa uma ferramenta em contêiner
