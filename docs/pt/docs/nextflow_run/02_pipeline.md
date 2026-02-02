# Parte 2: Executar pipelines reais

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na Parte 1 deste curso (Operações Básicas de Execução), começamos com um fluxo de trabalho de exemplo que tinha apenas recursos mínimos para manter a complexidade do código baixa.
Por exemplo, `1-hello.nf` usou um parâmetro de linha de comando (`--input`) para fornecer um único valor por vez.

No entanto, a maioria dos pipelines do mundo real usa recursos mais sofisticados para permitir o processamento eficiente de grandes quantidades de dados em escala e aplicar múltiplas etapas de processamento encadeadas por lógica às vezes complexa.

Nesta parte do treinamento, demonstramos recursos-chave de pipelines do mundo real experimentando versões expandidas do pipeline Hello World original.

## 1. Processando dados de entrada de um arquivo

Em um pipeline do mundo real, tipicamente queremos processar múltiplos pontos de dados (ou séries de dados) contidos em um ou mais arquivos de entrada.
E onde possível, queremos executar o processamento de dados independentes em paralelo, para encurtar o tempo gasto esperando pela análise.

Para demonstrar como o Nextflow faz isso, preparamos um arquivo CSV chamado `greetings.csv` que contém várias saudações de entrada, imitando o tipo de dados colunares que você pode querer processar em uma análise de dados real.
Note que os números não são significativos, eles estão lá apenas para fins ilustrativos.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Também escrevemos uma versão melhorada do fluxo de trabalho original, agora chamada `2a-inputs.nf`, que lerá o arquivo CSV, extrairá as saudações e escreverá cada uma delas em um arquivo separado.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Vamos executar o fluxo de trabalho primeiro, e depois olharemos o código Nextflow relevante.

### 1.1. Execute o fluxo de trabalho

Execute o seguinte comando no seu terminal.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Empolgante, isso parece indicar que '3 of 3' chamadas foram feitas para o processo, o que é encorajador, já que havia três linhas de dados no CSV que fornecemos como entrada.
Isso sugere que o processo `sayHello()` foi chamado três vezes, uma vez em cada linha de entrada.

### 1.2. Encontre as saídas publicadas no diretório `results`

Vamos olhar o diretório 'results' para ver se nosso fluxo de trabalho ainda está escrevendo uma cópia de nossas saídas lá.

??? abstract "Conteúdo do diretório"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Sim! Vemos um novo diretório chamado `2a-inputs` com três arquivos de saída com nomes diferentes, convenientemente.

Você pode abrir cada um deles para se satisfazer de que contêm a string de saudação apropriada.

??? abstract "Conteúdo dos arquivos"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Isso confirma que cada saudação no arquivo de entrada foi processada apropriadamente.

### 1.3. Encontre as saídas originais e os logs

Você pode ter notado que a saída do console acima se referiu a apenas um diretório de tarefa.
Isso significa que todas as três chamadas a `sayHello()` foram executadas dentro daquele único diretório de tarefa?

#### 1.3.1. Examine o diretório de tarefa dado no terminal

Vamos dar uma olhada dentro daquele diretório de tarefa `8e/0eb066`.

??? abstract "Conteúdo do diretório"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Encontramos apenas a saída correspondente a uma das saudações (bem como os arquivos acessórios se habilitarmos a exibição de arquivos ocultos).

Então o que está acontecendo aqui?

Por padrão, o sistema de logging ANSI escreve as informações de status para todas as chamadas ao mesmo processo na mesma linha.
Como resultado, ele nos mostrou apenas um dos três caminhos de diretório de tarefa (`8e/0eb066`) na saída do console.
Há outros dois que não estão listados lá.

#### 1.3.2. Faça o terminal mostrar mais detalhes

Podemos modificar o comportamento de logging para ver a lista completa de chamadas de processo adicionando o `-ansi-log false` ao comando da seguinte forma:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Saída do comando"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Desta vez vemos todas as três execuções de processo e seus subdiretórios de trabalho associados listados na saída.
Desabilitar o logging ANSI também impediu o Nextflow de usar cores na saída do terminal.

Note que a forma como o status é reportado é um pouco diferente entre os dois modos de logging.
No modo condensado, o Nextflow reporta se as chamadas foram completadas com sucesso ou não.
Neste modo expandido, ele apenas reporta que foram submetidas.

Isso confirma que o processo `sayHello()` é chamado três vezes, e um diretório de tarefa separado é criado para cada um.

Se olharmos dentro de cada um dos diretórios de tarefa listados lá, podemos verificar que cada um corresponde a uma das saudações.

??? abstract "Conteúdo do diretório"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Isso confirma que cada chamada de processo é executada em isolamento de todas as outras.
Isso tem muitas vantagens, incluindo evitar colisões se o processo produzir quaisquer arquivos intermediários com nomes não-únicos.

!!! tip "Dica"

    Para um fluxo de trabalho complexo, ou um grande número de entradas, ter a lista completa exibida no terminal pode ficar um pouco avassalador, então as pessoas normalmente não usam `-ansi-log false` no uso rotineiro.

### 1.4. Examine o código do fluxo de trabalho

Então esta versão do fluxo de trabalho é capaz de ler um arquivo CSV de entradas, processar as entradas separadamente e nomear as saídas de forma única.

Vamos dar uma olhada no que torna isso possível no código do fluxo de trabalho.

??? full-code "Arquivo de código completo"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Usa echo para imprimir 'Hello World!' em um arquivo
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Novamente, você não precisa memorizar sintaxe de código, mas é bom aprender a reconhecer componentes-chave do fluxo de trabalho que fornecem funcionalidade importante.

#### 1.4.1. Carregando os dados de entrada do CSV

Esta é a parte mais interessante: como mudamos de receber um único valor da linha de comando para receber um arquivo CSV, analisá-lo e processar as saudações individuais que ele contém?

No Nextflow, fazemos isso com um **canal**: uma construção projetada para lidar com entradas eficientemente e transportá-las de uma etapa para outra em fluxos de trabalho de múltiplas etapas, enquanto fornece paralelismo embutido e muitos benefícios adicionais.

Vamos analisar.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // cria um canal para entradas de um arquivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emite uma saudação
    sayHello(greeting_ch)
```

Este código cria um canal chamado `greeting_ch` que lê o arquivo CSV, analisa-o e extrai a primeira coluna de cada linha.
O resultado é um canal contendo `Hello`, `Bonjour` e `Holà`.

??? tip "Como isso funciona?"

    Aqui está o que essa linha significa em português simples:

    - `channel.fromPath` é uma **fábrica de canal** que cria um canal a partir de caminho(s) de arquivo
    - `(params.input)` especifica que o caminho do arquivo é fornecido por `--input` na linha de comando

    Em outras palavras, essa linha diz ao Nextflow: pegue o caminho do arquivo dado com `--input` e prepare-se para tratar seu conteúdo como dados de entrada.

    Então as próximas duas linhas aplicam **operadores** que fazem a análise real do arquivo e o carregamento dos dados na estrutura de dados apropriada:

    - `.splitCsv()` diz ao Nextflow para analisar o arquivo CSV em um array representando linhas e colunas
    - `.map { line -> line[0] }` diz ao Nextflow para pegar apenas o elemento na primeira coluna de cada linha

    Então na prática, começando do seguinte arquivo CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Transformamos isso em um array que se parece com isso:

    ```txt title="Conteúdo do array"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    E então pegamos o primeiro elemento de cada uma das três linhas e os carregamos em um canal Nextflow que agora contém: `Hello`, `Bonjour` e `Holà`.

    Se você quiser entender canais e operadores em profundidade, incluindo como escrevê-los você mesmo, veja [Hello Nextflow Parte 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Chame o processo em cada saudação

A seguir, na última linha do bloco `main:` do fluxo de trabalho, fornecemos o canal `greeting_ch` carregado como entrada para o processo `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // cria um canal para entradas de um arquivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emite uma saudação
    sayHello(greeting_ch)
```

Isso diz ao Nextflow para executar o processo individualmente em cada elemento no canal, _ou seja_, em cada saudação.
E porque o Nextflow é inteligente assim, ele executará essas chamadas de processo em paralelo se possível, dependendo da infraestrutura computacional disponível.

É assim que você pode alcançar processamento eficiente e escalável de muitos dados (muitas amostras, ou pontos de dados, seja qual for sua unidade de pesquisa) com comparativamente muito pouco código.

#### 1.4.3. Como as saídas são nomeadas

Finalmente, vale a pena dar uma olhada rápida no código do processo para ver como fazemos os arquivos de saída serem nomeados de forma única.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

Você vê que, comparado à versão deste processo em `1-hello.nf`, a declaração de saída e a parte relevante do comando mudaram para incluir o valor da saudação no nome do arquivo de saída.

Esta é uma forma de garantir que os nomes dos arquivos de saída não colidirão quando forem publicados no diretório de resultados comum.

E essa é a única mudança que tivemos que fazer dentro da declaração do processo!

### Conclusão

Você entende em um nível básico como canais e operadores nos permitem processar múltiplas entradas eficientemente.

### O que vem a seguir?

Descubra como fluxos de trabalho de múltiplas etapas são construídos e como eles operam.

---

## 2. Executando fluxos de trabalho de múltiplas etapas

A maioria dos fluxos de trabalho do mundo real envolve mais de uma etapa.
Vamos construir sobre o que acabamos de aprender sobre canais, e olhar como o Nextflow usa canais e operadores para conectar processos em um fluxo de trabalho de múltiplas etapas.

Para isso, fornecemos a você um fluxo de trabalho de exemplo que encadeia três etapas separadas e demonstra o seguinte:

1. Fazer dados fluírem de um processo para o próximo
2. Coletar saídas de múltiplas chamadas de processo em uma única chamada de processo

Especificamente, fizemos uma versão expandida do fluxo de trabalho chamada `2b-multistep.nf` que pega cada saudação de entrada, converte-a para maiúsculas, depois coleta todas as saudações em maiúsculas em um único arquivo de saída.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Como anteriormente, executaremos o fluxo de trabalho primeiro e depois olharemos o código para ver o que é novo.

### 2.1. Execute o fluxo de trabalho

Execute o seguinte comando no seu terminal:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Saída do comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Você vê que como prometido, múltiplas etapas foram executadas como parte do fluxo de trabalho; as duas primeiras (`sayHello` e `convertToUpper`) foram presumivelmente executadas em cada saudação individual, e a terceira (`collectGreetings`) terá sido executada apenas uma vez, nas saídas de todas as três chamadas `convertToUpper`.

### 2.2. Encontre as saídas

Vamos verificar que isso é de fato o que aconteceu olhando no diretório `results`.

??? abstract "Conteúdo do diretório"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Como você pode ver, temos um novo diretório chamado `2b-multistep`, e ele contém bem mais arquivos do que antes.
Alguns dos arquivos foram agrupados em um subdiretório chamado `intermediates`, enquanto dois arquivos estão localizados no nível superior.

Esses dois são os resultados finais do fluxo de trabalho de múltiplas etapas.
Tire um minuto para olhar os nomes dos arquivos e verificar seu conteúdo para confirmar que são o que você espera.

??? abstract "Conteúdo dos arquivos"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

O primeiro contém nossas três saudações, em maiúsculas e coletadas de volta em um único arquivo como prometido.
O segundo é um arquivo de relatório que resume algumas informações sobre a execução.

### 2.3. Examine o código

Vamos olhar o código e identificar os padrões-chave para fluxos de trabalho de múltiplas etapas.

??? full-code "Arquivo de código completo"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Usa echo para imprimir 'Hello World!' em um arquivo
    */
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

    /*
    * Usa uma ferramenta de substituição de texto para converter a saudação para maiúsculas
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Coleta saudações em maiúsculas em um único arquivo de saída
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Há muita coisa acontecendo ali, mas a diferença mais óbvia comparada à versão anterior do fluxo de trabalho é que agora há múltiplas definições de processo, e correspondentemente, várias chamadas de processo no bloco workflow.

Vamos dar uma olhada mais de perto e ver se conseguimos identificar as peças mais interessantes.

#### 2.3.1. Visualizando a estrutura do fluxo de trabalho

Se você está usando VSCode com a extensão Nextflow, você pode obter um diagrama útil de como os processos estão conectados clicando no pequeno link `DAG preview` exibido logo acima do bloco workflow em qualquer script Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Isso dá a você uma boa visão geral de como os processos estão conectados e o que eles produzem.

Você vê que além do processo `sayHello` original, agora também temos `convertToUpper` e `collectGreetings`, que correspondem aos nomes dos processos que vimos na saída do console.
As duas novas definições de processo são estruturadas da mesma forma que o processo `sayHello`, exceto que `collectGreetings` recebe um parâmetro de entrada adicional chamado `batch` e produz duas saídas.

Não entraremos no código de cada um em detalhes, mas se você está curioso, pode consultar os detalhes em [Parte 2 do Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Por agora, vamos investigar como os processos estão conectados uns aos outros.

#### 2.3.2. Como os processos estão conectados

A coisa realmente interessante a observar aqui é como as chamadas de processo estão encadeadas no bloco `main:` do fluxo de trabalho.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // cria um canal para entradas de um arquivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emite uma saudação
    sayHello(greeting_ch)
    // converte a saudação para maiúsculas
    convertToUpper(sayHello.out)
    // coleta todas as saudações em um arquivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Você pode ver que a primeira chamada de processo, `sayHello(greeting_ch)`, está inalterada.
Então a próxima chamada de processo, para `convertToUpper`, se refere à saída de `sayHello` como `sayHello.out`.

O padrão é simples: `processName.out` se refere ao canal de saída de um processo, que pode ser passado diretamente para o próximo processo.
É assim que transportamos dados de uma etapa para a próxima no Nextflow.

#### 2.3.3. Um processo pode receber múltiplas entradas

A terceira chamada de processo, para `collectGreetings`, é um pouco diferente.

```groovy title="2b-multistep.nf" linenums="77"
    // coleta todas as saudações em um arquivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Você vê que esta chamada recebe duas entradas, `convertToUpper.out.collect()` e `params.batch`.
Ignorando o bit `.collect()` por enquanto, podemos generalizar isso como `collectGreetings(input1, input2)`.

Isso corresponde às duas declarações de entrada no módulo do processo:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Quando o Nextflow analisa isso, ele atribuirá a primeira entrada na chamada a `path input_files`, e a segunda a `val batch_name`.

Então agora você sabe que um processo pode receber múltiplas entradas, e como a chamada se parece no bloco workflow.

Agora vamos dar uma olhada mais de perto naquela primeira entrada, `convertToUpper.out.collect()`.

#### 2.3.4. O que `collect()` faz na chamada `collectGreetings`

Para passar a saída de `sayHello` para `convertToUpper`, simplesmente nos referimos ao canal de saída de `sayHello` como `sayHello.out`. Mas para a próxima etapa, estamos vendo uma referência a `convertToUpper.out.collect()`.

O que é esse bit `collect()` e o que ele faz?

É um operador, é claro. Assim como os operadores `splitCsv` e `map` que encontramos anteriormente.
Desta vez o operador é chamado `collect`, e é aplicado ao canal de saída produzido por `convertToUpper`.

O operador `collect` é usado para coletar as saídas de múltiplas chamadas ao mesmo processo e empacotá-las em um único elemento de canal.

No contexto deste fluxo de trabalho, ele está pegando as três saudações em maiúsculas no canal `convertToUpper.out` --que são três itens de canal separados, e normalmente seriam tratados em chamadas separadas pelo próximo processo-- e empacotando-os em um único item.

Em termos mais práticos: se não aplicássemos `collect()` à saída de `convertToUpper()` antes de alimentá-la para `collectGreetings()`, o Nextflow simplesmente executaria `collectGreetings()` independentemente em cada saudação, o que não alcançaria nosso objetivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Em contraste, usar `collect()` nos permite pegar todas as saudações em maiúsculas separadas produzidas pela segunda etapa do fluxo de trabalho e alimentá-las todas juntas para uma única chamada na terceira etapa do pipeline.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

É assim que colocamos todas as saudações de volta no mesmo arquivo.

Há muitos outros [operadores](https://www.nextflow.io/docs/latest/reference/operator.html#operator-page) disponíveis para aplicar transformações ao conteúdo de canais entre chamadas de processo.

Isso dá aos desenvolvedores de pipeline muita flexibilidade para personalizar a lógica de fluxo de seu pipeline.
A desvantagem é que às vezes pode tornar mais difícil decifrar o que o pipeline está fazendo.

#### 2.3.5. Um parâmetro de entrada pode ter um valor padrão

Você pode ter notado que `collectGreetings` recebe uma segunda entrada, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // coleta todas as saudações em um arquivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Isso passa um parâmetro CLI chamado `--batch` para o fluxo de trabalho.
No entanto, quando executamos o fluxo de trabalho anteriormente, não especificamos um parâmetro `--batch`.

O que está acontecendo aí?
Dê uma olhada no bloco `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Há um valor padrão configurado no fluxo de trabalho, então não precisamos fornecê-lo.
Mas se fornecermos um na linha de comando, o valor que especificarmos será usado em vez do padrão.

Tente:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Saída do comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Você deve ver novas saídas finais nomeadas com seu nome de lote personalizado.

??? abstract "Conteúdo do diretório"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Este é um aspecto da configuração de entrada, que cobriremos com mais detalhes na Parte 3, mas por enquanto o importante é saber que parâmetros de entrada podem receber valores padrão.

#### 2.3.6. Um processo pode produzir múltiplas saídas

Na definição do processo `collectGreetings`, vemos as seguintes declarações de saída:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Que são então referenciadas pelo nome dado com `emit:` no bloco `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Isso torna fácil então passar saídas específicas individualmente para outros processos no fluxo de trabalho, em combinação com vários operadores.

#### 2.3.7. Saídas publicadas podem ser organizadas

No bloco `output`, usamos caminhos personalizados para agrupar resultados intermediários para tornar mais fácil destacar apenas as saídas finais do fluxo de trabalho.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Há formas mais sofisticadas de organizar saídas publicadas; abordaremos algumas na parte sobre configuração.

!!! tip "Quer aprender mais sobre construir fluxos de trabalho?"

    Para cobertura detalhada de construir fluxos de trabalho de múltiplas etapas, veja [Hello Nextflow Parte 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Conclusão

Você entende em um nível básico como fluxos de trabalho de múltiplas etapas são construídos usando canais e operadores e como eles operam.
Você também viu que processos podem receber múltiplas entradas e produzir múltiplas saídas, e que estas podem ser publicadas de forma estruturada.

### O que vem a seguir?

Aprenda como pipelines Nextflow podem ser modularizados para promover reutilização de código e manutenibilidade.

---

## 3. Executando pipelines modularizados

Até agora, todos os fluxos de trabalho que olhamos consistiram em um único arquivo de fluxo de trabalho contendo todo o código relevante.

No entanto, pipelines do mundo real tipicamente se beneficiam de serem _modularizados_, significando que o código é dividido em diferentes arquivos.
Isso pode tornar seu desenvolvimento e manutenção mais eficientes e sustentáveis.

Aqui vamos demonstrar a forma mais comum de modularidade de código no Nextflow, que é o uso de **módulos**.

No Nextflow, um **módulo** é uma única definição de processo que é encapsulada sozinha em um arquivo de código autônomo.
Para usar um módulo em um fluxo de trabalho, você apenas adiciona uma declaração de importação de uma linha ao seu arquivo de código de fluxo de trabalho; então você pode integrar o processo no fluxo de trabalho da mesma forma que normalmente faria.
Isso torna possível reutilizar definições de processo em múltiplos fluxos de trabalho sem produzir múltiplas cópias do código.

Até agora estávamos executando fluxos de trabalho que tinham todos os seus processos incluídos em um arquivo de código monolítico.
Agora vamos ver como fica quando os processos são armazenados em módulos individuais.

Preparamos novamente um fluxo de trabalho adequado para fins de demonstração, chamado `2c-modules.nf`, junto com um conjunto de módulos localizados no diretório `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Conteúdo do diretório"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Você vê que há quatro arquivos Nextflow, cada um nomeado após um dos processos.
Você pode ignorar o arquivo `cowpy.nf` por enquanto; chegaremos a ele mais tarde.

### 3.1. Examine o código

Desta vez vamos olhar o código primeiro.
Comece abrindo o arquivo de fluxo de trabalho `2c-modules.nf`.

??? full-code "Arquivo de código completo"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Inclui módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Você vê que a lógica do fluxo de trabalho é exatamente a mesma da versão anterior do fluxo de trabalho.
No entanto, o código do processo foi removido do arquivo de fluxo de trabalho, e em vez disso há declarações `include` apontando para arquivos separados em `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Inclui módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Abra um desses arquivos e você encontrará o código para o processo correspondente.

??? full-code "Arquivo de código completo"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usa echo para imprimir 'Hello World!' em um arquivo
    */
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

Como você pode ver, o código do processo não mudou; ele foi apenas copiado para um arquivo de módulo individual em vez de estar no arquivo de fluxo de trabalho principal.
O mesmo se aplica aos outros dois processos.

Então vamos ver como fica executar esta nova versão.

### 3.2. Execute o fluxo de trabalho

Execute este comando no seu terminal, com a flag `-resume`:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
    [95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
    [5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
    ```

Você notará que as execuções de processo foram todas cacheadas com sucesso, significando que o Nextflow reconheceu que já fez o trabalho solicitado, mesmo que o código tenha sido dividido e o arquivo de fluxo de trabalho principal tenha sido renomeado.

Nada disso importa para o Nextflow; o que importa é o script de job que é gerado uma vez que todo o código foi reunido e avaliado.

!!! tip "Dica"

    Também é possível encapsular uma seção de um fluxo de trabalho como um 'subworkflow' que pode ser importado em um pipeline maior, mas isso está fora do escopo deste curso.

    Você pode aprender mais sobre desenvolver fluxos de trabalho composíveis na Side Quest sobre [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Conclusão

Você sabe como processos podem ser armazenados em módulos autônomos para promover reutilização de código e melhorar a manutenibilidade.

### O que vem a seguir?

Aprenda a usar contêineres para gerenciar dependências de software.

---

## 4. Usando software em contêiner

Até agora os fluxos de trabalho que estamos usando como exemplos só precisavam executar operações de processamento de texto muito básicas usando ferramentas UNIX disponíveis em nosso ambiente.

No entanto, pipelines do mundo real tipicamente requerem ferramentas e pacotes especializados que não estão incluídos por padrão na maioria dos ambientes.
Normalmente, você precisaria instalar essas ferramentas, gerenciar suas dependências e resolver quaisquer conflitos.

Tudo isso é muito tedioso e irritante.
Uma forma muito melhor de abordar este problema é usar **contêineres**.

Um **contêiner** é uma unidade leve, autônoma e executável de software criada a partir de uma **imagem** de contêiner que inclui tudo necessário para executar uma aplicação incluindo código, bibliotecas de sistema e configurações.

!!! Tip "Dica"

    Ensinamos isso usando a tecnologia [Docker](https://www.docker.com/get-started/), mas o Nextflow suporta [várias outras tecnologias de contêiner](https://www.nextflow.io/docs/latest/container.html#) também.

### 4.1. Use um contêiner diretamente

Primeiro, vamos tentar interagir com um contêiner diretamente.
Isso ajudará a solidificar seu entendimento do que são contêineres antes de começarmos a usá-los no Nextflow.

#### 4.1.1. Baixe a imagem do contêiner

Para usar um contêiner, você normalmente baixa ou "puxa" uma imagem de contêiner de um registro de contêiner, e então executa a imagem de contêiner para criar uma instância de contêiner.

A sintaxe geral é a seguinte:

```bash title="Sintaxe"
docker pull '<container>'
```

- `docker pull` é a instrução para o sistema de contêiner para puxar uma imagem de contêiner de um repositório.
- `'<container>'` é o endereço URI da imagem de contêiner.

Como exemplo, vamos puxar uma imagem de contêiner que contém [cowpy](https://github.com/jeffbuttars/cowpy), uma implementação em python de uma ferramenta chamada `cowsay` que gera arte ASCII para exibir entradas de texto arbitrárias de forma divertida.

Há vários repositórios onde você pode encontrar contêineres publicados.
Usamos o serviço [Seqera Containers](https://seqera.io/containers/) para gerar esta imagem de contêiner Docker a partir do pacote Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Execute o comando de pull completo:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Saída do comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Isso diz ao sistema para baixar a imagem especificada.
Uma vez que o download esteja completo, você tem uma cópia local da imagem de contêiner.

#### 4.1.2. Inicie o contêiner

Contêineres podem ser executados como um comando único, mas você também pode usá-los interativamente, o que dá a você um prompt de shell dentro do contêiner e permite que você brinque com o comando.

A sintaxe geral é a seguinte:

```bash title="Sintaxe"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` é a instrução para o sistema de contêiner para iniciar uma instância de contêiner a partir de uma imagem de contêiner e executar um comando nela.
- `--rm` diz ao sistema para desligar a instância de contêiner após o comando ter sido concluído.

Completamente montado, o comando de execução do contêiner se parece com isto:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Execute esse comando, e você deve ver seu prompt mudar para algo como `(base) root@b645838b3314:/tmp#`, que indica que você agora está dentro do contêiner.

Você pode verificar isso executando `ls` para listar o conteúdo do diretório:

```bash
ls /
```

??? success "Saída do comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Você vê que o sistema de arquivos dentro do contêiner é diferente do sistema de arquivos no seu sistema host.

!!! Tip "Dica"

    Quando você executa um contêiner, ele é isolado do sistema host por padrão.
    Isso significa que o contêiner não pode acessar nenhum arquivo no sistema host a menos que você explicitamente permita fazê-lo especificando que quer montar um volume como parte do comando `docker run` usando a seguinte sintaxe:

    ```bash title="Sintaxe"
    -v <outside_path>:<inside_path>
    ```

    Isso efetivamente estabelece um túnel através da parede do contêiner que você pode usar para acessar essa parte do seu sistema de arquivos.

    Isso é coberto em mais detalhes na [Parte 5 do Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Execute a ferramenta `cowpy`

De dentro do contêiner, você pode executar o comando `cowpy` diretamente.

```bash
cowpy "Hello Containers"
```

??? success "Saída do comando"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Isso produz arte ASCII do personagem vaca padrão (ou 'cowacter') com um balão de fala contendo o texto que especificamos.

Agora que você testou o uso básico, pode tentar dar alguns parâmetros.
Por exemplo, a documentação da ferramenta diz que podemos definir o personagem com `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Saída do comando"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Desta vez a saída de arte ASCII mostra o pinguim do Linux, Tux, porque especificamos o parâmetro `-c tux`.

Como você está dentro do contêiner, pode executar o comando cowpy quantas vezes quiser, variando os parâmetros de entrada, sem se preocupar em instalar nenhuma biblioteca no seu próprio sistema.

??? tip "Outros personagens disponíveis"

    Use a flag '-c' para escolher um personagem diferente, incluindo:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Sinta-se livre para brincar com isso.
Quando terminar, saia do contêiner usando o comando `exit`:

```bash
exit
```

Você se encontrará de volta no seu shell normal.

### 4.2. Use um contêiner em um fluxo de trabalho

Quando executamos um pipeline, queremos ser capazes de dizer ao Nextflow qual contêiner usar em cada etapa, e importante, queremos que ele lide com todo aquele trabalho que acabamos de fazer: puxar o contêiner, iniciá-lo, executar o comando e derrubar o contêiner quando terminar.

Boa notícia: é exatamente isso que o Nextflow vai fazer por nós.
Nós só precisamos especificar um contêiner para cada processo.

Para demonstrar como isso funciona, fizemos outra versão do nosso fluxo de trabalho que executa `cowpy` no arquivo de saudações coletadas produzido na terceira etapa.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-cowpy.svg"
</figure>

Isso deve produzir um arquivo contendo a arte ASCII com as três saudações no balão de fala.

#### 4.2.1. Examine o código

O fluxo de trabalho é muito similar ao anterior, mais a etapa extra para executar `cowpy`.

??? full-code "Arquivo de código completo"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Inclui módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // gera arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Você vê que este fluxo de trabalho importa um processo `cowpy` de um arquivo de módulo, e o chama na saída da chamada `collectGreetings()`, mais um parâmetro de entrada chamado `params.character`.

```groovy title="2d-container.nf" linenums="25"
// gera arte ASCII com cowpy
cowpy(collectGreetings.out, params.character)
```

O processo `cowpy`, que envolve o comando cowpy para gerar arte ASCII, é definido no módulo `cowpy.nf`.

??? full-code "Arquivo de código completo"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
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

O processo `cowpy` requer duas entradas: o caminho para um arquivo de entrada contendo o texto para colocar no balão de fala (`input_file`), e um valor para a variável de personagem.

Importante, ele também inclui a linha `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, que aponta para o URI do contêiner que usamos anteriormente.

#### 4.2.2. Verifique se o Docker está habilitado na configuração

Vamos antecipar ligeiramente a Parte 3 deste curso de treinamento introduzindo o arquivo de configuração `nextflow.config`, que é uma das principais formas que o Nextflow oferece para configurar a execução de fluxo de trabalho.
Quando um arquivo chamado `nextflow.config` está presente no diretório atual, o Nextflow o carregará automaticamente e aplicará qualquer configuração que ele contenha.

Para isso, incluímos um arquivo `nextflow.config` com uma única linha de código que habilita o Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Esta configuração diz ao Nextflow para usar Docker para qualquer processo que especifique um contêiner compatível.

!!! tip "Dica"

    É tecnicamente possível habilitar a execução Docker da linha de comando, por execução, usando o parâmetro `-with-docker <container>` no seu comando.
    No entanto, isso só nos permite especificar um contêiner para todo o fluxo de trabalho, enquanto a abordagem que acabamos de mostrar permite especificar um contêiner diferente por processo.
    Este último é muito melhor para modularidade, manutenção de código e reprodutibilidade.

#### 4.2.3. Execute o fluxo de trabalho

Só para recapitular, isto é o que estamos prestes a executar:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Você acha que vai funcionar?

Vamos executar o fluxo de trabalho com a flag `-resume`, e especificar que queremos que o personagem seja o peru.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

As três primeiras etapas foram cacheadas já que as executamos antes, mas o processo `cowpy` é novo então esse realmente é executado.

Você pode encontrar a saída da etapa `cowpy` no diretório `results`.

??? abstract "Conteúdo do arquivo"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
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

Você vê que o personagem está dizendo todas as saudações, já que ele executou no arquivo de saudações em maiúsculas coletadas.

Mais ao ponto, conseguimos executar isso como parte do nosso pipeline sem ter que fazer uma instalação apropriada do cowpy e todas as suas dependências.
E agora podemos compartilhar o pipeline com colaboradores e fazer com que eles o executem em sua infraestrutura sem que precisem instalar nada também, além do Docker ou uma de suas alternativas (como Singularity/Apptainer) conforme mencionado acima.

#### 4.2.4. Inspecione como o Nextflow lançou a tarefa em contêiner

Como uma coda final para esta seção, vamos dar uma olhada no subdiretório de trabalho para uma das chamadas de processo `cowpy` para ter um pouco mais de visão sobre como o Nextflow funciona com contêineres por baixo dos panos.

Verifique a saída do seu comando `nextflow run` para encontrar o caminho para o subdiretório de trabalho do processo `cowpy`.
Olhando o que obtivemos para a execução mostrada acima, a linha de log do console para o processo `cowpy` começa com `[7f/caf718]`.
Isso corresponde ao seguinte caminho de diretório truncado: `work/7f/caf718`.

Naquele diretório, você encontrará o arquivo `.command.run` que contém todos os comandos que o Nextflow executou em seu nome durante a execução do pipeline.

??? abstract "Conteúdo do arquivo"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}

    ...

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    ...
    ```

Se você procurar por `nxf_launch` neste arquivo, você deve ver algo assim:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Este comando de lançamento mostra que o Nextflow está usando um comando `docker run` muito similar para lançar a chamada de processo como fizemos quando o executamos manualmente.
Ele também monta o subdiretório de trabalho correspondente no contêiner, define o diretório de trabalho dentro do contêiner de acordo, e executa nosso script bash modelado no arquivo `.command.sh`.

Isso confirma que todo o trabalho duro que tivemos que fazer manualmente na seção anterior agora é feito para nós pelo Nextflow!

### Conclusão

Você entende qual papel os contêineres desempenham no gerenciamento de versões de ferramentas de software e garantindo reprodutibilidade.

Mais geralmente, você tem uma compreensão básica de quais são os componentes principais de pipelines Nextflow do mundo real e como eles estão organizados.
Você conhece os fundamentos de como o Nextflow pode processar múltiplas entradas eficientemente, executar fluxos de trabalho compostos de múltiplas etapas conectadas, aproveitar componentes de código modulares e utilizar contêineres para maior reprodutibilidade e portabilidade.

### O que vem a seguir?

Faça outra pausa! Essa foi uma grande pilha de informações sobre como pipelines Nextflow funcionam.

Na última seção deste treinamento, vamos mergulhar mais profundamente no tópico de configuração.
Você aprenderá como configurar a execução do seu pipeline para se adequar à sua infraestrutura, bem como gerenciar a configuração de entradas e parâmetros.

---

## Quiz

<quiz>
Por que o Nextflow cria um diretório de tarefa separado para cada chamada de processo?
- [ ] Para melhorar a velocidade de execução
- [ ] Para reduzir o uso de memória
- [x] Para isolar execuções e evitar colisões entre saídas
- [ ] Para habilitar compressão paralela de arquivos

Saiba mais: [1.3. Encontre as saídas originais e os logs](#13-encontre-as-saidas-originais-e-os-logs)
</quiz>

<quiz>
O que a opção `-ansi-log false` faz ao executar um fluxo de trabalho?
- [ ] Desabilita toda saída do console
- [x] Remove cores da saída
- [x] Mostra todos os caminhos de diretório de tarefa em vez de condensá-los em uma linha
- [ ] Habilita modo de depuração verboso

Saiba mais: [1.3.2. Faça o terminal mostrar mais detalhes](#132-faca-o-terminal-mostrar-mais-detalhes)

Você também pode usar qualquer uma das seguintes variáveis de ambiente se preferir este estilo:

```bash
export NXF_ANSI_LOG=0
# ou
export NO_COLOR=1
```

</quiz>

<quiz>
No código `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, o que `#!groovy .map { line -> line[0] }` faz?
- [ ] Filtra linhas vazias
- [ ] Ordena as linhas alfabeticamente
- [x] Extrai a primeira coluna de cada linha CSV
- [ ] Conta o número de linhas

Saiba mais: [1.4.1. Carregando os dados de entrada do CSV](#141-carregando-os-dados-de-entrada-do-csv)
</quiz>

<quiz>
Por que é importante incluir o valor de entrada nos nomes de arquivos de saída (ex., `#!groovy "${greeting}-output.txt"`)?
- [ ] Para melhorar a velocidade de processamento
- [ ] Para habilitar funcionalidade de resume
- [x] Para evitar que arquivos de saída sobrescrevam uns aos outros ao processar múltiplas entradas
- [ ] Para tornar arquivos mais fáceis de comprimir

Saiba mais: [1.4.3. Como as saídas são nomeadas](#143-como-as-saidas-sao-nomeadas)
</quiz>

<quiz>
Qual é o propósito da declaração `include` em um fluxo de trabalho modularizado?
- [ ] Copiar código de processo para o arquivo de fluxo de trabalho
- [x] Importar uma definição de processo de um arquivo de módulo externo
- [ ] Incluir configurações de configuração
- [ ] Adicionar comentários de documentação

Saiba mais: [3. Executando pipelines modularizados](#3-executando-pipelines-modularizados)
</quiz>

<quiz>
Quando você modulariza um fluxo de trabalho e o executa com `-resume`, o que acontece?
- [ ] Caching é desabilitado para processos modulares
- [ ] Todas as tarefas devem ser re-executadas
- [x] Caching funciona normalmente baseado nos scripts de job gerados
- [ ] Apenas o arquivo de fluxo de trabalho principal é cacheado

Saiba mais: [3.2. Execute o fluxo de trabalho](#32-execute-o-fluxo-de-trabalho)
</quiz>

<quiz>
O que a diretiva `container` em uma definição de processo especifica?
- [ ] O diretório de trabalho para o processo
- [ ] A alocação máxima de memória
- [x] O URI da imagem de contêiner para usar na execução do processo
- [ ] O formato do arquivo de saída

Saiba mais: [4.2. Use um contêiner em um fluxo de trabalho](#42-use-um-conteiner-em-um-fluxo-de-trabalho)
</quiz>

<quiz>
No arquivo `.command.run`, o que a função `nxf_launch` contém?
- [ ] A informação de versão do Nextflow
- [ ] Os parâmetros do fluxo de trabalho
- [x] O comando `docker run` com montagens de volume e configurações de contêiner
- [ ] As declarações de entrada do processo

Saiba mais: [4.2.4. Inspecione como o Nextflow lançou a tarefa em contêiner](#424-inspecione-como-o-nextflow-lancou-a-tarefa-em-conteiner)
</quiz>

<quiz>
O que o Nextflow automaticamente lida ao executar um processo em contêiner? (Selecione todos que se aplicam)
- [x] Puxar a imagem do contêiner se necessário
- [x] Montar o diretório de trabalho no contêiner
- [x] Executar o script do processo dentro do contêiner
- [x] Limpar a instância do contêiner após a execução

Saiba mais: [4. Usando software em contêiner](#4-usando-software-em-conteiner)
</quiz>
