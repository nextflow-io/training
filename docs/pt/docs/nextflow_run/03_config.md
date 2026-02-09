# Parte 3: Configuração de execução

Esta seção explorará como gerenciar a configuração de um pipeline Nextflow para personalizar seu comportamento, adaptá-lo a diferentes ambientes e otimizar o uso de recursos _sem alterar uma única linha do código do fluxo de trabalho em si_.

Existem várias maneiras de fazer isso, que podem ser usadas em combinação e são interpretadas de acordo com a ordem de precedência descrita na documentação de [Configuração](https://nextflow.io/docs/latest/config.html).

Nesta parte do curso, vamos mostrar o mecanismo de arquivo de configuração mais simples e comum, o arquivo `nextflow.config`, que você já encontrou na seção sobre contêineres na Parte 2.

Vamos abordar componentes essenciais da configuração do Nextflow, como diretivas de processo, executores, perfis e arquivos de parâmetros.
Ao aprender a utilizar essas opções de configuração de forma eficaz, você pode aproveitar ao máximo a flexibilidade, escalabilidade e desempenho dos pipelines Nextflow.

Para exercitar esses elementos de configuração, vamos executar uma cópia nova do fluxo de trabalho que executamos pela última vez no final da Parte 2 deste curso de treinamento, renomeado como `3-main.nf`.

Se você não está familiarizado com o pipeline Hello ou precisa de uma revisão, consulte [esta página de informações](../info/hello_pipeline.md).

---

## 1. Gerenciar parâmetros de entrada do fluxo de trabalho

??? example "Cenário"

    Você baixou um pipeline e quer executá-lo repetidamente com os mesmos arquivos de entrada e configurações, mas não quer digitar todos os parâmetros toda vez.
    Ou talvez você esteja configurando o pipeline para um colega que não se sente confortável com argumentos de linha de comando.

Vamos começar com um aspecto da configuração que é simplesmente uma extensão do que temos trabalhado até agora: o gerenciamento de parâmetros de entrada.

Atualmente, nosso fluxo de trabalho está configurado para aceitar vários valores de parâmetros via linha de comando, declarados em um bloco `params` no próprio script do fluxo de trabalho.
Um deles tem um valor padrão definido como parte de sua declaração.

No entanto, você pode querer definir padrões para todos eles, ou substituir o padrão existente sem ter que especificar parâmetros na linha de comando ou modificar o arquivo de script original.

Existem várias maneiras de fazer isso; vamos mostrar três maneiras básicas que são muito comumente usadas.

### 1.1. Configurar valores em `nextflow.config`

Esta é a abordagem mais simples, embora seja possivelmente a menos flexível, já que o arquivo principal `nextflow.config` não é algo que você queira editar para cada execução.
Mas tem a vantagem de separar as preocupações de _declarar_ os parâmetros no fluxo de trabalho (que definitivamente pertence lá) versus fornecer _valores padrão_, que estão mais em casa em um arquivo de configuração.

Vamos fazer isso em duas etapas.

#### 1.1.1. Criar um bloco `params` no arquivo de configuração

Faça as seguintes alterações de código no arquivo `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Parâmetros do pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Note que não simplesmente copiamos o bloco `params` do fluxo de trabalho para o arquivo de configuração.
Para o parâmetro `batch` que já tinha um valor padrão declarado, a sintaxe é um pouco diferente.
No arquivo do fluxo de trabalho, essa é uma declaração tipada.
Na configuração, essas são atribuições de valor.

Tecnicamente, isso é suficiente para substituir os valores padrão ainda especificados no arquivo do fluxo de trabalho.
Você poderia modificar o valor padrão para `batch` e executar o fluxo de trabalho para se convencer de que o valor definido no arquivo de configuração substitui o definido no arquivo do fluxo de trabalho.

Mas no espírito de mover a configuração completamente para o arquivo de configuração, vamos remover esse valor padrão do arquivo do fluxo de trabalho inteiramente.

#### 1.1.2. Remover o valor padrão para `batch` no arquivo do fluxo de trabalho

Faça a seguinte alteração de código no arquivo de fluxo de trabalho `3-main.nf`:

=== "Depois"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Parâmetros do pipeline
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Parâmetros do pipeline
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Agora o arquivo do fluxo de trabalho em si não define nenhum valor padrão para esses parâmetros.

#### 1.1.3. Executar o pipeline

Vamos testar se funciona corretamente sem especificar nenhum parâmetro na linha de comando.

```bash
nextflow run 3-main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Isso ainda produz a mesma saída de antes.

A saída final em arte ASCII está no diretório `results/3-main/`, com o nome `cowpy-COLLECTED-batch-output.txt`, igual a antes.

??? abstract "Conteúdo do arquivo"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
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

Funcionalmente, essa mudança não alterou nada, mas conceitualmente é um pouco mais limpo ter os valores padrão definidos no arquivo de configuração.

### 1.2. Usar um arquivo de configuração específico para a execução

??? example "Cenário"

    Você quer experimentar com diferentes configurações sem modificar seu arquivo de configuração principal.

Você pode fazer isso criando um novo arquivo `nextflow.config` em um subdiretório que você usará como diretório de trabalho para seus experimentos.

#### 1.2.1. Criar o diretório de trabalho com uma configuração em branco

Vamos começar criando um novo diretório e entrando nele:

```bash
mkdir -p tux-run
cd tux-run
```

Então, crie um arquivo de configuração em branco nesse diretório:

```bash
touch nextflow.config
```

Isso produz um arquivo vazio.

#### 1.2.2. Configurar a configuração experimental

Agora abra o novo arquivo e adicione os parâmetros que você quer personalizar:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Note que o caminho para o arquivo de entrada deve refletir a estrutura de diretórios.

#### 1.2.3. Executar o pipeline

Agora podemos executar nosso pipeline de dentro do nosso novo diretório de trabalho.
Certifique-se de adaptar o caminho adequadamente!

```bash
nextflow run ../3-main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Isso criará um novo conjunto de diretórios em `tux-run/` incluindo `tux-run/work/` e `tux-run/results/`.

Nesta execução, o Nextflow combina o `nextflow.config` em nosso diretório atual com o `nextflow.config` no diretório raiz do pipeline, e assim substitui o caractere padrão (turkey) pelo caractere tux.

O arquivo de saída final deve conter o caractere tux dizendo as saudações.

??? abstract "Conteúdo do arquivo"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

É isso; agora você tem um espaço para experimentar sem modificar sua configuração 'normal'.

!!! warning

    Certifique-se de voltar ao diretório anterior antes de passar para a próxima seção!

    ```bash
    cd ..
    ```

Agora vamos ver outra maneira útil de definir valores de parâmetros.

### 1.3. Usar um arquivo de parâmetros

??? example "Cenário"

    Você precisa compartilhar parâmetros de execução exatos com um colaborador, ou registrá-los para uma publicação.

A abordagem de subdiretório funciona muito bem para experimentar, mas envolve um pouco de configuração e requer que você adapte os caminhos adequadamente.
Há uma abordagem mais simples para quando você quer executar seu pipeline com um conjunto específico de valores, ou permitir que outra pessoa faça isso com o mínimo de esforço.

O Nextflow nos permite especificar parâmetros via um [arquivo de parâmetros](https://nextflow.io/docs/latest/config.html#parameter-file) em formato YAML ou JSON, o que torna muito conveniente gerenciar e distribuir conjuntos alternativos de valores padrão, por exemplo, bem como valores de parâmetros específicos da execução.

#### 1.3.1. Examinar o arquivo de parâmetros de exemplo

Para demonstrar isso, fornecemos um arquivo de parâmetros de exemplo no diretório atual, chamado `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Este arquivo de parâmetros contém um par chave-valor para cada uma das entradas que queremos especificar.
Note o uso de dois pontos (`:`) em vez de sinais de igual (`=`) se você comparar a sintaxe com o arquivo de configuração.
O arquivo de configuração é escrito em Groovy, enquanto o arquivo de parâmetros é escrito em YAML.

!!! info

    Também fornecemos uma versão JSON do arquivo de parâmetros como exemplo, mas não vamos executá-la aqui.
    Sinta-se à vontade para experimentar essa por conta própria.

#### 1.3.2. Executar o pipeline

Para executar o fluxo de trabalho com este arquivo de parâmetros, simplesmente adicione `-params-file <nome_do_arquivo>` ao comando base.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

O arquivo de saída final deve conter o caractere stegosaurus dizendo as saudações.

??? abstract "Conteúdo do arquivo"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Usar um arquivo de parâmetros pode parecer exagero quando você tem apenas alguns parâmetros para especificar, mas alguns pipelines esperam dezenas de parâmetros.
Nesses casos, usar um arquivo de parâmetros nos permitirá fornecer valores de parâmetros em tempo de execução sem ter que digitar linhas de comando massivas e sem modificar o script do fluxo de trabalho.

Também torna mais fácil distribuir conjuntos de parâmetros para colaboradores, ou como informação de suporte para uma publicação, por exemplo.
Isso torna seu trabalho mais reproduzível por outros.

### Conclusão

Você sabe como aproveitar as principais opções de configuração para gerenciar entradas do fluxo de trabalho.

### O que vem a seguir?

Aprenda como gerenciar onde e como as saídas do seu fluxo de trabalho são publicadas.

---

## 2. Gerenciar saídas do fluxo de trabalho

??? example "Cenário"

    Seu pipeline publica saídas em um diretório codificado, mas você quer organizar os resultados por projeto ou nome de experimento sem editar o código do fluxo de trabalho toda vez.

O fluxo de trabalho que herdamos usa caminhos para declarações de saída no nível do fluxo de trabalho, o que não é terrivelmente flexível e envolve muita repetição.

Vamos ver algumas maneiras comuns de configurar isso para ser mais flexível.

### 2.1. Personalizar o nome do diretório `outputDir`

Cada versão do fluxo de trabalho que executamos até agora publicou suas saídas em um subdiretório diferente codificado nas definições de saída.

Mudamos onde esse subdiretório estava na Parte 1 usando a flag CLI `-output-dir`, mas isso ainda é apenas uma string estática.
Vamos em vez disso configurar isso em um arquivo de configuração, onde podemos definir caminhos dinâmicos mais complexos.
Poderíamos criar um parâmetro totalmente novo para isso, mas vamos usar o parâmetro `batch` já que ele está bem ali.

#### 2.1.1. Definir um valor para `outputDir` no arquivo de configuração

O caminho que o Nextflow usa para publicar saídas é controlado pela opção `outputDir`.
Para alterar o caminho para todas as saídas, você pode definir um valor para esta opção no arquivo de configuração `nextflow.config`.

Adicione o seguinte código ao arquivo `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Parâmetros do pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Configurações de saída
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Parâmetros do pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Isso substituirá o caminho padrão integrado, `results/`, por `results_config/` mais o valor do parâmetro `batch` como subdiretório.

Lembre-se de que você também pode definir esta opção a partir da linha de comando usando o parâmetro `-output-dir` em seu comando (`-o` para abreviar), mas então você não poderia usar o valor do parâmetro `batch`.
Usar a flag CLI sobrescreverá `outputDir` na configuração se estiver definido.

#### 2.1.2. Remover a parte repetida do caminho codificado

Ainda temos um subdiretório codificado nas opções de saída, então vamos nos livrar disso agora.

Faça as seguintes alterações de código no arquivo do fluxo de trabalho:

=== "Depois"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

Também poderíamos ter apenas adicionado `${params.batch}` a cada caminho em vez de modificar o padrão `outputDir`, mas isso é mais conciso.

#### 2.1.3. Executar o pipeline

Vamos testar se funciona corretamente, definindo o nome do lote como `outdir` a partir da linha de comando.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Isso ainda produz a mesma saída de antes, exceto que desta vez encontramos nossas saídas em `results_config/outdir/`.

??? abstract "Conteúdo do diretório"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Você pode combinar esta abordagem com definições de caminho personalizadas para construir qualquer hierarquia de diretórios que desejar.

### 2.2. Organizar saídas por processo

Uma maneira popular de organizar ainda mais as saídas é fazê-lo por processo, _ou seja_, criar subdiretórios para cada processo executado no pipeline.

#### 2.2.1. Substituir os caminhos de saída por uma referência aos nomes dos processos

Tudo o que você precisa fazer é referenciar o nome do processo como `<processo>.name` na declaração do caminho de saída.

Faça as seguintes alterações no arquivo do fluxo de trabalho:

=== "Depois"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Isso remove os elementos codificados restantes da configuração do caminho de saída.

#### 2.2.2. Executar o pipeline

Vamos testar se funciona corretamente, definindo o nome do lote como `pnames` a partir da linha de comando.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Isso ainda produz a mesma saída de antes, exceto que desta vez encontramos nossas saídas em `results_config/pnames/`, e elas estão agrupadas por processo.

??? abstract "Conteúdo do diretório"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note

    Note que aqui apagamos a distinção entre `intermediates` versus saídas finais estando no nível superior.
    Você pode misturar e combinar essas abordagens e até incluir múltiplas variáveis, por exemplo, definindo o caminho da primeira saída como `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Definir o modo de publicação no nível do fluxo de trabalho

Finalmente, no espírito de reduzir a quantidade de código repetitivo, podemos substituir as declarações `mode` por saída com uma única linha na configuração.

#### 2.3.1. Adicionar `workflow.output.mode` ao arquivo de configuração

Adicione o seguinte código ao arquivo `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Configurações de saída
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Configurações de saída
    */
    outputDir = "results_config/${params.batch}"
    ```

Assim como a opção `outputDir`, dar a `workflow.output.mode` um valor no arquivo de configuração seria suficiente para substituir o que está definido no arquivo do fluxo de trabalho, mas vamos remover o código desnecessário de qualquer forma.

#### 2.3.2. Remover o modo de saída do arquivo do fluxo de trabalho

Faça as seguintes alterações no arquivo do fluxo de trabalho:

=== "Depois"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

Isso é mais conciso, não é?

#### 2.3.3. Executar o pipeline

Vamos testar se funciona corretamente, definindo o nome do lote como `outmode` a partir da linha de comando.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Isso ainda produz a mesma saída de antes, exceto que desta vez encontramos nossas saídas em `results_config/outmode/`.
Todas ainda são cópias adequadas, não links simbólicos.

??? abstract "Conteúdo do diretório"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

A principal razão pela qual você ainda pode querer usar a maneira por saída de definir o modo é se você quiser misturar e combinar dentro do mesmo fluxo de trabalho, _ou seja_, ter algumas saídas copiadas e algumas com links simbólicos.

Existem muitas outras opções que você pode personalizar dessa maneira, mas esperamos que isso lhe dê uma noção da gama de opções e como utilizá-las efetivamente para atender às suas preferências.

### Conclusão

Você sabe como controlar a nomenclatura e estrutura dos diretórios onde suas saídas são publicadas, bem como o modo de publicação de saída do fluxo de trabalho.

### O que vem a seguir?

Aprenda como adaptar a configuração do seu fluxo de trabalho ao seu ambiente de computação, começando com a tecnologia de empacotamento de software.

---

## 3. Selecionar uma tecnologia de empacotamento de software

Até agora estivemos olhando para elementos de configuração que controlam como as entradas entram e onde as saídas saem. Agora é hora de focar mais especificamente em adaptar a configuração do seu fluxo de trabalho ao seu ambiente de computação.

O primeiro passo nesse caminho é especificar de onde virão os pacotes de software que serão executados em cada etapa.
Eles já estão instalados no ambiente de computação local?
Precisamos recuperar imagens e executá-las via um sistema de contêineres?
Ou precisamos recuperar pacotes Conda e construir um ambiente Conda local?

Na primeira parte deste curso de treinamento (Partes 1-4) apenas usamos software instalado localmente em nosso fluxo de trabalho.
Então na Parte 5, introduzimos contêineres Docker e o arquivo `nextflow.config`, que usamos para habilitar o uso de contêineres Docker.

Agora vamos ver como podemos configurar uma opção alternativa de empacotamento de software via o arquivo `nextflow.config`.

### 3.1. Desabilitar Docker e habilitar Conda no arquivo de configuração

??? example "Cenário"

    Você está movendo seu pipeline para um cluster HPC onde o Docker não é permitido por razões de segurança.
    O cluster suporta Singularity e Conda, então você precisa mudar sua configuração adequadamente.

Como observado anteriormente, o Nextflow suporta múltiplas tecnologias de contêineres incluindo Singularity (que é mais amplamente usado em HPC), bem como gerenciadores de pacotes de software como Conda.

Podemos alterar nosso arquivo de configuração para usar Conda em vez de Docker.
Para fazer isso, vamos mudar o valor de `docker.enabled` para `false`, e adicionar uma diretiva habilitando o uso de Conda:

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Isso permitirá que o Nextflow crie e utilize ambientes Conda para processos que têm pacotes Conda especificados.
O que significa que agora precisamos adicionar um desses ao nosso processo `cowpy`!

### 3.2. Especificar um pacote Conda na definição do processo

Já recuperamos o URI para um pacote Conda contendo a ferramenta `cowpy`: `conda-forge::cowpy==1.1.5`

Agora adicionamos o URI à definição do processo `cowpy` usando a diretiva `conda`:

=== "Depois"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Para ser claro, não estamos _substituindo_ a diretiva `docker`, estamos _adicionando_ uma opção alternativa.

!!! tip

    Existem algumas maneiras diferentes de obter o URI para um determinado pacote conda.
    Recomendamos usar a consulta de pesquisa [Seqera Containers](https://seqera.io/containers/), que lhe dará um URI que você pode copiar e colar, mesmo que não esteja planejando criar um contêiner a partir dele.

### 3.3. Executar o fluxo de trabalho para verificar que pode usar Conda

Vamos experimentar.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Saída do comando"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Isso deve funcionar sem problemas e produzir as mesmas saídas de antes em `results_config/conda`.

Nos bastidores, o Nextflow recuperou os pacotes Conda e criou o ambiente, o que normalmente dá um pouco de trabalho; então é bom que não tenhamos que fazer nada disso nós mesmos!

!!! info

    Isso executa rapidamente porque o pacote `cowpy` é bastante pequeno, mas se você estiver trabalhando com pacotes grandes, pode demorar um pouco mais do que o normal na primeira vez, e você pode ver a saída do console ficar 'presa' por um minuto ou mais antes de completar.
    Isso é normal e se deve ao trabalho extra que o Nextflow faz na primeira vez que você usa um novo pacote.

Do nosso ponto de vista, parece que funciona exatamente da mesma forma que executar com Docker, embora no backend a mecânica seja um pouco diferente.

Isso significa que estamos prontos para executar com ambientes Conda se necessário.

??? info "Misturando e combinando Docker e Conda"

    Como essas diretivas são atribuídas por processo, é possível 'misturar e combinar', _ou seja_, configurar alguns dos processos em seu fluxo de trabalho para executar com Docker e outros com Conda, por exemplo, se a infraestrutura de computação que você está usando suportar ambos.
    Nesse caso, você habilitaria tanto Docker quanto Conda em seu arquivo de configuração.
    Se ambos estiverem disponíveis para um determinado processo, o Nextflow priorizará contêineres.

    E como observado anteriormente, o Nextflow suporta múltiplas outras tecnologias de empacotamento de software e contêineres, então você não está limitado a apenas essas duas.

### Conclusão

Você sabe como configurar qual pacote de software cada processo deve usar, e como alternar entre tecnologias.

### O que vem a seguir?

Aprenda como alterar a plataforma de execução usada pelo Nextflow para realmente fazer o trabalho.

---

## 4. Selecionar uma plataforma de execução

??? example "Cenário"

    Você tem desenvolvido e testado seu pipeline em seu laptop, mas agora precisa executá-lo em milhares de amostras.
    Sua instituição tem um cluster HPC com um agendador Slurm que você gostaria de usar em vez disso.

Até agora, temos executado nosso pipeline com o executor local.
Isso executa cada tarefa na máquina em que o Nextflow está sendo executado.
Quando o Nextflow começa, ele olha para as CPUs e memória disponíveis.
Se os recursos das tarefas prontas para executar excederem os recursos disponíveis, o Nextflow reterá as últimas tarefas da execução até que uma ou mais das tarefas anteriores tenham terminado, liberando os recursos necessários.

O executor local é conveniente e eficiente, mas é limitado a essa única máquina. Para cargas de trabalho muito grandes, você pode descobrir que sua máquina local é um gargalo, seja porque você tem uma única tarefa que requer mais recursos do que você tem disponível, ou porque você tem tantas tarefas que esperar que uma única máquina as execute levaria muito tempo.

O Nextflow suporta [muitos backends de execução diferentes](https://nextflow.io/docs/latest/executor.html), incluindo agendadores HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor e outros), bem como backends de execução em nuvem (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e mais).

### 4.1. Direcionando para um backend diferente

A escolha do executor é definida por uma diretiva de processo chamada `executor`.
Por padrão, está definido como `local`, então a seguinte configuração está implícita:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Para definir o executor para direcionar a um backend diferente, você simplesmente especificaria o executor que deseja usando sintaxe similar à descrita acima para alocações de recursos (veja [Executors](https://nextflow.io/docs/latest/executor.html) para todas as opções).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Não podemos realmente testar isso no ambiente de treinamento porque ele não está configurado para se conectar a um HPC.

### 4.2. Lidando com sintaxe específica do backend para parâmetros de execução

A maioria das plataformas de computação de alto desempenho permite (e às vezes requer) que você especifique certos parâmetros, como solicitações e limitações de alocação de recursos (por exemplo, número de CPUs e memória) e nome da fila de trabalhos a usar.

Infelizmente, cada um desses sistemas usa diferentes tecnologias, sintaxes e configurações para definir como um trabalho deve ser definido e submetido ao agendador relevante.

??? abstract "Exemplos"

    Por exemplo, o mesmo trabalho requerendo 8 CPUs e 4GB de RAM para ser executado na fila "my-science-work" precisa ser expresso das seguintes maneiras diferentes dependendo do backend.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Felizmente, o Nextflow simplifica tudo isso.
Ele fornece uma sintaxe padronizada para que você possa especificar as propriedades relevantes como `cpus`, `memory` e `queue` apenas uma vez (veja [Diretivas de processo](https://nextflow.io/docs/latest/reference/process.html#process-directives) para todas as opções disponíveis).
Então, em tempo de execução, o Nextflow usará essas configurações para gerar os scripts específicos do backend apropriados com base na configuração do executor.

Vamos cobrir essa sintaxe padronizada na próxima seção.

### Conclusão

Você agora sabe como alterar o executor para usar diferentes tipos de infraestrutura de computação.

### O que vem a seguir?

Aprenda como avaliar e expressar alocações e limitações de recursos no Nextflow.

---

## 5. Controlar alocações de recursos de computação

??? example "Cenário"

    Seu pipeline continua falhando no cluster porque as tarefas estão sendo eliminadas por exceder os limites de memória.
    Ou talvez você esteja sendo cobrado por recursos que não está usando e queira otimizar custos.

A maioria das plataformas de computação de alto desempenho permite (e às vezes requer) que você especifique certos parâmetros de alocação de recursos, como número de CPUs e memória.

Por padrão, o Nextflow usará uma única CPU e 2GB de memória para cada processo.
As diretivas de processo correspondentes são chamadas `cpus` e `memory`, então a seguinte configuração está implícita:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Você pode modificar esses valores, seja para todos os processos ou para processos nomeados específicos, usando diretivas de processo adicionais em seu arquivo de configuração.
O Nextflow os traduzirá nas instruções apropriadas para o executor escolhido.

Mas como você sabe quais valores usar?

### 5.1. Executar o fluxo de trabalho para gerar um relatório de utilização de recursos

??? example "Cenário"

    Você não sabe quanta memória ou CPU seus processos precisam e quer evitar desperdiçar recursos ou ter trabalhos eliminados.

Se você não sabe antecipadamente quanta CPU e memória seus processos provavelmente precisarão, você pode fazer algum perfil de recursos, o que significa que você executa o fluxo de trabalho com algumas alocações padrão, registra quanto cada processo usou e, a partir daí, estima como ajustar as alocações base.

Convenientemente, o Nextflow inclui ferramentas integradas para fazer isso, e gerará alegremente um relatório para você mediante solicitação.

Para fazer isso, adicione `-with-report <nome_do_arquivo>.html` à sua linha de comando.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

O relatório é um arquivo html, que você pode baixar e abrir em seu navegador. Você também pode clicar com o botão direito nele no explorador de arquivos à esquerda e clicar em `Show preview` para visualizá-lo no ambiente de treinamento.

Reserve alguns minutos para examinar o relatório e ver se você consegue identificar algumas oportunidades para ajustar recursos.
Certifique-se de clicar nas abas que mostram os resultados de utilização como uma porcentagem do que foi alocado.

Veja [Relatórios](https://nextflow.io/docs/latest/reports.html) para documentação sobre todos os recursos disponíveis.

### 5.2. Definir alocações de recursos para todos os processos

O perfil mostra que os processos em nosso fluxo de trabalho de treinamento são muito leves, então vamos reduzir a alocação de memória padrão para 1GB por processo.

Adicione o seguinte ao seu arquivo `nextflow.config`, antes da seção de parâmetros do pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Configurações de processo
*/
process {
    memory = 1.GB
}
```

Isso ajudará a reduzir a quantidade de computação que consumimos.

### 5.3. Definir alocações de recursos para um processo específico

Ao mesmo tempo, vamos fingir que o processo `cowpy` requer mais recursos do que os outros, apenas para que possamos demonstrar como ajustar alocações para um processo individual.

=== "Depois"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Configurações de processo
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Configurações de processo
    */
    process {
        memory = 1.GB
    }
    ```

Com esta configuração, todos os processos solicitarão 1GB de memória e uma única CPU (o padrão implícito), exceto o processo `cowpy`, que solicitará 2GB e 2 CPUs.

!!! info

    Se você tiver uma máquina com poucas CPUs e alocar um número alto por processo, você pode ver chamadas de processo sendo enfileiradas uma atrás da outra.
    Isso ocorre porque o Nextflow garante que não solicitemos mais CPUs do que estão disponíveis.

### 5.4. Executar o fluxo de trabalho com a configuração atualizada

Vamos experimentar isso, fornecendo um nome de arquivo diferente para o relatório de perfil para que possamos comparar o desempenho antes e depois das alterações de configuração.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Você provavelmente não notará nenhuma diferença real, já que esta é uma carga de trabalho tão pequena, mas esta é a abordagem que você usaria para analisar o desempenho e os requisitos de recursos de um fluxo de trabalho do mundo real.

É muito útil quando seus processos têm diferentes requisitos de recursos. Isso permite que você dimensione corretamente as alocações de recursos que você configura para cada processo com base em dados reais, não em suposições.

!!! tip

    Isso é apenas uma pequena amostra do que você pode fazer para otimizar seu uso de recursos.
    O próprio Nextflow tem uma [lógica de repetição dinâmica](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) muito interessante integrada para repetir trabalhos que falham devido a limitações de recursos.
    Além disso, a Seqera Platform oferece ferramentas orientadas por IA para otimizar suas alocações de recursos automaticamente também.

### 5.5. Adicionar limites de recursos

Dependendo de qual executor de computação e infraestrutura de computação você está usando, pode haver algumas restrições sobre o que você pode (ou deve) alocar.
Por exemplo, seu cluster pode exigir que você permaneça dentro de certos limites.

Você pode usar a diretiva `resourceLimits` para definir as limitações relevantes. A sintaxe se parece com isso quando está sozinha em um bloco de processo:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

O Nextflow traduzirá esses valores nas instruções apropriadas dependendo do executor que você especificou.

Não vamos executar isso, já que não temos acesso à infraestrutura relevante no ambiente de treinamento.
No entanto, se você tentasse executar o fluxo de trabalho com alocações de recursos que excedem esses limites, então procurar o comando `sbatch` no arquivo de script `.command.run`, você veria que as solicitações que realmente são enviadas ao executor são limitadas aos valores especificados por `resourceLimits`.

??? info "Configurações de referência institucionais"

    O projeto nf-core compilou uma [coleção de arquivos de configuração](https://nf-co.re/configs/) compartilhados por várias instituições ao redor do mundo, cobrindo uma ampla gama de executores HPC e em nuvem.

    Essas configurações compartilhadas são valiosas tanto para pessoas que trabalham lá e podem, portanto, apenas utilizar a configuração de sua instituição pronta para uso, quanto como modelo para pessoas que estão procurando desenvolver uma configuração para sua própria infraestrutura.

### Conclusão

Você sabe como gerar um relatório de perfil para avaliar a utilização de recursos e como modificar alocações de recursos para todos os processos e/ou para processos individuais, bem como definir limitações de recursos para execução em HPC.

### O que vem a seguir?

Aprenda como configurar perfis de configuração predefinidos e alternar entre eles em tempo de execução.

---

## 6. Usar perfis para alternar entre configurações predefinidas

??? example "Cenário"

    Você alterna regularmente entre executar pipelines em seu laptop para desenvolvimento e no HPC de sua instituição para execuções de produção.
    Você está cansado de alterar manualmente as configurações toda vez que muda de ambiente.

Mostramos várias maneiras de personalizar a configuração do seu pipeline dependendo do projeto em que você está trabalhando ou do ambiente de computação que está usando.

Você pode querer alternar entre configurações alternativas dependendo de qual infraestrutura de computação você está usando. Por exemplo, você pode querer desenvolver e executar testes em pequena escala localmente em seu laptop, depois executar cargas de trabalho em escala completa em HPC ou nuvem.

O Nextflow permite que você configure qualquer número de [**perfis**](https://nextflow.io/docs/latest/config.html#profiles) que descrevem diferentes configurações, que você pode então selecionar em tempo de execução usando um argumento de linha de comando, em vez de ter que modificar o próprio arquivo de configuração.

### 6.1. Criar perfis para alternar entre desenvolvimento local e execução em HPC

Vamos configurar dois perfis alternativos; um para executar cargas em pequena escala em um computador regular, onde usaremos contêineres Docker, e um para executar em um HPC universitário com um agendador Slurm, onde usaremos pacotes Conda.

#### 6.1.1. Configurar os perfis

Adicione o seguinte ao seu arquivo `nextflow.config`, após a seção de parâmetros do pipeline, mas antes das configurações de saída:

```groovy title="nextflow.config" linenums="24"
/*
* Perfis
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

Você vê que para o HPC universitário, também estamos especificando limitações de recursos.

#### 6.1.2. Executar o fluxo de trabalho com um perfil

Para especificar um perfil em nossa linha de comando Nextflow, usamos o argumento `-profile`.

Vamos tentar executar o fluxo de trabalho com a configuração `my_laptop`.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Como você pode ver, isso nos permite alternar entre configurações muito convenientemente em tempo de execução.

!!! warning

    O perfil `univ_hpc` não será executado corretamente no ambiente de treinamento, já que não temos acesso a um agendador Slurm.

Se no futuro encontrarmos outros elementos de configuração que estão sempre co-ocorrendo com estes, podemos simplesmente adicioná-los ao(s) perfil(is) correspondente(s).
Também podemos criar perfis adicionais se houver outros elementos de configuração que queremos agrupar.

### 6.2. Criar um perfil de parâmetros de teste

??? example "Cenário"

    Você quer que outros possam experimentar seu pipeline rapidamente sem reunir seus próprios dados de entrada.

Os perfis não são apenas para configuração de infraestrutura.
Também podemos usá-los para definir valores padrão para parâmetros do fluxo de trabalho, para facilitar que outros experimentem o fluxo de trabalho sem ter que reunir valores de entrada apropriados eles mesmos.
Você pode considerar isso uma alternativa ao uso de um arquivo de parâmetros.

#### 6.2.1. Configurar o perfil

A sintaxe para expressar valores padrão neste contexto se parece com isso, para um perfil que nomeamos `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Se adicionarmos um perfil de teste para nosso fluxo de trabalho, o bloco `profiles` se torna:

```groovy title="nextflow.config" linenums="24"
/*
* Perfis
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Assim como para perfis de configuração técnica, você pode configurar múltiplos perfis diferentes especificando parâmetros sob qualquer nome arbitrário que desejar.

#### 6.2.2. Executar o fluxo de trabalho localmente com o perfil de teste

Convenientemente, os perfis não são mutuamente exclusivos, então podemos especificar múltiplos perfis em nossa linha de comando usando a seguinte sintaxe `-profile <perfil1>,<perfil2>` (para qualquer número de perfis).

Se você combinar perfis que definem valores para os mesmos elementos de configuração e são descritos no mesmo arquivo de configuração, o Nextflow resolverá o conflito usando qualquer valor que leu por último (_ou seja_, o que vem depois no arquivo).
Se as configurações conflitantes forem definidas em diferentes fontes de configuração, a [ordem de precedência](https://www.nextflow.io/docs/latest/config.html) padrão se aplica.

Vamos tentar adicionar o perfil de teste ao nosso comando anterior:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Isso usará Docker onde possível e produzirá saídas em `results_config/test`, e desta vez o caractere é a dupla cômica `dragonandcow`.

??? abstract "Conteúdo do arquivo"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Isso significa que, desde que distribuamos quaisquer arquivos de dados de teste com o código do fluxo de trabalho, qualquer pessoa pode experimentar rapidamente o fluxo de trabalho sem ter que fornecer suas próprias entradas via linha de comando ou arquivo de parâmetros.

!!! tip

    Podemos apontar para URLs para arquivos maiores que são armazenados externamente.
    O Nextflow os baixará automaticamente, desde que haja uma conexão aberta.

    Para mais detalhes, veja a Missão Secundária [Trabalhando com Arquivos](../side_quests/working_with_files.md)

### 6.3. Usar `nextflow config` para ver a configuração resolvida

Como observado acima, às vezes o mesmo parâmetro pode ser definido com valores diferentes em perfis que você deseja combinar.
E mais geralmente, existem numerosos lugares onde elementos de configuração podem ser armazenados, e às vezes as mesmas propriedades podem ser definidas com valores diferentes em lugares diferentes.

O Nextflow aplica uma [ordem de precedência](https://nextflow.io/docs/latest/config.html#configuration-file) definida para resolver quaisquer conflitos, mas isso pode ser complicado de determinar você mesmo.
E mesmo que nada esteja em conflito, pode ser tedioso procurar todos os lugares possíveis onde as coisas poderiam ser configuradas.

Felizmente, o Nextflow inclui uma ferramenta utilitária conveniente chamada `config` que pode automatizar todo esse processo para você.

A ferramenta `config` explorará todo o conteúdo em seu diretório de trabalho atual, coletará quaisquer arquivos de configuração e produzirá a configuração totalmente resolvida que o Nextflow usaria para executar o fluxo de trabalho.
Isso permite que você descubra quais configurações serão usadas sem ter que iniciar nada.

#### 6.3.1. Resolver a configuração padrão

Execute este comando para resolver a configuração que seria aplicada por padrão.

```bash
nextflow config
```

??? success "Saída do comando"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Isso mostra a configuração base que você obtém se não especificar nada extra na linha de comando.

#### 6.3.2. Resolver a configuração com configurações específicas ativadas

Se você fornecer parâmetros de linha de comando, por exemplo, habilitando um ou mais perfis ou carregando um arquivo de parâmetros, o comando adicionalmente levará isso em conta.

```bash
nextflow config -profile my_laptop,test
```

??? success "Saída do comando"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Isso se torna especialmente útil para projetos complexos que envolvem múltiplas camadas de configuração.

### Conclusão

Você sabe como usar perfis para selecionar uma configuração predefinida em tempo de execução com o mínimo de complicação.
De forma mais geral, você sabe como configurar suas execuções de fluxo de trabalho para se adequar a diferentes plataformas de computação e melhorar a reprodutibilidade de suas análises.

### O que vem a seguir?

Aprenda como executar pipelines diretamente de repositórios remotos como GitHub.

---

## 7. Executar pipelines de repositórios remotos

??? example "Cenário"

    Você quer executar um pipeline bem estabelecido como aqueles do nf-core sem ter que baixar e gerenciar o código você mesmo.

Até agora temos executado scripts de fluxo de trabalho localizados no diretório atual.
Na prática, você frequentemente vai querer executar pipelines armazenados em repositórios remotos, como GitHub.

O Nextflow torna isso simples: você pode executar qualquer pipeline diretamente de uma URL de repositório Git sem baixá-lo manualmente primeiro.

### 7.1. Executar um pipeline do GitHub

A sintaxe básica para executar um pipeline remoto é `nextflow run <repositório>`, onde `<repositório>` pode ser um caminho de repositório GitHub como `nextflow-io/hello`, uma URL completa, ou um caminho para GitLab, Bitbucket ou outros serviços de hospedagem Git.

Tente executar o pipeline de demonstração oficial "hello" do Nextflow:

```bash
nextflow run nextflow-io/hello
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Na primeira vez que você executa um pipeline remoto, o Nextflow o baixa e o armazena em cache localmente.
Execuções subsequentes usam a versão em cache, a menos que você solicite explicitamente uma atualização.

### 7.2. Especificar uma versão para reprodutibilidade

Por padrão, o Nextflow executa a versão mais recente do branch padrão.
Você pode especificar uma versão particular (tag), branch ou commit usando a flag `-r`:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Especificar versões exatas é essencial para reprodutibilidade.

### Conclusão

Você sabe como executar pipelines diretamente do GitHub e outros repositórios remotos, e como especificar versões para reprodutibilidade.

### O que vem a seguir?

Dê um grande tapinha nas costas!
Você sabe tudo o que precisa saber para começar a executar e gerenciar pipelines Nextflow.

Isso conclui este curso, mas se você está ansioso para continuar aprendendo, temos duas recomendações principais:

- Se você quer se aprofundar no desenvolvimento de seus próprios pipelines, dê uma olhada em [Hello Nextflow](../hello_nextflow/index.md), um curso para iniciantes que cobre a mesma progressão geral que este, mas entra em muito mais detalhes sobre canais e operadores.
- Se você gostaria de continuar aprendendo como executar pipelines Nextflow sem se aprofundar no código, dê uma olhada na primeira parte de [Hello nf-core](../hello_nf-core/index.md), que introduz as ferramentas para encontrar e executar pipelines do projeto [nf-core](https://nf-co.re/), extremamente popular.

Divirta-se!

---

## Quiz

<quiz>
Quando valores de parâmetros são definidos tanto no arquivo do fluxo de trabalho quanto em `nextflow.config`, qual tem precedência?
- [ ] O valor do arquivo do fluxo de trabalho
- [x] O valor do arquivo de configuração
- [ ] O primeiro valor encontrado
- [ ] Causa um erro

Saiba mais: [1.1. Configurar valores em `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Qual é a diferença de sintaxe entre definir um padrão de parâmetro em um arquivo de fluxo de trabalho vs. um arquivo de configuração?
- [ ] Eles usam a mesma sintaxe
- [x] O fluxo de trabalho usa declaração tipada (`#!groovy param: Type = value`), a configuração usa atribuição (`#!groovy param = value`)
- [ ] A configuração usa declaração tipada, o fluxo de trabalho usa atribuição
- [ ] Apenas arquivos de configuração podem definir valores padrão

Saiba mais: [1.1. Configurar valores em `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Como você especifica um arquivo de parâmetros ao executar um fluxo de trabalho?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Saiba mais: [1.3. Usar um arquivo de parâmetros](#13-use-a-parameter-file)
</quiz>

<quiz>
O que a opção de configuração `outputDir` controla?
- [ ] A localização do diretório de trabalho
- [x] O caminho base onde as saídas do fluxo de trabalho são publicadas
- [ ] O diretório para arquivos de log
- [ ] A localização dos arquivos de módulo

Saiba mais: [2.1. Personalizar o nome do diretório outputDir](#21-customize-the-outputdir-directory-name)
</quiz>

<quiz>
Como você referencia um nome de processo dinamicamente na configuração do caminho de saída?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

Saiba mais: [2.2. Organizar saídas por processo](#22-organize-outputs-by-process)
</quiz>

<quiz>
Se tanto Docker quanto Conda estiverem habilitados e um processo tiver ambas as diretivas, qual é priorizado?
- [x] Docker (contêineres)
- [ ] Conda
- [ ] O primeiro definido no processo
- [ ] Causa um erro

Saiba mais: [3. Selecionar uma tecnologia de empacotamento de software](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Qual é o executor padrão no Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Saiba mais: [4. Selecionar uma plataforma de execução](#4-select-an-execution-platform)
</quiz>

<quiz>
Qual comando gera um relatório de utilização de recursos?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Saiba mais: [5.1. Executar o fluxo de trabalho para gerar um relatório de utilização de recursos](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
Como você define requisitos de recursos para um processo específico chamado `cowpy` no arquivo de configuração?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Saiba mais: [5.3. Definir alocações de recursos para um processo específico](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
O que a diretiva `resourceLimits` faz?
- [ ] Define requisitos mínimos de recursos
- [ ] Aloca recursos para processos
- [x] Limita os recursos máximos que podem ser solicitados
- [ ] Monitora o uso de recursos em tempo real

Saiba mais: [5.5. Adicionar limites de recursos](#55-add-resource-limits)
</quiz>

<quiz>
Como você especifica múltiplos perfis em um único comando?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Saiba mais: [6. Usar perfis para alternar entre configurações predefinidas](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
Qual comando mostra a configuração totalmente resolvida que o Nextflow usaria?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Saiba mais: [6.3. Usar `nextflow config` para ver a configuração resolvida](#63-use-nextflow-config-to-see-the-resolved-configuration)
</quiz>

<quiz>
Para que os perfis podem ser usados? (Selecione todas as opções aplicáveis)
- [x] Definir configurações específicas de infraestrutura (executores, contêineres)
- [x] Definir limites de recursos para diferentes ambientes
- [x] Fornecer parâmetros de teste para facilitar o teste do fluxo de trabalho
- [ ] Definir novos processos

Saiba mais: [6. Usar perfis para alternar entre configurações predefinidas](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
