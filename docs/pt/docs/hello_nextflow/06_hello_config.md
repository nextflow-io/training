# Parte 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=pt" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja [a playlist completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) no canal do YouTube do Nextflow.

:green_book: A transcrição do vídeo está disponível [aqui](./transcripts/06_hello_config.md).
///

Esta seção explorará como configurar e gerenciar a configuração do seu pipeline Nextflow para que você possa personalizar seu comportamento, adaptá-lo a diferentes ambientes e otimizar o uso de recursos _sem alterar uma única linha do código do fluxo de trabalho_.

Existem várias maneiras de fazer isso, que podem ser usadas em combinação e são interpretadas de acordo com a [ordem de precedência](https://nextflow.io/docs/latest/config.html) descrita na documentação de configuração.

Nesta parte do curso, vamos mostrar o mecanismo de arquivo de configuração mais simples e comum, o arquivo [`nextflow.config`](https://nextflow.io/docs/latest/config.html), que você já encontrou na Parte 5: Hello Containers.

Vamos abordar componentes essenciais da configuração do Nextflow, como diretivas de processos, executores, perfis e arquivos de parâmetros.
Ao aprender a utilizar essas opções de configuração de forma eficaz, você pode melhorar a flexibilidade, escalabilidade e desempenho dos seus pipelines.

??? info "Como começar a partir desta seção"

    Esta seção do curso assume que você completou as Partes 1-5 do curso [Hello Nextflow](./index.md) e tem um pipeline completo funcionando.

    Se você está começando o curso a partir deste ponto, precisará copiar o diretório `modules` e o arquivo `nextflow.config` das soluções:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    O arquivo `nextflow.config` contém a linha `docker.enabled = true` que habilita o uso de contêineres Docker.

    Se você não está familiarizado com o pipeline Hello ou precisa relembrar, veja [esta página de informações](../info/hello_pipeline.md).

---

## 0. Aquecimento: Execute `hello-config.nf`

Vamos usar o script do fluxo de trabalho `hello-config.nf` como ponto de partida.
Ele é equivalente ao script produzido ao trabalhar na Parte 5 deste curso de treinamento, exceto que mudamos os destinos de saída:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Apenas para garantir que tudo está funcionando, execute o script uma vez antes de fazer qualquer alteração:

```bash
nextflow run hello-config.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Como anteriormente, você encontrará os arquivos de saída no diretório especificado no bloco `output` (`results/hello_config/`).

??? abstract "Conteúdo do diretório"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

A saída final em arte ASCII está no diretório `results/hello_config/`, com o nome `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Se isso funcionou para você, está pronto para aprender como configurar seus pipelines.

---

## 1. Gerencie parâmetros de entrada do fluxo de trabalho

Vamos começar com um aspecto da configuração que é simplesmente uma extensão do que temos trabalhado até agora: o gerenciamento de parâmetros de entrada.

Atualmente, nosso fluxo de trabalho está configurado para aceitar vários valores de parâmetros via linha de comando, com valores padrão definidos em um bloco `params` no próprio script do fluxo de trabalho.
No entanto, você pode querer substituir esses padrões sem ter que especificar parâmetros na linha de comando ou modificar o arquivo de script original.

Existem várias maneiras de fazer isso; vamos mostrar três formas básicas que são muito comumente usadas.

### 1.1. Mova os valores padrão para o `nextflow.config`

Esta é a abordagem mais simples, embora seja possivelmente a menos flexível, já que o arquivo `nextflow.config` principal não é algo que você quer estar editando para cada execução.
Mas tem a vantagem de separar as preocupações de _declarar_ os parâmetros no fluxo de trabalho (que definitivamente pertence lá) versus fornecer _valores padrão_, que estão mais em casa em um arquivo de configuração.

Vamos fazer isso em dois passos.

#### 1.1.1. Crie um bloco `params` no arquivo de configuração

Faça as seguintes alterações de código no arquivo `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
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
A sintaxe é um pouco diferente.
No arquivo de fluxo de trabalho, essas são declarações tipadas.
Na configuração, essas são atribuições de valores.

Tecnicamente, isso é suficiente para substituir os valores padrão ainda especificados no arquivo de fluxo de trabalho.
Você poderia modificar o caractere, por exemplo, e executar o fluxo de trabalho para se convencer de que o valor definido no arquivo de configuração substitui o definido no arquivo de fluxo de trabalho.

Mas no espírito de mover a configuração completamente para o arquivo de configuração, vamos remover esses valores do arquivo de fluxo de trabalho inteiramente.

#### 1.1.2. Remova os valores do bloco `params` no arquivo de fluxo de trabalho

Faça as seguintes alterações de código no arquivo de fluxo de trabalho `hello-config.nf`:

=== "Depois"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Antes"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Agora o arquivo de fluxo de trabalho em si não define nenhum valor padrão para esses parâmetros.

#### 1.1.3. Execute o pipeline

Vamos testar se funciona corretamente.

```bash
nextflow run hello-config.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Isso ainda produz a mesma saída de antes.

A saída final em arte ASCII está no diretório `results/hello_config/`, com o nome `cowpy-COLLECTED-batch-output.txt`, igual a antes.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

### 1.2. Use um arquivo de configuração específico para a execução

Isso é ótimo, mas às vezes você pode querer executar alguns experimentos temporários com diferentes valores padrão sem mexer no arquivo de configuração principal.
Você pode fazer isso criando um novo arquivo `nextflow.config` em um subdiretório que você usará como diretório de trabalho para seus experimentos.

#### 1.2.1. Crie o diretório de trabalho com uma configuração em branco

Vamos começar criando um novo diretório e entrando nele:

```bash
mkdir -p tux-run
cd tux-run
```

Em seguida, crie um arquivo de configuração em branco nesse diretório:

```bash
touch nextflow.config
```

Isso produz um arquivo vazio.

#### 1.2.2. Configure a configuração experimental

Agora abra o novo arquivo e adicione os parâmetros que você deseja personalizar:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Note que o caminho para o arquivo de entrada deve refletir a estrutura de diretórios.

#### 1.2.3. Execute o pipeline

Agora podemos executar nosso pipeline de dentro do nosso novo diretório de trabalho.
Certifique-se de adaptar o caminho adequadamente!

```bash
nextflow run ../hello-config.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Isso criará um novo conjunto de diretórios em `tux-run/` incluindo `tux-run/work/` e `tux-run/results/`.

Nesta execução, o Nextflow combina o `nextflow.config` no nosso diretório atual com o `nextflow.config` no diretório raiz do pipeline, e assim substitui o caractere padrão (turkey) pelo caractere tux.

O arquivo de saída final deve conter o caractere tux dizendo as saudações.

??? abstract "Conteúdo do arquivo"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

### 1.3. Use um arquivo de parâmetros

A abordagem de subdiretório funciona muito bem para experimentar, mas envolve um pouco de configuração e requer que você adapte os caminhos adequadamente.
Existe uma abordagem mais simples para quando você quer executar seu pipeline com um conjunto específico de valores, ou permitir que outra pessoa faça isso com o mínimo de esforço.

O Nextflow nos permite especificar parâmetros via um [arquivo de parâmetros](https://nextflow.io/docs/latest/config.html#params-file) no formato YAML ou JSON, o que torna muito conveniente gerenciar e distribuir conjuntos alternativos de valores padrão, por exemplo, assim como valores de parâmetros específicos da execução.

#### 1.3.1. Examine o arquivo de parâmetros de exemplo

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

    Também fornecemos uma versão JSON do arquivo de parâmetros como exemplo, mas não vamos executar com ela aqui.
    Sinta-se livre para tentar essa por conta própria.

#### 1.3.2. Execute o pipeline

Para executar o fluxo de trabalho com este arquivo de parâmetros, simplesmente adicione `-params-file <filename>` ao comando base.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

O arquivo de saída final deve conter o caractere stegosaurus dizendo as saudações.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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
Nesses casos, usar um arquivo de parâmetros nos permitirá fornecer valores de parâmetros no tempo de execução sem ter que digitar linhas de comando enormes e sem modificar o script do fluxo de trabalho.

Também torna mais fácil distribuir conjuntos de parâmetros para colaboradores, ou como informação de apoio para uma publicação, por exemplo.
Isso torna seu trabalho mais reproduzível por outros.

### Conclusão

Você sabe como aproveitar as principais opções de configuração para gerenciar entradas do fluxo de trabalho.

### O que vem a seguir?

Aprenda como gerenciar onde e como as saídas do seu fluxo de trabalho são publicadas.

---

## 2. Gerencie saídas do fluxo de trabalho

Até agora temos codificado todos os caminhos para declarações de saída no nível do fluxo de trabalho, e como notamos quando começamos a adicionar múltiplas saídas, pode haver um pouco de repetição envolvida.

Vamos ver algumas maneiras comuns de configurar isso para ser mais flexível.

### 2.1. Personalize o diretório de saída com `-output-dir`

Quando estamos controlando como nossas saídas 'publicadas' são organizadas, temos duas prioridades distintas:

- O diretório de saída de nível superior
- Como os arquivos são organizados dentro deste diretório

Temos usado o diretório de nível superior padrão até agora: `results`.
Vamos começar personalizando isso, usando a opção CLI `-output-dir`.

#### 2.1.1. Execute o pipeline com `-output-dir`

A opção `-output-dir` (abreviação: `-o`) substitui o diretório de saída padrão (`results/`) para todas as saídas do fluxo de trabalho.
Esta é a maneira recomendada de controlar o caminho raiz onde as saídas são publicadas.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

Isso publica saídas em `custom-outdir-cli/` em vez de `results/`:

??? abstract "Conteúdo do diretório"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Note que ainda temos o subdiretório `hello_config` das declarações `path` no bloco de saída.
Vamos limpar isso.

#### 2.1.2. Remova caminhos codificados do bloco de saída

O prefixo `hello_config/` foi codificado em capítulos anteriores, mas como agora estamos aprendendo a configurar caminhos de saída de forma flexível, podemos remover essa codificação.
Para saídas que não precisam de um subdiretório, podemos definir a diretiva `path` como uma string vazia, ou removê-la completamente.

Faça as seguintes alterações de código no arquivo de fluxo de trabalho:

=== "Depois"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Execute o pipeline novamente:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Agora as saídas são publicadas diretamente em `custom-outdir-cli-2/`, sem o subdiretório `hello_config`:

??? abstract "Conteúdo do diretório"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip

    A opção `-output-dir` é usada para controlar _onde_ as saídas vão, enquanto a diretiva `path` no bloco de saída controla a _estrutura de subdiretórios_.

### 2.2. Caminhos de saída dinâmicos

Além de mudar o diretório de saída via CLI, também podemos definir um valor padrão personalizado no arquivo de configuração usando `outputDir`.
Isso nos permite definir o caminho do diretório dinamicamente - não apenas usando strings estáticas.

#### 2.2.1. Defina `outputDir` no arquivo de configuração

Adicione o seguinte código ao arquivo `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Isso define o diretório de saída como `custom-outdir-config/` mais o valor do parâmetro `batch` como subdiretório.
Agora você pode mudar o local de saída definindo o parâmetro `--batch`:

```bash
nextflow run hello-config.nf --batch my_run
```

Isso publica saídas em `custom-outdir-config/my_run/`.

!!! note

    A opção CLI `-output-dir` tem precedência sobre a configuração `outputDir`.
    Se estiver definida, a opção de configuração será ignorada completamente.

#### 2.2.2. Subdiretórios com nomes de lote e processo

Também podemos definir declarações de `path` de saída de subdiretórios dinamicamente, por saída.

Por exemplo, podemos organizar nossas saídas por processo referenciando `<processo>.name` na declaração do caminho de saída:

=== "Depois"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Podemos ir além e compor caminhos de subdiretórios mais complexos.

Na edição acima, apagamos a distinção entre `intermediates` versus saídas finais estando no nível superior.
Vamos colocar isso de volta, e também colocar os arquivos em um subdiretório `params.batch`.

!!! tip

    Incluir `params.batch` no `path` do bloco de saída, em vez do `outputDir` de configuração, significa que não será sobrescrito com `-output-dir` no CLI.

Primeiro, atualize o arquivo de configuração para remover `${params.batch}` do `outputDir` (já que estamos movendo-o para as declarações de caminho):

=== "Depois"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Então, faça as seguintes alterações no arquivo de fluxo de trabalho:

=== "Depois"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. Execute o pipeline

Vamos ver como isso funciona na prática, definindo tanto `-output-dir` (ou `-o` para abreviar) como `custom-outdir-config-2` quanto o nome do lote como `rep2` a partir da linha de comando:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

Isso publica saídas em `custom-outdir-config-2/rep2/`, com o caminho base especificado _e_ o subdiretório do nome do lote _e_ resultados agrupados por processo:

??? abstract "Conteúdo do diretório"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Defina o modo de publicação no nível do fluxo de trabalho

Finalmente, no espírito de reduzir a quantidade de código repetitivo, podemos substituir as declarações `mode` por saída com uma única linha na configuração.

#### 2.3.1. Adicione `workflow.output.mode` ao arquivo de configuração

Adicione o seguinte código ao arquivo `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    ```

Definir `workflow.output.mode` no arquivo de configuração é suficiente para substituir o que está definido no arquivo de fluxo de trabalho, mas vamos remover o código desnecessário de qualquer forma.

#### 2.3.2. Remova o modo de saída do arquivo de fluxo de trabalho

Faça as seguintes alterações no arquivo de fluxo de trabalho:

=== "Depois"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

Isso é mais conciso, não é?

#### 2.3.3. Execute o pipeline

Vamos testar se funciona corretamente:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

Isso publica saídas em `config-output-mode/`, e elas ainda são todas cópias adequadas, não symlinks.

??? abstract "Conteúdo do diretório"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

A principal razão pela qual você ainda pode querer usar a maneira por saída de definir o modo é se você quiser misturar e combinar dentro do mesmo fluxo de trabalho, _ou seja_, ter algumas saídas sendo copiadas e algumas sendo symlinkadas.

Existem muitas outras opções que você pode personalizar dessa maneira, mas esperamos que isso lhe dê uma noção da gama de opções e como utilizá-las efetivamente para atender às suas preferências.

### Conclusão

Você sabe como controlar a nomenclatura e estrutura dos diretórios onde suas saídas são publicadas, bem como o modo de publicação de saída do fluxo de trabalho.

### O que vem a seguir?

Aprenda como adaptar a configuração do seu fluxo de trabalho ao seu ambiente de computação, começando com a tecnologia de empacotamento de software.

---

## 3. Selecione uma tecnologia de empacotamento de software

Até agora temos visto elementos de configuração que controlam como as entradas entram e onde as saídas saem. Agora é hora de focar mais especificamente em adaptar a configuração do seu fluxo de trabalho ao seu ambiente de computação.

O primeiro passo nesse caminho é especificar de onde virão os pacotes de software que serão executados em cada etapa.
Eles já estão instalados no ambiente de computação local?
Precisamos recuperar imagens e executá-las via um sistema de contêineres?
Ou precisamos recuperar pacotes Conda e construir um ambiente Conda local?

Na primeira parte deste curso de treinamento (Partes 1-4) apenas usamos software instalado localmente em nosso fluxo de trabalho.
Então na Parte 5, introduzimos contêineres Docker e o arquivo `nextflow.config`, que usamos para habilitar o uso de contêineres Docker.

Agora vamos ver como podemos configurar uma opção alternativa de empacotamento de software via o arquivo `nextflow.config`.

### 3.1. Desabilite o Docker e habilite o Conda no arquivo de configuração

Vamos fingir que estamos trabalhando em um cluster HPC e o administrador não permite o uso do Docker por razões de segurança.
Felizmente para nós, o Nextflow suporta múltiplas outras tecnologias de contêineres, incluindo Singularity (que é mais amplamente usado em HPC), e gerenciadores de pacotes de software como Conda.

Podemos mudar nosso arquivo de configuração para usar [Conda](https://nextflow.io/docs/latest/conda.html) em vez de Docker.
Para fazer isso, vamos mudar o valor de `docker.enabled` para `false`, e adicionar uma diretiva habilitando o uso do Conda:

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

### 3.2. Especifique um pacote Conda na definição do processo

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
    Recomendamos usar a consulta de busca [Seqera Containers](https://seqera.io/containers/), que lhe dará um URI que você pode copiar e colar, mesmo que você não esteja planejando criar um contêiner a partir dele.

### 3.3. Execute o fluxo de trabalho para verificar que ele pode usar Conda

Vamos experimentar.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Saída do comando"

    ```console title="Saída"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

Isso deve funcionar sem problemas e produzir as mesmas saídas de antes em `custom-outdir-config/conda`.

Nos bastidores, o Nextflow recuperou os pacotes Conda e criou o ambiente, o que normalmente requer um pouco de trabalho; então é bom que não tenhamos que fazer nada disso nós mesmos!

!!! note

    Isso executa rapidamente porque o pacote `cowpy` é bastante pequeno, mas se você estiver trabalhando com pacotes grandes, pode levar um pouco mais de tempo do que o normal na primeira vez, e você pode ver a saída do console ficar 'presa' por um minuto ou mais antes de completar.
    Isso é normal e se deve ao trabalho extra que o Nextflow faz na primeira vez que você usa um novo pacote.

Do nosso ponto de vista, parece que funciona exatamente da mesma forma que executar com Docker, embora no backend a mecânica seja um pouco diferente.

Isso significa que estamos todos prontos para executar com ambientes Conda se necessário.

??? info "Misturando e combinando Docker e Conda"

    Como essas diretivas são atribuídas por processo, é possível 'misturar e combinar', _ou seja_, configurar alguns dos processos em seu fluxo de trabalho para executar com Docker e outros com Conda, por exemplo, se a infraestrutura de computação que você está usando suporta ambos.
    Nesse caso, você habilitaria tanto Docker quanto Conda no seu arquivo de configuração.
    Se ambos estiverem disponíveis para um determinado processo, o Nextflow priorizará contêineres.

    E como observado anteriormente, o Nextflow suporta múltiplas outras tecnologias de empacotamento de software e contêineres, então você não está limitado a apenas essas duas.

### Conclusão

Você sabe como configurar qual pacote de software cada processo deve usar, e como alternar entre tecnologias.

### O que vem a seguir?

Aprenda como mudar a plataforma de execução usada pelo Nextflow para realmente fazer o trabalho.

---

## 4. Selecione uma plataforma de execução

Até agora, temos executado nosso pipeline com o executor local.
Isso executa cada tarefa na máquina em que o Nextflow está sendo executado.
Quando o Nextflow começa, ele verifica as CPUs e memória disponíveis.
Se os recursos das tarefas prontas para executar excedem os recursos disponíveis, o Nextflow manterá as últimas tarefas de volta da execução até que uma ou mais das tarefas anteriores tenham terminado, liberando os recursos necessários.

O executor local é conveniente e eficiente, mas é limitado àquela única máquina. Para cargas de trabalho muito grandes, você pode descobrir que sua máquina local é um gargalo, seja porque você tem uma única tarefa que requer mais recursos do que você tem disponíveis, ou porque você tem tantas tarefas que esperar por uma única máquina para executá-las levaria muito tempo.

O Nextflow suporta [muitos executores diferentes](https://nextflow.io/docs/latest/executor.html), incluindo agendadores HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor e outros), bem como backends de execução em nuvem (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e mais).

### 4.1. Direcionando um backend diferente

A escolha do executor é definida por uma diretiva de processo chamada `executor`.
Por padrão, é definido como `local`, então a seguinte configuração está implícita:

```groovy title="Configuração integrada"
process {
    executor = 'local'
}
```

Para definir o executor para direcionar um backend diferente, você simplesmente especificaria o executor que deseja usando sintaxe similar à descrita acima para alocações de recursos (veja a [documentação de executores](https://nextflow.io/docs/latest/executor.html) para todas as opções).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Não podemos realmente testar isso no ambiente de treinamento porque não está configurado para se conectar a um HPC.

### 4.2. Lidando com sintaxe específica do backend para parâmetros de execução

A maioria das plataformas de computação de alto desempenho permite (e às vezes exige) que você especifique certos parâmetros, como solicitações e limitações de alocação de recursos (por exemplo, número de CPUs e memória) e nome da fila de trabalhos a ser usada.

Infelizmente, cada um desses sistemas usa tecnologias, sintaxes e configurações diferentes para definir como um trabalho deve ser definido e submetido ao agendador relevante.

??? abstract "Exemplos"

    Por exemplo, o mesmo trabalho requerendo 8 CPUs e 4GB de RAM para ser executado na fila "my-science-work" precisa ser expresso das seguintes maneiras diferentes dependendo do backend.

    ```bash title="Config para SLURM / submeter usando sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config para PBS / submeter usando qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config para SGE / submeter usando qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Felizmente, o Nextflow simplifica tudo isso.
Ele fornece uma sintaxe padronizada para que você possa especificar as propriedades relevantes como [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) e [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (veja [diretivas de processo](https://nextflow.io/docs/latest/reference/process.html#process-directives) para outras propriedades) apenas uma vez.
Então, no tempo de execução, o Nextflow usará essas configurações para gerar os scripts específicos do backend apropriados com base na configuração do executor.

Vamos cobrir essa sintaxe padronizada na próxima seção.

### Conclusão

Você agora sabe como mudar o executor para usar diferentes tipos de infraestrutura de computação.

### O que vem a seguir?

Aprenda como avaliar e expressar alocações e limitações de recursos no Nextflow.

---

## 5. Controle alocações de recursos de computação

A maioria das plataformas de computação de alto desempenho permite (e às vezes exige) que você especifique certos parâmetros de alocação de recursos, como número de CPUs e memória.

Por padrão, o Nextflow usará uma única CPU e 2GB de memória para cada processo.
As diretivas de processo correspondentes são chamadas `cpus` e `memory`, então a seguinte configuração está implícita:

```groovy title="Configuração integrada" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Você pode modificar esses valores, seja para todos os processos ou para processos nomeados específicos, usando diretivas de processo adicionais no seu arquivo de configuração.
O Nextflow os traduzirá nas instruções apropriadas para o executor escolhido.

Mas como você sabe quais valores usar?

### 5.1. Execute o fluxo de trabalho para gerar um relatório de utilização de recursos

Se você não sabe de antemão quanta CPU e memória seus processos provavelmente precisarão, você pode fazer algum perfil de recursos, o que significa que você executa o fluxo de trabalho com algumas alocações padrão, registra quanto cada processo usou e, a partir daí, estima como ajustar as alocações base.

Convenientemente, o Nextflow inclui ferramentas integradas para fazer isso, e irá gerar alegremente um relatório para você mediante solicitação.

Para fazer isso, adicione `-with-report <filename>.html` à sua linha de comando.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

O relatório é um arquivo html, que você pode baixar e abrir no seu navegador. Você também pode clicar com o botão direito nele no explorador de arquivos à esquerda e clicar em `Show preview` para visualizá-lo no ambiente de treinamento.

Reserve alguns minutos para examinar o relatório e ver se você consegue identificar algumas oportunidades para ajustar recursos.
Certifique-se de clicar nas abas que mostram os resultados de utilização como uma porcentagem do que foi alocado.

Veja [Relatórios](https://nextflow.io/docs/latest/reports.html) para documentação sobre todos os recursos disponíveis.

### 5.2. Defina alocações de recursos para todos os processos

O perfil mostra que os processos em nosso fluxo de trabalho de treinamento são muito leves, então vamos reduzir a alocação de memória padrão para 1GB por processo.

Adicione o seguinte ao seu arquivo `nextflow.config`, antes da seção de parâmetros do pipeline:

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    conda.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Isso ajudará a reduzir a quantidade de computação que consumimos.

### 5.3. Defina alocações de recursos para um processo específico

Ao mesmo tempo, vamos fingir que o processo `cowpy` requer mais recursos do que os outros, apenas para que possamos demonstrar como ajustar alocações para um processo individual.

=== "Depois"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
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
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Com esta configuração, todos os processos solicitarão 1GB de memória e uma única CPU (o padrão implícito), exceto o processo `cowpy`, que solicitará 2GB e 2 CPUs.

!!! tip

    Se você tiver uma máquina com poucas CPUs e alocar um número alto por processo, poderá ver chamadas de processo sendo enfileiradas uma atrás da outra.
    Isso ocorre porque o Nextflow garante que não solicitemos mais CPUs do que as disponíveis.

### 5.4. Execute o fluxo de trabalho com a configuração atualizada

Vamos experimentar isso, fornecendo um nome de arquivo diferente para o relatório de perfil para que possamos comparar o desempenho antes e depois das mudanças de configuração.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Você provavelmente não notará nenhuma diferença real, pois esta é uma carga de trabalho tão pequena, mas esta é a abordagem que você usaria para analisar o desempenho e os requisitos de recursos de um fluxo de trabalho do mundo real.

É muito útil quando seus processos têm requisitos de recursos diferentes. Ele o capacita a dimensionar corretamente as alocações de recursos que você configura para cada processo com base em dados reais, não em suposições.

!!! tip

    Este é apenas um pequeno aperitivo do que você pode fazer para otimizar seu uso de recursos.
    O próprio Nextflow tem uma [lógica de repetição dinâmica](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) realmente interessante embutida para repetir trabalhos que falham devido a limitações de recursos.
    Além disso, a Seqera Platform oferece ferramentas orientadas por IA para otimizar suas alocações de recursos automaticamente também.

### 5.5. Adicione limites de recursos

Dependendo de qual executor de computação e infraestrutura de computação você está usando, pode haver algumas restrições sobre o que você pode (ou deve) alocar.
Por exemplo, seu cluster pode exigir que você permaneça dentro de certos limites.

Você pode usar a diretiva `resourceLimits` para definir as limitações relevantes. A sintaxe se parece com isso quando está sozinha em um bloco de processo:

```groovy title="Exemplo de sintaxe"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

O Nextflow traduzirá esses valores nas instruções apropriadas dependendo do executor que você especificou.

Não vamos executar isso, pois não temos acesso à infraestrutura relevante no ambiente de treinamento.
No entanto, se você tentasse executar o fluxo de trabalho com alocações de recursos que excedem esses limites, então procurasse o comando `sbatch` no arquivo de script `.command.run`, você veria que as solicitações que realmente são enviadas ao executor são limitadas aos valores especificados por `resourceLimits`.

??? info "Configurações de referência institucionais"

    O projeto nf-core compilou uma [coleção de arquivos de configuração](https://nf-co.re/configs/) compartilhados por várias instituições ao redor do mundo, cobrindo uma ampla gama de executores HPC e nuvem.

    Essas configurações compartilhadas são valiosas tanto para pessoas que trabalham lá e, portanto, podem simplesmente utilizar a configuração de sua instituição pronta para uso, quanto como modelo para pessoas que estão procurando desenvolver uma configuração para sua própria infraestrutura.

### Conclusão

Você sabe como gerar um relatório de perfil para avaliar a utilização de recursos e como modificar alocações de recursos para todos os processos e/ou para processos individuais, bem como definir limitações de recursos para executar em HPC.

### O que vem a seguir?

Aprenda como configurar perfis de configuração predefinidos e alternar entre eles no tempo de execução.

---

## 6. Use perfis para alternar entre configurações predefinidas

Mostramos a você várias maneiras de personalizar a configuração do seu pipeline dependendo do projeto em que você está trabalhando ou do ambiente de computação que está usando.

Você pode querer alternar entre configurações alternativas dependendo de qual infraestrutura de computação está usando. Por exemplo, você pode querer desenvolver e executar testes em pequena escala localmente no seu laptop, depois executar cargas de trabalho em escala completa em HPC ou nuvem.

O Nextflow permite que você configure qualquer número de [perfis](https://nextflow.io/docs/latest/config.html#config-profiles) que descrevem diferentes configurações, que você pode então selecionar no tempo de execução usando um argumento de linha de comando, em vez de ter que modificar o arquivo de configuração em si.

### 6.1. Crie perfis para alternar entre desenvolvimento local e execução em HPC

Vamos configurar dois perfis alternativos; um para executar cargas em pequena escala em um computador normal, onde usaremos contêineres Docker, e outro para executar em um HPC universitário com um agendador Slurm, onde usaremos pacotes Conda.

#### 6.1.1. Configure os perfis

Adicione o seguinte ao seu arquivo `nextflow.config`, após a seção de parâmetros do pipeline, mas antes das configurações de saída:

=== "Depois"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Profiles
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

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Você vê que para o HPC universitário, também estamos especificando limitações de recursos.

#### 6.1.2. Execute o fluxo de trabalho com um perfil

Para especificar um perfil na nossa linha de comando do Nextflow, usamos o argumento `-profile`.

Vamos tentar executar o fluxo de trabalho com a configuração `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Como você pode ver, isso nos permite alternar entre configurações muito convenientemente no tempo de execução.

!!! warning

    O perfil `univ_hpc` não será executado corretamente no ambiente de treinamento, pois não temos acesso a um agendador Slurm.

Se no futuro encontrarmos outros elementos de configuração que estão sempre co-ocorrendo com esses, podemos simplesmente adicioná-los ao(s) perfil(is) correspondente(s).
Também podemos criar perfis adicionais se houver outros elementos de configuração que queremos agrupar.

### 6.2. Crie um perfil de parâmetros de teste

Perfis não são apenas para configuração de infraestrutura.
Também podemos usá-los para definir valores padrão para parâmetros de fluxo de trabalho, para facilitar que outros experimentem o fluxo de trabalho sem ter que reunir valores de entrada apropriados por conta própria.
Você pode considerar isso uma alternativa ao uso de um arquivo de parâmetros.

#### 6.2.1. Configure o perfil

A sintaxe para expressar valores padrão neste contexto se parece com isso, para um perfil que nomeamos `test`:

```groovy title="Exemplo de sintaxe"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Se adicionarmos um perfil de teste para nosso fluxo de trabalho, o bloco `profiles` se torna:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* Profiles
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

Assim como para perfis de configuração técnica, você pode configurar vários perfis diferentes especificando parâmetros sob qualquer nome arbitrário que desejar.

#### 6.2.2. Execute o fluxo de trabalho localmente com o perfil de teste

Convenientemente, perfis não são mutuamente exclusivos, então podemos especificar múltiplos perfis em nossa linha de comando usando a seguinte sintaxe `-profile <profile1>,<profile2>` (para qualquer número de perfis).

Se você combinar perfis que definem valores para os mesmos elementos de configuração e estão descritos no mesmo arquivo de configuração, o Nextflow resolverá o conflito usando qualquer valor que ele leu por último (_ou seja_, o que vem depois no arquivo).
Se as configurações conflitantes são definidas em diferentes fontes de configuração, a [ordem de precedência](https://nextflow.io/docs/latest/config.html) padrão se aplica.

Vamos tentar adicionar o perfil de teste ao nosso comando anterior:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

Isso usará Docker onde possível e produzirá saídas em `custom-outdir-config/test`, e desta vez o caractere é a dupla cômica `dragonandcow`.

??? abstract "Conteúdo do arquivo"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

Isso significa que, desde que distribuamos quaisquer arquivos de dados de teste com o código do fluxo de trabalho, qualquer pessoa pode rapidamente experimentar o fluxo de trabalho sem ter que fornecer suas próprias entradas via linha de comando ou arquivo de parâmetros.

!!! tip

    Podemos apontar para URLs para arquivos maiores que são armazenados externamente.
    O Nextflow os baixará automaticamente desde que haja uma conexão aberta.

    Para mais detalhes, veja a Side Quest [Trabalhando com Arquivos](../side_quests/working_with_files.md)

### 6.3. Use `nextflow config` para ver a configuração resolvida

Como observado acima, às vezes o mesmo parâmetro pode ser definido com valores diferentes em perfis que você deseja combinar.
E de forma mais geral, existem numerosos lugares onde elementos de configuração podem ser armazenados, e às vezes as mesmas propriedades podem ser definidas com valores diferentes em lugares diferentes.

O Nextflow aplica uma [ordem de precedência](https://nextflow.io/docs/latest/config.html) definida para resolver quaisquer conflitos, mas isso pode ser difícil de determinar por conta própria.
E mesmo que nada esteja conflitando, pode ser tedioso procurar em todos os lugares possíveis onde as coisas poderiam estar configuradas.

Felizmente, o Nextflow inclui uma ferramenta de utilitário conveniente chamada `config` que pode automatizar todo esse processo para você.

A ferramenta `config` explorará todo o conteúdo no seu diretório de trabalho atual, coletará quaisquer arquivos de configuração e produzirá a configuração totalmente resolvida que o Nextflow usaria para executar o fluxo de trabalho.
Isso permite que você descubra quais configurações serão usadas sem ter que lançar nada.

#### 6.3.1. Resolva a configuração padrão

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Isso mostra a configuração base que você obtém se não especificar nada extra na linha de comando.

#### 6.3.2. Resolva a configuração com configurações específicas ativadas

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Isso fica especialmente útil para projetos complexos que envolvem múltiplas camadas de configuração.

### Conclusão

Você sabe como usar perfis para selecionar uma configuração predefinida no tempo de execução com o mínimo de esforço.
De forma mais geral, você sabe como configurar suas execuções de fluxo de trabalho para se adequar a diferentes plataformas de computação e melhorar a reprodutibilidade de suas análises.

### O que vem a seguir?

Comemore e dê um grande tapinha nas costas! Você completou seu primeiro curso de desenvolvedor Nextflow.

Vá para o [resumo final do curso](./next_steps.md) para revisar o que você aprendeu e descobrir o que vem a seguir.

---

## Quiz

<quiz>
Qual é o nome do arquivo de configuração que o Nextflow carrega automaticamente?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
O que tem precedência quando o mesmo parâmetro é definido tanto no arquivo de configuração quanto na linha de comando?
- [ ] O valor do arquivo de configuração
- [x] O valor da linha de comando
- [ ] O primeiro valor encontrado
- [ ] Nenhum; causa um erro

Saiba mais: [1.1. Mova os valores padrão para o `nextflow.config`](#11-mova-os-valores-padrão-para-o-nextflowconfig)
</quiz>

<quiz>
Você pode ter tanto Docker quanto Conda habilitados na mesma configuração?
- [x] Sim, o Nextflow pode usar ambos dependendo das diretivas de processo
- [ ] Não, apenas um pode ser habilitado por vez
- [ ] Sim, mas apenas em perfis
- [ ] Não, eles são mutuamente exclusivos
</quiz>

<quiz>
Se tanto Docker quanto Conda estão habilitados e um processo tem ambas as diretivas, qual é priorizado?
- [x] Docker (contêineres)
- [ ] Conda
- [ ] O primeiro definido
- [ ] Causa um erro

Saiba mais: [3. Selecione uma tecnologia de empacotamento de software](#3-selecione-uma-tecnologia-de-empacotamento-de-software)
</quiz>

<quiz>
Qual é a alocação de memória padrão para processos Nextflow?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Sem limite
</quiz>

<quiz>
Como você define requisitos de recursos para um processo específico no arquivo de configuração?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Saiba mais: [5.3. Defina alocações de recursos para um processo específico](#53-defina-alocações-de-recursos-para-um-processo-específico)
</quiz>

<quiz>
Qual opção de linha de comando gera um relatório de utilização de recursos?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Saiba mais: [5.1. Execute o fluxo de trabalho para gerar um relatório de utilização de recursos](#51-execute-o-fluxo-de-trabalho-para-gerar-um-relatório-de-utilização-de-recursos)
</quiz>

<quiz>
O que a diretiva `resourceLimits` faz?
- [ ] Define requisitos mínimos de recursos
- [ ] Aloca recursos aos processos
- [x] Limita os recursos máximos que podem ser solicitados
- [ ] Monitora o uso de recursos

Saiba mais: [5.5. Adicione limites de recursos](#55-adicione-limites-de-recursos)
</quiz>

<quiz>
Qual é o executor padrão no Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Saiba mais: [4. Selecione uma plataforma de execução](#4-selecione-uma-plataforma-de-execução)
</quiz>

<quiz>
Como você especifica um arquivo de parâmetros ao executar o Nextflow?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Saiba mais: [1.3. Use um arquivo de parâmetros](#13-use-um-arquivo-de-parâmetros)
</quiz>

<quiz>
Para que os perfis podem ser usados? (Selecione todas as opções que se aplicam)
- [x] Definir configurações específicas de infraestrutura
- [x] Definir limites de recursos para diferentes ambientes
- [x] Fornecer parâmetros de teste
- [ ] Definir novos processos

Saiba mais: [6. Use perfis para alternar entre configurações predefinidas](#6-use-perfis-para-alternar-entre-configurações-predefinidas)
</quiz>

<quiz>
Como você especifica múltiplos perfis em um único comando?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Saiba mais: [6. Use perfis para alternar entre configurações predefinidas](#6-use-perfis-para-alternar-entre-configurações-predefinidas)
</quiz>
