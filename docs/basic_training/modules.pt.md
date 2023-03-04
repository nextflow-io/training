---
title: Módulos
description: Material de treinamento básico do Nextflow
---

# Modularização

A definição de bibliotecas modulares simplifica a escrita de pipelines complexos de análise de dados, além tornar o reuso de processos mais fácil.

Ao usar o exemplo `hello.nf` da seção de introdução, nós converteremos os processos do pipeline em módulos e, em seguida, executaremos estes processos dentro do escopo do workflow de diferentes formas.

## Módulos

A DSL2 do Nextflow permite a definição de scripts de módulos autônomos que podem ser incluídos e compartilhados em vários fluxos de trabalho. Cada módulo pode conter sua própria definição de `process` ou `workflow`.

### Importando módulos

Os componentes definidos no script do módulo podem ser importados para outros scripts do Nextflow usando a instrução `include`. Isso permite que você armazene esses componentes em arquivos separados para que possam ser reutilizados em vários fluxos de trabalho.

Usando o exemplo `hello.nf`, podemos fazer isso:

-   Criando um arquivo chamado `modules.nf` no diretório de nível superior.
-   Recortando e colando as duas definições de processo para `SPLITLETTERS` e `CONVERTTOUPPER` em `modules.nf`.
-   Removendo as definições `process` no script `hello.nf`.
-   Importando os processos de `modules.nf` dentro do script `hello.nf` em qualquer lugar acima da definição de `workflow`:

```groovy linenums="1"
include { SPLITLETTERS } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'
```

!!! note

    Em geral, você deve usar caminhos relativos para definir a localização dos scripts do módulo usando o prefixo `./`.

!!! exercise

    Crie um arquivo `modules.nf` com os processos previamente definidos de `hello.nf`. Em seguida, remova esses processos de `hello.nf` e adicione as definições `include` mostradas acima.

    ??? solution

        O script `hello.nf` deve ser similar a este:

        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        params.greeting  = 'Hello world!'
        greeting_ch = Channel.of(params.greeting)

        include { SPLITLETTERS   } from './modules.nf'
        include { CONVERTTOUPPER } from './modules.nf'

        workflow {
            letters_ch = SPLITLETTERS(greeting_ch)
            results_ch = CONVERTTOUPPER(letters_ch.flatten())
            results_ch.view{ it }
        }
        ```

        Você deve ter o seguinte código em `./modules.nf`:

        ```groovy linenums="1"
        process SPLITLETTERS {

            input:
            val x

            output:
            path 'chunk_*'

            """
            printf '$x' | split -b 6 - chunk_
            """
        }

        process CONVERTTOUPPER {

            input:
            path y

            output:
            stdout

            """
            cat $y | tr '[a-z]' '[A-Z]'
            """
        }
        ```

        Agora nós modularizamos os processos, o que faz com que o código seja reutilizável.

### Importações múltiplas

Se um script de módulo Nextflow contiver várias definições de `process`, elas também podem ser importadas usando uma única instrução `include`, conforme mostrado no exemplo abaixo:

```groovy linenums="1"
include { SPLITLETTERS; CONVERTTOUPPER } from './modules.nf'
```

### Apelidos dos módulos

Ao incluir um componente de um módulo, é possível especificar um apelido para os processos usando a declaração `as`. Isso permite a inclusão e a invocação do mesmo componente várias vezes usando diferentes apelidos:

```groovy linenums="1"
#!/usr/bin/env nextflow

params.greeting = 'Hello world!'
greeting_ch = Channel.of(params.greeting)

include { SPLITLETTERS as SPLITLETTERS_one } from './modules.nf'
include { SPLITLETTERS as SPLITLETTERS_two } from './modules.nf'

include { CONVERTTOUPPER as CONVERTTOUPPER_one } from './modules.nf'
include { CONVERTTOUPPER as CONVERTTOUPPER_two } from './modules.nf'

workflow {
    letters_ch1 = SPLITLETTERS_one(greeting_ch)
    results_ch1 = CONVERTTOUPPER_one(letters_ch1.flatten())
    results_ch1.view{ it }

    letters_ch2 = SPLITLETTERS_two(greeting_ch)
    results_ch2 = CONVERTTOUPPER_two(letters_ch2.flatten())
    results_ch2.view{ it }
}
```

!!! exercise

    Salve o trecho anterior como `hello.2.nf`, e preveja a saída na tela.

    ??? solution

        A saída de `hello.2.nf` deve ser semelhante a essa:

        ```console title="Output"
        N E X T F L O W  ~  version 22.04.3
        Launching `hello.2.nf` [goofy_goldstine] DSL2 - revision: 449cf82eaf
        executor >  local (6)
        [e1/5e6523] process > SPLITLETTERS_one (1)   [100%] 1 of 1 ✔
        [14/b77deb] process > CONVERTTOUPPER_one (1) [100%] 2 of 2 ✔
        [c0/115bd6] process > SPLITLETTERS_two (1)   [100%] 1 of 1 ✔
        [09/f9072d] process > CONVERTTOUPPER_two (2) [100%] 2 of 2 ✔
        WORLD!
        HELLO
        WORLD!
        HELLO
        ```

!!! tip

    Você pode armazenar cada processo em arquivos separados em subpastas separadas ou combinados em um arquivo grande (ambos são válidos).
    Você pode encontrar exemplos disso em repositórios públicos, como o [tutorial Seqera RNA-Seq](https://github.com/seqeralabs/rnaseq-nf/tree/master/modules) ou em pipelines do nf-core, como [nf-core/rnaseq](https://github.com/nf-core/rnaseq/tree/master/modules/nf-core/modules).

### Definição de saída

O Nextflow permite o uso de definições de saída alternativas em fluxos de trabalho para simplificar seu código.

No exemplo básico anterior (`hello.nf`), definimos os nomes dos canais para especificar a entrada para o próximo processo:

```groovy linenums="1"
workflow  {
    greeting_ch = Channel.of(params.greeting)
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}
```

!!! note

    Nós movemos o `greeting_ch` para o escopo do `workflow` para este exercício.

Também podemos definir explicitamente a saída de um canal para outro usando o atributo `.out`, removendo completamente as definições de canal:

```groovy linenums="1" hl_lines="3-5"
workflow  {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.view()
}
```

Se um processo define dois ou mais canais de saída, cada canal pode ser acessado indexando o atributo `.out`, por exemplo, `.out[0]`, `.out[1]`, etc. Em nosso exemplo, temos apenas a saída `[0]'th`:

```groovy linenums="1" hl_lines="5"
workflow  {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out[0].view()
}
```

Alternativamente, a definição de `output` do processo permite o uso da instrução `emit` para definir um identificador nomeado que pode ser usado para referenciar o canal no escopo externo.

Por exemplo, tente adicionar a instrução `emit` no processo `convertToUpper` em seu arquivo `modules.nf`:

```groovy linenums="1" title="modules.nf"
process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout emit: upper

    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}
```

Em seguida, altere o escopo do workflow em `hello.nf` para chamar essa saída nomeada específica (observe o `.upper` adicionado):

```groovy linenums="1" title="hello.nf"
workflow {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}
```

### Usando saídas canalizadas

Outra maneira de lidar com as saídas no escopo do workflow é usar pipes `|`.

!!! exercício

     Tente alterar o script do fluxo de trabalho para o trecho abaixo:

    ```groovy linenums="1"
    workflow {
        Channel.of(params.greeting) | SPLITLETTERS | flatten() | CONVERTTOUPPER | view
    }
    ```

     Aqui usamos um [pipe](https://www.nextflow.io/docs/latest/dsl2.html#pipes) que passa a saída de um processo como um canal para o próximo processo.

## Definição de workflow

O escopo `workflow` permite a definição de componentes que definem a invocação de um ou mais processos ou operadores:

```groovy linenums="1"
#!/usr/bin/env nextflow

params.greeting = 'Hello world!'

include { SPLITLETTERS } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'


workflow my_pipeline {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}

workflow {
    my_pipeline()
}
```

Por exemplo, o trecho acima define um `workflow` chamado `my_pipeline`, que pode ser chamado por meio de outra definição de `workflow`.

!!! note

    Certifique-se de que seu arquivo `modules.nf` é o que contém o `emit` no processo `CONVERTTOUPPER`.

!!! warning

    Um componente de um workflow pode acessar qualquer variável ou parâmetro definido no escopo externo. No exemplo em execução, também podemos acessar `params.greeting` diretamente na definição de `workflow`.

### Entradas de workflow

Um componente `workflow` pode declarar um ou mais canais de entrada usando a instrução `take`. Por exemplo:

```groovy linenums="1"
#!/usr/bin/env nextflow

params.greeting = 'Hello world!'

include { SPLITLETTERS } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'

workflow my_pipeline {
    take:
    greeting

    main:
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}
```

!!! observação

    Quando a instrução `take` é usada, a definição `workflow` precisa ser declarada dentro do bloco `main`.

A entrada para o `workflow` pode então ser especificada como um argumento:

```groovy linenums="1"
workflow {
    my_pipeline(Channel.of(params.greeting))
}
```

### Saídas do workflow

Um `workflow` pode declarar um ou mais canais de saída usando a instrução `emit`. Por exemplo:

```groovy linenums="1"
workflow my_pipeline {
    take:
    greeting

    main:
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())

    emit:
    CONVERTTOUPPER.out.upper
}

workflow {
    my_pipeline(Channel.of(params.greeting))
    my_pipeline.out.view()
}
```

Como resultado, podemos usar a notação `my_pipeline.out` para acessar as saídas de `my_pipeline` na chamada `workflow`.

```groovy linenums="1" hl_lines="10 15"
workflow my_pipeline {
    take:
    greeting

    main:
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())

    emit:
    my_data = CONVERTTOUPPER.out.upper
}

workflow {
    my_pipeline(Channel.of(params.greeting))
    my_pipeline.out.my_data.view()
}
```

O resultado do snippet acima pode ser acessado usando `my_pipeline.out.my_data`.

### Chamando workflows nomeados

Dentro de um script `main.nf` (chamado `hello.nf` em nosso exemplo), também podemos ter vários fluxos de trabalho. Nesse caso, podemos chamar um fluxo de trabalho específico ao executar o código. Para isso, usamos a chamada de ponto de entrada `-entry <workflow_name>`.

O trecho a seguir tem dois fluxos de trabalho nomeados (`my_pipeline_one` e `my_pipeline_two`):

```groovy linenums="1"
#!/usr/bin/env nextflow

params.greeting = 'Hello world!'

include { SPLITLETTERS as SPLITLETTERS_one } from './modules.nf'
include { SPLITLETTERS as SPLITLETTERS_two } from './modules.nf'

include { CONVERTTOUPPER as CONVERTTOUPPER_one } from './modules.nf'
include { CONVERTTOUPPER as CONVERTTOUPPER_two } from './modules.nf'


workflow my_pipeline_one {
    letters_ch1 = SPLITLETTERS_one(params.greeting)
    results_ch1 = CONVERTTOUPPER_one(letters_ch1.flatten())
    results_ch1.view{ it }
}

workflow my_pipeline_two {
    letters_ch2 = SPLITLETTERS_two(params.greeting)
    results_ch2 = CONVERTTOUPPER_two(letters_ch2.flatten())
    results_ch2.view{ it }
}

workflow {
    my_pipeline_one(Channel.of(params.greeting))
    my_pipeline_two(Channel.of(params.greeting))
}
```

Você pode escolher qual pipeline é executado usando o sinalizador `entry`:

```bash
nextflow run hello.2.nf -entry my_pipeline_one
```

### Escopos de parâmetros

Um script de módulo pode definir um ou mais parâmetros ou funções personalizadas usando a mesma sintaxe de qualquer outro script Nextflow. Usando os exemplos mínimos abaixo:

```groovy linenums="1" title="Script do módulo (<code>./modules.nf</code>)"
params.foo = 'Hello'
params.bar = 'world!'

def SAYHELLO() {
    println "$params.foo $params.bar"
}

```

```groovy linenums="1" title="Script principal (<code>./hello.nf</code>)"
#!/usr/bin/env nextflow

params.foo = 'Hola'
params.bar = 'mundo!'

include { SAYHELLO } from './modules.nf'

workflow {
    SAYHELLO()
}
```

A execução de `hello.nf` deve imprimir:

```console
Hola mundo!
```

Como destacado acima, o script imprimirá `Hola mundo!` em vez de `Hello world!` porque os parâmetros são herdados do contexto de inclusão.

!!! info

    Para evitar que sejam ignorados, os parâmetros do pipeline devem ser definidos no início do script antes de qualquer declaração de inclusão.

A opção `addParams` pode ser usada para estender os parâmetros do módulo sem afetar o escopo externo. Por exemplo:

```groovy linenums="1"
#!/usr/bin/env nextflow

params.foo = 'Hola'
params.bar = 'mundo!'

include { SAYHELLO } from './modules.nf' addParams(foo: 'Olá')

workflow {
    SAYHELLO()
}
```

A execução do script principal acima deve imprimir:

```console
Olá mundo!
```

## Notas de migração DSL2

Para visualizar um resumo das alterações introduzidas quando o Nextflow migrou de DSL1 para DSL2, consulte as [notas de migração DSL2](https://www.nextflow.io/docs/latest/dsl2.html#dsl2-migration-notes) na documentação oficial do Nextflow.
