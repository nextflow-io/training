# Fluxos de Trabalho de Fluxos de Trabalho

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ao desenvolver um pipeline, você frequentemente se encontra criando sequências similares de processos para diferentes tipos de dados ou etapas de análise. Você pode acabar copiando e colando essas sequências de processos, gerando código duplicado que é difícil de manter; ou pode criar um fluxo de trabalho massivo que é difícil de entender e modificar.

Um dos recursos mais poderosos do Nextflow é sua capacidade de compor pipelines complexos a partir de módulos de fluxo de trabalho menores e reutilizáveis. Essa abordagem modular torna os pipelines mais fáceis de desenvolver, testar e manter.

### Objetivos de aprendizado

Nesta side quest, vamos explorar como desenvolver módulos de fluxo de trabalho que podem ser testados e usados separadamente, compor esses módulos em um pipeline maior e gerenciar o fluxo de dados entre os módulos.

Ao final desta side quest, você será capaz de:

- Dividir pipelines complexos em unidades lógicas e reutilizáveis
- Testar cada módulo de fluxo de trabalho de forma independente
- Combinar fluxos de trabalho para criar novos pipelines
- Compartilhar módulos de fluxo de trabalho comuns entre diferentes pipelines
- Tornar seu código mais fácil de manter e entender

Essas habilidades vão ajudá-lo a construir pipelines complexos mantendo uma estrutura de código limpa e de fácil manutenção.

### Pré-requisitos

Antes de embarcar nesta side quest, você deve:

- Ter concluído o tutorial [Hello Nextflow](../../hello_nextflow/index.md) ou um curso equivalente para iniciantes.
- Estar confortável com os conceitos e mecanismos básicos do Nextflow (processos, canais, operadores, módulos)

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se ainda não tiver feito isso, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Acesse o diretório do projeto

Vamos acessar o diretório onde estão os arquivos deste tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

O editor abre com o diretório do projeto em foco.

#### Revise os materiais

Você encontrará um diretório `modules` com definições de processos, um diretório `workflows` com dois scripts de fluxo de trabalho pré-escritos e um arquivo `main.nf` que você irá atualizar progressivamente:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

O diretório `modules/` contém as definições individuais de processos, e o diretório `workflows/` contém os dois scripts de fluxo de trabalho pré-escritos com os quais você trabalhará nesta side quest.

#### Revise a tarefa

Seu desafio é montar esses módulos em dois fluxos de trabalho separados que serão compostos em um fluxo de trabalho principal:

- Um `GREETING_WORKFLOW` que valida nomes, cria saudações e adiciona timestamps
- Um `TRANSFORM_WORKFLOW` que converte texto para maiúsculas e o inverte

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar de cabeça?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho corretamente
- [ ] Entendo a tarefa

Se você conseguir marcar todas as caixas, pode começar.

---

## 1. Adicionar o fluxo de trabalho de saudação ao pipeline

O fluxo de trabalho de saudação valida nomes e gera saudações com timestamp.

### 1.1. Revisar e executar o fluxo de trabalho de saudação

Abra `workflows/greeting.nf` e examine o código:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Encadeia os processos: validar -> criar saudação -> adicionar timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Este é um fluxo de trabalho completo e independente com a mesma estrutura que você viu no tutorial 'Hello Nextflow'.
Ele define os nomes de entrada diretamente no código, encadeia três processos e publica duas saídas.

Execute-o para verificar que tudo funciona:

```bash
nextflow run workflows/greeting.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Para torná-lo composável com outros fluxos de trabalho, algumas coisas precisam mudar.

### 1.2. Tornar o fluxo de trabalho composável

Para tornar um fluxo de trabalho composável, quatro coisas precisam mudar:
o fluxo de trabalho recebe um nome, as entradas são movidas para um bloco `take:`, as saídas são movidas para um bloco `emit:`,
e os blocos independentes `publish:`/`output {}` são removidos (eles pertencem ao entry workflow).

Vamos percorrer essas mudanças uma a uma.

#### 1.2.1. Nomear o fluxo de trabalho

Dê um nome ao fluxo de trabalho para que ele possa ser importado por um fluxo de trabalho pai.

=== "Depois"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Antes"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Com um nome, o fluxo de trabalho pode ser importado em outros scripts.

#### 1.2.2. Declarar entradas com `take:`

Substitua a declaração de canal direta no código por um bloco `take:` que declara quais entradas o fluxo de trabalho espera.
O bloco `take:` vai antes de `main:`, e a linha `names_ch = channel.of(...)` é removida.

=== "Depois"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Canal de entrada com nomes

        main:
        // Encadeia os processos: validar -> criar saudação -> adicionar timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Antes"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Encadeia os processos: validar -> criar saudação -> adicionar timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

O bloco `take:` declara o canal apenas pelo nome — os detalhes do que entra nele serão definidos pelo fluxo de trabalho pai.

#### 1.2.3. Declarar saídas com `emit:`

Substitua a seção `publish:` e remova o bloco `output {}`, substituindo-os por um bloco `emit:` que nomeia as saídas.

=== "Depois"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Saudações originais
        timestamped = timestamped_ch // Saudações com timestamp
    }
    ```

=== "Antes"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

O bloco `emit:` expõe saídas nomeadas que os fluxos de trabalho pai podem acessar via `GREETING_WORKFLOW.out.greetings` e `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Verificar o resultado e testá-lo

Após as três mudanças, o arquivo completo deve ficar assim:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Canal de entrada com nomes

    main:
    // Encadeia os processos: validar -> criar saudação -> adicionar timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Saudações originais
    timestamped = timestamped_ch // Saudações com timestamp
}
```

Agora tente executá-lo diretamente:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Isso introduz um conceito importante: o **entry workflow**.
O Nextflow usa um bloco `workflow {}` sem nome como ponto de entrada quando você executa um script diretamente.
`GREETING_WORKFLOW` tem um nome, então o Nextflow não sabe como executá-lo por conta própria.

Isso é intencional — fluxos de trabalho composáveis são projetados para serem chamados a partir de um entry workflow, não executados diretamente.
A solução é um entry workflow em `main.nf` que importa e chama o `GREETING_WORKFLOW`.

### 1.3. Atualizar e testar o fluxo de trabalho principal

Agora vamos atualizar o fluxo de trabalho principal para chamar o fluxo de trabalho de saudação.

#### 1.3.1. Incluir o fluxo de trabalho de saudação e chamá-lo

Adicione a instrução `include`, atualize o corpo do fluxo de trabalho para chamar `GREETING_WORKFLOW` e substitua o placeholder `channel.empty()` em `publish:`:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Executa o fluxo de trabalho de saudação
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

O entry workflow permanece sem nome para que o Nextflow o use como ponto de entrada do pipeline.

#### 1.3.2. Atualizar o bloco de saída

Adicione uma diretiva `path` para direcionar as saudações publicadas para um subdiretório `greetings/`:

=== "Depois"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Executar o fluxo de trabalho

Execute o fluxo de trabalho para testar que funciona:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    ```

??? abstract "Conteúdo do diretório"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Os arquivos de saudação são publicados em `results/greetings/`.
O fluxo de trabalho principal chama `GREETING_WORKFLOW` e conecta sua saída diretamente à seção `publish:`.

### Conclusão

Nesta seção, você aprendeu vários conceitos importantes:

- **Fluxos de Trabalho Nomeados**: Criar um fluxo de trabalho nomeado (`GREETING_WORKFLOW`) que pode ser importado e reutilizado
- **Interfaces de Fluxo de Trabalho**: Definir entradas claras com `take:` e saídas com `emit:` para criar um fluxo de trabalho composável
- **Entry Points**: Entender que o Nextflow precisa de um entry workflow sem nome para executar um script
- **Composição de Fluxos de Trabalho**: Importar e usar um fluxo de trabalho nomeado dentro de outro fluxo de trabalho
- **Namespaces de Fluxo de Trabalho**: Acessar saídas do fluxo de trabalho usando o namespace `.out` (`GREETING_WORKFLOW.out.greetings`)

Agora você tem um fluxo de trabalho de saudação funcional que:

- Recebe um canal de nomes como entrada
- Valida cada nome
- Cria uma saudação para cada nome válido
- Adiciona timestamps às saudações
- Expõe tanto as saudações originais quanto as com timestamp como saídas

Essa abordagem modular permite que você teste o fluxo de trabalho de saudação de forma independente ou o use como componente em pipelines maiores.

---

## 2. Adicionar o fluxo de trabalho de transformação ao pipeline

O fluxo de trabalho de transformação aplica transformações de texto às saudações com timestamp.

### 2.1. Revisar e executar o fluxo de trabalho

Abra `workflows/transform.nf` e examine o código:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Aplica as transformações em sequência
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Este fluxo de trabalho independente lê os arquivos de saudação com timestamp do diretório `results/` produzido por `greeting.nf`, converte-os para maiúsculas e depois inverte o texto.

Execute-o para verificar que funciona com os resultados de saudação da seção 1.1:

```bash
nextflow run workflows/transform.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

Para torná-lo composável com o `GREETING_WORKFLOW`, as mesmas três mudanças da seção 1.2 se aplicam.

### 2.2. Torná-lo composável

Aplique as mesmas três mudanças da seção 1.2: nomeie o fluxo de trabalho, substitua a entrada direta no código por `take:`, e substitua `publish:`/`output {}` por `emit:`.

O arquivo finalizado deve ficar assim:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Canal de entrada com mensagens

    main:
    // Aplica as transformações em sequência
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Saudações em maiúsculas
    reversed = reversed_ch // Saudações em maiúsculas invertidas
}
```

O fluxo de trabalho de transformação agora é composável e está pronto para ser importado no fluxo de trabalho principal.

### 2.3. Atualizar e testar o fluxo de trabalho principal

Agora vamos atualizar o fluxo de trabalho principal para chamar o fluxo de trabalho de transformação.

#### 2.3.1. Incluir o fluxo de trabalho de transformação e chamá-lo

Adicione a instrução include, uma chamada ao `TRANSFORM_WORKFLOW` encadeada nas saudações com timestamp e as duas novas entradas em `publish:`:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Executa o fluxo de trabalho de saudação
        GREETING_WORKFLOW(names)

        // Executa o fluxo de trabalho de transformação
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Executa o fluxo de trabalho de saudação
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Isso executará o fluxo de trabalho de transformação nas saudações com timestamp.

#### 2.3.2. Atualizar o bloco de saída

Adicione as entradas `upper` e `reversed` ao bloco `output {}`, cada uma com uma diretiva `path` para seu subdiretório:

=== "Depois"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Isso publicará as saídas finais nos diretórios apropriados.

#### 2.3.3. Executar o pipeline completo

Execute o pipeline para testar que tudo funciona:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Conteúdo do diretório"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

O pipeline está funcionando de ponta a ponta: a saudação foi convertida para maiúsculas e invertida.

### Conclusão

Agora você deve ter um pipeline completo que:

- Processa nomes através do fluxo de trabalho de saudação
- Alimenta as saudações com timestamp no fluxo de trabalho de transformação
- Produz versões em maiúsculas e invertidas das saudações

---

## Resumo

Nesta side quest, exploramos o poderoso conceito de composição de fluxos de trabalho no Nextflow, que nos permite construir pipelines complexos a partir de componentes menores e reutilizáveis.

Essa abordagem modular oferece várias vantagens em relação a pipelines monolíticos:

- Cada fluxo de trabalho pode ser desenvolvido, testado e depurado de forma independente
- Os fluxos de trabalho podem ser reutilizados em diferentes pipelines
- A estrutura geral do pipeline se torna mais legível e fácil de manter
- Mudanças em um fluxo de trabalho não afetam necessariamente os outros, desde que as interfaces permaneçam consistentes
- Os entry points podem ser configurados para executar diferentes partes do seu pipeline conforme necessário

_É importante notar, no entanto, que embora chamar fluxos de trabalho seja um pouco parecido com chamar processos, não é exatamente a mesma coisa. Você não pode, por exemplo, executar um fluxo de trabalho N vezes chamando-o com um canal de tamanho N — você precisaria passar um canal de tamanho N para o fluxo de trabalho e iterar internamente._

Aplicar essas técnicas no seu próprio trabalho permitirá que você construa pipelines Nextflow mais sofisticados, capazes de lidar com tarefas complexas de processamento de dados, mantendo-os fáceis de manter e escaláveis.

### Padrões principais

1.  **Estrutura do fluxo de trabalho**: Definimos entradas e saídas claras para cada fluxo de trabalho usando a sintaxe `take:` e `emit:`, criando interfaces bem definidas entre os componentes, e envolvemos a lógica do fluxo de trabalho dentro do bloco `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Os canais de entrada são declarados aqui
            input_ch

        main:
            // A lógica do fluxo de trabalho vai aqui
            // É aqui que os processos são chamados e os canais são manipulados
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Os canais de saída são declarados aqui
            output_ch = result_ch
    }
    ```

2.  **Importações de fluxo de trabalho:** Construímos dois módulos de fluxo de trabalho independentes e os importamos em um pipeline principal com instruções `include`.

    - Incluir um único fluxo de trabalho

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Incluir múltiplos fluxos de trabalho

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Incluir com alias para evitar conflitos de nomes

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry points**: O Nextflow requer um entry workflow sem nome para saber por onde começar a execução. Esse entry workflow chama seus fluxos de trabalho nomeados.

    - Fluxo de trabalho sem nome (entry point)

    ```groovy
    workflow {
        // Este é o entry point quando o script é executado
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Fluxo de trabalho nomeado (chamado a partir do entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Deve ser chamado a partir do entry workflow
    }
    ```

4.  **Gerenciamento do fluxo de dados:** Aprendemos como acessar as saídas do fluxo de trabalho usando a notação de namespace (`WORKFLOW_NAME.out.channel_name`) e passá-las para outros fluxos de trabalho.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Recursos adicionais

- [Documentação de Workflow do Nextflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Referência de Operadores de Canal](https://www.nextflow.io/docs/latest/operator.html)
- [Documentação de Estratégia de Erros](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## O que vem a seguir?

Retorne ao [menu de Side Quests](../index.md) ou clique no botão no canto inferior direito da página para avançar para o próximo tópico da lista.
