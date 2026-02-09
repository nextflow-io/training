# Workflows de Workflows

Quando você está desenvolvendo um pipeline, frequentemente se vê criando sequências similares de processos para diferentes tipos de dados ou etapas de análise. Você pode acabar copiando e colando essas sequências de processos, levando a código duplicado que é difícil de manter; ou pode criar um fluxo de trabalho massivo que é difícil de entender e modificar.

Uma das características mais poderosas do Nextflow é sua capacidade de compor pipelines complexos a partir de módulos de fluxo de trabalho menores e reutilizáveis. Essa abordagem modular torna os pipelines mais fáceis de desenvolver, testar e manter.

### Objetivos de aprendizado

Nesta missão secundária, vamos explorar como desenvolver módulos de fluxo de trabalho que podem ser testados e usados separadamente, compor esses módulos em um pipeline maior e gerenciar o fluxo de dados entre módulos.

Ao final desta missão secundária, você será capaz de:

- Dividir pipelines complexos em unidades lógicas e reutilizáveis
- Testar cada módulo de fluxo de trabalho independentemente
- Misturar e combinar fluxos de trabalho para criar novos pipelines
- Compartilhar módulos de fluxo de trabalho comuns entre diferentes pipelines
- Tornar seu código mais sustentável e fácil de entender

Essas habilidades ajudarão você a construir pipelines complexos mantendo uma estrutura de código limpa e sustentável.

### Pré-requisitos

Antes de assumir esta missão secundária você deve:

- Ter completado o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores, módulos)

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mova-se para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/workflows_of_workflows
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um diretório `modules` contendo várias definições de processos que se baseiam no que você aprendeu em 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### Revise a tarefa

Seu desafio é montar esses módulos em dois fluxos de trabalho separados que então comporemos em um fluxo de trabalho principal:

- Um `GREETING_WORKFLOW` que valida nomes, cria saudações e adiciona timestamps
- Um `TRANSFORM_WORKFLOW` que converte texto para maiúsculas e o reverte

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Eu configurei meu diretório de trabalho apropriadamente
- [ ] Eu entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Crie o Greeting Workflow

Vamos começar criando um fluxo de trabalho que valida nomes e gera saudações com timestamp.

### 1.1. Crie a estrutura do fluxo de trabalho

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Adicione o código do primeiro (sub)fluxo de trabalho

Adicione este código a `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Este é um fluxo de trabalho completo, com uma estrutura similar às que você viu no tutorial 'Hello Nextflow', que podemos testar independentemente. Vamos tentar isso agora:

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

Isso funciona como esperado, mas para torná-lo componível há algumas coisas que precisamos mudar.

### 1.3. Torne o fluxo de trabalho componível

Fluxos de trabalho componíveis têm algumas diferenças dos que você viu no tutorial 'Hello Nextflow':

- O bloco workflow precisa ser nomeado
- Entradas são declaradas usando a palavra-chave `take:`
- O conteúdo do fluxo de trabalho é colocado dentro do bloco `main:`
- Saídas são declaradas usando a palavra-chave `emit:`

Vamos atualizar o fluxo de trabalho de saudação para corresponder a esta estrutura. Mude o código para o seguinte:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Input channel with names

    main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
}
```

Você pode ver que o fluxo de trabalho agora está nomeado e tem um bloco `take:` e `emit:`, e essas são as conexões que usaremos para compor um fluxo de trabalho de nível superior.
O conteúdo do fluxo de trabalho também está colocado dentro do bloco `main:`. Note também que removemos a declaração do canal de entrada `names_ch`, pois agora ele é passado como um argumento para o fluxo de trabalho.

Vamos testar o fluxo de trabalho novamente para ver se funciona como esperado:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Isso informa sobre outro novo conceito, um 'fluxo de trabalho de entrada'. O fluxo de trabalho de entrada é o fluxo de trabalho que é chamado quando você executa um script Nextflow. Por padrão, o Nextflow usará um fluxo de trabalho sem nome como o fluxo de trabalho de entrada, quando presente, e é isso que você tem feito até agora, com blocos workflow começando assim:

```groovy title="hello.nf" linenums="1"
workflow {
```

Mas nosso fluxo de trabalho de saudação não tem um fluxo de trabalho sem nome, em vez disso temos um fluxo de trabalho nomeado:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

É por isso que o Nextflow lançou um erro e não fez o que queríamos.

Não adicionamos a sintaxe `take:`/`emit:` para que pudéssemos chamar o fluxo de trabalho diretamente - fizemos isso para que pudéssemos compô-lo com outros fluxos de trabalho. A solução é criar um script principal com um fluxo de trabalho de entrada sem nome que importa e chama nosso fluxo de trabalho nomeado.

### 1.4. Crie e teste o fluxo de trabalho principal

Agora vamos criar um fluxo de trabalho principal que importa e usa o fluxo de trabalho `greeting`.

Crie `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Note que nossa entrada de fluxo de trabalho neste arquivo não tem nome, e isso é porque vamos usá-la como um fluxo de trabalho de entrada.

Execute isso e veja a saída:

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
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Funciona! Nós envolvemos o fluxo de trabalho de saudação nomeado em um fluxo de trabalho principal com um bloco `workflow` de entrada sem nome. O fluxo de trabalho principal está usando o fluxo de trabalho `GREETING_WORKFLOW` quase (não exatamente) como um processo, e passando o canal `names` como um argumento.

### Conclusão

Nesta seção, você aprendeu vários conceitos importantes:

- **Workflows Nomeados**: Criar um fluxo de trabalho nomeado (`GREETING_WORKFLOW`) que pode ser importado e reutilizado
- **Interfaces de Workflow**: Definir entradas claras com `take:` e saídas com `emit:` para criar um fluxo de trabalho componível
- **Pontos de Entrada**: Entender que o Nextflow precisa de um fluxo de trabalho de entrada sem nome para executar um script
- **Composição de Workflow**: Importar e usar um fluxo de trabalho nomeado dentro de outro fluxo de trabalho
- **Namespaces de Workflow**: Acessar saídas de fluxo de trabalho usando o namespace `.out` (`GREETING_WORKFLOW.out.greetings`)

Você agora tem um fluxo de trabalho de saudação funcionando que:

- Recebe um canal de nomes como entrada
- Valida cada nome
- Cria uma saudação para cada nome válido
- Adiciona timestamps às saudações
- Expõe tanto as saudações originais quanto as com timestamp como saídas

Essa abordagem modular permite que você teste o fluxo de trabalho de saudação independentemente ou o use como um componente em pipelines maiores.

---

## 2. Adicione o Transform Workflow

Agora vamos criar um fluxo de trabalho que aplica transformações de texto às saudações.

### 2.1. Crie o arquivo do fluxo de trabalho

```bash
touch workflows/transform.nf
```

### 2.2. Adicione o código do fluxo de trabalho

Adicione este código a `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Input channel with messages

    main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase greetings
        reversed = reversed_ch  // Reversed uppercase greetings
}
```

Não vamos repetir a explicação da sintaxe componível aqui, mas note que o fluxo de trabalho nomeado é novamente declarado com um bloco `take:` e `emit:`, e o conteúdo do fluxo de trabalho é colocado dentro do bloco `main:`.

### 2.3. Atualize o fluxo de trabalho principal

Atualize `main.nf` para usar ambos os fluxos de trabalho:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Execute o pipeline completo:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Se você der uma olhada em um desses arquivos revertidos, verá que é a versão em maiúsculas da saudação revertida:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Conclusão

Você agora deve ter um pipeline completo que:

- Processa nomes através do fluxo de trabalho de saudação
- Alimenta as saudações com timestamp no fluxo de trabalho de transformação
- Produz versões tanto em maiúsculas quanto revertidas das saudações

---

## Resumo

Nesta missão secundária, exploramos o poderoso conceito de composição de fluxo de trabalho no Nextflow, que nos permite construir pipelines complexos a partir de componentes menores e reutilizáveis.

Essa abordagem modular oferece várias vantagens sobre pipelines monolíticos:

- Cada fluxo de trabalho pode ser desenvolvido, testado e depurado independentemente
- Fluxos de trabalho podem ser reutilizados em diferentes pipelines
- A estrutura geral do pipeline se torna mais legível e sustentável
- Mudanças em um fluxo de trabalho não necessariamente afetam outros se as interfaces permanecerem consistentes
- Pontos de entrada podem ser configurados para executar diferentes partes do seu pipeline conforme necessário

_É importante notar, no entanto, que embora chamar fluxos de trabalho seja um pouco como chamar processos, não é realmente a mesma coisa. Você não pode, por exemplo, executar um fluxo de trabalho N vezes chamando-o com um canal de tamanho N - você precisaria passar um canal de tamanho N para o fluxo de trabalho e iterar internamente._

Aplicar essas técnicas em seu próprio trabalho permitirá que você construa pipelines Nextflow mais sofisticados que podem lidar com tarefas complexas de bioinformática enquanto permanecem sustentáveis e escaláveis.

### Padrões principais

1.  **Estrutura de workflow**: Definimos entradas e saídas claras para cada fluxo de trabalho usando a sintaxe `take:` e `emit:`, criando interfaces bem definidas entre componentes, e envolvemos a lógica do fluxo de trabalho dentro do bloco `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **Importações de workflow:** Construímos dois módulos de fluxo de trabalho independentes e os importamos em um pipeline principal com declarações include.

    - Incluir um único fluxo de trabalho

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Incluir múltiplos fluxos de trabalho

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Incluir com alias para evitar conflitos de nome

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Pontos de entrada**: O Nextflow requer um fluxo de trabalho de entrada sem nome para saber onde começar a execução. Este fluxo de trabalho de entrada chama seus fluxos de trabalho nomeados.

    - Workflow sem nome (ponto de entrada)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow nomeado (chamado do fluxo de trabalho de entrada)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **Gerenciando fluxo de dados:** Aprendemos como acessar saídas de fluxo de trabalho usando a notação de namespace (`WORKFLOW_NAME.out.channel_name`) e passá-las para outros fluxos de trabalho.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Recursos adicionais

- [Documentação de Workflow do Nextflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Referência de Operadores de Canal](https://www.nextflow.io/docs/latest/operator.html)
- [Documentação de Error Strategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## O que vem a seguir?

Retorne ao [menu de Missões Secundárias](./index.md) ou clique no botão no canto inferior direito da página para passar para o próximo tópico da lista.
