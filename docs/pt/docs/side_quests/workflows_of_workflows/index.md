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

- Ter concluído o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou um curso equivalente para iniciantes.
- Estar confortável com os conceitos e mecanismos básicos do Nextflow (processos, canais, operadores, módulos)

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se ainda não tiver feito isso, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

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

## 1. Criar o Fluxo de Trabalho de Saudação

Vamos começar criando um fluxo de trabalho que valida nomes e gera saudações com timestamp.

### 1.1. Criar a estrutura do fluxo de trabalho

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Adicionar o código do primeiro (sub)fluxo de trabalho

Adicione este código em `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Encadeia os processos: validar -> criar saudação -> adicionar timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Este é um fluxo de trabalho completo, com uma estrutura similar às que você viu no tutorial 'Hello Nextflow', que podemos testar de forma independente. Vamos fazer isso agora:

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

Isso funciona como esperado, mas para torná-lo composável há algumas coisas que precisamos mudar.

### 1.3. Tornar o fluxo de trabalho composável

Fluxos de trabalho composáveis têm algumas diferenças em relação aos que você viu no tutorial 'Hello Nextflow':

- O bloco workflow precisa ter um nome
- As entradas são declaradas usando a palavra-chave `take:`
- O conteúdo do fluxo de trabalho é colocado dentro do bloco `main:`
- As saídas são declaradas usando a palavra-chave `emit:`

Vamos atualizar o fluxo de trabalho de saudação para corresponder a essa estrutura. Altere o código para o seguinte:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Canal de entrada com nomes

    main:
        // Encadeia os processos: validar -> criar saudação -> adicionar timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Saudações originais
        timestamped = timestamped_ch  // Saudações com timestamp
}
```

Você pode ver que o fluxo de trabalho agora tem um nome e possui um bloco `take:` e `emit:`, e essas são as conexões que usaremos para compor um fluxo de trabalho de nível superior.
O conteúdo do fluxo de trabalho também é colocado dentro do bloco `main:`. Note também que removemos a declaração do canal de entrada `names_ch`, pois ele agora é passado como argumento para o fluxo de trabalho.

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

Isso indica um novo conceito: o 'entry workflow' (fluxo de trabalho de entrada). O entry workflow é o fluxo de trabalho que é chamado quando você executa um script Nextflow. Por padrão, o Nextflow usará um fluxo de trabalho sem nome como entry workflow, quando presente, e é isso que você tem feito até agora, com blocos workflow começando assim:

```groovy title="hello.nf" linenums="1"
workflow {
```

Mas nosso fluxo de trabalho de saudação não tem um fluxo de trabalho sem nome — em vez disso, temos um fluxo de trabalho nomeado:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

É por isso que o Nextflow gerou um erro e não fez o que queríamos.

Não adicionamos a sintaxe `take:`/`emit:` para poder chamar o fluxo de trabalho diretamente — fizemos isso para poder compô-lo com outros fluxos de trabalho. A solução é criar um script principal com um entry workflow sem nome que importa e chama nosso fluxo de trabalho nomeado.

### 1.4. Criar e testar o fluxo de trabalho principal

Agora vamos criar um fluxo de trabalho principal que importa e usa o fluxo de trabalho `greeting`.

Crie o `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Note que o entry workflow neste arquivo não tem nome, e isso ocorre porque vamos usá-lo como entry workflow.

Execute e veja a saída:

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

Funcionou! Envolvemos o fluxo de trabalho de saudação nomeado em um fluxo de trabalho principal com um bloco `workflow` de entrada sem nome. O fluxo de trabalho principal está usando o `GREETING_WORKFLOW` quase (mas não exatamente) como um processo, e passando o canal `names` como argumento.

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

## 2. Adicionar o Fluxo de Trabalho de Transformação

Agora vamos criar um fluxo de trabalho que aplica transformações de texto às saudações.

### 2.1. Criar o arquivo do fluxo de trabalho

```bash
touch workflows/transform.nf
```

### 2.2. Adicionar o código do fluxo de trabalho

Adicione este código em `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Canal de entrada com mensagens

    main:
        // Aplica as transformações em sequência
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Saudações em maiúsculas
        reversed = reversed_ch  // Saudações em maiúsculas invertidas
}
```

Não vamos repetir a explicação da sintaxe composável aqui, mas note que o fluxo de trabalho nomeado é novamente declarado com um bloco `take:` e `emit:`, e o conteúdo do fluxo de trabalho é colocado dentro do bloco `main:`.

### 2.3. Atualizar o fluxo de trabalho principal

Atualize o `main.nf` para usar ambos os fluxos de trabalho:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Executa o fluxo de trabalho de saudação
    GREETING_WORKFLOW(names)

    // Executa o fluxo de trabalho de transformação
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Exibe os resultados
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

Se você abrir um desses arquivos invertidos, verá que é a versão em maiúsculas da saudação invertida:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

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

Aplicar essas técnicas no seu próprio trabalho permitirá que você construa pipelines Nextflow mais sofisticados, capazes de lidar com tarefas complexas de bioinformática, mantendo-os fáceis de manter e escaláveis.

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

2.  **Importações de fluxo de trabalho:** Construímos dois módulos de fluxo de trabalho independentes e os importamos em um pipeline principal com instruções include.

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

Retorne ao [menu de Side Quests](../) ou clique no botão no canto inferior direito da página para avançar para o próximo tópico da lista.
