# Parte 2: Reescrever Hello para nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta segunda parte do curso de treinamento Hello nf-core, mostramos como criar uma versão compatível com nf-core do pipeline produzido no curso para iniciantes [Hello Nextflow](../hello_nextflow/index.md).

Você deve ter notado na primeira seção do treinamento que os pipelines nf-core seguem uma estrutura bastante elaborada com muitos arquivos acessórios.
Criar tudo isso do zero seria muito tedioso, então a comunidade nf-core desenvolveu ferramentas para fazer isso a partir de um template, para inicializar o processo.

Vamos mostrar como usar essas ferramentas para criar uma estrutura de pipeline e então adaptar o código de pipeline 'regular' existente para a estrutura nf-core.

Se você não está familiarizado com o pipeline Hello ou precisa relembrar, consulte [esta página de informações](../info/hello_pipeline.md).

---

## 1. Criar um novo projeto de pipeline

Primeiro, criamos a estrutura para o novo pipeline.

!!! note "Nota"

    Certifique-se de estar no diretório `hello-nf-core` no seu terminal.

### 1.1. Executar a ferramenta de criação de pipeline baseada em template

Vamos começar criando um novo pipeline com o comando `nf-core pipelines create`.
Isso criará uma nova estrutura de pipeline usando o template base nf-core, customizado com um nome de pipeline, descrição e autor.

```bash
nf-core pipelines create
```

Executar este comando abrirá uma Interface de Usuário de Texto (TUI) para criação de pipeline:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Esta TUI pedirá que você forneça informações básicas sobre seu pipeline e oferecerá uma escolha de recursos para incluir ou excluir na estrutura do seu pipeline.

- Na tela de boas-vindas, clique em **Let's go!**.
- Na tela `Choose pipeline type`, clique em **Custom**.
- Insira os detalhes do seu pipeline como a seguir (substituindo `< SEU NOME >` pelo seu próprio nome), então clique em **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < SEU NOME >
```

- Na tela Template features, defina `Toggle all features` como **off**, então **habilite** seletivamente os seguintes. Verifique suas seleções e clique em **Continue**.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- Na tela `Final details`, clique em **Finish**. Aguarde o pipeline ser criado, então clique em **Continue**.
- Na tela Create GitHub repository, clique em **Finish without creating a repo**. Isso exibirá instruções para criar um repositório GitHub posteriormente. Ignore-as e clique em **Close**.

Quando a TUI fechar, você deverá ver a seguinte saída no console.

??? success "Saída do comando"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

Não há confirmação explícita na saída do console de que a criação do pipeline funcionou, mas você deverá ver um novo diretório chamado `core-hello`.

Visualize o conteúdo do novo diretório para ver quanto trabalho você economizou usando o template.

```bash
tree core-hello
```

??? abstract "Conteúdo do diretório"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

São muitos arquivos!

Esperamos que você reconheça muitos deles como os mesmos que encontramos quando exploramos a estrutura do pipeline `nf-core/demo`.
Mas não se preocupe se ainda estiver se sentindo um pouco perdido; vamos percorrer as partes importantes juntos no decorrer deste treinamento.

!!! note "Nota"

    Uma diferença importante em comparação com o pipeline `nf-core/demo` que examinamos na primeira parte deste treinamento é que não há diretório `modules`.
    Isso ocorre porque não optamos por incluir nenhum dos módulos padrão nf-core.

### 1.2. Testar que a estrutura é funcional

Acredite ou não, mesmo que você ainda não tenha adicionado nenhum módulo para fazer trabalho real, a estrutura do pipeline pode realmente ser executada usando o perfil de teste, da mesma forma que executamos o pipeline `nf-core/demo`.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Isso mostra que toda a configuração básica está no lugar.
Então onde estão as saídas? Existem algumas?

Na verdade, um novo diretório de resultados chamado `core-hello-results` foi criado contendo os relatórios de execução padrão:

```bash
tree core-hello-results
```

??? abstract "Conteúdo do diretório"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Você pode dar uma olhada nos relatórios para ver o que foi executado, e a resposta é: nada!

![relatório de timeline de execução vazio](./img/execution_timeline_empty.png)

Vamos dar uma olhada no que realmente está no código.

### 1.3. Examinar o fluxo de trabalho placeholder

Se você olhar dentro do arquivo `main.nf`, verá que ele importa um fluxo de trabalho chamado `HELLO` de `workflows/hello`.

Isso é equivalente ao fluxo de trabalho `workflows/demo.nf` que encontramos na Parte 1, e serve como um fluxo de trabalho placeholder para nosso fluxo de trabalho de interesse, com alguma funcionalidade nf-core já implementada.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Agrupar e salvar versões de software
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Comparado a um fluxo de trabalho Nextflow básico como o desenvolvido em [Hello Nextflow](../hello_nextflow/index.md), você notará algumas coisas novas aqui (linhas destacadas acima):

- O bloco workflow tem um nome
- As entradas do fluxo de trabalho são declaradas usando a palavra-chave `take:` e a construção do canal é movida para o fluxo de trabalho pai
- O conteúdo do fluxo de trabalho é colocado dentro de um bloco `main:`
- As saídas são declaradas usando a palavra-chave `emit:`

Esses são recursos opcionais do Nextflow que tornam o fluxo de trabalho **componível**, o que significa que ele pode ser chamado de dentro de outro fluxo de trabalho.

!!! note "Fluxos de trabalho componíveis em profundidade"

    A [Side Quest Workflows of Workflows](../side_quests/workflows_of_workflows.md) explora a composição de fluxos de trabalho em muito mais profundidade, incluindo como compor múltiplos fluxos de trabalho juntos e gerenciar fluxos de dados complexos entre eles. Estamos introduzindo a componibilidade aqui porque é um requisito fundamental da arquitetura de template nf-core, que usa fluxos de trabalho aninhados para organizar a inicialização do pipeline, o fluxo de trabalho de análise principal e tarefas de conclusão em componentes separados e reutilizáveis.

Vamos precisar conectar a lógica relevante do nosso fluxo de trabalho de interesse nessa estrutura.
O primeiro passo para isso é tornar nosso fluxo de trabalho original componível.

### Conclusão

Agora você sabe como criar uma estrutura de pipeline usando ferramentas nf-core.

### O que vem a seguir?

Aprenda como tornar um fluxo de trabalho simples componível como prelúdio para torná-lo compatível com nf-core.

---

## 2. Tornar o fluxo de trabalho Hello Nextflow original componível

Agora é hora de trabalhar na integração do nosso fluxo de trabalho na estrutura nf-core.
Como lembrete, estamos trabalhando com o fluxo de trabalho apresentado no nosso curso de treinamento [Hello Nextflow](../hello_nextflow/index.md).

!!! tip "Dica"

    Se você não está familiarizado com esse pipeline ou precisa relembrar, consulte [O pipeline Hello](../info/hello_pipeline.md).

Fornecemos uma cópia limpa e totalmente funcional do fluxo de trabalho Hello Nextflow completo no diretório `original-hello` junto com seus módulos e o arquivo CSV padrão que ele espera usar como entrada.

```bash
tree original-hello/
```

??? abstract "Conteúdo do diretório"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

Sinta-se à vontade para executá-lo para se convencer de que funciona:

```bash
nextflow run original-hello/hello.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Vamos abrir o arquivo de fluxo de trabalho `hello.nf` para inspecionar o código, que é mostrado completo abaixo (sem contar os processos, que estão em módulos):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Inclui módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // cria um canal para entradas de um arquivo CSV
  greeting_ch = channel.fromPath(params.greeting)
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
}
```

Como você pode ver, este fluxo de trabalho foi escrito como um fluxo de trabalho simples sem nome que pode ser executado sozinho.
Para torná-lo executável de dentro de um fluxo de trabalho pai como o template nf-core requer, precisamos torná-lo **componível**.

Vamos percorrer as mudanças necessárias uma a uma.

### 2.1. Nomear o fluxo de trabalho

Primeiro, vamos dar um nome ao fluxo de trabalho para que possamos nos referir a ele de um fluxo de trabalho pai.

=== "Depois"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

As mesmas convenções se aplicam aos nomes de fluxos de trabalho e aos nomes de módulos.

### 2.2. Substituir construção de canal por `take`

Agora, substitua a construção do canal por uma simples declaração `take` declarando as entradas esperadas.

=== "Depois"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // canal de saudações
        greeting_ch
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Isso deixa os detalhes de como as entradas são fornecidas para o fluxo de trabalho pai.

Já que estamos nisso, também podemos comentar a linha `params.greeting = 'greetings.csv'`

=== "Depois"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Nota"

    Se você tiver a extensão do servidor de linguagem Nextflow instalada, o verificador de sintaxe iluminará seu código com rabiscos vermelhos.
    Isso ocorre porque se você colocar uma declaração `take:`, também precisa ter um `main:`.

    Adicionaremos isso no próximo passo.

### 2.3. Prefaciar operações do fluxo de trabalho com declaração `main`

Em seguida, adicione uma declaração `main` antes das demais operações chamadas no corpo do fluxo de trabalho.

=== "Depois"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // emite uma saudação
        sayHello(greeting_ch)

        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // gera arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emite uma saudação
        sayHello(greeting_ch)

        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // gera arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Isso basicamente diz 'isto é o que este fluxo de trabalho _faz_'.

### 2.4. Adicionar declaração `emit`

Finalmente, adicione uma declaração `emit` declarando quais são as saídas finais do fluxo de trabalho.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Esta é uma adição completamente nova ao código em comparação com o fluxo de trabalho original.

### 2.5. Recapitulação das mudanças concluídas

Se você fez todas as mudanças conforme descrito, seu fluxo de trabalho deve agora se parecer com isto:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Inclui módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // canal de saudações
    greeting_ch

    main:

    // emite uma saudação
    sayHello(greeting_ch)

    // converte a saudação para maiúsculas
    convertToUpper(sayHello.out)

    // coleta todas as saudações em um arquivo
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // gera arte ASCII das saudações com cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Isso descreve tudo que o Nextflow precisa EXCETO o que alimentar no canal de entrada.
Isso será definido no fluxo de trabalho pai, também chamado de fluxo de trabalho **entrypoint**.

### 2.6. Fazer um fluxo de trabalho entrypoint de teste

Antes de integrar nosso fluxo de trabalho componível na estrutura nf-core complexa, vamos verificar se funciona corretamente.
Podemos fazer um fluxo de trabalho entrypoint de teste simples para testar o fluxo de trabalho componível isoladamente.

Crie um arquivo em branco chamado `main.nf` no mesmo diretório `original-hello`.

```bash
touch original-hello/main.nf
```

Copie o seguinte código para o arquivo `main.nf`.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// importar o código do fluxo de trabalho do arquivo hello.nf
include { HELLO } from './hello.nf'

// declarar parâmetro de entrada
params.greeting = 'greetings.csv'

workflow {
  // criar um canal para entradas de um arquivo CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // chamar o fluxo de trabalho importado no canal de saudações
  HELLO(greeting_ch)

  // visualizar as saídas emitidas pelo fluxo de trabalho
  HELLO.out.view { output -> "Output: $output" }
}
```

Há duas observações importantes a fazer aqui:

- A sintaxe para chamar o fluxo de trabalho importado é essencialmente a mesma que a sintaxe para chamar módulos.
- Tudo que está relacionado a trazer as entradas para o fluxo de trabalho (parâmetro de entrada e construção de canal) agora é declarado neste fluxo de trabalho pai.

!!! note "Nota"

    Nomear o arquivo de fluxo de trabalho entrypoint `main.nf` é uma convenção, não um requisito.

    Se você seguir esta convenção, pode omitir a especificação do nome do arquivo de fluxo de trabalho no seu comando `nextflow run`.
    O Nextflow procurará automaticamente por um arquivo chamado `main.nf` no diretório de execução.

    No entanto, você pode nomear o arquivo de fluxo de trabalho entrypoint de outra forma se preferir.
    Nesse caso, certifique-se de especificar o nome do arquivo de fluxo de trabalho no seu comando `nextflow run`.

### 2.7. Testar que o fluxo de trabalho executa

Finalmente temos todas as peças que precisamos para verificar que o fluxo de trabalho componível funciona.
Vamos executá-lo!

```bash
nextflow run ./original-hello
```

Aqui você vê a vantagem de usar a convenção de nomenclatura `main.nf`.
Se tivéssemos nomeado o fluxo de trabalho entrypoint `algo_diferente.nf`, teríamos que fazer `nextflow run original-hello/algo_diferente.nf`.

Se você fez todas as mudanças corretamente, isso deve executar até a conclusão.

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Isso significa que atualizamos com sucesso nosso fluxo de trabalho HELLO para ser componível.

### Conclusão

Você sabe como tornar um fluxo de trabalho componível dando-lhe um nome e adicionando declarações `take`, `main` e `emit`, e como chamá-lo de um fluxo de trabalho entrypoint.

### O que vem a seguir?

Aprenda como enxertar um fluxo de trabalho componível básico na estrutura nf-core.

---

## 3. Ajustar a lógica do fluxo de trabalho atualizado no fluxo de trabalho placeholder

Agora que verificamos que nosso fluxo de trabalho componível funciona corretamente, vamos retornar à estrutura do pipeline nf-core que criamos na seção 1.
Queremos integrar o fluxo de trabalho componível que acabamos de desenvolver na estrutura do template nf-core, então o resultado final deve se parecer com isto.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Então como fazemos isso acontecer? Vamos dar uma olhada no conteúdo atual do fluxo de trabalho `HELLO` em `core-hello/workflows/hello.nf` (a estrutura nf-core).

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Agrupar e salvar versões de software
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

No geral, este código faz muito pouco além de algumas tarefas administrativas relacionadas à captura da versão de qualquer ferramenta de software que seja executada no pipeline.

Precisamos adicionar o código relevante da versão componível do fluxo de trabalho original que desenvolvemos na seção 2.

Vamos abordar isso nos seguintes estágios:

1. Copiar os módulos e configurar importações de módulos
2. Deixar a declaração `take` como está
3. Adicionar a lógica do fluxo de trabalho ao bloco `main`
4. Atualizar o bloco `emit`

!!! note "Nota"

    Vamos ignorar a captura de versão nesta primeira passagem e veremos como conectar isso em uma parte posterior deste treinamento.

### 3.1. Copiar os módulos e configurar importações de módulos

Os quatro processos do nosso fluxo de trabalho Hello Nextflow são armazenados como módulos em `original-hello/modules/`.
Precisamos copiar esses módulos para a estrutura do projeto nf-core (em `core-hello/modules/local/`) e adicionar declarações de importação ao arquivo de fluxo de trabalho nf-core.

Primeiro vamos copiar os arquivos de módulo de `original-hello/` para `core-hello/`:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Você deve agora ver o diretório de módulos listado em `core-hello/`.

```bash
tree core-hello/modules
```

??? abstract "Conteúdo do diretório"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Agora vamos configurar as declarações de importação de módulos.

Estas eram as declarações de importação no fluxo de trabalho `original-hello/hello.nf`:

```groovy title="original-hello/hello.nf" linenums="9"
// Inclui módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Abra o arquivo `core-hello/workflows/hello.nf` e transponha essas declarações de importação para ele conforme mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Mais duas observações interessantes aqui:

- Adaptamos a formatação das declarações de importação para seguir a convenção de estilo nf-core.
- Atualizamos os caminhos relativos para os módulos para refletir que eles agora estão armazenados em um nível diferente de aninhamento.

### 3.2. Deixar a declaração `take` como está

O projeto nf-core tem muita funcionalidade pré-construída em torno do conceito de samplesheet, que é tipicamente um arquivo CSV contendo dados em colunas.
Como isso é essencialmente o que nosso arquivo `greetings.csv` é, vamos manter a declaração `take` atual como está, e simplesmente atualizar o nome do canal de entrada no próximo passo.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

O tratamento de entrada será feito antes deste fluxo de trabalho (não neste arquivo de código).

### 3.3. Adicionar a lógica do fluxo de trabalho ao bloco `main`

Agora que nossos módulos estão disponíveis para o fluxo de trabalho, podemos conectar a lógica do fluxo de trabalho no bloco `main`.

Como lembrete, este é o código relevante no fluxo de trabalho original, que não mudou muito quando o tornamos componível (apenas adicionamos a linha `main:`):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // emite uma saudação
    sayHello(greeting_ch)

    // converte a saudação para maiúsculas
    convertToUpper(sayHello.out)

    // coleta todas as saudações em um arquivo
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // gera arte ASCII das saudações com cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

Precisamos copiar o código que vem depois de `main:` para a nova versão do fluxo de trabalho.

Já existe algum código lá que tem a ver com capturar as versões das ferramentas que são executadas pelo fluxo de trabalho. Vamos deixar isso em paz por enquanto (lidaremos com as versões de ferramentas mais tarde).
Manteremos a inicialização `ch_versions = channel.empty()` no topo, depois inseriremos nossa lógica de fluxo de trabalho, mantendo o código de coleta de versões no final.
Esta ordenação faz sentido porque em um pipeline real, os processos emitiriam informações de versão que seriam adicionadas ao canal `ch_versions` conforme o fluxo de trabalho executa.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // emitir uma saudação
        sayHello(greeting_ch)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // coletar todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Agrupar e salvar versões de software
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

        //
        // Agrupar e salvar versões de software
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

Você notará que também adicionamos uma linha em branco antes de `main:` para tornar o código mais legível.

Isso parece ótimo, mas ainda precisamos atualizar o nome do canal que estamos passando para o processo `sayHello()` de `greeting_ch` para `ch_samplesheet` conforme mostrado abaixo, para corresponder ao que está escrito sob a palavra-chave `take:`.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir uma saudação (atualizado para usar a convenção nf-core para samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emite uma saudação
        sayHello(greeting_ch)
    ```

Agora a lógica do fluxo de trabalho está corretamente conectada.

### 3.4. Atualizar o bloco `emit`

Finalmente, precisamos atualizar o bloco `emit` para incluir a declaração das saídas finais do fluxo de trabalho.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Isso conclui as modificações que precisamos fazer no próprio fluxo de trabalho HELLO.
Neste ponto, alcançamos a estrutura geral de código que nos propusemos a implementar.

### Conclusão

Você sabe como ajustar as peças principais de um fluxo de trabalho componível em um fluxo de trabalho placeholder nf-core.

### O que vem a seguir?

Aprenda como adaptar o tratamento de entradas na estrutura do pipeline nf-core.

---

## 4. Adaptar o tratamento de entrada

Agora que integramos com sucesso nossa lógica de fluxo de trabalho na estrutura nf-core, precisamos abordar mais uma peça crítica: garantir que nossos dados de entrada sejam processados corretamente.
O template nf-core vem com tratamento de entrada sofisticado projetado para conjuntos de dados de genômica complexos, então precisamos adaptá-lo para funcionar com nosso arquivo `greetings.csv` mais simples.

### 4.1. Identificar onde as entradas são tratadas

O primeiro passo é descobrir onde o tratamento de entrada é feito.

Você pode se lembrar que quando reescrevemos o fluxo de trabalho Hello Nextflow para ser componível, movemos a declaração de parâmetro de entrada um nível acima, no fluxo de trabalho entrypoint `main.nf`.
Então vamos dar uma olhada no fluxo de trabalho entrypoint `main.nf` de nível superior que foi criado como parte da estrutura do pipeline:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

O projeto nf-core faz uso intenso de subfluxos de trabalho aninhados, então esta parte pode ser um pouco confusa na primeira abordagem.

O que importa aqui é que existem dois fluxos de trabalho definidos:

- `CORE_HELLO` é um wrapper fino para executar o fluxo de trabalho HELLO que acabamos de terminar de adaptar em `core-hello/workflows/hello.nf`.
- Um fluxo de trabalho sem nome que chama `CORE_HELLO` bem como dois outros subfluxos de trabalho, `PIPELINE_INITIALISATION` e `PIPELINE_COMPLETION`.

Aqui está um diagrama de como eles se relacionam:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Importante, não podemos encontrar nenhum código construindo um canal de entrada neste nível, apenas referências a uma samplesheet fornecida via parâmetro `--input`.

Um pouco de investigação revela que o tratamento de entrada é feito pelo subfluxo de trabalho `PIPELINE_INITIALISATION`, apropriadamente, que é importado de `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`.

Se abrirmos esse arquivo e rolarmos para baixo, chegamos a este pedaço de código:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Cria canal a partir do arquivo de entrada fornecido através de params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

Esta é a fábrica de canais que analisa a samplesheet e a passa adiante em uma forma que está pronta para ser consumida pelo fluxo de trabalho HELLO.

!!! note "Nota"

    A sintaxe acima é um pouco diferente do que usamos anteriormente, mas basicamente isto:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    é equivalente a isto:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Este código envolve algumas etapas de análise e validação que são altamente específicas para a samplesheet de exemplo incluída no template de pipeline nf-core, que no momento da escrita é muito específica de domínio e não adequada para nosso projeto de pipeline simples.

### 4.2. Substituir o código de canal de entrada do template

A boa notícia é que as necessidades do nosso pipeline são muito mais simples, então podemos substituir tudo isso pelo código de construção de canal que desenvolvemos no fluxo de trabalho Hello Nextflow original.

Como lembrete, isto é como a construção de canal se parecia (como visto no diretório de soluções):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // criar um canal para entradas de um arquivo CSV
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Então só precisamos conectar isso no fluxo de trabalho de inicialização, com pequenas mudanças: atualizamos o nome do canal de `greeting_ch` para `ch_samplesheet`, e o nome do parâmetro de `params.greeting` para `params.input` (veja a linha destacada).

=== "Depois"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Cria canal a partir do arquivo de entrada fornecido através de params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Cria canal a partir do arquivo de entrada fornecido através de params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Isso completa as mudanças que precisamos fazer para que o processamento de entrada funcione.

Na sua forma atual, isso não nos permitirá aproveitar as capacidades integradas do nf-core para validação de schema, mas podemos adicionar isso mais tarde.
Por enquanto, estamos focados em mantê-lo o mais simples possível para chegar a algo que possamos executar com sucesso em dados de teste.

### 4.3. Atualizar o perfil de teste

Falando em dados de teste e parâmetros, vamos atualizar o perfil de teste para este pipeline para usar a mini-samplesheet `greetings.csv` em vez da samplesheet de exemplo fornecida no template.

Em `core-hello/conf`, encontramos dois perfis de teste do template: `test.config` e `test_full.config`, que são destinados a testar uma pequena amostra de dados e uma de tamanho completo.
Dado o propósito do nosso pipeline, não há realmente sentido em configurar um perfil de teste de tamanho completo, então sinta-se à vontade para ignorar ou excluir `test_full.config`.
Vamos focar em configurar `test.config` para executar em nosso arquivo `greetings.csv` com alguns parâmetros padrão.

#### 4.3.1. Copiar o arquivo `greetings.csv`

Primeiro precisamos copiar o arquivo `greetings.csv` para um local apropriado em nosso projeto de pipeline.
Tipicamente pequenos arquivos de teste são armazenados no diretório `assets`, então vamos copiar o arquivo do nosso diretório de trabalho.

```bash
cp greetings.csv core-hello/assets/.
```

Agora o arquivo `greetings.csv` está pronto para ser usado como entrada de teste.

#### 4.3.2. Atualizar o arquivo `test.config`

Agora podemos atualizar o arquivo `test.config` da seguinte forma:

=== "Depois"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Dados de entrada
        input  = "${projectDir}/assets/greetings.csv"

        // Outros parâmetros
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Antes"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Dados de entrada
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Pontos-chave:

- **Usando `${projectDir}`**: Esta é uma variável implícita do Nextflow que aponta para o diretório onde o script de fluxo de trabalho principal está localizado (a raiz do pipeline). Usá-la garante que o caminho funcione independentemente de onde o pipeline seja executado.
- **Caminhos absolutos**: Ao usar `${projectDir}`, criamos um caminho absoluto, o que é importante para dados de teste que são enviados com o pipeline.
- **Localização de dados de teste**: pipelines nf-core tipicamente armazenam dados de teste no diretório `assets/` dentro do repositório do pipeline para pequenos arquivos de teste, ou referenciam conjuntos de dados de teste externos para arquivos maiores.

E já que estamos nisso, vamos apertar os limites de recursos padrão para garantir que isso seja executado em máquinas muito básicas (como as VMs mínimas no GitHub Codespaces):

=== "Depois"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Antes"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Isso completa as modificações de código que precisamos fazer.

### 4.4. Executar o pipeline com o perfil de teste

Isso foi muito, mas finalmente podemos tentar executar o pipeline!
Note que temos que adicionar `--validate_params false` à linha de comando porque ainda não configuramos a validação (isso virá mais tarde).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Se você fez todas as modificações corretamente, deve executar até a conclusão.

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Como você pode ver, isso produziu o resumo típico nf-core no início graças ao subfluxo de trabalho de inicialização, e as linhas para cada módulo agora mostram os nomes completos PIPELINE:WORKFLOW:módulo.

### 4.5. Encontrar as saídas do pipeline

A questão agora é: onde estão as saídas do pipeline?
E a resposta é bastante interessante: agora há dois lugares diferentes para procurar os resultados.

Como você pode se lembrar anteriormente, nossa primeira execução do fluxo de trabalho recém-criado produziu um diretório chamado `core-hello-results/` que continha vários relatórios de execução e metadados.

```bash
tree core-hello-results
```

??? abstract "Conteúdo do diretório"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

Você vê que obtivemos outro conjunto de relatórios de execução além dos que obtivemos da primeira execução, quando o fluxo de trabalho ainda era apenas um placeholder.
Desta vez você vê todas as tarefas que foram executadas como esperado.

![relatório de timeline de execução para o pipeline Hello](./img/execution_timeline_hello.png)

!!! note "Nota"

    Mais uma vez as tarefas não foram executadas em paralelo porque estamos executando em uma máquina minimalista no GitHub Codespaces.
    Para ver essas executadas em paralelo, tente aumentar a alocação de CPU do seu codespace e os limites de recursos na configuração de teste.

Isso é ótimo, mas nossos resultados reais do pipeline não estão lá!

Aqui está o que aconteceu: não mudamos nada nos próprios módulos, então as saídas tratadas pelas diretivas `publishDir` em nível de módulo ainda vão para um diretório `results` conforme especificado no pipeline original.

```bash
tree results
```

??? abstract "Conteúdo do diretório"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, ali estão eles, misturados com as saídas de execuções anteriores do pipeline Hello original.

Se quisermos que eles sejam organizados de forma organizada como as saídas do pipeline demo foram, precisaremos mudar como configuramos as saídas para serem publicadas.
Mostraremos como fazer isso mais tarde neste curso de treinamento.

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

E aí está! Pode parecer muito trabalho para obter o mesmo resultado que o pipeline original, mas você obtém todos aqueles relatórios lindos gerados automaticamente, e agora tem uma base sólida para aproveitar recursos adicionais do nf-core, incluindo validação de entrada e alguns recursos legais de tratamento de metadados que cobriremos em uma seção posterior.

---

### Conclusão

Você sabe como converter um pipeline Nextflow regular em um pipeline no estilo nf-core usando o template nf-core.
Como parte disso, você aprendeu como tornar um fluxo de trabalho componível, e como identificar os elementos do template nf-core que mais comumente precisam ser adaptados ao desenvolver um pipeline customizado no estilo nf-core.

### O que vem a seguir?

Faça uma pausa, isso foi trabalho duro! Quando estiver pronto, prossiga para [Parte 3: Usar um módulo nf-core](./03_use_module.md) para aprender como aproveitar módulos mantidos pela comunidade do repositório nf-core/modules.
