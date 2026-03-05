# Parte 2: Processamento de amostra única

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na Parte 1, você testou os comandos {TOOL_A} e {TOOL_B} manualmente em seus respectivos contêineres.
Agora vamos encapsular esses mesmos comandos em um fluxo de trabalho Nextflow.

## Tarefa

Nesta parte do curso, vamos desenvolver um fluxo de trabalho que faz o seguinte:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Isso replica as etapas da Parte 1, onde você executou esses comandos manualmente em seus contêineres.

Como ponto de partida, fornecemos um arquivo de fluxo de trabalho, `{DOMAIN_DIR}.nf`, que descreve as partes principais do fluxo de trabalho, bem como dois arquivos de módulo, {TOOL_A_MODULE}.nf e {TOOL_B_MODULE}.nf, que descrevem a estrutura dos módulos.
Esses arquivos não são funcionais; seu propósito é apenas servir como estruturas para você preencher com as partes interessantes do código.

## Plano de aula

Para tornar o processo de desenvolvimento mais educativo, dividimos isso em {N} etapas:

1. **Escrever um fluxo de trabalho de estágio único que executa {TOOL_A_ACTION}.**
   Isso cobre a criação de um módulo, sua importação e chamada em um fluxo de trabalho.
2. **Adicionar um segundo processo para executar {TOOL_B_ACTION}.**
   Isso introduz o encadeamento de saídas de processos para entradas e o manuseio de arquivos acessórios.
3. **Adaptar o fluxo de trabalho para executar em um lote de amostras.**
   Isso cobre a execução paralela e introduz tuplas para manter arquivos associados juntos.
4. **Fazer o fluxo de trabalho aceitar uma planilha de amostras contendo um lote de arquivos de entrada.**
   Isso demonstra um padrão comum para fornecer entradas em massa.

Cada etapa se concentra em um aspecto específico do desenvolvimento de fluxo de trabalho.

---

## 1. Escrever um fluxo de trabalho de estágio único que executa {TOOL_A_ACTION}

Esta primeira etapa se concentra no básico: carregar {PRIMARY_INPUT_TYPE} e {TOOL_A_OUTPUT_DESCRIPTION}.

Lembre-se do comando `{TOOL_A_COMMAND_NAME}` da [Parte 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

O comando recebe {INPUT_DESCRIPTION} e produz {OUTPUT_DESCRIPTION}.
O URI do contêiner era `{TOOL_A_CONTAINER_URI}`.

Vamos pegar essas informações e encapsulá-las no Nextflow em três estágios:

1. Configurar a entrada
2. Escrever o processo e chamá-lo no fluxo de trabalho
3. Configurar o manuseio da saída

### 1.1. Configurar a entrada

Precisamos declarar um parâmetro de entrada, criar um perfil de teste para fornecer um valor padrão conveniente e criar um canal de entrada.

#### 1.1.1. Adicionar uma declaração de parâmetro de entrada

No arquivo principal do fluxo de trabalho `{DOMAIN_DIR}.nf`, na seção `Pipeline parameters`, declare um parâmetro CLI chamado `{PRIMARY_PARAM_NAME}`.

=== "Depois"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Antes"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Isso configura o parâmetro CLI, mas não queremos digitar o caminho do arquivo toda vez que executarmos o fluxo de trabalho durante o desenvolvimento.
Existem várias opções para fornecer um valor padrão; aqui usamos um perfil de teste.

#### 1.1.2. Criar um perfil de teste com um valor padrão em `nextflow.config`

Um perfil de teste fornece valores padrão convenientes para experimentar um fluxo de trabalho sem especificar entradas na linha de comando.
Esta é uma convenção comum no ecossistema Nextflow (veja [Hello Config](../../hello_nextflow/06_hello_config.md) para mais detalhes).

Adicione um bloco `profiles` ao `nextflow.config` com um perfil `test` que define o parâmetro `{PRIMARY_PARAM_NAME}` para um dos {PRIMARY_INPUT_TYPE}s de teste.

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aqui, estamos usando `${projectDir}`, uma variável integrada do Nextflow que aponta para o diretório onde o script do fluxo de trabalho está localizado.
Isso facilita a referência a arquivos de dados e outros recursos sem codificar caminhos absolutos.

#### 1.1.3. Configurar o canal de entrada

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Escrever o módulo {TOOL_A_NAME}

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Importar e chamar o módulo no fluxo de trabalho

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Executar o fluxo de trabalho

Neste ponto, temos um fluxo de trabalho de uma etapa que deve estar totalmente funcional.

Podemos executá-lo com `-profile test` para usar o valor padrão configurado no perfil de teste e evitar ter que escrever o caminho na linha de comando.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusão

Você sabe como criar um módulo contendo um processo, importá-lo para um fluxo de trabalho, chamá-lo com um canal de entrada e publicar os resultados.

### O que vem a seguir?

Adicionar um segundo processo para encadear etapas de análise adicionais.

---

## 2. Adicionar um segundo processo para executar {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Lembre-se do comando `{TOOL_B_COMMAND_NAME}` da [Parte 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Escrever o módulo {TOOL_B_NAME}

{MODULE_INSTRUCTIONS}

### 2.2. Importar e chamar o módulo no fluxo de trabalho

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Executar o fluxo de trabalho

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Saída do comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusão

Você sabe como encadear saídas de processos para entradas e manusear arquivos acessórios no fluxo de trabalho.

### O que vem a seguir?

Escalar o fluxo de trabalho para processar múltiplas amostras em paralelo.

---

## 3. Adaptar o fluxo de trabalho para executar em um lote de amostras

Até agora, o fluxo de trabalho processa uma única amostra.
Para lidar com múltiplas amostras, precisamos modificar como as entradas são fornecidas e aproveitar o paradigma de fluxo de dados do Nextflow para paralelizar a execução.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Executar o fluxo de trabalho em múltiplas amostras

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Saída do comando"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Conclusão

Você sabe como aproveitar o paradigma de fluxo de dados do Nextflow para paralelizar o processamento por amostra em múltiplas amostras de entrada.

### O que vem a seguir?

Na [Parte 3](03_multi_sample.md), você adicionará agregação de múltiplas amostras para combinar resultados de todas as amostras.
