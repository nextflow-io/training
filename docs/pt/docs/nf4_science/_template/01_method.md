# Parte 1: Visão geral do método e testes manuais

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![Visão geral do pipeline](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Métodos

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Antes de mergulharmos na escrita de qualquer código de fluxo de trabalho, vamos testar os comandos manualmente em alguns dados de teste.

### Conjunto de dados

Fornecemos os seguintes dados e recursos relacionados:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Software

As principais ferramentas envolvidas são [{TOOL_A}]({TOOL_A_URL}) e [{TOOL_B}]({TOOL_B_URL}).

Essas ferramentas não estão instaladas no ambiente GitHub Codespaces, então vamos usá-las via contêineres (veja [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

    Certifique-se de estar no diretório `nf4-science/{DOMAIN_DIR}` para que a última parte do caminho mostrada quando você digitar `pwd` seja `{DOMAIN_DIR}`.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

Nesta seção testamos os comandos que compõem a abordagem de processamento de amostra única.
Estes são os comandos que vamos encapsular em um fluxo de trabalho Nextflow na Parte 2 deste curso.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Começamos testando os comandos em apenas uma amostra.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Baixar o contêiner

Execute o comando `docker pull` para baixar a imagem do contêiner:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Saída do comando"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Iniciar o contêiner interativamente

Inicie o contêiner e monte o diretório `data` para que as ferramentas possam acessar os arquivos de entrada:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

Seu prompt muda para indicar que você está dentro do contêiner.

#### 1.1.3. Executar o comando

```bash
{TOOL_A_COMMAND}
```

??? success "Saída do comando"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Baixar o contêiner

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Saída do comando"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Iniciar o contêiner interativamente

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Executar o comando

```bash
{TOOL_B_COMMAND}
```

??? success "Saída do comando"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal.
Isso conclui o teste de processamento de amostra única.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

Nesta seção testamos os comandos adicionais necessários para o processamento de múltiplas amostras.
Estes são os comandos que vamos encapsular em um fluxo de trabalho Nextflow na Parte 3 deste curso.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Conclusão

Você sabe como testar os comandos {TOOL_A} e {TOOL_B} em seus respectivos contêineres, incluindo como {MULTI_SAMPLE_SUMMARY}.

### O que vem a seguir?

Aprenda como encapsular esses mesmos comandos em fluxos de trabalho que usam contêineres para executar o trabalho.
