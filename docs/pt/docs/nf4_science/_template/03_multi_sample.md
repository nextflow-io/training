# Parte 3: Agregação de múltiplas amostras

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na Parte 2, você construiu um pipeline de processamento por amostra que tratava cada amostra de forma independente.
Agora vamos estendê-lo para implementar {AGGREGATION_METHOD} de múltiplas amostras, conforme abordado na [Parte 1](01_method.md).

## Tarefa

Nesta parte do curso, vamos estender o fluxo de trabalho para fazer o seguinte:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Esta parte se baseia diretamente no fluxo de trabalho produzido pela Parte 2.

??? info "Como começar a partir desta seção"

    Esta seção do curso pressupõe que você completou a [Parte 2: Processamento de amostra única](./02_single_sample.md) e tem um pipeline `{DOMAIN_DIR}.nf` funcionando.

    Se você não completou a Parte 2 ou quer começar do zero para esta parte, você pode usar a solução da Parte 2 como ponto de partida.
    Execute estes comandos de dentro do diretório `nf4-science/{DOMAIN_DIR}/`:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Isso fornece um fluxo de trabalho completo de processamento de amostra única.
    Você pode testar se ele executa com sucesso executando o seguinte comando:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Plano da lição

Dividimos isso em duas etapas:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Isso cobre a atualização de comandos e saídas de processos.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Isso introduz o operador `collect()` {AND_OTHER_CONCEPTS}.

!!! note "Nota"

     Certifique-se de estar no diretório de trabalho correto:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Lembre-se do comando modificado da [Parte 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Depois"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Antes"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Execute o fluxo de trabalho para verificar a modificação

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Saída do comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusão

Você sabe como modificar comandos e saídas de processos para adaptar o comportamento do fluxo de trabalho.

### O que vem a seguir?

Adicione a etapa de agregação de múltiplas amostras.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Escreva o módulo de agregação

{MODULE_INSTRUCTIONS}

### 2.2. Colete as saídas por amostra e alimente-as no processo de agregação

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Execute o fluxo de trabalho completo

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Saída do comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Conclusão

Você tem um pipeline completo que processa amostras individualmente e agrega resultados de todas as amostras.
Você sabe como usar operadores de canal como `collect()` para agregar saídas por amostra para análise de múltiplas amostras.

### O que vem a seguir?

Parabéns por completar este curso! Vá para o [resumo do curso](next_steps.md) para revisar o que você aprendeu e explorar os próximos passos.
