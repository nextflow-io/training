# Part 3: Agregació multi-mostra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la Part 2, vau construir un pipeline de processament per mostra que gestionava cada mostra de manera independent.
Ara l'ampliarem per implementar {AGGREGATION_METHOD} multi-mostra, tal com es va tractar a la [Part 1](01_method.md).

## Assignació

En aquesta part del curs, ampliarem el workflow per fer el següent:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Aquesta part es basa directament en el workflow produït a la Part 2.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat la [Part 2: Processament d'una sola mostra](./02_single_sample.md) i teniu un pipeline `{DOMAIN_DIR}.nf` funcional.

    Si no vau completar la Part 2 o voleu començar de nou per a aquesta part, podeu utilitzar la solució de la Part 2 com a punt de partida.
    Executeu aquestes comandes des de dins del directori `nf4-science/{DOMAIN_DIR}/`:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Això us proporciona un workflow complet de processament d'una sola mostra.
    Podeu comprovar que s'executa correctament executant la comanda següent:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Pla de la lliçó

Ho hem dividit en dos passos:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Això cobreix l'actualització de comandes i sortides de processos.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Això introdueix l'operador `collect()` {AND_OTHER_CONCEPTS}.

!!! note "Nota"

     Assegureu-vos que esteu al directori de treball correcte:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Recordeu la comanda modificada de la [Part 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Després"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Abans"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Executeu el workflow per verificar la modificació

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusió

Sabeu com modificar comandes i sortides de processos per adaptar el comportament del workflow.

### Què segueix?

Afegiu el pas d'agregació multi-mostra.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Escriviu el mòdul d'agregació

{MODULE_INSTRUCTIONS}

### 2.2. Recolliu les sortides per mostra i alimenteu-les al procés d'agregació

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Executeu el workflow complet

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Conclusió

Teniu un pipeline complet que processa mostres individualment i agrega resultats de totes les mostres.
Sabeu com utilitzar operadors de canal com `collect()` per agregar sortides per mostra per a l'anàlisi multi-mostra.

### Què segueix?

Felicitats per completar aquest curs! Aneu al [resum del curs](next_steps.md) per revisar el que heu après i explorar els passos següents.
