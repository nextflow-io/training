# Parte 3: Aggregazione multi-campione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella Parte 2, avete costruito una pipeline di elaborazione per-campione che gestiva ogni campione in modo indipendente.
Ora la estenderemo per implementare l'{AGGREGATION_METHOD} multi-campione, come trattato nella [Parte 1](01_method.md).

## Compito

In questa parte del corso, estenderemo il flusso di lavoro per fare quanto segue:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Questa parte si basa direttamente sul flusso di lavoro prodotto dalla Parte 2.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato la [Parte 2: Elaborazione di un singolo campione](./02_single_sample.md) e che abbiate una pipeline `{DOMAIN_DIR}.nf` funzionante.

    Se non avete completato la Parte 2 o volete ripartire da zero per questa parte, potete usare la soluzione della Parte 2 come punto di partenza.
    Eseguite questi comandi dall'interno della directory `nf4-science/{DOMAIN_DIR}/`:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Questo vi fornisce un flusso di lavoro completo di elaborazione per singolo campione.
    Potete verificare che funzioni correttamente eseguendo il seguente comando:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Piano della lezione

Abbiamo suddiviso questo in due passaggi:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Questo copre l'aggiornamento dei comandi e degli output dei processi.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Questo introduce l'operatore `collect()` {AND_OTHER_CONCEPTS}.

!!! note "Nota"

     Assicuratevi di essere nella directory di lavoro corretta:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Ricordiamo il comando modificato dalla [Parte 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Dopo"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Prima"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Eseguire il flusso di lavoro per verificare la modifica

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Output del comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Takeaway

Sapete come modificare i comandi e gli output dei processi per adattare il comportamento del flusso di lavoro.

### Cosa c'è dopo?

Aggiungere il passaggio di aggregazione multi-campione.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Scrivere il modulo di aggregazione

{MODULE_INSTRUCTIONS}

### 2.2. Raccogliere gli output per-campione e fornirli al processo di aggregazione

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Eseguire il flusso di lavoro completo

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Output del comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Takeaway

Avete una pipeline completa che elabora i campioni individualmente e aggrega i risultati di tutti i campioni.
Sapete come usare operatori di canale come `collect()` per aggregare gli output per-campione per l'analisi multi-campione.

### Cosa c'è dopo?

Congratulazioni per aver completato questo corso! Andate al [riepilogo del corso](next_steps.md) per rivedere ciò che avete imparato ed esplorare i prossimi passi.
