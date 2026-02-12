# Parte 2: Elaborazione di un singolo campione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella Parte 1, avete testato manualmente i comandi {TOOL_A} e {TOOL_B} nei rispettivi container.
Ora andremo a integrare quegli stessi comandi in un flusso di lavoro Nextflow.

## Obiettivo

In questa parte del corso, svilupperemo un flusso di lavoro che esegue le seguenti operazioni:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Questo replica i passaggi della Parte 1, dove avete eseguito questi comandi manualmente nei loro container.

Come punto di partenza, vi forniamo un file di flusso di lavoro, `{DOMAIN_DIR}.nf`, che delinea le parti principali del flusso di lavoro, insieme a due file modulo, {TOOL_A_MODULE}.nf e {TOOL_B_MODULE}.nf, che delineano la struttura dei moduli.
Questi file non sono funzionali; il loro scopo è solo quello di servire come strutture da riempire con le parti interessanti del codice.

## Piano della lezione

Per rendere il processo di sviluppo più educativo, abbiamo suddiviso il tutto in {N} passaggi:

1. **Scrivere un flusso di lavoro a singolo stadio che esegue {TOOL_A_ACTION}.**
   Questo copre la creazione di un modulo, la sua importazione e la sua chiamata in un flusso di lavoro.
2. **Aggiungere un secondo processo per eseguire {TOOL_B_ACTION}.**
   Questo introduce il concatenamento degli output dei processi agli input e la gestione dei file accessori.
3. **Adattare il flusso di lavoro per eseguirlo su un batch di campioni.**
   Questo copre l'esecuzione parallela e introduce le tuple per mantenere insieme i file associati.
4. **Far accettare al flusso di lavoro un foglio campioni contenente un batch di file di input.**
   Questo dimostra un pattern comune per fornire input in blocco.

Ogni passaggio si concentra su un aspetto specifico dello sviluppo del flusso di lavoro.

---

## 1. Scrivere un flusso di lavoro a singolo stadio che esegue {TOOL_A_ACTION}

Questo primo passaggio si concentra sulle basi: caricare {PRIMARY_INPUT_TYPE} e {TOOL_A_OUTPUT_DESCRIPTION}.

Ricordate il comando `{TOOL_A_COMMAND_NAME}` dalla [Parte 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

Il comando prende {INPUT_DESCRIPTION} e produce {OUTPUT_DESCRIPTION}.
L'URI del container era `{TOOL_A_CONTAINER_URI}`.

Andremo a integrare queste informazioni in Nextflow in tre fasi:

1. Configurare l'input
2. Scrivere il processo e chiamarlo nel flusso di lavoro
3. Configurare la gestione dell'output

### 1.1. Configurare l'input

Dobbiamo dichiarare un parametro di input, creare un profilo di test per fornire un valore predefinito conveniente e creare un canale di input.

#### 1.1.1. Aggiungere una dichiarazione di parametro di input

Nel file principale del flusso di lavoro `{DOMAIN_DIR}.nf`, sotto la sezione `Pipeline parameters`, dichiarate un parametro CLI chiamato `{PRIMARY_PARAM_NAME}`.

=== "Dopo"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Prima"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Questo configura il parametro CLI, ma non vogliamo digitare il percorso del file ogni volta che eseguiamo il flusso di lavoro durante lo sviluppo.
Ci sono diverse opzioni per fornire un valore predefinito; qui usiamo un profilo di test.

#### 1.1.2. Creare un profilo di test con un valore predefinito in `nextflow.config`

Un profilo di test fornisce valori predefiniti convenienti per provare un flusso di lavoro senza specificare input dalla riga di comando.
Questa è una convenzione comune nell'ecosistema Nextflow (vedi [Hello Config](../../hello_nextflow/06_hello_config.md) per maggiori dettagli).

Aggiungete un blocco `profiles` a `nextflow.config` con un profilo `test` che imposta il parametro `{PRIMARY_PARAM_NAME}` a uno dei {PRIMARY_INPUT_TYPE} di test.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Qui stiamo usando `${projectDir}`, una variabile integrata di Nextflow che punta alla directory dove si trova lo script del flusso di lavoro.
Questo rende facile fare riferimento a file di dati e altre risorse senza codificare percorsi assoluti.

#### 1.1.3. Configurare il canale di input

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Scrivere il modulo {TOOL_A_NAME}

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Importare e chiamare il modulo nel flusso di lavoro

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Eseguire il flusso di lavoro

A questo punto, abbiamo un flusso di lavoro a un solo passaggio che dovrebbe essere completamente funzionale.

Possiamo eseguirlo con `-profile test` per usare il valore predefinito configurato nel profilo di test ed evitare di dover scrivere il percorso sulla riga di comando.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Takeaway

Sapete come creare un modulo contenente un processo, importarlo in un flusso di lavoro, chiamarlo con un canale di input e pubblicare i risultati.

### Cosa c'è dopo?

Aggiungere un secondo processo per concatenare ulteriori passaggi di analisi.

---

## 2. Aggiungere un secondo processo per eseguire {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Ricordate il comando `{TOOL_B_COMMAND_NAME}` dalla [Parte 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Scrivere il modulo {TOOL_B_NAME}

{MODULE_INSTRUCTIONS}

### 2.2. Importare e chiamare il modulo nel flusso di lavoro

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Eseguire il flusso di lavoro

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Output del comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Takeaway

Sapete come concatenare gli output dei processi agli input e gestire i file accessori nel flusso di lavoro.

### Cosa c'è dopo?

Scalare il flusso di lavoro per elaborare più campioni in parallelo.

---

## 3. Adattare il flusso di lavoro per eseguirlo su un batch di campioni

Finora, il flusso di lavoro elabora un singolo campione.
Per gestire più campioni, dobbiamo modificare il modo in cui vengono forniti gli input e sfruttare il paradigma dataflow di Nextflow per parallelizzare l'esecuzione.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Eseguire il flusso di lavoro su più campioni

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Output del comando"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Takeaway

Sapete come sfruttare il paradigma dataflow di Nextflow per parallelizzare l'elaborazione per campione su più campioni di input.

### Cosa c'è dopo?

Nella [Parte 3](03_multi_sample.md), aggiungerete l'aggregazione multi-campione per combinare i risultati di tutti i campioni.
