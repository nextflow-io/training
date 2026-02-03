# Debug dei Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Il debug è una competenza critica che può farvi risparmiare ore di frustrazione e aiutarvi a diventare sviluppatori Nextflow più efficaci. Durante la vostra carriera, specialmente quando state iniziando, incontrerete bug durante la costruzione e manutenzione dei vostri workflow. Imparare approcci sistematici al debug vi aiuterà a identificare e risolvere i problemi rapidamente.

### Obiettivi di apprendimento

In questa missione secondaria, esploreremo **tecniche sistematiche di debug** per i workflow Nextflow:

- **Debug degli errori di sintassi**: Utilizzo efficace delle funzionalità IDE e dei messaggi di errore di Nextflow
- **Debug dei canali**: Diagnosi dei problemi di flusso dati e dei problemi di struttura dei canali
- **Debug dei processi**: Investigazione dei fallimenti di esecuzione e dei problemi di risorse
- **Strumenti di debug integrati**: Utilizzo della modalità preview, stub running e directory di lavoro di Nextflow
- **Approcci sistematici**: Una metodologia in quattro fasi per un debug efficiente

Alla fine, avrete una metodologia di debug robusta che trasforma frustranti messaggi di errore in chiare roadmap per le soluzioni.

### Prerequisiti

Prima di intraprendere questa missione secondaria, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso per principianti equivalente.
- Essere a vostro agio nell'uso di concetti e meccanismi base di Nextflow (processi, canali, operatori)

**Opzionale:** Raccomandiamo di completare prima la missione secondaria [Funzionalità IDE per lo Sviluppo Nextflow](./ide_features.md).
Quella copre in modo completo le funzionalità IDE che supportano il debug (evidenziazione sintassi, rilevamento errori, ecc.), che useremo intensamente qui.

---

## 0. Iniziare

#### Aprire il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/debugging
```

Potete impostare VSCode per focalizzarsi su questa directory:

```bash
code .
```

#### Rivedere i materiali

Troverete un insieme di workflow di esempio con vari tipi di bug che useremo per la pratica:

??? abstract "Contenuto della directory"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Questi file rappresentano scenari comuni di debug che incontrerete nello sviluppo reale.

#### Rivedere l'assegnazione

La vostra sfida è eseguire ogni workflow, identificare gli errori e correggerli.

Per ogni workflow con bug:

1. **Eseguire il workflow** e osservare l'errore
2. **Analizzare il messaggio di errore**: cosa vi sta dicendo Nextflow?
3. **Localizzare il problema** nel codice usando gli indizi forniti
4. **Correggere il bug** e verificare che la vostra soluzione funzioni
5. **Ripristinare il file** prima di passare alla sezione successiva (usare `git checkout <filename>`)

Gli esercizi procedono da semplici errori di sintassi a problemi di runtime più sottili.
Le soluzioni sono discusse inline, ma provate a risolvere ognuno da soli prima di leggere avanti.

#### Checklist di preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro appropriatamente
- [ ] Comprendo l'assegnazione

Se potete spuntare tutte le caselle, siete pronti per iniziare.

---

## 1. Errori di Sintassi

Gli errori di sintassi sono il tipo più comune di errore che incontrerete quando scrivete codice Nextflow. Avvengono quando il codice non si conforma alle regole di sintassi attese del DSL Nextflow. Questi errori impediscono al vostro workflow di eseguire del tutto, quindi è importante imparare come identificarli e correggerli rapidamente.

### 1.1. Parentesi graffe mancanti

Uno degli errori di sintassi più comuni, e talvolta uno dei più complessi da debuggare è quello delle **parentesi mancanti o non corrispondenti**.

Iniziamo con un esempio pratico.

#### Eseguire il pipeline

```bash
nextflow run bad_syntax.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Elementi chiave dei messaggi di errore di sintassi:**

- **File e posizione**: Mostra quale file e quale riga/colonna contengono l'errore (`bad_syntax.nf:24:1`)
- **Descrizione dell'errore**: Spiega cosa ha trovato il parser che non si aspettava (`Unexpected input: '<EOF>'`)
- **Indicatore EOF**: Il messaggio `<EOF>` (End Of File) indica che il parser ha raggiunto la fine del file mentre ancora si aspettava più contenuto - un segno classico di parentesi graffe non chiuse

#### Controllare il codice

Ora, esaminiamo `bad_syntax.nf` per capire cosa sta causando l'errore:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Missing closing brace for the process

workflow {

    // Crea canale di input
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Chiama il processo con il canale di input
    PROCESS_FILES(input_ch)
}
```

Per gli scopi di questo esempio abbiamo lasciato un commento per mostrarvi dove si trova l'errore. L'estensione Nextflow per VSCode dovrebbe anche darvi alcuni suggerimenti su cosa potrebbe essere sbagliato, mettendo la parentesi graffa non corrispondente in rosso e evidenziando la fine prematura del file:

![Bad syntax](img/bad_syntax.png)

**Strategia di debug per errori di parentesi:**

1. Usare la corrispondenza delle parentesi di VS Code (posizionare il cursore vicino a una parentesi)
2. Controllare il pannello Problemi per messaggi relativi alle parentesi
3. Assicurarsi che ogni `{` di apertura abbia una `}` di chiusura corrispondente

#### Correggere il codice

Sostituire il commento con la parentesi graffa di chiusura mancante:

=== "Dopo"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Add the missing closing brace

    workflow {

        // Crea canale di input
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chiama il processo con il canale di input
        PROCESS_FILES(input_ch)
    }
    ```

=== "Prima"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Missing closing brace for the process

    workflow {

        // Crea canale di input
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chiama il processo con il canale di input
        PROCESS_FILES(input_ch)
    }
    ```

#### Eseguire il pipeline

Ora eseguite nuovamente il workflow per confermare che funzioni:

```bash
nextflow run bad_syntax.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Uso di parole chiave o direttive di processo errate

Un altro errore di sintassi comune è una **definizione di processo non valida**. Questo può accadere se dimenticate di definire blocchi richiesti o usate direttive errate nella definizione del processo.

#### Eseguire il pipeline

```bash
nextflow run invalid_process.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Controllare il codice

L'errore indica una "Invalid process definition" e mostra il contesto intorno al problema. Guardando le righe 3-7, possiamo vedere `inputs:` alla riga 4, che è il problema. Esaminiamo `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Crea canale di input
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Chiama il processo con il canale di input
    PROCESS_FILES(input_ch)
}
```

Guardando la riga 4 nel contesto dell'errore, possiamo individuare il problema: stiamo usando `inputs` invece della corretta direttiva `input`. L'estensione Nextflow per VSCode segnalerà anche questo:

![Invalid process message](img/invalid_process_message.png)

#### Correggere il codice

Sostituire la parola chiave errata con quella corretta facendo riferimento [alla documentazione](https://www.nextflow.io/docs/latest/process.html#):

=== "Dopo"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crea canale di input
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chiama il processo con il canale di input
        PROCESS_FILES(input_ch)
    }
    ```

=== "Prima"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crea canale di input
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chiama il processo con il canale di input
        PROCESS_FILES(input_ch)
    }
    ```

#### Eseguire il pipeline

Ora eseguite nuovamente il workflow per confermare che funzioni:

```bash
nextflow run invalid_process.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Uso di nomi di variabili errati

I nomi delle variabili che usate nei vostri blocchi script devono essere validi, derivati o dagli input o da codice groovy inserito prima dello script. Ma quando state gestendo la complessità all'inizio dello sviluppo del pipeline, è facile fare errori nella denominazione delle variabili, e Nextflow ve lo farà sapere rapidamente.

#### Eseguire il pipeline

```bash
nextflow run no_such_var.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

L'errore viene catturato al momento della compilazione e punta direttamente alla variabile non definita alla riga 17, con un cursore che indica esattamente dove si trova il problema.

#### Controllare il codice

Esaminiamo `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Il messaggio di errore indica che la variabile non è riconosciuta nel template dello script, e eccola lì - dovreste essere in grado di vedere `${undefined_var}` usato nel blocco script, ma non definito altrove.

#### Correggere il codice

Se ricevete un errore 'No such variable', potete correggerlo definendo la variabile (correggendo i nomi delle variabili di input o modificando il codice groovy prima dello script), o rimuovendola dal blocco script se non è necessaria:

=== "Dopo"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Prima"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Eseguire il pipeline

Ora eseguite nuovamente il workflow per confermare che funzioni:

```bash
nextflow run no_such_var.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Uso errato delle variabili Bash

Iniziando con Nextflow, può essere difficile capire la differenza tra variabili Nextflow (Groovy) e Bash. Questo può generare un'altra forma di errore di variabile errata che appare quando si cerca di usare variabili nel contenuto Bash del blocco script.

#### Eseguire il pipeline

```bash
nextflow run bad_bash_var.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Controllare il codice

L'errore punta alla riga 13 dove viene usato `${prefix}`. Esaminiamo `bad_bash_var.nf` per vedere cosa sta causando il problema:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

In questo esempio, stiamo definendo la variabile `prefix` in Bash, ma in un processo Nextflow la sintassi `$` che abbiamo usato per riferirci ad essa (`${prefix}`) viene interpretata come una variabile Groovy, non Bash. La variabile non esiste nel contesto Groovy, quindi otteniamo un errore 'no such variable'.

#### Correggere il codice

Se volete usare una variabile Bash, dovete escludere il simbolo del dollaro in questo modo:

=== "Dopo"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Prima"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
        """
    }
    ```

Questo dice a Nextflow di interpretarla come una variabile Bash.

#### Eseguire il pipeline

Ora eseguite nuovamente il workflow per confermare che funzioni:

```bash
nextflow run bad_bash_var.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Variabili Groovy vs Bash"

    Per semplici manipolazioni di variabili come concatenazione di stringhe o operazioni di prefisso/suffisso, è solitamente più leggibile usare variabili Groovy nella sezione script piuttosto che variabili Bash nel blocco script:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Questo approccio evita la necessità di escludere i simboli del dollaro e rende il codice più facile da leggere e mantenere.

### 1.5. Istruzioni fuori dal blocco Workflow

L'estensione Nextflow per VSCode evidenzia problemi con la struttura del codice che causeranno errori. Un esempio comune è la definizione di canali fuori dal blocco `workflow {}` - questo ora è imposto come errore di sintassi.

#### Eseguire il pipeline

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Il messaggio di errore indica chiaramente il problema: le istruzioni (come le definizioni di canali) non possono essere mischiate con le dichiarazioni di script fuori da un blocco workflow o process.

#### Controllare il codice

Esaminiamo `badpractice_syntax.nf` per vedere cosa sta causando l'errore:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

Anche l'estensione VSCode evidenzierà la variabile `input_ch` come definita fuori dal blocco workflow:

![Non-lethal syntax error](img/nonlethal.png)

#### Correggere il codice

Spostare la definizione del canale all'interno del blocco workflow:

=== "Dopo"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "Prima"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Eseguire il pipeline

Eseguite nuovamente il workflow per confermare che la correzione funzioni:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Mantenete i vostri canali di input definiti all'interno del blocco workflow, e in generale seguite qualsiasi altra raccomandazione che l'estensione fa.

### Takeaway

Potete identificare e correggere sistematicamente gli errori di sintassi usando i messaggi di errore di Nextflow e gli indicatori visivi dell'IDE. Gli errori di sintassi comuni includono parentesi graffe mancanti, parole chiave di processo errate, variabili non definite e uso improprio di variabili Bash vs. Nextflow. L'estensione VSCode aiuta a catturare molti di questi prima del runtime. Con queste competenze di debug della sintassi nel vostro toolkit, sarete in grado di risolvere rapidamente gli errori di sintassi più comuni di Nextflow e passare ad affrontare problemi di runtime più complessi.

### Prossimi passi?

Imparate a debuggare errori di struttura dei canali più complessi che si verificano anche quando la sintassi è corretta.

---

## 2. Errori di Struttura dei Canali

Gli errori di struttura dei canali sono più sottili degli errori di sintassi perché il codice è sintatticamente corretto, ma le forme dei dati non corrispondono a ciò che i processi si aspettano. Nextflow cercherà di eseguire il pipeline, ma potrebbe trovare che il numero di input non corrisponde a ciò che si aspetta e fallire. Questi errori tipicamente appaiono solo al runtime e richiedono una comprensione dei dati che fluiscono attraverso il vostro workflow.

!!! tip "Debug dei Canali con `.view()`"

    Durante questa sezione, ricordate che potete usare l'operatore `.view()` per ispezionare il contenuto del canale in qualsiasi punto del vostro workflow. Questo è uno degli strumenti di debug più potenti per capire i problemi di struttura dei canali. Esploreremo questa tecnica in dettaglio nella sezione 2.4, ma sentitevi liberi di usarla mentre lavorate attraverso gli esempi.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Numero Errato di Canali di Input

Questo errore si verifica quando passate un numero diverso di canali rispetto a quelli che un processo si aspetta.

#### Eseguire il pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Controllare il codice

Il messaggio di errore afferma chiaramente che la chiamata si aspettava 1 argomento ma ne ha ricevuti 2, e punta alla riga 23. Esaminiamo `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Dovreste vedere la chiamata a `PROCESS_FILES` non corrispondente, che fornisce più canali di input quando il processo ne definisce solo uno. L'estensione VSCode sottolineerà anche la chiamata al processo in rosso e fornirà un messaggio diagnostico quando ci passate sopra con il mouse:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Correggere il codice

Per questo esempio specifico, il processo si aspetta un singolo canale e non richiede il secondo canale, quindi possiamo correggerlo passando solo il canale `samples_ch`:

=== "Dopo"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Prima"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Eseguire il pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Più comunemente di questo esempio, potreste aggiungere input aggiuntivi a un processo e dimenticare di aggiornare la chiamata al workflow di conseguenza, il che può portare a questo tipo di errore. Fortunatamente, questo è uno degli errori più facili da capire e correggere, poiché il messaggio di errore è abbastanza chiaro sulla mancata corrispondenza.

### 2.2. Esaurimento del Canale (Il Processo Viene Eseguito Meno Volte del Previsto)

Alcuni errori di struttura dei canali sono molto più sottili e non producono alcun errore. Probabilmente il più comune di questi riflette una sfida che i nuovi utenti Nextflow affrontano nel capire che i canali queue possono esaurirsi e rimanere senza elementi, il che significa che il workflow finisce prematuramente.

#### Eseguire il pipeline

```bash
nextflow run exhausted.nf
```

??? success "Output del comando"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Questo workflow si completa senza errori, ma elabora solo un singolo campione!

#### Controllare il codice

Esaminiamo `exhausted.nf` per vedere se è corretto:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variables in Groovy code before the script
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

Il processo viene eseguito solo una volta invece di tre volte perché il canale `reference_ch` è un canale queue che si esaurisce dopo la prima esecuzione del processo. Quando un canale si esaurisce, l'intero processo si ferma, anche se altri canali hanno ancora elementi.

Questo è un pattern comune in cui avete un singolo file di riferimento che deve essere riutilizzato attraverso più campioni. La soluzione è convertire il canale di riferimento in un canale value che può essere riutilizzato indefinitamente.

#### Correggere il codice

Ci sono un paio di modi per affrontare questo a seconda di quanti file sono coinvolti.

**Opzione 1**: Avete un singolo file di riferimento che state riutilizzando molto. Potete semplicemente creare un tipo di canale value, che può essere usato più e più volte. Ci sono tre modi per farlo:

**1a** Usare `channel.value()`:

```groovy title="exhausted.nf (corretto - Opzione 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Usare l'[operatore](https://www.nextflow.io/docs/latest/reference/operator.html#first) `first()`:

```groovy title="exhausted.nf (corretto - Opzione 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Usare l'[operatore](https://www.nextflow.io/docs/latest/reference/operator.html#collect) `collect()`:

```groovy title="exhausted.nf (corretto - Opzione 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opzione 2**: In scenari più complessi, forse dove avete più file di riferimento per tutti i campioni nel canale dei campioni, potete usare l'operatore `combine` per creare un nuovo canale che combina i due canali in tuple:

```groovy title="exhausted.nf (corretto - Opzione 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

L'operatore `.combine()` genera un prodotto cartesiano dei due canali, quindi ogni elemento in `reference_ch` sarà accoppiato con ogni elemento in `input_ch`. Questo permette al processo di eseguire per ogni campione mentre si usa ancora il riferimento.

Questo richiede che l'input del processo sia regolato. Nel nostro esempio, l'inizio della definizione del processo dovrebbe essere regolato come segue:

```groovy title="exhausted.nf (corretto - Opzione 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Questo approccio potrebbe non essere adatto in tutte le situazioni.

#### Eseguire il pipeline

Provate una delle correzioni sopra ed eseguite nuovamente il workflow:

```bash
nextflow run exhausted.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Dovreste ora vedere tutti e tre i campioni elaborati invece di uno solo.

### 2.3. Struttura del Contenuto del Canale Errata

Quando i workflow raggiungono un certo livello di complessità, può essere un po' difficile tenere traccia delle strutture interne di ogni canale, e le persone generano comunemente mancate corrispondenze tra ciò che il processo si aspetta e ciò che il canale contiene effettivamente. Questo è più sottile del problema di cui abbiamo discusso in precedenza, dove il numero di canali era errato. In questo caso, potete avere il numero corretto di canali di input, ma la struttura interna di uno o più di quei canali non corrisponde a ciò che il processo si aspetta.

#### Eseguire il pipeline

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Output del comando"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Controllare il codice

Le parentesi quadre nel messaggio di errore forniscono l'indizio qui - il processo sta trattando la tupla come un singolo valore, che non è ciò che vogliamo. Esaminiamo `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Potete vedere che stiamo generando un canale composto da tuple: `['sample1', 'file1.txt']`, ma il processo si aspetta un singolo valore, `val sample_name`. Il comando eseguito mostra che il processo sta cercando di creare un file chiamato `[sample3, file3.txt]_output.txt`, che non è l'output previsto.

#### Correggere il codice

Per correggere questo, se il processo richiede entrambi gli input potremmo regolare il processo per accettare una tupla:

=== "Opzione 1: Accettare tupla nel processo"

    === "Dopo"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Prima"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Opzione 2: Estrarre primo elemento"

    === "Dopo"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "Prima"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Eseguire il pipeline

Scegliete una delle soluzioni ed eseguite nuovamente il workflow:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Tecniche di Debug dei Canali

#### Uso di `.view()` per l'Ispezione dei Canali

Lo strumento di debug più potente per i canali è l'operatore `.view()`. Con `.view()`, potete capire la forma dei vostri canali in tutte le fasi per aiutare con il debug.

#### Eseguire il pipeline

Eseguite `bad_channel_shape_viewed.nf` per vedere questo in azione:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Controllare il codice

Esaminiamo `bad_channel_shape_viewed.nf` per vedere come viene usato `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### Correggere il codice

Per evitarvi di usare eccessivamente operazioni `.view()` in futuro per capire il contenuto dei canali, è consigliabile aggiungere alcuni commenti per aiutare:

```groovy title="bad_channel_shape_viewed.nf (con commenti)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Questo diventerà più importante man mano che i vostri workflow crescono in complessità e la struttura dei canali diventa più opaca.

#### Eseguire il pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### Takeaway

Molti errori di struttura dei canali possono essere creati con sintassi Nextflow valida. Potete debuggare gli errori di struttura dei canali comprendendo il flusso dei dati, usando operatori `.view()` per l'ispezione e riconoscendo pattern di errore come parentesi quadre che indicano strutture di tuple inaspettate.

### Prossimi passi?

Imparate sugli errori creati dalle definizioni dei processi.

---

## 3. Errori di Struttura dei Processi

La maggior parte degli errori che incontrerete relativi ai processi sarà legata a errori che avete fatto nella formazione del comando, o a problemi relativi al software sottostante. Detto questo, similmente ai problemi dei canali sopra, potete fare errori nella definizione del processo che non si qualificano come errori di sintassi, ma che causeranno errori al runtime.

### 3.1. File di Output Mancanti

Un errore comune quando si scrivono processi è fare qualcosa che genera una mancata corrispondenza tra ciò che il processo si aspetta e ciò che viene generato.

#### Eseguire il pipeline

```bash
nextflow run missing_output.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Controllare il codice

Il messaggio di errore indica che il processo si aspettava di produrre un file di output chiamato `sample3.txt`, ma lo script crea effettivamente `sample3_output.txt`. Esaminiamo la definizione del processo in `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

Dovreste vedere che c'è una mancata corrispondenza tra il nome del file di output nel blocco `output:` e quello usato nello script. Questa mancata corrispondenza causa il fallimento del processo. Se incontrate questo tipo di errore, tornate indietro e verificate che gli output corrispondano tra la definizione del vostro processo e il vostro blocco output.

Se il problema ancora non è chiaro, controllate la directory di lavoro stessa per identificare i file di output effettivamente creati:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Per questo esempio questo ci evidenzierebbe che un suffisso `_output` viene incorporato nel nome del file di output, contrariamente alla nostra definizione `output:`.

#### Correggere il codice

Correggete la mancata corrispondenza rendendo il nome del file di output coerente:

=== "Dopo"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Prima"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### Eseguire il pipeline

```bash
nextflow run missing_output.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Software mancante

Un'altra classe di errori si verifica a causa di errori nel provisioning del software. `missing_software.nf` è un workflow sintatticamente valido, ma dipende da del software esterno per fornire il comando `cowpy` che usa.

#### Eseguire il pipeline

```bash
nextflow run missing_software.nf
```

??? failure "Output del comando"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Il processo non ha accesso al comando che stiamo specificando. A volte questo è perché uno script è presente nella directory `bin` del workflow, ma non è stato reso eseguibile. Altre volte è perché il software non è installato nel container o nell'ambiente dove il workflow sta eseguendo.

#### Controllare il codice

Fate attenzione a quel codice di uscita `127` - vi dice esattamente il problema. Esaminiamo `missing_software.nf`:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Correggere il codice

Siamo stati un po' disonesti qui, e in realtà non c'è niente di sbagliato con il codice. Dobbiamo solo specificare la configurazione necessaria per eseguire il processo in modo tale che abbia accesso al comando in questione. In questo caso il processo ha una definizione container, quindi tutto ciò che dobbiamo fare è eseguire il workflow con Docker abilitato.

#### Eseguire il pipeline

Abbiamo configurato un profilo Docker per voi in `nextflow.config`, quindi potete eseguire il workflow con:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Per saperne di più su come Nextflow usa i container, vedere [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Configurazione delle risorse errata

Nell'uso in produzione, configurerete le risorse sui vostri processi. Per esempio `memory` definisce la quantità massima di memoria disponibile per il vostro processo, e se il processo supera quella, il vostro scheduler tipicamente ucciderà il processo e restituirà un codice di uscita `137`. Non possiamo dimostrarlo qui perché stiamo usando l'executor `local`, ma possiamo
