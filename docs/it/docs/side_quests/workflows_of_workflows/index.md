# Flussi di Lavoro di Flussi di Lavoro

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Quando si sviluppa una pipeline, spesso ci si ritrova a creare sequenze simili di processi per diversi tipi di dati o fasi di analisi. Si potrebbe finire per copiare e incollare queste sequenze di processi, portando a codice duplicato difficile da mantenere; oppure si potrebbe creare un unico flusso di lavoro massiccio difficile da comprendere e modificare.

Una delle funzionalità più potenti di Nextflow è la sua capacità di comporre pipeline complesse da moduli di flusso di lavoro più piccoli e riutilizzabili. Questo approccio modulare rende le pipeline più facili da sviluppare, testare e mantenere.

### Obiettivi di apprendimento

In questa side quest, esploreremo come sviluppare moduli di flusso di lavoro che possono essere testati e utilizzati separatamente, come comporre questi moduli in una pipeline più grande e come gestire il flusso di dati tra i moduli.

Al termine di questa side quest, sarete in grado di:

- Suddividere pipeline complesse in unità logiche e riutilizzabili
- Testare ogni modulo di flusso di lavoro in modo indipendente
- Combinare flussi di lavoro per creare nuove pipeline
- Condividere moduli di flusso di lavoro comuni tra diverse pipeline
- Rendere il codice più manutenibile e facile da comprendere

Queste competenze vi aiuteranno a costruire pipeline complesse mantenendo una struttura del codice pulita e manutenibile.

### Prerequisiti

Prima di affrontare questa side quest dovreste:

- Aver completato il tutorial [Hello Nextflow](../../hello_nextflow/index.md) o un corso equivalente per principianti.
- Essere a proprio agio con i concetti e i meccanismi di base di Nextflow (processi, canali, operatori, moduli)

---

## 0. Iniziamo

#### Aprite il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto nella sezione [Configurazione dell'Ambiente](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostatevi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

L'editor si apre con la directory del progetto in primo piano.

#### Esaminate i materiali

Troverete una directory `modules` con le definizioni dei processi, una directory `workflows` con due script di flusso di lavoro già scritti, e un file `main.nf` che aggiornerete progressivamente:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

La directory `modules/` contiene le definizioni dei singoli processi, e la directory `workflows/` contiene i due script di flusso di lavoro già scritti con cui lavorerete in questa side quest.

#### Esaminate il compito

La vostra sfida è assemblare questi moduli in due flussi di lavoro separati che comporremo poi in un flusso di lavoro principale:

- Un `GREETING_WORKFLOW` che valida i nomi, crea saluti e aggiunge timestamp
- Un `TRANSFORM_WORKFLOW` che converte il testo in maiuscolo e lo inverte

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista di controllo per la preparazione

Pensate di essere pronti a tuffarvi?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro
- [ ] Comprendo il compito assegnato

Se potete spuntare tutte le caselle, siete pronti a partire.

---

## 1. Aggiungere il greeting workflow alla pipeline

Il greeting workflow valida i nomi e genera saluti con timestamp.

### 1.1. Esaminare ed eseguire il greeting workflow

Aprite `workflows/greeting.nf` e date un'occhiata al codice:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Concatena i processi: valida -> crea saluto -> aggiungi timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Questo è un flusso di lavoro completo e autonomo con la stessa struttura vista nel tutorial 'Hello Nextflow'.
I nomi di input sono hardcoded, tre processi vengono concatenati e due output vengono pubblicati.

Eseguitelo per verificare che tutto funzioni:

```bash
nextflow run workflows/greeting.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Per renderlo componibile con altri flussi di lavoro, alcune cose devono cambiare.

### 1.2. Rendere il flusso di lavoro componibile

Per rendere un flusso di lavoro componibile, quattro cose devono cambiare:
il flusso di lavoro riceve un nome, gli input si spostano in un blocco `take:`, gli output si spostano in un blocco `emit:`,
e i blocchi standalone `publish:`/`output {}` vengono rimossi (appartengono all'entry workflow).

Esaminiamo queste modifiche una per una.

#### 1.2.1. Nominare il flusso di lavoro

Assegnate un nome al flusso di lavoro in modo che possa essere importato da un flusso di lavoro padre.

=== "Dopo"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Prima"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Con un nome, il flusso di lavoro può essere importato in altri script.

#### 1.2.2. Dichiarare gli input con `take:`

Sostituite la dichiarazione hardcoded del canale con un blocco `take:` che dichiara quali input il flusso di lavoro si aspetta.
Il blocco `take:` va prima di `main:`, e la riga `names_ch = channel.of(...)` viene rimossa.

=== "Dopo"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Canale di input con i nomi

        main:
        // Concatena i processi: valida -> crea saluto -> aggiungi timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Prima"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Concatena i processi: valida -> crea saluto -> aggiungi timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

Il blocco `take:` dichiara il canale solo per nome — i dettagli di cosa vi entra saranno definiti dal flusso di lavoro padre.

#### 1.2.3. Dichiarare gli output con `emit:`

Sostituite la sezione `publish:` e rimuovete il blocco `output {}`, sostituendoli con un blocco `emit:` che nomina gli output.

=== "Dopo"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Saluti originali
        timestamped = timestamped_ch // Saluti con timestamp
    }
    ```

=== "Prima"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

Il blocco `emit:` espone output nominati a cui i flussi di lavoro padre possono accedere tramite `GREETING_WORKFLOW.out.greetings` e `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Verificare il risultato e testarlo

Dopo tutte e tre le modifiche, il file completo dovrebbe apparire così:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Canale di input con i nomi

    main:
    // Concatena i processi: valida -> crea saluto -> aggiungi timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Saluti originali
    timestamped = timestamped_ch // Saluti con timestamp
}
```

Ora provate ad eseguirlo direttamente:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Questo introduce un concetto chiave: l'**entry workflow**.
Nextflow utilizza un blocco `workflow {}` senza nome come punto di ingresso quando si esegue uno script direttamente.
`GREETING_WORKFLOW` è nominato, quindi Nextflow non sa come eseguirlo da solo.

Questo è intenzionale — i flussi di lavoro componibili sono progettati per essere chiamati da un entry workflow, non eseguiti direttamente.
La soluzione è un entry workflow in `main.nf` che importa e chiama `GREETING_WORKFLOW`.

### 1.3. Aggiornare e testare il flusso di lavoro principale

Aggiorniamo ora il flusso di lavoro principale per chiamare il greeting workflow.

#### 1.3.1. Includere il greeting workflow e chiamarlo

Aggiungete l'istruzione `include`, aggiornate il corpo del flusso di lavoro per chiamare `GREETING_WORKFLOW` e sostituite il segnaposto `channel.empty()` in `publish:`:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Esegui il greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

L'entry workflow rimane senza nome in modo che Nextflow lo utilizzi come punto di ingresso della pipeline.

#### 1.3.2. Aggiornare il blocco output

Aggiungete una direttiva `path` per instradare i saluti pubblicati in una sottodirectory `greetings/`:

=== "Dopo"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Eseguire il flusso di lavoro

Eseguite il flusso di lavoro per verificare che funzioni:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    ```

??? abstract "Contenuto della directory"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "Contenuto del file"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

I file dei saluti vengono pubblicati in `results/greetings/`.
Il flusso di lavoro principale chiama `GREETING_WORKFLOW` e collega il suo output direttamente alla sezione `publish:`.

### Takeaway

In questa sezione, avete appreso diversi concetti importanti:

- **Flussi di lavoro nominati**: Creare un flusso di lavoro nominato (`GREETING_WORKFLOW`) che può essere importato e riutilizzato
- **Interfacce del flusso di lavoro**: Definire input chiari con `take:` e output con `emit:` per creare un flusso di lavoro componibile
- **Entry point**: Capire che Nextflow ha bisogno di un entry workflow senza nome per eseguire uno script
- **Composizione del flusso di lavoro**: Importare e utilizzare un flusso di lavoro nominato all'interno di un altro flusso di lavoro
- **Namespace del flusso di lavoro**: Accedere agli output del flusso di lavoro usando il namespace `.out` (`GREETING_WORKFLOW.out.greetings`)

Avete ora un greeting workflow funzionante che:

- Riceve un canale di nomi come input
- Valida ogni nome
- Crea un saluto per ogni nome valido
- Aggiunge timestamp ai saluti
- Espone sia i saluti originali che quelli con timestamp come output

Questo approccio modulare vi permette di testare il greeting workflow in modo indipendente o di utilizzarlo come componente in pipeline più grandi.

---

## 2. Aggiungere il transform workflow alla pipeline

Il transform workflow applica trasformazioni di testo ai saluti con timestamp.

### 2.1. Esaminare ed eseguire il flusso di lavoro

Aprite `workflows/transform.nf` e date un'occhiata al codice:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Applica le trasformazioni in sequenza
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Questo flusso di lavoro autonomo legge i file dei saluti con timestamp dalla directory `results/` prodotta da `greeting.nf`, li converte in maiuscolo e poi inverte il testo.

Eseguitelo per verificare che funzioni con i risultati del greeting dalla sezione 1.1:

```bash
nextflow run workflows/transform.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

Per renderlo componibile con `GREETING_WORKFLOW`, si applicano le stesse tre modifiche della sezione 1.2.

### 2.2. Renderlo componibile

Applicate le stesse tre modifiche della sezione 1.2: nominate il flusso di lavoro, sostituite l'input hardcoded con `take:`, e sostituite `publish:`/`output {}` con `emit:`.

Il file completato dovrebbe apparire così:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Canale di input con i messaggi

    main:
    // Applica le trasformazioni in sequenza
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Saluti in maiuscolo
    reversed = reversed_ch // Saluti in maiuscolo invertiti
}
```

Il transform workflow è ora componibile e pronto per essere importato nel flusso di lavoro principale.

### 2.3. Aggiornare e testare il flusso di lavoro principale

Aggiorniamo ora il flusso di lavoro principale per chiamare il transform workflow.

#### 2.3.1. Includere il transform workflow e chiamarlo

Aggiungete l'istruzione include, una chiamata a `TRANSFORM_WORKFLOW` collegata ai saluti con timestamp, e le due nuove voci `publish:`:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Esegui il greeting workflow
        GREETING_WORKFLOW(names)

        // Esegui il transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Esegui il greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Questo eseguirà il transform workflow sui saluti con timestamp.

#### 2.3.2. Aggiornare il blocco output

Aggiungete le voci `upper` e `reversed` al blocco `output {}`, ciascuna con una direttiva `path` per la sua sottodirectory:

=== "Dopo"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Questo pubblicherà gli output finali nelle directory appropriate.

#### 2.3.3. Eseguire la pipeline completa

Eseguite la pipeline per verificare che tutto funzioni:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Contenuto della directory"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "Contenuto del file"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

La pipeline funziona dall'inizio alla fine: il saluto è stato convertito in maiuscolo e invertito.

### Takeaway

Dovreste ora avere una pipeline completa che:

- Elabora i nomi attraverso il greeting workflow
- Invia i saluti con timestamp al transform workflow
- Produce versioni sia in maiuscolo che invertite dei saluti

---

## Riepilogo

In questa side quest, abbiamo esplorato il potente concetto di composizione di flussi di lavoro in Nextflow, che ci permette di costruire pipeline complesse da componenti più piccoli e riutilizzabili.

Questo approccio modulare offre diversi vantaggi rispetto alle pipeline monolitiche:

- Ogni flusso di lavoro può essere sviluppato, testato e debuggato in modo indipendente
- I flussi di lavoro possono essere riutilizzati in diverse pipeline
- La struttura complessiva della pipeline diventa più leggibile e manutenibile
- Le modifiche a un flusso di lavoro non influenzano necessariamente gli altri se le interfacce rimangono coerenti
- Gli entry point possono essere configurati per eseguire diverse parti della pipeline secondo le necessità

È importante notare tuttavia che, sebbene chiamare flussi di lavoro sia un po' come chiamare processi, non è esattamente la stessa cosa. Non è possibile, ad esempio, eseguire un flusso di lavoro N volte chiamandolo con un canale di dimensione N — sarebbe necessario passare un canale di dimensione N al flusso di lavoro e iterare internamente.

Applicare queste tecniche nel vostro lavoro vi permetterà di costruire pipeline Nextflow più sofisticate, in grado di gestire attività di elaborazione dati complesse rimanendo manutenibili e scalabili.

### Pattern chiave

1.  **Struttura del flusso di lavoro**: Abbiamo definito input e output chiari per ogni flusso di lavoro usando la sintassi `take:` e `emit:`, creando interfacce ben definite tra i componenti, e racchiuso la logica del flusso di lavoro all'interno del blocco `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // I canali di input sono dichiarati qui
            input_ch

        main:
            // La logica del flusso di lavoro va qui
            // Qui vengono chiamati i processi e manipolati i canali
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // I canali di output sono dichiarati qui
            output_ch = result_ch
    }
    ```

2.  **Import del flusso di lavoro:** Abbiamo costruito due moduli di flusso di lavoro indipendenti e li abbiamo importati in una pipeline principale con istruzioni `include`.

    - Includere un singolo flusso di lavoro

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Includere più flussi di lavoro

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Includere con alias per evitare conflitti di nomi

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry point**: Nextflow richiede un entry workflow senza nome per sapere da dove iniziare l'esecuzione. Questo entry workflow chiama i vostri flussi di lavoro nominati.

    - Flusso di lavoro senza nome (entry point)

    ```groovy
    workflow {
        // Questo è il punto di ingresso quando lo script viene eseguito
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Flusso di lavoro nominato (chiamato dall'entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Deve essere chiamato dall'entry workflow
    }
    ```

4.  **Gestione del flusso di dati:** Abbiamo imparato come accedere agli output del flusso di lavoro usando la notazione namespace (`WORKFLOW_NAME.out.channel_name`) e passarli ad altri flussi di lavoro.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Risorse aggiuntive

- [Documentazione Nextflow sui Workflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Riferimento agli Operatori dei Canali](https://www.nextflow.io/docs/latest/operator.html)
- [Documentazione sull'Error Strategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Cosa c'è dopo?

Tornate al [menu delle Side Quest](../index.md) oppure cliccate sul pulsante in basso a destra della pagina per passare all'argomento successivo nell'elenco.
