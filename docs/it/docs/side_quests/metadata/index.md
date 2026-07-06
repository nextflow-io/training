# Metadati e meta map

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In qualsiasi analisi scientifica, raramente lavoriamo solo con i file di dati grezzi.
Ogni file porta con sé informazioni aggiuntive: cosa contiene, da dove proviene e cosa lo rende speciale.
Queste informazioni extra sono ciò che chiamiamo metadati.

I metadati sono dati che descrivono altri dati.
I metadati tracciano dettagli importanti sui file e sulle condizioni sperimentali, e aiutano ad adattare le analisi alle caratteristiche uniche di ogni dataset.

Pensateci come a un catalogo di biblioteca: mentre i libri contengono il contenuto effettivo (dati grezzi), le schede catalografiche forniscono informazioni essenziali su ogni libro — quando è stato pubblicato, chi lo ha scritto, dove trovarlo (metadati).
Nei pipeline Nextflow, i metadati possono essere utilizzati per:

- Tracciare informazioni specifiche dei file durante il flusso di lavoro
- Configurare i processi in base alle caratteristiche dei file
- Raggruppare file correlati per analisi congiunte

### Obiettivi di apprendimento

In questa side quest, esploreremo come gestire i metadati nei flussi di lavoro.
Partendo da un semplice foglio dati (spesso chiamato samplesheet in bioinformatica) contenente informazioni di base sui file, imparerete come:

- Leggere e analizzare i metadati dei file da file CSV
- Capire perché l'interfaccia "meta map + file di dati" è una convenzione ampiamente utilizzata
- Aggiungere nuovi campi di metadati durante l'esecuzione del flusso di lavoro
- Utilizzare i metadati per personalizzare il comportamento dei processi e organizzare gli output

Queste competenze vi aiuteranno a costruire pipeline più robusti e flessibili, in grado di gestire relazioni complesse tra file e requisiti di elaborazione.

### Prerequisiti

Prima di affrontare questa side quest, dovreste:

- Aver completato il tutorial [Hello Nextflow](../../hello_nextflow/index.md) o un corso equivalente per principianti.
- Essere a proprio agio con i concetti e i meccanismi di base di Nextflow (processi, canali, operatori)

---

## 0. Iniziamo

#### Aprite il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto nella sezione [Configurazione dell'ambiente](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostatevi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/metadata
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

L'editor si apre con la directory del progetto in primo piano.

#### Esaminate i materiali

Troverete un file principale del flusso di lavoro e una directory `data` contenente un foglio dati e alcuni file di dati.

??? abstract "Contenuto della directory"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

Il flusso di lavoro nel file `main.nf` è uno stub che espanderete gradualmente in un flusso di lavoro completamente funzionante.

Il foglio dati elenca i percorsi dei file di dati e alcuni metadati associati, organizzati in 3 colonne:

- `id`: autoesplicativo, un ID assegnato al file
- `character`: il nome di un personaggio, che useremo in seguito per disegnare diverse creature
- `data`: percorsi ai file `.txt` che contengono saluti in diverse lingue

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Ogni file di dati contiene del testo di saluto in una delle cinque lingue (fr: francese, de: tedesco, es: spagnolo, it: italiano, en: inglese).

Useremo uno strumento chiamato [`COWPY`](https://github.com/jeffbuttars/cowpy) per generare ASCII art di ogni personaggio che pronuncia il saluto registrato.

??? info "Cosa fa `COWPY`?"

    `COWPY` è uno strumento da riga di comando che genera ASCII art per visualizzare input di testo arbitrari in modo divertente.
    È un'implementazione Python del classico strumento [cowsay](https://en.wikipedia.org/wiki/Cowsay) di Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Facoltativamente, potete selezionare un personaggio (o 'cowacter') da usare al posto della mucca predefinita.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Inoltre, useremo uno strumento di analisi linguistica chiamato `langid` per identificare la lingua parlata da ogni personaggio e organizzare di conseguenza gli output del pipeline.

#### Esaminate il compito

La vostra sfida è scrivere un flusso di lavoro Nextflow che:

1. **Generi ASCII art** di ogni personaggio
2. **Organizzi** gli output per famiglia linguistica (lingue germaniche vs lingue romanze)

Questo rappresenta un tipico schema di flusso di lavoro in cui i metadati specifici dei file guidano le decisioni di elaborazione; esattamente il tipo di problema che le meta map risolvono elegantemente.

#### Lista di controllo per la preparazione

Pensate di essere pronti a tuffarvi?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro
- [ ] Comprendo il compito assegnato

Se potete spuntare tutte le caselle, siete pronti a partire.

---

## 1. Opzioni di base per caricare e usare i metadati

Aprite il file del flusso di lavoro `main.nf` per esaminare lo stub del flusso di lavoro che vi forniamo come punto di partenza.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

L'operatore [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) legge ogni riga del file come elemento nel canale.
Questo è lo stesso approccio che usiamo per caricare dati CSV in Hello Nextflow, il nostro corso per principianti.
Date un'occhiata a [questa sezione](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) se avete bisogno di un ripasso su come funziona.

Con `header: true`, la prima riga viene trattata come intestazione delle colonne, quindi ogni elemento diventa una mappa di coppie chiave-valore indicizzate per nome di colonna.

Notate che poiché non stiamo ancora eseguendo alcun processo sui dati, i blocchi `publish` e `output` sono solo stub.

### 1.1. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per vedere come è strutturato il contenuto del canale una volta che tutto è caricato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Come potete vedere, l'operatore ha costruito una mappa di coppie chiave-valore per ogni riga del file CSV, con le intestazioni delle colonne come chiavi per i valori corrispondenti.

Ogni voce della mappa corrisponde a una colonna nel nostro foglio dati:

- `id`
- `character`
- `recording`

Questo rende facile accedere a campi specifici da ogni riga.
Ad esempio, potremmo accedere all'ID del file con `id` o al percorso del file txt con `recording`.

??? info "(Opzionale) Ulteriori informazioni sulle map Groovy"

    In Groovy, il linguaggio di programmazione su cui è costruito Nextflow, una map è una struttura dati chiave-valore simile ai dizionari in Python, agli oggetti in JavaScript o agli hash in Ruby.

    Ecco uno script eseguibile che mostra come definire una map e accedere al suo contenuto in pratica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crea una semplice map
    def my_map = [id:'sampleA', character:'squirrel']

    // Stampa l'intera map
    println "map: ${my_map}"

    // Accede ai singoli valori usando la notazione con punto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Anche se non ha un blocco `workflow` vero e proprio, Nextflow può eseguirlo come se fosse un flusso di lavoro:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Ed ecco cosa potete aspettarvi di vedere nell'output:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Selezionare un campo specifico con `map`

Useremo l'operatore `map` per iterare su ogni elemento in un canale ed estrarre solo il campo `character`, a cui possiamo accedere per nome usando la notazione con punto.

#### 1.2.1. Aggiungere l'operazione map

Per accedere alla colonna `character`, aggiungete l'operazione `map` prima dell'operazione `.view()` come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Questo modo di accedere a un campo specifico è spiegato in maggior dettaglio in [questa sezione](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) di Hello Nextflow, se avete bisogno di un ripasso su come funziona.

#### 1.2.2. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per verificare che possiate visualizzare i nomi dei personaggi estratti.

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Questo mostra che siamo in grado di accedere ai valori dalla colonna `character` per ogni riga.

Ora facciamo qualcosa con questi dati: usiamo i campi `character` e `recording` insieme per generare ASCII art usando `COWPY`.

### 1.3. Emettere sotto-canali con `multiMap`

Vi forniamo un modulo di processo pre-scritto per `COWPY`, quindi prima dovete esaminare i requisiti di input del processo.

Potete aprire il file per vedere come appare il processo:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Genera ASCII art con cowpy
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Come potete vedere, il processo riceve due input separati: un file di registrazione e un nome di personaggio.
Abbiamo valori per entrambi, ma sono attualmente raggruppati all'interno di ogni elemento nel canale.

Un modo per estrarre più campi in canali separati è l'operatore [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap), che divide un canale in più sotto-canali con nome in una singola operazione.

#### 1.3.1. Aggiungere l'operazione multiMap

Sostituite l'operazione `map` con `multiMap`:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

Il blocco `multiMap` definisce due sotto-canali con nome (`file` e `character`) da ogni riga, a cui possiamo accedere come `ch_datasheet.file` e `ch_datasheet.character`.

#### 1.3.2. Chiamare COWPY sui sotto-canali

Ora includete il processo `COWPY` e passategli ogni sotto-canale come argomento separato:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Questo ci permette di passare i due campi separatamente come richiede `COWPY`.

#### 1.3.3. Configurare la pubblicazione dell'output

Infine, aggiungete l'output di `COWPY` al blocco `publish:`:

=== "Dopo"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Questo ci permetterà di visualizzare facilmente gli output prodotti dal flusso di lavoro.

#### 1.3.4. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per verificare che `COWPY` venga eseguito sugli input che abbiamo fornito:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Come potete vedere, `COWPY` è stato eseguito su ogni file usando il personaggio corretto per ognuno.

??? abstract "Contenuti della directory dei risultati"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Contenuto di results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Questo approccio funziona, ma ha un limite: abbiamo dovuto dividere il canale in due sotto-canali separati.
Se volessimo passare più campi al processo, dovremmo dividerli in ancora più sotto-canali.
Questo potrebbe diventare fastidioso e disordinato.

Buone notizie: esiste un modo più semplice per farlo.

### 1.4. Raggruppare tutto come un singolo input al processo

Invece di dividere i campi in canali separati, possiamo aggiornare il processo per ricevere tutti gli input come una singola tupla, il che semplifica la chiamata al processo.

#### 1.4.1. Aggiornare il processo COWPY

Aggiornate `COWPY` per accettare una tupla corrispondente ai tre elementi in ogni riga:

=== "Dopo"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera ASCII art con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Prima"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Genera ASCII art con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
        """
    }
    ```

Ora il processo riceve un solo input contenente tutto ciò che potremmo volergli passare.

#### 1.4.2. Usare `map()` per creare la tupla di input

Dobbiamo comunque usare un'operazione di mapping per enumerare gli elementi che vogliamo passare nella tupla al processo:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Potreste chiedervi perché non possiamo semplicemente passare l'intera map Groovy proveniente da `splitCsv` così com'è.
Il motivo è che dobbiamo dire esplicitamente a Nextflow che il file di registrazione deve essere gestito come un path (cioè deve essere messo in staging correttamente).
Questo avviene a livello dell'interfaccia di input di `COWPY`, dove l'elemento `recording` è esplicitamente designato come `path`.

#### 1.4.3. Aggiornare la chiamata al processo

Infine, sostituiamo i due input separati nella chiamata al processo con la singola tupla che abbiamo appena creato:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Questo semplifica un po' la chiamata al processo.

#### 1.4.4. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per verificare che `COWPY` possa ancora elaborare i dati correttamente:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

L'output è lo stesso sette file `cowpy-*.txt` di prima, ora prodotti con una chiamata più semplice a `COWPY`.

??? abstract "Contenuti della directory dei risultati"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Contenuto di results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Questo è un leggero miglioramento rispetto all'approccio con `multiMap`.
Ma abbiamo comunque dovuto spacchettare la map Groovy originale per creare la tupla di input, e c'è un forte accoppiamento tra il processo e il foglio dati: la definizione di input di `COWPY` ora fa riferimento direttamente ai nomi delle colonne `id`, `character` e `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Se un collaboratore usa un foglio dati strutturato diversamente — con colonne aggiuntive, o colonne in un ordine diverso — questo processo non funzionerà senza modifiche.
Questo rende il processo fragile, perché la sua struttura di input è legata alla composizione esatta del foglio dati.

Per risolvere questo problema, abbiamo bisogno di un modo per passare tutti i metadati come un pacchetto senza codificare la sua struttura esatta nell'interfaccia del processo.

### 1.5. Usare un'interfaccia meta map + file

La soluzione è separare due aspetti distinti nel canale: i **metadati su un campione** e il **file di dati** stesso.
Raggruppando tutti i metadati in una singola map — la "meta map" — otteniamo una tupla consistente a due elementi indipendentemente da quante colonne di metadati contiene il foglio dati:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Aggiungere o rimuovere colonne dal foglio dati cambia ciò che è dentro `meta`, ma la forma della tupla `[meta, file]` rimane costante.
I processi che accettano questa struttura non hanno bisogno di sapere o preoccuparsi di quanti campi di metadati esistono.

#### 1.5.1. Riorganizzare il contenuto della tupla in una meta map

Ristrutturiamo l'operazione `map` per produrre una tupla `[meta, file]`:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Aggiorneremo nel passo successivo

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Noterete che abbiamo anche aggiunto un'istruzione `view()`, commentato la chiamata a `COWPY` e sostituito `COWPY.out` con `channel.empty()` perché la definizione di input del processo non corrisponde ancora alla nuova struttura.

#### 1.5.2. Eseguire il flusso di lavoro per ispezionare il contenuto riorganizzato

Eseguiamo il flusso di lavoro per vedere la nuova forma del canale:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Ogni elemento nel canale è ora una tupla a due elementi: prima la meta map, poi il file.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Se in seguito aggiungiamo una colonna `language` al foglio dati, sarà disponibile come `meta.language` senza richiedere alcuna modifica alla definizione di input del processo.

#### 1.5.3. Aggiornare il processo `COWPY` per usare la meta map

Aggiornate `COWPY` per accettare la struttura tupla `[meta, file]`:

=== "Dopo"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera ASCII art con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Prima"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera ASCII art con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

All'interno del blocco script, `meta.character` accede al campo `character` dalla meta map.
Qualsiasi campo nella meta map è accessibile nello stesso modo.

#### 1.5.4. Aggiornare la chiamata al processo

Ripristinate la chiamata a `COWPY` e collegate il suo output per la pubblicazione:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Aggiorneremo nel passo successivo

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

Abbiamo anche ripristinato la pubblicazione dell'output.

#### 1.5.5. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per verificare che tutto funzioni:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

La directory dei risultati ora contiene i file ASCII art.

??? abstract "Contenuto della directory"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Contenuto di results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Il processo ora riceve tutti i metadati come pacchetto tramite `meta`, usa ciò di cui ha bisogno (`meta.character`) e ignora il resto.

Questa è l'interfaccia standard usata da tutti i moduli [nf-core](https://nf-co.re/).
Lo schema `tuple val(meta), path(file)` appare in modo consistente in tutta la libreria di moduli nf-core, motivo per cui i flussi di lavoro che adottano questa convenzione possono integrare moduli nf-core con il minimo attrito.

### Takeaway

In questa sezione, avete imparato:

- **Come leggere i fogli dati:** Usare `splitCsv` per analizzare file CSV con informazioni di intestazione
- **Perché esiste la convenzione della meta map:** Separare i metadati dai file di dati in tuple `[meta, file]` mantiene la struttura del canale stabile man mano che il foglio dati evolve
- **Come usare i campi della meta map all'interno di un processo:** Qualsiasi campo nella meta map è accessibile tramite notazione con punto nel blocco script

---

## 2. Manipolazioni aggiuntive dei metadati

Ora che l'interfaccia della meta map è in place, possiamo arricchirla man mano che i dati scorrono attraverso il pipeline.

Useremo uno strumento chiamato [`langid`](https://github.com/saffsd/langid.py) per identificare la lingua in ogni file di registrazione.
Dato un frammento di testo, produce una previsione della lingua e un punteggio di probabilità su `stdout`.

### 2.1. Aggiungere un passo di identificazione della lingua

Vi forniamo un modulo di processo pre-scritto chiamato `IDENTIFY_LANGUAGE` che racchiude lo strumento `langid`.

Aprite il file del modulo per esaminarne il codice:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
// Usa langid per prevedere la lingua di ogni file di input
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

La definizione di input usa la stessa struttura `tuple val(meta), path(file)` che abbiamo appena costruito nella Sezione 1, quindi `ch_datasheet` può alimentare direttamente questo processo senza alcun adattamento.

L'output aggiunge `stdout` come terzo elemento: questo cattura la previsione della lingua che `langid` stampa sulla console.
Il comando `sed` rimuove il punteggio di probabilità e il carattere di nuova riga finale, lasciando solo il codice lingua a due lettere.

#### 2.1.1. Aggiungere una chiamata a `IDENTIFY_LANGUAGE`

Includete il modulo del processo `IDENTIFY_LANGUAGE` e chiamatelo sul canale del foglio dati:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

L'output principale di questo processo è solo una stringa, quindi non ci sono file di output da pubblicare.
Usiamo invece `IDENTIFY_LANGUAGE.out.view()` per visualizzare i risultati dell'operazione.

#### 2.1.2. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per produrre l'identificazione della lingua, usando `-resume` per evitare di rieseguire le attività `COWPY`:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Ora abbiamo una previsione della lingua per ogni file nel dataset.

Notate che la tupla di output è composta da `[meta, file, lang_id]`, il che significa che la meta map e il file vengono trasportati insieme al nuovo risultato.

!!! note "Nota"

    Questo schema di mantenere la meta map associata ai risultati rende più facile unire i risultati tra canali in seguito.
    Non potete fare affidamento sull'ordine degli elementi nei canali per associare i dati correttamente.
    Dovete invece usare le chiavi.
    Le meta map forniscono una struttura ideale per questo scopo.

    Questo caso d'uso è esplorato in dettaglio nella side quest [Splitting & Grouping](../splitting_and_grouping/index.md).

### 2.2. Arricchire i metadati con gli output dei processi

La previsione della lingua è essa stessa un metadato sul contenuto del file.
Invece di mantenerla come elemento separato, incorporiamola nella meta map.

#### 2.2.1. Creare una meta map nuova e ampliata

Possiamo creare una nuova meta map per sostituire quella originale usando l'operatore Groovy `+`:

=== "Dopo"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Il cuore di questa operazione è `#!groovy meta + [lang: lang_id]`.

Questo codice crea essenzialmente una map temporanea con una singola coppia chiave-valore contenente il codice lingua (`[lang: lang_id]`), poi usa l'operatore Groovy `+` per combinarla con la map `meta` originale contenente i metadati preesistenti, producendo una meta map nuova e ampliata.

Per una spiegazione più dettagliata, vedete il riquadro qui sotto.

??? info "Creazione della nuova meta map usando l'operatore `+`"

    **Prima di tutto, dovete sapere che possiamo unire il contenuto di due map usando l'operatore Groovy `+`.**

    Supponiamo di avere le seguenti map:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Possiamo unirle così:

    ```groovy
    new_map = map1 + map2
    ```

    Il contenuto di `new_map` sarà:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Ottimo!

    **Ma cosa succede se dovete aggiungere un campo che non fa già parte di una map?**

    Supponiamo di ricominciare da `map1`, ma la previsione della lingua non è nella sua propria map (non esiste `map2`).
    Invece, è contenuta in una variabile chiamata `lang_id`, e sapete di voler memorizzare il suo valore (`'fr'`) con la chiave `lang`.

    Potete effettivamente fare quanto segue:

    ```groovy
    new_map = map1 + [lang: lang_id]
    ```

    Qui, `[lang: lang_id]` crea una nuova map senza nome al volo, e `map1 + ` unisce `map1` con la nuova map senza nome, producendo lo stesso contenuto di `new_map` come prima.

    Elegante, vero?

    **Ora trasponiamo questo nel contesto di un'operazione `channel.map()` di Nextflow.**

    Il codice diventa:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Questo fa quanto segue:

    - `#!groovy map1, lang_id ->` prende i due elementi nella tupla
    - `#!groovy map1 + [lang: lang_id]` crea la nuova map come descritto sopra

    L'output è una singola map senza nome con lo stesso contenuto di `new_map` nel nostro esempio precedente.
    Quindi abbiamo effettivamente trasformato:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    in:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Speriamo che possiate vedere che se cambiamo `map1` in `meta`, è fondamentalmente tutto ciò di cui abbiamo bisogno per aggiungere la previsione della lingua alla nostra meta map nel nostro flusso di lavoro.

    Tranne per una cosa!

    Nel caso del nostro flusso di lavoro, **dobbiamo anche tenere conto della presenza dell'oggetto `file` nella tupla**, che è composta da `meta, file, lang_id`.

    Quindi il codice qui diventerebbe:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Se avete difficoltà a capire perché il `file` sembra spostarsi nell'operazione `map`, immaginate che invece di `#!groovy [meta + [lang: lang_id], file]`, quella riga reciti `[new_map, file]`.
    Questo dovrebbe rendere più chiaro che stiamo semplicemente lasciando il `file` nella sua posizione originale in seconda posizione nella tupla. Abbiamo semplicemente preso il valore `new_info` e lo abbiamo incorporato nella mappa che si trova in prima posizione.

    **E questo ci riporta alla struttura del canale `tuple val(meta), path(file)`!**

#### 2.2.2. Eseguire il flusso di lavoro

Una volta che siete sicuri di capire cosa fa questo codice, eseguiamo il flusso di lavoro per vedere se funziona:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Sì, funziona!
Abbiamo riorganizzato ordinatamente l'output del processo da `meta, file, lang_id` in modo che `lang_id` sia ora una delle chiavi nella meta map, e le tuple del canale si adattano nuovamente al modello `meta, file`.

!!! tip "Rimuovere chiavi da una meta map"

    Potete rimuovere una chiave da una meta map usando il metodo Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)), che restituisce una nuova map contenente solo le chiavi specificate:

    ```groovy
    meta.subMap(['id', 'character'])  // restituisce una map con solo 'id' e 'character'
    ```

    Questo è utile quando un processo o modulo a valle non ha bisogno di tutti i campi che si sono accumulati nella meta map.

### 2.3. Assegnare un gruppo linguistico usando i condizionali

Con la previsione della lingua nella meta map, possiamo derivarne ulteriori metadati.
Le lingue nel nostro dataset appartengono a due famiglie: germaniche (inglese, tedesco) e romanze (francese, spagnolo, italiano).
Aggiungere un campo `lang_group` renderà quella classificazione disponibile a valle.

#### 2.3.1. Aggiungere un'operazione `map` con la logica condizionale

Useremo una seconda operazione `map` con logica condizionale per assegnare la famiglia linguistica:

```groovy
.map { meta, file ->

    // la logica condizionale che definisce lang_group va qui

    [meta + [lang_group: lang_group], file]
}
```

Ecco la logica da applicare:

- Iniziare con `lang_group = 'unknown'` come valore predefinito.
- Se `meta.lang` è `'de'` o `'en'`, impostare `lang_group` su `'germanic'`.
- Altrimenti se `meta.lang` è in `['fr', 'es', 'it']`, impostare `lang_group` su `'romance'`.

!!! tip "Suggerimento"

    Potete accedere al valore di `lang` all'interno dell'operazione map con `meta.lang`.

Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        ch_languages.view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Punti chiave:

- `def lang_group = "unknown"` inizializza la variabile con un valore predefinito sicuro.
- La struttura `if / else if` gestisce le due famiglie linguistiche; tutto il resto rimane `'unknown'`.
- `#!groovy .set { ch_languages }` assegna un nome al canale risultante per l'uso nel passo successivo.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per verificare che funzioni:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

La meta map ora contiene quattro campi: `id`, `character`, `lang` e `lang_group`.
La struttura del canale è ancora `[meta, file]`.

### 2.4. Usare i metadati per nominare e organizzare gli output

Con `lang` e `lang_group` ora disponibili nella meta map, possiamo usarli per aggiungere un codice lingua ai nomi dei file di output e organizzarli in sottodirectory per famiglia linguistica.

Questo richiede tre modifiche: aggiornare il processo `COWPY` per rinominare il suo output e includere `meta` in ciò che emette, aggiornare la chiamata a `COWPY` per eseguirlo su `ch_languages`, e aggiornare il blocco output per specificare il percorso della sottodirectory.

#### 2.4.1. Aggiornare il processo `COWPY`

Rinominate il file di output usando il codice lingua dalla meta map, e aggiungete `meta` all'output in modo che il blocco output possa accedere a `lang_group` per il routing nelle sottodirectory:

=== "Dopo"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Prima"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Questo mostra come possiamo sfruttare altri campi dei metadati per personalizzare il comportamento di un processo, senza dover modificare affatto la definizione di input.

#### 2.4.2. Aggiornare la chiamata a `COWPY` per eseguirlo su `ch_languages`

Sostituite `COWPY(ch_datasheet)` con `COWPY(ch_languages)`:

=== "Dopo"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

Rimuoviamo anche la riga `ch_languages.view()` poiché non abbiamo più bisogno di ispezionare il contenuto del canale.

#### 2.4.3. Aggiornare il blocco output

Aggiungete una closure `path` al blocco `output {}` per instradare ogni file nella sottodirectory del suo gruppo linguistico:

=== "Dopo"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Questo mostra come possiamo usare i metadati per organizzare gli output con grande flessibilità.

#### 2.4.4. Eseguire il pipeline completo

Eliminate i risultati dell'esecuzione precedente ed eseguiamo il pipeline completo:

```bash
rm -r results
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

La directory dei risultati è ora organizzata per famiglia linguistica, con ogni file nominato in base alla lingua rilevata:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

La closure `path` nel blocco `output {}` riceve ogni tupla `[meta, file]` e restituisce `meta.lang_group` come nome della sottodirectory.
Il nome del file stesso proviene da ciò che il processo produce (`#!groovy "${meta.lang}-${input_file}"`).
Entrambi i pezzi di metadati (codice lingua e gruppo linguistico) provengono dalla meta map arricchita costruita in questa sezione.

### Takeaway

In questa sezione, avete imparato:

- **Come arricchire la meta map con gli output dei processi:** Aggiungere nuove chiavi con `#!groovy meta + [key: value]` mantiene intatta la struttura del canale `[meta, file]` arricchendo i metadati.
- **Come derivare metadati dai metadati:** La logica condizionale all'interno di un'operazione `map` può calcolare nuovi campi da quelli esistenti.
- **Come usare i metadati per l'organizzazione degli output:** La closure `path` nel blocco `output {}` può leggere dalla meta map per instradare i file nelle sottodirectory.

---

## 3. Considerazioni sulla robustezza

Quando i valori dei metadati guidano il comportamento dei processi, dati mancanti o incompleti possono causare problemi difficili da diagnosticare.
Ecco cosa aspettarsi e come gestirlo.

### 3.1. Cosa succede quando manca un campo di metadati obbligatorio

Il valore `character` è necessario affinché il processo `COWPY` produca un risultato valido.
La modalità di errore dipende dal fatto che la colonna esista nel foglio dati ma sia vuota, o sia assente del tutto.

#### 3.1.1. La colonna esiste ma un valore è vuoto

Supponiamo che una voce nel foglio dati abbia il campo `character` vuoto:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La chiave `character` viene creata per tutte le voci quando il foglio dati viene analizzato, ma `meta.character` per `sampleA` sarà una stringa vuota.
Quando Nextflow sostituisce `#!groovy ${meta.character}` nel comando, lo strumento `COWPY` riceve un argomento vuoto per `-c` e fallisce:

??? failure "Output del comando"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Il messaggio di errore (`expected one argument`) indica il flag `-c` vuoto.
Controllare il file `.command.sh` nella directory di lavoro conferma che il comando è stato eseguito con un valore vuoto.

#### 3.1.2. La colonna non esiste nel foglio dati

Se la colonna `character` è assente del tutto:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La chiave `character` non viene mai creata nella meta map.
Quando lo script del processo valuta `#!groovy ${meta.character}`, la chiave mancante restituisce `null`, e Nextflow sostituisce letteralmente la stringa `null` nel comando:

??? failure "Output del comando"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Il `cowpy -c null` nel comando eseguito è l'indizio diagnostico.

### 3.2. Strategie per gestire i metadati mancanti

Esistono due approcci complementari per rendere i flussi di lavoro più robusti contro i metadati mancanti.

**1. Validazione dell'input**

La soluzione più affidabile è validare il foglio dati prima che inizi qualsiasi elaborazione, in modo che i problemi vengano rilevati tempestivamente con un messaggio di errore chiaro piuttosto che emergere come un errore criptico del processo a metà esecuzione.
La formazione [Hello nf-core](../../hello_nf-core/05_input_validation.md) spiega come aggiungere la validazione dell'input usando il plugin nf-schema. <!-- TODO (future) pending a proper Validation side quest -->

**2. Input espliciti del processo per i valori obbligatori**

Se volete che l'interfaccia del processo comunichi che un particolare valore è obbligatorio, considerate di estrarlo dalla meta map come input esplicito:

=== "Definizione del processo"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Chiamata nel flusso di lavoro"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Questo approccio rende `character` una parte visibile e obbligatoria del contratto del processo.
Chiunque legga il modulo può vedere immediatamente che deve essere fornito un valore per il personaggio.
Se il campo è assente, il flusso di lavoro fallisce chiaramente a livello del canale prima ancora che il processo venga eseguito.

Questo evidenzia un utile principio di progettazione:

**Usate la meta map per informazioni opzionali o descrittive; estraete i valori obbligatori come input espliciti.**

La meta map mantiene le strutture dei canali pulite e stabili, ma per i valori genuinamente richiesti da un processo, esporli come input con nome migliora la chiarezza e rende il modulo più facile da usare correttamente in altri contesti.

### Takeaway

In questa sezione, avete visto:

- **Come si manifestano i metadati mancanti:** Un campo vuoto produce un argomento vuoto; un campo assente produce `null` sostituito letteralmente nel comando.
- **Due strategie complementari:** La validazione dell'input per rilevare i problemi tempestivamente, e gli input espliciti del processo per comunicare chiaramente i requisiti.

---

## Riepilogo

In questa side quest, avete esplorato come lavorare efficacemente con i metadati nei flussi di lavoro Nextflow.

Lo schema della tupla "meta map + file di dati" è una convenzione fondamentale in Nextflow, che offre diversi vantaggi rispetto al passaggio dei metadati come valori individuali:

- La struttura del canale rimane stabile man mano che il foglio dati evolve
- Il comportamento dei processi può essere personalizzato per campione senza codificare i nomi dei campi
- I metadati sono disponibili durante tutto il pipeline per nominare, raggruppare e organizzare gli output
- I moduli scritti per questa interfaccia sono intercambiabili, inclusi i moduli nf-core

### Schemi chiave

1.  **Lettura e strutturazione dei metadati:** Analizzare un foglio dati CSV e creare una meta map.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Espansione dei metadati durante il flusso di lavoro:** Aggiungere nuove chiavi dagli output dei processi o dalla logica derivata.

    ```groovy
    // Da un output di processo
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // Da logica condizionale
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Usare i metadati all'interno di un processo:** Accedere a qualsiasi campo tramite notazione con punto nel blocco script.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organizzare gli output per valore dei metadati:** Usare la closure `path` nel blocco `output {}`.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Risorse aggiuntive

- [operatore map](https://www.nextflow.io/docs/latest/operator.html#map)
- [operatore multiMap](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [qualificatore di output stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Cosa c'è dopo?

Tornate al [menu delle Side Quest](../index.md) o cliccate il pulsante in basso a destra della pagina per passare all'argomento successivo nell'elenco.
