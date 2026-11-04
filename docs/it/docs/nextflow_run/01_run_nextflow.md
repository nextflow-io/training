# Parte 1: Eseguire Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa parte, introduciamo i concetti fondamentali dell'esecuzione di pipeline Nextflow.
Iniziamo con un semplice flusso di lavoro Hello World, per poi passare a una pipeline completa multi-step che elabora più input in parallelo usando container.

---

## 1. Hello World

Il flusso di lavoro `1-hello.nf` riceve un saluto tramite un argomento da riga di comando e lo scrive in un file.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Lanciare il flusso di lavoro

Eseguite il seguente comando nel vostro terminale.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Output del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

La riga chiave nell'output è la riga di stato del processo:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Questo ci dice che il processo `sayHello` è stato eseguito con successo una volta.
Il prefisso `[6d/740edd]` è un percorso abbreviato alla directory di lavoro dell'attività — ne parleremo più avanti.
Il blocco `Outputs:` che segue elenca ogni file pubblicato dalla pipeline, etichettato secondo il blocco `output` descritto nella sezione [1.4](#14-optional-code-walkthrough) qui sotto.

### 1.2. Trovare l'output

Questo flusso di lavoro è configurato per pubblicare il proprio output in una directory `results`.
Dopo l'esecuzione, dovreste trovare l'output lì:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Aprite il file per verificare che contenga `Hello World!`.

### 1.3. Esplorare la directory `work/`

Dietro le quinte, Nextflow crea una directory di attività unica per ogni chiamata di processo all'interno di una directory chiamata `work/`.
L'hash mostrato nell'output della console (`[6d/740edd]`) è il percorso verso quella directory.

```bash
ls work/6d/740edd*
```

Al suo interno troverete il file di output insieme a diversi file di log nascosti:

- **`.command.sh`**: il comando esatto eseguito da Nextflow
- **`.command.out`** / **`.command.err`**: stdout e stderr del processo
- **`.command.log`**: output combinato del log
- **`.exitcode`**: il codice di uscita del processo

Il file `.command.sh` è particolarmente utile durante il debug — mostra esattamente cosa è stato eseguito.

### 1.4. Opzionale: Analisi del codice

Capire il codice non è essenziale se volete solo eseguire pipeline, ma se siete curiosi, vale la pena darci un'occhiata.

??? optional "Cliccate per esplorare il codice associato a questo esercizio"

    Apriamo `1-hello.nf` e osserviamo i suoi componenti principali.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Parametri della pipeline
     */
    params {
        input: String
    }

    workflow {

        main:
        // emette un saluto
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Vediamo quanto segue:

    - un'istruzione `include` che punta a un modulo `process`
    - un blocco `params` che definisce i parametri della pipeline
    - un blocco `workflow` che descrive il lavoro da svolgere
    - un blocco `output` che descrive cosa fare con gli output

    Esaminiamo ciascuno di essi.

    ### Il modulo `process`

    L'istruzione `include` dice a Nextflow di caricare qualcosa chiamato `sayHello` da un file di codice separato.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    In quel file troviamo la definizione di un processo chiamato `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

    Un **process** definisce un singolo step della pipeline.
    Dichiara i propri input, output e lo script da eseguire.
    Il qualificatore `val` indica che l'input è un valore semplice (stringa, numero, ecc.).
    Il qualificatore `path` indica che l'output è un percorso di file.

    È possibile scrivere la definizione del processo nel file principale del flusso di lavoro, ma tenerli in file di modulo separati li rende riutilizzabili: lo stesso modulo può essere importato da più script di flusso di lavoro.

    ### Il blocco `params`

    Il blocco `params` dichiara i parametri da riga di comando accettati dal flusso di lavoro:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Ogni parametro dichiarato qui diventa disponibile da riga di comando con un doppio trattino (`--input`).
    I tipi supportati includono `String`, `Integer`, `Float`, `Boolean` e `Path`.

    !!! tip "Suggerimento"

        I parametri del flusso di lavoro usano sempre due trattini (`--input`) per distinguerli dai flag CLI propri di Nextflow, che usano un solo trattino (ad es. `-resume`).

    ### Il blocco `workflow`

    Il blocco **workflow** definisce la logica del flusso di dati: quali processi eseguire e in quale ordine.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // emette un saluto
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Qui viene chiamato un solo processo, quindi è molto semplice; vedremo esempi più realistici in seguito.

    La sezione `main:` chiama il processo `sayHello` con il valore di `--input`.
    La sezione `publish:` elenca quali output devono essere copiati nella directory dei risultati.

    ### Il blocco `output`

    Il blocco `output` in fondo al file specifica il percorso di destinazione e la modalità di copia.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Ogni voce con nome corrisponde a un'etichetta `publish:` nel flusso di lavoro e la mappa a una sottodirectory sotto `results/`.

### Takeaway

Sapete come eseguire una pipeline Nextflow e trovare i suoi output, e sapete che il lavoro viene eseguito nelle directory di attività sotto `work/`.

### Cosa c'è dopo?

Scopriamo come Nextflow gestisce più input in modo efficiente.

---

## 2. Elaborare più input

Le pipeline del mondo reale tipicamente elaborano molti dati, non solo uno.
Il flusso di lavoro `2-inputs.nf` legge da un file CSV ed esegue `sayHello` una volta per ogni riga, in parallelo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Eseguiamo prima il flusso di lavoro, poi vedremo quale meccanismo usa Nextflow per gestire questi input multipli.

### 2.1. Eseguire il flusso di lavoro

Eseguite il seguente comando nel vostro terminale.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Output del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

Il `3 of 3` ci dice che il processo `sayHello` è stato chiamato tre volte, una per ogni riga del CSV.

Nella directory `results`, dovreste ora vedere tre file di output, uno per ogni saluto:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Aprite uno qualsiasi dei file di output per verificare che ciascuno contenga un saluto.

L'output condensato mostrato sopra mostra una singola riga di riepilogo per `sayHello`, ma Nextflow ha effettivamente avviato tre esecuzioni di attività separate, una per ogni riga del CSV, eseguendole in parallelo non appena la macchina aveva le risorse necessarie.

Proprio come la singola attività che avete esplorato nella sezione [1.3](#13-explore-the-work-directory), ognuna di queste tre esecuzioni ottiene la propria directory di attività sotto `work/`, completamente isolata dalle altre:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Ogni `.command.sh` contiene solo il comando per quel singolo saluto:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Questo isolamento è ciò che rende sicura l'esecuzione parallela: tre attività in esecuzione contemporaneamente non condividono mai una directory di lavoro, quindi nulla di ciò che scrive un'attività può collidere con o sovrascrivere ciò che scrive un'altra attività, anche se producono file con lo stesso nome.
È anche per questo che `-resume` (trattato di seguito) può memorizzare nella cache e riutilizzare le singole attività in modo indipendente: gli input, gli output e i log di ogni attività risiedono interamente nella propria directory, senza nulla di condiviso tra le attività che potrebbe andare fuori sincronia.

### 2.2. Eseguire il flusso di lavoro di nuovo con `-ansi-log false`

Per impostazione predefinita, Nextflow condensa l'output in una singola riga di riepilogo per processo.
Per vedere ogni chiamata di processo elencata individualmente, aggiungete `-ansi-log false`:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Questo mostra tutte e tre le chiamate di processo e la sottodirectory di lavoro unica creata per ciascuna.

### 2.3. Usare `-resume` per saltare il lavoro già completato

Ora passate al file di input esteso, che aggiunge altri due saluti, e aggiungete `-resume` alla riga di comando:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Output del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow ha eseguito solo i due nuovi input.
I tre saluti elaborati nell'esecuzione precedente sono stati memorizzati nella cache e riutilizzati automaticamente.

Questo funziona anche per saltare l'esecuzione di processi per step già completati con successo in una pipeline multi-step.
Ad esempio, se l'esecuzione di una pipeline è stata interrotta da un errore di sistema, o se avete aggiunto nuovi step a una pipeline in sviluppo.

La funzionalità `-resume` è particolarmente preziosa nelle pipeline lunghe, dove il recupero da un errore può far risparmiare tempo e risorse critiche.

### 2.4. Opzionale: Analisi del codice

Capire il codice non è essenziale se volete solo eseguire pipeline, ma se siete curiosi, vale la pena darci un'occhiata.

??? optional "Cliccate per esplorare il codice associato a questo esercizio"

    La modifica chiave in `2-inputs.nf` si trova nella sezione `main:` del flusso di lavoro:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
        sayHello(greeting_ch)
    ```

    Quello che vedete qui si chiama **canale**: un costrutto a coda che gestisce i dati di input in modo da rendere semplice parallelizzare le operazioni.

    - `channel.fromPath(params.input)` crea un canale dal percorso del file fornito con `--input`
    - `.splitCsv()` analizza il CSV in righe
    - `#!groovy .map { line -> line[0] }` estrae la prima colonna da ogni riga

    Il risultato è un canale contenente `Hello`, `Bonjour` e `Hola`.
    Quando viene passato a `sayHello(greeting_ch)`, Nextflow chiama automaticamente il processo una volta per ogni elemento, eseguendoli in parallelo quando le risorse lo consentono.

### Takeaway

Sapete come elaborare più input da un file CSV in parallelo, e come usare `-resume` per evitare di ripetere il lavoro già completato.

### Cosa c'è dopo?

Scopriamo come una pipeline multi-step completa collega i processi tra loro usando i canali, e come usare i container per gestire gli strumenti di analisi e le loro dipendenze.

---

## 3. Eseguire una pipeline multi-step

Finora avete eseguito un singolo processo, poi lo avete eseguito più volte in parallelo su un insieme di input.
Le pipeline reali di solito vanno oltre: collegano diversi processi tra loro, passando l'output di uno come input del successivo, e spesso si affidano a più di un software lungo il percorso.
Il flusso di lavoro `main.nf` mette insieme entrambe queste caratteristiche in una pipeline completa.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Ogni saluto di input scorre attraverso tutti e quattro gli step: `sayHello` lo scrive in un file, `convertToUpper` converte il testo in maiuscolo, `collectGreetings` unisce tutti i risultati in un unico file, e `cowpy` genera ASCII art dall'output unito usando uno strumento containerizzato.
Nextflow collega questi step con i canali: l'output di un processo diventa l'input del successivo, così l'intera catena viene eseguita automaticamente man mano che i dati diventano disponibili, senza che dobbiate orchestrare ogni step manualmente.

Si noti che questo flusso di lavoro usa i moduli: ogni processo è definito nel proprio file sotto `modules/`, e `main.nf` li importa con istruzioni `include` invece di definirli inline.
Questo rende ogni processo riutilizzabile in più flussi di lavoro senza duplicare il codice. Per saperne di più, consultate la sezione di analisi del codice più avanti.

### 3.1. Eseguire il flusso di lavoro

Eseguite il seguente comando nel vostro terminale.

```bash
nextflow run main.nf --input data/greetings.csv
```

Il parametro `character` ha come valore predefinito `turkey` in `nextflow.config`, quindi l'ASCII art usa un tacchino a meno che non lo sovrascriviate (provate ad aggiungere `--character tux`).

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Quattro processi sono stati eseguiti, ma non lo stesso numero di volte.
`sayHello` e `convertToUpper` sono stati eseguiti ciascuno una volta per input (3 of 3): ogni saluto deve essere scritto e convertito in maiuscolo singolarmente.
`collectGreetings` e `cowpy` sono stati eseguiti ciascuno una sola volta (1 of 1): unire i saluti e generare l'ASCII art ha senso solo quando tutti i risultati individuali sono pronti.
Questa forma di fan-out-poi-fan-in, diverse attività parallele che confluiscono in un numero minore di attività a valle, è comune nelle pipeline reali.

Nextflow non aspetta che un intero step sia completato prima di avviare il successivo.
Non appena un output di `sayHello` è pronto, la corrispondente attività `convertToUpper` può iniziare, quindi le attività di processi diversi vengono eseguite in modo concorrente piuttosto che in batch rigidi.
`collectGreetings` e `cowpy` devono invece aspettare, poiché ciascuno dipende dalla disponibilità di tutti i risultati a monte.

La directory `results` riflette quel fan-in, più ciò che l'autore della pipeline ha scelto di pubblicare e dove: ricordate il blocco `output` dall'analisi del codice nella sezione 1.4, che è ciò che definisce questa struttura.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

La directory di primo livello prende il nome dal parametro `batch`, che ha come valore predefinito `batch`; lo vedrete cambiare negli esercizi successivi.

Controllate `cowpy-COLLECTED-batch-output.txt` per il file di ASCII art.

??? abstract "Contenuto del file"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
     ---------
      \                                  ,+*^^*+___+++_
       \                           ,*^^^^              )
        \                       _+*                     ^**+_
         \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                 {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
               {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
               U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
             (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
               (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                 ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                         ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Proprio come nella sezione [2.1](#21-run-the-workflow), ognuna di queste 8 esecuzioni di attività, su tutti e quattro i processi, ottiene la propria directory sotto `work/`, completamente isolata dalle altre.
`collectGreetings` è un buon esempio di perché questo è importante: dipende dagli output di tutte e tre le attività `convertToUpper`, che risiedono in tre directory di attività diverse, quindi Nextflow crea symlink a quei file all'interno della propria directory di `collectGreetings` invece di fargli leggere direttamente dalle directory delle attività a monte:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Ogni attività vede solo i file specifici di cui ha bisogno, indipendentemente da dove provengono, e mai il contenuto interno della directory di un'altra attività.
In tutta una pipeline, lo stesso isolamento che avete visto con un singolo processo nella sezione [2.1](#21-run-the-workflow) è ciò che permette a Nextflow di eseguire ogni attività di ogni processo in modo concorrente e sicuro.

!!! note "Nota"

    Lo step `cowpy` viene eseguito all'interno di un container Docker invece di affidarsi al software installato localmente.
    Un container raggruppa un'applicazione insieme a tutto ciò di cui ha bisogno per essere eseguita, così non dovete installare e gestire le dipendenze voi stessi, e la pipeline si comporta allo stesso modo su qualsiasi macchina in grado di eseguire il container.
    Nextflow supporta anche Conda come alternativa ai container; consultate la [Parte 2](./02_configure_pipeline.md) per sapere come passare da uno all'altro.

### 3.2. Opzionale: Analisi del codice

Capire il codice non è essenziale se volete solo eseguire pipeline, ma se siete curiosi, vale la pena darci un'occhiata.

??? optional "Cliccate per esplorare il codice associato a questo esercizio"

    ### Come i dati scorrono da uno step al successivo

    Ogni processo passa il proprio canale di output al successivo:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    Il pattern `processName.out` fa riferimento al canale di output di un processo.

    L'operatore `.collect()` raccoglie tutti gli output individuali di `convertToUpper` in un singolo elemento del canale prima di passarli a `collectGreetings`.

    ### Usare i moduli di processo

    `main.nf` non definisce direttamente alcun codice di processo.
    Invece, importa ogni processo dal proprio file sotto `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Ogni file di modulo contiene una singola definizione di processo, strutturata allo stesso modo del modulo `sayHello` nella sezione [1.4](#14-optional-code-walkthrough).
    Tenere i processi in file separati li rende riutilizzabili in più flussi di lavoro senza duplicare il codice.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Usare software containerizzato

    Il processo `cowpy` viene eseguito all'interno di un container Docker specificato nel suo file di modulo:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

    Nextflow scarica automaticamente l'immagine, esegue lo script all'interno del container e fa pulizia al termine.
    Docker è abilitato per questo progetto in `nextflow.config`:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Questa singola riga abilita Docker per qualsiasi processo nella pipeline che abbia un container specificato.

### Takeaway

Avete eseguito una pipeline multi-step completa che elabora più input in parallelo usando uno strumento containerizzato.

### Cosa c'è dopo?

Passate alla [Parte 2](./02_configure_pipeline.md), dove imparerete come configurare il comportamento della pipeline usando `nextflow.config`.

---

## Riepilogo

In questa parte avete imparato a:

- Eseguire un flusso di lavoro Nextflow e trovare i suoi output
- Esplorare la directory `work/` e i suoi file di log
- Elaborare più input da un file CSV in parallelo
- Usare `-resume` per saltare il lavoro già completato quando si aggiungono nuovi input
- Eseguire una pipeline multi-step che usa uno strumento containerizzato
