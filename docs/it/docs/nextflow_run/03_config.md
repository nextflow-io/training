# Parte 3: Configurazione dell'esecuzione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Questa sezione esplorerà come gestire la configurazione di una pipeline Nextflow per personalizzarne il comportamento, adattarla a diversi ambienti e ottimizzare l'uso delle risorse _senza alterare una singola riga del codice del workflow stesso_.

Ci sono molteplici modi per farlo, che possono essere usati in combinazione e vengono interpretati secondo l'ordine di precedenza descritto nella documentazione [Configuration](https://nextflow.io/docs/latest/config.html).

In questa parte del corso, vi mostreremo il meccanismo di file di configurazione più semplice e comune, il file `nextflow.config`, che avete già incontrato nella sezione sui container nella Parte 2.

Esamineremo i componenti essenziali della configurazione di Nextflow come le direttive dei process, gli executor, i profili e i file di parametri.
Imparando a utilizzare efficacemente queste opzioni di configurazione, potete sfruttare appieno la flessibilità, la scalabilità e le prestazioni delle pipeline Nextflow.

Per esercitare questi elementi di configurazione, eseguiremo una copia fresca del workflow che abbiamo eseguito per ultimo alla fine della Parte 2 di questo corso di formazione, rinominato `3-main.nf`.

Se non avete familiarità con la pipeline Hello o potreste aver bisogno di un promemoria, consultate [questa pagina informativa](../info/hello_pipeline.md).

---

## 1. Gestire i parametri di input del workflow

??? example "Scenario"

    Avete scaricato una pipeline e volete eseguirla ripetutamente con gli stessi file di input e impostazioni, ma non volete digitare tutti i parametri ogni volta.
    O forse state configurando la pipeline per un collega che non è a suo agio con gli argomenti da riga di comando.

Inizieremo con un aspetto della configurazione che è semplicemente un'estensione di ciò con cui abbiamo lavorato finora: la gestione dei parametri di input.

Attualmente, il nostro workflow è configurato per accettare diversi valori di parametri tramite la riga di comando, dichiarati in un blocco `params` nello script del workflow stesso.
Uno ha un valore predefinito impostato come parte della sua dichiarazione.

Tuttavia, potreste voler impostare valori predefiniti per tutti loro, o sovrascrivere il valore predefinito esistente senza dover specificare parametri sulla riga di comando o modificare il file di script originale.

Ci sono molteplici modi per farlo; vi mostreremo tre modi base che sono molto comunemente usati.

### 1.1. Imposta i valori in `nextflow.config`

Questo è l'approccio più semplice, sebbene sia probabilmente il meno flessibile poiché il file `nextflow.config` principale non è qualcosa che volete modificare per ogni esecuzione.
Ma ha il vantaggio di separare le preoccupazioni della _dichiarazione_ dei parametri nel workflow (che decisamente appartiene lì) rispetto alla fornitura dei _valori predefiniti_, che sono più a loro agio in un file di configurazione.

Facciamolo in due step.

#### 1.1.1. Crea un blocco `params` nel file di configurazione

Fate le seguenti modifiche al codice nel file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Notate che non abbiamo semplicemente copiato il blocco `params` dal workflow al file di configurazione.
Per il parametro `batch` che aveva già un valore predefinito dichiarato, la sintassi è un po' diversa.
Nel file del workflow, quella è una dichiarazione tipizzata.
Nella configurazione, quelle sono assegnazioni di valori.

Tecnicamente, questo è sufficiente per sovrascrivere i valori predefiniti ancora specificati nel file del workflow.
Potreste modificare il valore predefinito per `batch` ed eseguire il workflow per verificare che il valore impostato nel file di configurazione sovrascrive quello impostato nel file del workflow.

Ma nello spirito di spostare la configurazione completamente al file di configurazione, rimuoviamo quel valore predefinito dal file del workflow completamente.

#### 1.1.2. Rimuovi il valore predefinito per `batch` nel file del workflow

Fate la seguente modifica al codice nel file del workflow `3-main.nf`:

=== "Dopo"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Prima"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Ora il file del workflow stesso non imposta alcun valore predefinito per questi parametri.

#### 1.1.3. Esegui la pipeline

Testiamo che funzioni correttamente senza specificare alcun parametro nella riga di comando.

```bash
nextflow run 3-main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima.

L'output finale dell'arte ASCII è nella directory `results/3-main/`, sotto il nome `cowpy-COLLECTED-batch-output.txt`, come prima.

??? abstract "Contenuto del file"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
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

Funzionalmente, questo spostamento non ha cambiato nulla, ma concettualmente è un po' più pulito avere i valori predefiniti impostati nel file di configurazione.

### 1.2. Usa un file di configurazione specifico per l'esecuzione

??? example "Scenario"

    Volete sperimentare con diverse impostazioni senza modificare il vostro file di configurazione principale.

Potete farlo creando un nuovo file `nextflow.config` in una sottodirectory che userete come directory di lavoro per i vostri esperimenti.

#### 1.2.1. Crea la directory di lavoro con una configurazione vuota

Iniziamo creando una nuova directory e spostandoci al suo interno:

```bash
mkdir -p tux-run
cd tux-run
```

Poi, create un file di configurazione vuoto in quella directory:

```bash
touch nextflow.config
```

Questo produce un file vuoto.

#### 1.2.2. Imposta la configurazione sperimentale

Ora aprite il nuovo file e aggiungete i parametri che volete personalizzare:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Notate che il percorso del file di input deve riflettere la struttura delle directory.

#### 1.2.3. Esegui la pipeline

Ora possiamo eseguire la nostra pipeline dall'interno della nostra nuova directory di lavoro.
Assicuratevi di adattare il percorso di conseguenza!

```bash
nextflow run ../3-main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Questo creerà un nuovo set di directory sotto `tux-run/` incluse `tux-run/work/` e `tux-run/results/`.

In questa esecuzione, Nextflow combina il `nextflow.config` nella nostra directory corrente con il `nextflow.config` nella directory root della pipeline, e quindi sovrascrive il personaggio predefinito (turkey) con il personaggio tux.

Il file di output finale dovrebbe contenere il personaggio tux che dice i saluti.

??? abstract "Contenuto del file"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

Ecco fatto; ora avete uno spazio per sperimentare senza modificare la vostra configurazione 'normale'.

!!! warning "Avviso"

    Assicuratevi di tornare alla directory precedente prima di passare alla prossima sezione!

    ```bash
    cd ..
    ```

Ora guardiamo un altro modo utile per impostare i valori dei parametri.

### 1.3. Usa un file di parametri

??? example "Scenario"

    Dovete condividere i parametri esatti dell'esecuzione con un collaboratore, o registrarli per una pubblicazione.

L'approccio con la sottodirectory funziona benissimo per sperimentare, ma comporta un po' di setup e richiede che adattiate i percorsi di conseguenza.
C'è un approccio più semplice per quando volete eseguire la vostra pipeline con un set specifico di valori, o permettere a qualcun altro di farlo con il minimo sforzo.

Nextflow ci permette di specificare parametri tramite un [file di parametri](https://nextflow.io/docs/latest/config.html#parameter-file) in formato YAML o JSON, il che lo rende molto conveniente per gestire e distribuire set alternativi di valori predefiniti, per esempio, così come valori di parametri specifici per l'esecuzione.

#### 1.3.1. Esamina il file di parametri di esempio

Per dimostrare questo, forniamo un file di parametri di esempio nella directory corrente, chiamato `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Questo file di parametri contiene una coppia chiave-valore per ciascuno degli input che vogliamo specificare.
Notate l'uso dei due punti (`:`) invece dei segni di uguale (`=`) se confrontate la sintassi con il file di configurazione.
Il file config è scritto in Groovy, mentre il file di parametri è scritto in YAML.

!!! info "Informazione"

    Forniamo anche una versione JSON del file di parametri come esempio ma non la eseguiremo qui.
    Sentitevi liberi di provare quella da soli.

#### 1.3.2. Esegui la pipeline

Per eseguire il workflow con questo file di parametri, aggiungete semplicemente `-params-file <nomefile>` al comando base.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

Il file di output finale dovrebbe contenere il personaggio stegosaurus che dice i saluti.

??? abstract "Contenuto del file"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Usare un file di parametri potrebbe sembrare eccessivo quando avete solo pochi parametri da specificare, ma alcune pipeline si aspettano dozzine di parametri.
In quei casi, usare un file di parametri vi permetterà di fornire valori di parametri a runtime senza dover digitare linee di comando massive e senza modificare lo script del workflow.

Rende anche più facile distribuire set di parametri ai collaboratori, o come informazione di supporto per una pubblicazione, per esempio.
Questo rende il vostro lavoro più riproducibile da altri.

### Riepilogo

Sapete come sfruttare le opzioni di configurazione chiave per gestire gli input del workflow.

### Cosa c'è dopo?

Imparate come gestire dove e come vengono pubblicati gli output del vostro workflow.

---

## 2. Gestire gli output del workflow

??? example "Scenario"

    La vostra pipeline pubblica gli output in una directory hardcoded, ma volete organizzare i risultati per progetto o nome dell'esperimento senza modificare il codice del workflow ogni volta.

Il workflow che abbiamo ereditato usa percorsi per le dichiarazioni di output a livello di workflow, il che non è terribilmente flessibile e comporta molta ripetizione.

Guardiamo alcuni modi comuni in cui potreste configurare questo per essere più flessibile.

### 2.1. Personalizza il nome della directory `outputDir`

Ogni versione del workflow che abbiamo eseguito finora ha pubblicato i suoi output in una sottodirectory diversa hardcoded nelle definizioni di output.

Abbiamo cambiato dove quella sottodirectory si trovava nella Parte 1 usando il flag CLI `-output-dir`, ma quella è ancora solo una stringa statica.
Configuriamo invece questo in un file di config, dove possiamo definire percorsi dinamici più complessi.
Potremmo creare un parametro completamente nuovo per questo, ma usiamo il parametro `batch` dato che è già lì.

#### 2.1.1. Imposta un valore per `outputDir` nel file di configurazione

Il percorso che Nextflow usa per pubblicare gli output è controllato dall'opzione `outputDir`.
Per cambiare il percorso per tutti gli output, potete impostare un valore per questa opzione nel file di configurazione `nextflow.config`.

Aggiungete il seguente codice al file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Questo sostituirà il percorso predefinito integrato, `results/`, con `results_config/` più il valore del parametro `batch` come sottodirectory.

Ricordate che potete anche impostare questa opzione dalla riga di comando usando il parametro `-output-dir` nel vostro comando (`-o` in forma breve), ma allora non potreste usare il valore del parametro `batch`.
Usare il flag CLI sovrascriverà `outputDir` nel config se è impostato.

#### 2.1.2. Rimuovi la parte ripetuta del percorso hardcoded

Abbiamo ancora una sottodirectory hardcoded nelle opzioni di output, quindi togliamola ora.

Fate le seguenti modifiche al codice nel file del workflow:

=== "Dopo"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Prima"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

Avremmo anche potuto semplicemente aggiungere `${params.batch}` a ogni percorso invece di modificare il predefinito di `outputDir`, ma questo è più conciso.

#### 2.1.3. Esegui la pipeline

Testiamo che funzioni correttamente, impostando il nome del batch a `outdir` dalla riga di comando.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results_config/outdir/`.

??? abstract "Contenuti della directory"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Potete combinare questo approccio con definizioni di percorsi personalizzate per costruire qualsiasi gerarchia di directory vi piaccia.

### 2.2. Organizza gli output per process

Un modo popolare per organizzare ulteriormente gli output è farlo per process, _cioè_ creare sottodirectory per ogni process eseguito nella pipeline.

#### 2.2.1. Sostituisci i percorsi di output con un riferimento ai nomi dei process

Tutto ciò che dovete fare è fare riferimento al nome del process come `<process>.name` nella dichiarazione del percorso di output.

Fate le seguenti modifiche nel file del workflow:

=== "Dopo"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Prima"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Questo rimuove gli elementi hardcoded rimanenti dalla configurazione del percorso di output.

#### 2.2.2. Esegui la pipeline

Testiamo che funzioni correttamente, impostando il nome del batch a `pnames` dalla riga di comando.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results_config/pnames/`, e sono raggruppati per process.

??? abstract "Contenuti della directory"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note "Nota"

    Notate che qui abbiamo eliminato la distinzione tra `intermediates` rispetto agli output finali al livello superiore.
    Potete mescolare e abbinare questi approcci e includere anche variabili multiple, per esempio impostando il percorso del primo output come `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Imposta la modalità di pubblicazione a livello di workflow

Infine, nello spirito di ridurre la quantità di codice ripetitivo, possiamo sostituire le dichiarazioni `mode` per-output con una singola riga nella configurazione.

#### 2.3.1. Aggiungi `workflow.output.mode` al file di configurazione

Aggiungete il seguente codice al file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

Proprio come l'opzione `outputDir`, dare a `workflow.output.mode` un valore nel file di configurazione sarebbe sufficiente per sovrascrivere ciò che è impostato nel file del workflow, ma rimuoviamo comunque il codice non necessario.

#### 2.3.2. Rimuovi la modalità di output dal file del workflow

Fate le seguenti modifiche nel file del workflow:

=== "Dopo"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Prima"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

È più conciso, vero?

#### 2.3.3. Esegui la pipeline

Testiamo che funzioni correttamente, impostando il nome del batch a `outmode` dalla riga di comando.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results_config/outmode/`.
Sono ancora tutte copie vere, non symlink.

??? abstract "Contenuti della directory"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Il motivo principale per cui potreste ancora voler usare il modo per-output di impostare la modalità è se volete mescolare e abbinare all'interno dello stesso workflow, _cioè_ avere alcuni output copiati e alcuni linkati simbolicamente.

Ci sono molte altre opzioni che potete personalizzare in questo modo, ma speriamo che questo vi dia un'idea della gamma di opzioni e di come utilizzarle efficacemente per adattarle alle vostre preferenze.

### Riepilogo

Sapete come controllare la denominazione e la struttura delle directory dove vengono pubblicati i vostri output, così come la modalità di pubblicazione dell'output del workflow.

### Cosa c'è dopo?

Imparate come adattare la configurazione del vostro workflow al vostro ambiente di calcolo, partendo dalla tecnologia di pacchettizzazione software.

---

## 3. Seleziona una tecnologia di pacchettizzazione software

Finora abbiamo guardato elementi di configurazione che controllano come gli input entrano e dove gli output escono. Ora è il momento di concentrarsi più specificamente sull'adattare la configurazione del vostro workflow al vostro ambiente di calcolo.

Il primo passo su quel percorso è specificare da dove verranno i pacchetti software che verranno eseguiti in ogni step.
Sono già installati nell'ambiente di calcolo locale?
Dobbiamo recuperare immagini ed eseguirle tramite un sistema container?
O dobbiamo recuperare pacchetti Conda e costruire un ambiente Conda locale?

Nella prima parte di questo corso di formazione (Parti 1-4) abbiamo semplicemente usato software installato localmente nel nostro workflow.
Poi nella Parte 5, abbiamo introdotto i container Docker e il file `nextflow.config`, che abbiamo usato per abilitare l'uso dei container Docker.

Ora vediamo come possiamo configurare un'opzione alternativa di pacchettizzazione software tramite il file `nextflow.config`.

### 3.1. Disabilita Docker e abilita Conda nel file di config

??? example "Scenario"

    State spostando la vostra pipeline su un cluster HPC dove Docker non è permesso per motivi di sicurezza.
    Il cluster supporta Singularity e Conda, quindi dovete cambiare la vostra configurazione di conseguenza.

Come notato precedentemente, Nextflow supporta molteplici tecnologie container incluso Singularity (che è più ampiamente usato su HPC), così come gestori di pacchetti software come Conda.

Possiamo cambiare il nostro file di configurazione per usare Conda invece di Docker.
Per farlo, cambiamo il valore di `docker.enabled` a `false`, e aggiungiamo una direttiva che abilita l'uso di Conda:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Questo permetterà a Nextflow di creare e utilizzare ambienti Conda per i process che hanno pacchetti Conda specificati.
Il che significa che ora dobbiamo aggiungerne uno al nostro process `cowpy`!

### 3.2. Specifica un pacchetto Conda nella definizione del process

Abbiamo già recuperato l'URI per un pacchetto Conda contenente lo strumento `cowpy`: `conda-forge::cowpy==1.1.5`

Ora aggiungiamo l'URI alla definizione del process `cowpy` usando la direttiva `conda`:

=== "Dopo"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Prima"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Per essere chiari, non stiamo _sostituendo_ la direttiva `docker`, stiamo _aggiungendo_ un'opzione alternativa.

!!! tip "Suggerimento"

    Ci sono alcuni modi diversi per ottenere l'URI per un dato pacchetto conda.
    Raccomandiamo di usare la query di ricerca di [Seqera Containers](https://seqera.io/containers/), che vi darà un URI che potete copiare e incollare, anche se non state pianificando di creare un container da esso.

### 3.3. Esegui il workflow per verificare che possa usare Conda

Proviamolo.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Output del comando"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Questo dovrebbe funzionare senza problemi e produrre gli stessi output di prima sotto `results_config/conda`.

Dietro le quinte, Nextflow ha recuperato i pacchetti Conda e creato l'ambiente, che normalmente richiede un po' di lavoro; quindi è bello che non dobbiamo fare niente di tutto ciò noi stessi!

!!! info "Informazione"

    Questo viene eseguito velocemente perché il pacchetto `cowpy` è abbastanza piccolo, ma se state lavorando con pacchetti grandi, potrebbe richiedere un po' più di tempo del solito la prima volta, e potreste vedere l'output della console rimanere 'bloccato' per un minuto circa prima di completarsi.
    Questo è normale ed è dovuto al lavoro extra che Nextflow fa la prima volta che usate un nuovo pacchetto.

Dal nostro punto di vista, sembra che funzioni esattamente come l'esecuzione con Docker, anche se nel backend la meccanica è un po' diversa.

Questo significa che siamo pronti per eseguire con ambienti Conda se necessario.

??? info "Mescolare e abbinare Docker e Conda"

    Dato che queste direttive sono assegnate per process, è possibile 'mescolare e abbinare', _cioè_ configurare alcuni dei process nel vostro workflow per essere eseguiti con Docker e altri con Conda, per esempio, se l'infrastruttura di calcolo che state usando supporta entrambi.
    In quel caso, abilitereste sia Docker che Conda nel vostro file di configurazione.
    Se entrambi sono disponibili per un dato process, Nextflow darà priorità ai container.

    E come notato prima, Nextflow supporta molteplici altre tecnologie di pacchettizzazione software e container, quindi non siete limitati a solo quelle due.

### Riepilogo

Sapete come configurare quale pacchetto software ogni process dovrebbe usare, e come passare da una tecnologia all'altra.

### Cosa c'è dopo?

Imparate come cambiare la piattaforma di esecuzione usata da Nextflow per fare effettivamente il lavoro.

---

## 4. Seleziona una piattaforma di esecuzione

??? example "Scenario"

    Avete sviluppato e testato la vostra pipeline sul vostro laptop, ma ora dovete eseguirla su migliaia di campioni.
    La vostra istituzione ha un cluster HPC con uno scheduler Slurm che vorreste usare invece.

Finora, abbiamo eseguito la nostra pipeline con l'executor locale.
Questo esegue ogni attività sulla macchina su cui Nextflow è in esecuzione.
Quando Nextflow inizia, guarda le CPU e la memoria disponibili.
Se le risorse delle attività pronte per l'esecuzione superano le risorse disponibili, Nextflow tratterrà le ultime attività dall'esecuzione fino a quando una o più delle attività precedenti sono terminate, liberando le risorse necessarie.

L'executor locale è conveniente ed efficiente, ma è limitato a quella singola macchina. Per carichi di lavoro molto grandi, potreste scoprire che la vostra macchina locale è un collo di bottiglia, o perché avete una singola attività che richiede più risorse di quelle disponibili, o perché avete così tante attività che aspettare che una singola macchina le esegua richiederebbe troppo tempo.

Nextflow supporta [molti backend di esecuzione diversi](https://nextflow.io/docs/latest/executor.html), inclusi scheduler HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor e altri) così come backend di esecuzione cloud come (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e altri).

### 4.1. Puntare a un backend diverso

La scelta dell'executor è impostata da una direttiva del process chiamata `executor`.
Per impostazione predefinita è impostata su `local`, quindi la seguente configurazione è implicita:

```groovy title="Configurazione integrata"
process {
    executor = 'local'
}
```

Per impostare l'executor per puntare a un backend diverso, specifichereste semplicemente l'executor che volete usando una sintassi simile a quella descritta sopra per le allocazioni di risorse (vedi [Executors](https://nextflow.io/docs/latest/executor.html) per tutte le opzioni).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Avviso"

    Non possiamo effettivamente testare questo nell'ambiente di formazione perché non è configurato per connettersi a un HPC.

### 4.2. Gestire la sintassi specifica del backend per i parametri di esecuzione

La maggior parte delle piattaforme di calcolo ad alte prestazioni permettono (e a volte richiedono) di specificare certi parametri come richieste e limitazioni di allocazione delle risorse (per es. numero di CPU e memoria) e nome della coda di job da usare.

Sfortunatamente, ciascuno di questi sistemi usa tecnologie, sintassi e configurazioni diverse per definire come un job dovrebbe essere definito e sottomesso allo scheduler rilevante.

??? abstract "Esempi"

    Per esempio, lo stesso job che richiede 8 CPU e 4GB di RAM per essere eseguito sulla coda "my-science-work" deve essere espresso in modi diversi a seconda del backend.

    ```bash title="Config per SLURM / sottometti usando sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config per PBS / sottometti usando qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config per SGE / sottometti usando qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Fortunatamente, Nextflow semplifica tutto questo.
Fornisce una sintassi standardizzata in modo che possiate specificare le proprietà rilevanti come `cpus`, `memory` e `queue` solo una volta (vedi [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) per tutte le opzioni disponibili).
Poi, a runtime, Nextflow userà quelle impostazioni per generare gli script specifici del backend appropriati basati sull'impostazione dell'executor.

Copriremo quella sintassi standardizzata nella prossima sezione.

### Riepilogo

Ora sapete come cambiare l'executor per usare diversi tipi di infrastruttura di calcolo.

### Cosa c'è dopo?

Imparate come valutare ed esprimere allocazioni e limitazioni delle risorse in Nextflow.

---

## 5. Controlla le allocazioni delle risorse di calcolo

??? example "Scenario"

    La vostra pipeline continua a fallire sul cluster perché le attività vengono terminate per aver superato i limiti di memoria.
    O forse vi vengono addebitate risorse che non state usando e volete ottimizzare i costi.

La maggior parte delle piattaforme di calcolo ad alte prestazioni permettono (e a volte richiedono) di specificare certi parametri di allocazione delle risorse come il numero di CPU e la memoria.

Per impostazione predefinita, Nextflow userà una singola CPU e 2GB di memoria per ogni process.
Le corrispondenti direttive del process sono chiamate `cpus` e `memory`, quindi la seguente configurazione è implicita:

```groovy title="Configurazione integrata" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Potete modificare questi valori, sia per tutti i process che per process specifici nominati, usando direttive di process aggiuntive nel vostro file di configurazione.
Nextflow le tradurrà nelle istruzioni appropriate per l'executor scelto.

Ma come sapete quali valori usare?

### 5.1. Esegui il workflow per generare un report di utilizzo delle risorse

??? example "Scenario"

    Non sapete quanta memoria o CPU i vostri process hanno bisogno e volete evitare di sprecare risorse o avere job terminati.

Se non sapete in anticipo quanta CPU e memoria i vostri process probabilmente avranno bisogno, potete fare del profiling delle risorse, il che significa eseguire il workflow con alcune allocazioni predefinite, registrare quanto ogni process ha usato, e da lì, stimare come regolare le allocazioni base.

Convenientemente, Nextflow include strumenti integrati per fare questo, e genererà felicemente un report per voi su richiesta.

Per farlo, aggiungete `-with-report <nomefile>.html` alla vostra riga di comando.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Il report è un file html, che potete scaricare e aprire nel vostro browser. Potete anche fare clic destro su di esso nell'esploratore file a sinistra e cliccare su `Show preview` per visualizzarlo nell'ambiente di formazione.

Prendetevi qualche minuto per guardare il report e vedere se riuscite a identificare alcune opportunità per regolare le risorse.
Assicuratevi di cliccare sulle schede che mostrano i risultati di utilizzo come percentuale di ciò che è stato allocato.

Vedi [Reports](https://nextflow.io/docs/latest/reports.html) per la documentazione su tutte le funzionalità disponibili.

### 5.2. Imposta le allocazioni delle risorse per tutti i process

Il profiling mostra che i process nel nostro workflow di formazione sono molto leggeri, quindi riduciamo l'allocazione di memoria predefinita a 1GB per process.

Aggiungete quanto segue al vostro file `nextflow.config`, prima della sezione dei parametri della pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

Questo aiuterà a ridurre la quantità di calcolo che consumiamo.

### 5.3. Imposta le allocazioni delle risorse per un process specifico

Allo stesso tempo, faremo finta che il process `cowpy` richieda più risorse degli altri, solo per poter dimostrare come regolare le allocazioni per un process individuale.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Con questa configurazione, tutti i process richiederanno 1GB di memoria e una singola CPU (il predefinito implicito), tranne il process `cowpy`, che richiederà 2GB e 2 CPU.

!!! info "Informazione"

    Se avete una macchina con poche CPU e allocate un numero elevato per process, potreste vedere chiamate ai process mettersi in coda l'una dietro l'altra.
    Questo perché Nextflow assicura che non richiediamo più CPU di quelle disponibili.

### 5.4. Esegui il workflow con la configurazione aggiornata

Proviamolo, fornendo un nome file diverso per il report di profiling così possiamo confrontare le prestazioni prima e dopo i cambiamenti di configurazione.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Probabilmente non noterete alcuna differenza reale dato che questo è un carico di lavoro così piccolo, ma questo è l'approccio che usereste per analizzare le prestazioni e i requisiti di risorse di un workflow del mondo reale.

È molto utile quando i vostri process hanno requisiti di risorse diversi. Vi permette di dimensionare correttamente le allocazioni delle risorse che impostate per ogni process basandovi su dati reali, non su supposizioni.

!!! tip "Suggerimento"

    Questo è solo un piccolo assaggio di ciò che potete fare per ottimizzare il vostro uso delle risorse.
    Nextflow stesso ha una logica di [retry dinamico](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) davvero interessante integrata per ritentare i job che falliscono a causa di limitazioni delle risorse.
    Inoltre, Seqera Platform offre anche strumenti basati sull'IA per ottimizzare automaticamente le vostre allocazioni delle risorse.

### 5.5. Aggiungi limiti alle risorse

A seconda di quale executor di calcolo e infrastruttura di calcolo state usando, potrebbero esserci alcuni vincoli su ciò che potete (o dovete) allocare.
Per esempio, il vostro cluster potrebbe richiedere di rimanere entro certi limiti.

Potete usare la direttiva `resourceLimits` per impostare le limitazioni rilevanti. La sintassi appare così quando è da sola in un blocco process:

```groovy title="Esempio di sintassi"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow tradurrà questi valori nelle istruzioni appropriate a seconda dell'executor che avete specificato.

Non eseguiremo questo, dato che non abbiamo accesso a un'infrastruttura rilevante nell'ambiente di formazione.
Tuttavia, se provaste a eseguire il workflow con allocazioni di risorse che superano questi limiti, poi cercaste il comando `sbatch` nel file di script `.command.run`, vedreste che le richieste che vengono effettivamente inviate all'executor sono limitate ai valori specificati da `resourceLimits`.

??? info "Configurazioni di riferimento istituzionali"

    Il progetto nf-core ha compilato una [collezione di file di configurazione](https://nf-co.re/configs/) condivisi da varie istituzioni in tutto il mondo, coprendo una vasta gamma di executor HPC e cloud.

    Quei config condivisi sono preziosi sia per le persone che lavorano lì e possono quindi semplicemente utilizzare la configurazione della loro istituzione out of the box, sia come modello per le persone che cercano di sviluppare una configurazione per la propria infrastruttura.

### Riepilogo

Sapete come generare un report di profiling per valutare l'utilizzo delle risorse e come modificare le allocazioni delle risorse per tutti i process e/o per process individuali, oltre a impostare limitazioni delle risorse per l'esecuzione su HPC.

### Cosa c'è dopo?

Imparate come impostare profili di configurazione preimpostati e passare da uno all'altro a runtime.

---

## 6. Usa i profili per passare da una configurazione preimpostata all'altra

??? example "Scenario"

    Passate regolarmente dall'esecuzione di pipeline sul vostro laptop per lo sviluppo all'HPC della vostra istituzione per le esecuzioni in produzione.
    Siete stanchi di cambiare manualmente le impostazioni di configurazione ogni volta che cambiate ambiente.

Vi abbiamo mostrato diversi modi in cui potete personalizzare la configurazione della vostra pipeline a seconda del progetto su cui state lavorando o dell'ambiente di calcolo che state usando.

Potreste voler passare da un'impostazione alternativa all'altra a seconda di quale infrastruttura di calcolo state usando. Per esempio, potreste voler sviluppare e eseguire test su piccola scala localmente sul vostro laptop, poi eseguire carichi di lavoro su scala completa su HPC o cloud.

Nextflow vi permette di impostare un numero qualsiasi di [**profili**](https://nextflow.io/docs/latest/config.html#profiles) che descrivono diverse configurazioni, che potete poi selezionare a runtime usando un argomento da riga di comando, piuttosto che dover modificare il file di configurazione stesso.

### 6.1. Crea profili per passare dallo sviluppo locale all'esecuzione su HPC

Impostiamo due profili alternativi; uno per eseguire carichi su piccola scala su un computer normale, dove useremo container Docker, e uno per l'esecuzione su un HPC universitario con uno scheduler Slurm, dove useremo pacchetti Conda.

#### 6.1.1. Imposta i profili

Aggiungete quanto segue al vostro file `nextflow.config`, dopo la sezione dei parametri della pipeline ma prima delle impostazioni di output:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

Vedete che per l'HPC universitario, stiamo anche specificando limitazioni delle risorse.

#### 6.1.2. Esegui il workflow con un profilo

Per specificare un profilo nella nostra riga di comando di Nextflow, usiamo l'argomento `-profile`.

Proviamo a eseguire il workflow con la configurazione `my_laptop`.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Come potete vedere, questo ci permette di alternare tra configurazioni molto comodamente a runtime.

!!! warning "Avviso"

    Il profilo `univ_hpc` non funzionerà correttamente nell'ambiente di formazione dato che non abbiamo accesso a uno scheduler Slurm.

Se in futuro troviamo altri elementi di configurazione che co-occorrono sempre con questi, possiamo semplicemente aggiungerli al/i profilo/i corrispondente/i.
Possiamo anche creare profili aggiuntivi se ci sono altri elementi di configurazione che vogliamo raggruppare insieme.

### 6.2. Crea un profilo di parametri di test

??? example "Scenario"

    Volete che altri possano provare rapidamente la vostra pipeline senza dover raccogliere i propri dati di input.

I profili non sono solo per la configurazione dell'infrastruttura.
Possiamo anche usarli per impostare valori predefiniti per i parametri del workflow, per rendere più facile per gli altri provare il workflow senza dover raccogliere valori di input appropriati loro stessi.
Potete considerare questo un'alternativa all'uso di un file di parametri.

#### 6.2.1. Imposta il profilo

La sintassi per esprimere valori predefiniti in questo contesto appare così, per un profilo che chiamiamo `test`:

```groovy title="Esempio di sintassi"
    test {
        params.<parametro1>
        params.<parametro2>
        ...
    }
```

Se aggiungiamo un profilo test per il nostro workflow, il blocco `profiles` diventa:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Proprio come per i profili di configurazione tecnica, potete impostare molteplici profili diversi che specificano parametri sotto qualsiasi nome arbitrario vi piaccia.

#### 6.2.2. Esegui il workflow localmente con il profilo test

Convenientemente, i profili non sono mutuamente esclusivi, quindi possiamo specificare profili multipli nella nostra riga di comando usando la seguente sintassi `-profile <profilo1>,<profilo2>` (per qualsiasi numero di profili).

Se combinate profili che impostano valori per gli stessi elementi di configurazione e sono descritti nello stesso file di configurazione, Nextflow risolverà il conflitto usando qualunque valore abbia letto per ultimo (_cioè_ qualunque cosa venga dopo nel file).
Se le impostazioni in conflitto sono impostate in diverse fonti di configurazione, si applica l'[ordine di precedenza](https://www.nextflow.io/docs/latest/config.html#configuration-file) predefinito.

Proviamo ad aggiungere il profilo test al nostro comando precedente:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Questo userà Docker dove possibile e produrrà output sotto `results_config/test`, e questa volta il personaggio è il duo comico `dragonandcow`.

??? abstract "Contenuto del file"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Questo significa che finché distribuiamo tutti i file di dati di test con il codice del workflow, chiunque può provare rapidamente il workflow senza dover fornire i propri input tramite la riga di comando o un file di parametri.

!!! tip "Suggerimento"

    Possiamo puntare a URL per file più grandi che sono memorizzati esternamente.
    Nextflow li scaricherà automaticamente finché c'è una connessione aperta.

    Per maggiori dettagli, vedi la Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Usa `nextflow config` per vedere la configurazione risolta

Come notato sopra, a volte lo stesso parametro può essere impostato a valori diversi in profili che volete combinare.
E più in generale, ci sono numerosi posti dove gli elementi di configurazione possono essere memorizzati, e a volte le stesse proprietà possono essere impostate a valori diversi in posti diversi.

Nextflow applica un [ordine di precedenza](https://www.nextflow.io/docs/latest/config.html#configuration-file) stabilito per risolvere qualsiasi conflitto, ma può essere difficile determinarlo da soli.
E anche se nulla è in conflitto, può essere noioso cercare tutti i possibili posti dove le cose potrebbero essere configurate.

Fortunatamente, Nextflow include un conveniente strumento utility chiamato `config` che può automatizzare tutto quel processo per voi.

Lo strumento `config` esplorerà tutti i contenuti nella vostra directory di lavoro corrente, raccoglierà tutti i file di configurazione, e produrrà la configurazione completamente risolta che Nextflow userebbe per eseguire il workflow.
Questo vi permette di scoprire quali impostazioni verranno usate senza dover lanciare nulla.

#### 6.3.1. Risolvi la configurazione predefinita

Eseguite questo comando per risolvere la configurazione che verrebbe applicata per impostazione predefinita.

```bash
nextflow config
```

??? success "Output del comando"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Questo vi mostra la configurazione base che ottenete se non specificate nulla di extra nella riga di comando.

#### 6.3.2. Risolvi la configurazione con impostazioni specifiche attivate

Se fornite parametri da riga di comando, es. abilitando uno o più profili o caricando un file di parametri, il comando terrà conto anche di quelli.

```bash
nextflow config -profile my_laptop,test
```

??? success "Output del comando"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Questo diventa particolarmente utile per progetti complessi che coinvolgono molteplici livelli di configurazione.

### Riepilogo

Sapete come usare i profili per selezionare una configurazione preimpostata a runtime con il minimo fastidio.
Più in generale, sapete come configurare le esecuzioni del vostro workflow per adattarle a diverse piattaforme di calcolo e migliorare la riproducibilità delle vostre analisi.

### Cosa c'è dopo?

Imparate come eseguire pipeline direttamente da repository remoti come GitHub.

---

## 7. Esegui pipeline da repository remoti

??? example "Scenario"

    Volete eseguire una pipeline ben consolidata come quelle di nf-core senza dover scaricare e gestire il codice voi stessi.

Finora abbiamo eseguito script di workflow situati nella directory corrente.
In pratica, vorrete spesso eseguire pipeline memorizzate in repository remoti, come GitHub.

Nextflow rende questo semplice: potete eseguire qualsiasi pipeline direttamente da un URL di un repository Git senza scaricarlo manualmente prima.

### 7.1. Esegui una pipeline da GitHub

La sintassi base per eseguire una pipeline remota è `nextflow run <repository>`, dove `<repository>` può essere un percorso di repository GitHub come `nextflow-io/hello`, un URL completo, o un percorso a GitLab, Bitbucket, o altri servizi di hosting Git.

Provate a eseguire la pipeline demo ufficiale di Nextflow "hello":

```bash
nextflow run nextflow-io/hello
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

La prima volta che eseguite una pipeline remota, Nextflow la scarica e la mette in cache localmente.
Le esecuzioni successive usano la versione in cache a meno che non richiediate esplicitamente un aggiornamento.

### 7.2. Specifica una versione per la riproducibilità

Per impostazione predefinita, Nextflow esegue l'ultima versione dal branch predefinito.
Potete specificare una particolare versione (tag), branch, o commit usando il flag `-r`:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Specificare versioni esatte è essenziale per la riproducibilità.

### Riepilogo

Sapete come eseguire pipeline direttamente da GitHub e altri repository remoti, e come specificare versioni per la riproducibilità.

### Cosa c'è dopo?

Datevi una bella pacca sulla spalla!
Sapete tutto ciò che dovete sapere per iniziare a eseguire e gestire pipeline Nextflow.

Questo conclude questo corso, ma se siete desiderosi di continuare a imparare, abbiamo due raccomandazioni principali:

- Se volete approfondire lo sviluppo delle vostre pipeline, date un'occhiata a [Hello Nextflow](../hello_nextflow/index.md), un corso per principianti che copre la stessa progressione generale di questo ma va molto più in dettaglio su channel e operatori.
- Se vorreste continuare a imparare come eseguire pipeline Nextflow senza andare più in profondità nel codice, date un'occhiata alla prima parte di [Hello nf-core](../hello_nf-core/index.md), che introduce gli strumenti per trovare e eseguire pipeline dall'estremamente popolare progetto [nf-core](https://nf-co.re/).

Buon divertimento!

---

## Quiz

<quiz>
Quando i valori dei parametri sono impostati sia nel file del workflow che in `nextflow.config`, quale ha la precedenza?
- [ ] Il valore del file del workflow
- [x] Il valore del file di configurazione
- [ ] Il primo valore incontrato
- [ ] Causa un errore

Approfondisci: [1.1. Imposta i valori in `nextflow.config`](#11-imposta-i-valori-in-nextflowconfig)
</quiz>

<quiz>
Qual è la differenza di sintassi tra l'impostazione di un valore predefinito di un parametro in un file workflow vs. un file config?
- [ ] Usano la stessa sintassi
- [x] Il workflow usa dichiarazione tipizzata (`#!groovy param: Type = value`), il config usa assegnazione (`#!groovy param = value`)
- [ ] Il config usa dichiarazione tipizzata, il workflow usa assegnazione
- [ ] Solo i file config possono impostare valori predefiniti

Approfondisci: [1.1. Imposta i valori in `nextflow.config`](#11-imposta-i-valori-in-nextflowconfig)
</quiz>

<quiz>
Come specificate un file di parametri quando eseguite un workflow?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Approfondisci: [1.3. Usa un file di parametri](#13-usa-un-file-di-parametri)
</quiz>

<quiz>
Cosa controlla l'opzione di configurazione `outputDir`?
- [ ] La posizione della directory di lavoro
- [x] Il percorso base dove vengono pubblicati gli output del workflow
- [ ] La directory per i file di log
- [ ] La posizione dei file modulo

Approfondisci: [2.1. Personalizza il nome della directory outputDir](#21-personalizza-il-nome-della-directory-outputdir)
</quiz>

<quiz>
Come si fa riferimento dinamicamente al nome di un process nella configurazione del percorso di output?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

Approfondisci: [2.2. Organizza gli output per process](#22-organizza-gli-output-per-process)
</quiz>

<quiz>
Se sia Docker che Conda sono abilitati e un process ha entrambe le direttive, quale ha la priorità?
- [x] Docker (container)
- [ ] Conda
- [ ] Il primo definito nel process
- [ ] Causa un errore

Approfondisci: [3. Seleziona una tecnologia di pacchettizzazione software](#3-seleziona-una-tecnologia-di-pacchettizzazione-software)
</quiz>

<quiz>
Qual è l'executor predefinito in Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Approfondisci: [4. Seleziona una piattaforma di esecuzione](#4-seleziona-una-piattaforma-di-esecuzione)
</quiz>

<quiz>
Quale comando genera un report di utilizzo delle risorse?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Approfondisci: [5.1. Esegui il workflow per generare un report di utilizzo delle risorse](#51-esegui-il-workflow-per-generare-un-report-di-utilizzo-delle-risorse)
</quiz>

<quiz>
Come si impostano i requisiti di risorse per un process specifico chiamato `cowpy` nel file config?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Approfondisci: [5.3. Imposta le allocazioni delle risorse per un process specifico](#53-imposta-le-allocazioni-delle-risorse-per-un-process-specifico)
</quiz>

<quiz>
Cosa fa la direttiva `resourceLimits`?
- [ ] Imposta i requisiti minimi di risorse
- [ ] Alloca risorse ai process
- [x] Limita le risorse massime che possono essere richieste
- [ ] Monitora l'utilizzo delle risorse in tempo reale

Approfondisci: [5.5. Aggiungi limiti alle risorse](#55-aggiungi-limiti-alle-risorse)
</quiz>

<quiz>
Come si specificano profili multipli in un singolo comando?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Approfondisci: [6. Usa i profili per passare da una configurazione preimpostata all'altra](#6-usa-i-profili-per-passare-da-una-configurazione-preimpostata-allaltra)
</quiz>

<quiz>
Quale comando mostra la configurazione completamente risolta che Nextflow userebbe?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Approfondisci: [6.3. Usa `nextflow config` per vedere la configurazione risolta](#63-usa-nextflow-config-per-vedere-la-configurazione-risolta)
</quiz>

<quiz>
Per cosa possono essere usati i profili? (Selezionate tutte le risposte corrette)
- [x] Definire impostazioni specifiche dell'infrastruttura (executor, container)
- [x] Impostare limiti di risorse per diversi ambienti
- [x] Fornire parametri di test per testare facilmente il workflow
- [ ] Definire nuovi process

Approfondisci: [6. Usa i profili per passare da una configurazione preimpostata all'altra](#6-usa-i-profili-per-passare-da-una-configurazione-preimpostata-allaltra)
</quiz>
