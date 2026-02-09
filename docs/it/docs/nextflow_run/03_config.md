# Parte 3: Configurazione dell'esecuzione

Questa sezione esplorerà come gestire la configurazione di una pipeline Nextflow per personalizzarne il comportamento, adattarla a diversi ambienti e ottimizzare l'uso delle risorse _senza modificare una sola riga del codice del workflow stesso_.

Ci sono diversi modi per farlo, che possono essere utilizzati in combinazione e vengono interpretati secondo l'ordine di precedenza descritto nella documentazione sulla [Configurazione](https://nextflow.io/docs/latest/config.html).

In questa parte del corso, vi mostreremo il meccanismo di configurazione più semplice e comune, il file `nextflow.config`, che avete già incontrato nella sezione sui container nella Parte 2.

Esamineremo componenti essenziali della configurazione di Nextflow come le direttive di processo, gli executor, i profili e i file di parametri.
Imparando a utilizzare efficacemente queste opzioni di configurazione, potrete sfruttare appieno la flessibilità, la scalabilità e le prestazioni delle pipeline Nextflow.

Per esercitarci con questi elementi di configurazione, eseguiremo una copia aggiornata del workflow che abbiamo eseguito alla fine della Parte 2 di questo corso di formazione, rinominato `3-main.nf`.

Se non conoscete la pipeline Hello o avete bisogno di un ripasso, consultate [questa pagina informativa](../info/hello_pipeline.md).

---

## 1. Gestire i parametri di input del workflow

??? example "Scenario"

    Avete scaricato una pipeline e volete eseguirla ripetutamente con gli stessi file di input e impostazioni, ma non volete digitare tutti i parametri ogni volta.
    O forse state configurando la pipeline per un collega che non ha dimestichezza con gli argomenti da riga di comando.

Inizieremo con un aspetto della configurazione che è semplicemente un'estensione di ciò con cui abbiamo lavorato finora: la gestione dei parametri di input.

Attualmente, il nostro workflow è configurato per accettare diversi valori di parametri tramite la riga di comando, dichiarati in un blocco `params` nello script del workflow stesso.
Uno ha un valore predefinito impostato come parte della sua dichiarazione.

Tuttavia, potreste voler impostare valori predefiniti per tutti, o sovrascrivere il valore predefinito esistente senza dover specificare parametri dalla riga di comando o modificare il file di script originale.

Ci sono diversi modi per farlo; vi mostreremo tre modi di base molto comunemente utilizzati.

### 1.1. Impostare i valori in `nextflow.config`

Questo è l'approccio più semplice, anche se forse il meno flessibile poiché il file principale `nextflow.config` non è qualcosa che vorreste modificare per ogni esecuzione.
Ma ha il vantaggio di separare le preoccupazioni di _dichiarare_ i parametri nel workflow (che appartiene sicuramente lì) rispetto a fornire _valori predefiniti_, che sono più a loro agio in un file di configurazione.

Facciamolo in due passaggi.

#### 1.1.1. Creare un blocco `params` nel file di configurazione

Apportate le seguenti modifiche al codice nel file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Parametri della pipeline
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
Per il parametro `batch` che aveva già un valore predefinito dichiarato, la sintassi è leggermente diversa.
Nel file del workflow, quella è una dichiarazione tipizzata.
Nella configurazione, sono assegnazioni di valori.

Tecnicamente, questo è sufficiente per sovrascrivere i valori predefiniti ancora specificati nel file del workflow.
Potreste modificare il valore predefinito per `batch` ed eseguire il workflow per verificare che il valore impostato nel file di configurazione sovrascriva quello impostato nel file del workflow.

Ma nello spirito di spostare completamente la configurazione nel file di configurazione, rimuoviamo del tutto quel valore predefinito dal file del workflow.

#### 1.1.2. Rimuovere il valore predefinito per `batch` nel file del workflow

Apportate la seguente modifica al codice nel file del workflow `3-main.nf`:

=== "Dopo"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Parametri della pipeline
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
    * Parametri della pipeline
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Ora il file del workflow stesso non imposta alcun valore predefinito per questi parametri.

#### 1.1.3. Eseguire la pipeline

Verifichiamo che funzioni correttamente senza specificare alcun parametro nella riga di comando.

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

L'output finale in ASCII art si trova nella directory `results/3-main/`, con il nome `cowpy-COLLECTED-batch-output.txt`, come prima.

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

Funzionalmente, questa modifica non ha cambiato nulla, ma concettualmente è un po' più pulito avere i valori predefiniti impostati nel file di configurazione.

### 1.2. Utilizzare un file di configurazione specifico per l'esecuzione

??? example "Scenario"

    Volete sperimentare con diverse impostazioni senza modificare il vostro file di configurazione principale.

Potete farlo creando un nuovo file `nextflow.config` in una sottodirectory che userete come directory di lavoro per i vostri esperimenti.

#### 1.2.1. Creare la directory di lavoro con una configurazione vuota

Iniziamo creando una nuova directory e spostandoci al suo interno:

```bash
mkdir -p tux-run
cd tux-run
```

Quindi, create un file di configurazione vuoto in quella directory:

```bash
touch nextflow.config
```

Questo produce un file vuoto.

#### 1.2.2. Configurare la configurazione sperimentale

Ora aprite il nuovo file e aggiungete i parametri che volete personalizzare:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Notate che il percorso del file di input deve riflettere la struttura delle directory.

#### 1.2.3. Eseguire la pipeline

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

Questo creerà un nuovo insieme di directory sotto `tux-run/` incluse `tux-run/work/` e `tux-run/results/`.

In questa esecuzione, Nextflow combina il `nextflow.config` nella nostra directory corrente con il `nextflow.config` nella directory radice della pipeline, e quindi sovrascrive il carattere predefinito (turkey) con il carattere tux.

Il file di output finale dovrebbe contenere il carattere tux che dice i saluti.

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

!!! warning

    Assicuratevi di tornare alla directory precedente prima di passare alla sezione successiva!

    ```bash
    cd ..
    ```

Ora vediamo un altro modo utile per impostare i valori dei parametri.

### 1.3. Utilizzare un file di parametri

??? example "Scenario"

    Dovete condividere parametri di esecuzione esatti con un collaboratore, o registrarli per una pubblicazione.

L'approccio della sottodirectory funziona benissimo per sperimentare, ma richiede un po' di configurazione e richiede di adattare i percorsi di conseguenza.
C'è un approccio più semplice per quando volete eseguire la vostra pipeline con un insieme specifico di valori, o consentire a qualcun altro di farlo con il minimo sforzo.

Nextflow ci consente di specificare parametri tramite un [file di parametri](https://nextflow.io/docs/latest/config.html#parameter-file) in formato YAML o JSON, il che rende molto conveniente gestire e distribuire insiemi alternativi di valori predefiniti, ad esempio, così come valori di parametri specifici per l'esecuzione.

#### 1.3.1. Esaminare il file di parametri di esempio

Per dimostrarlo, forniamo un file di parametri di esempio nella directory corrente, chiamato `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Questo file di parametri contiene una coppia chiave-valore per ciascuno degli input che vogliamo specificare.
Notate l'uso dei due punti (`:`) invece dei segni di uguale (`=`) se confrontate la sintassi con il file di configurazione.
Il file di configurazione è scritto in Groovy, mentre il file di parametri è scritto in YAML.

!!! info

    Forniamo anche una versione JSON del file di parametri come esempio ma non la eseguiremo qui.
    Sentitevi liberi di provarla da soli.

#### 1.3.2. Eseguire la pipeline

Per eseguire il workflow con questo file di parametri, aggiungete semplicemente `-params-file <filename>` al comando base.

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

Il file di output finale dovrebbe contenere il carattere stegosauro che dice i saluti.

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

Utilizzare un file di parametri può sembrare eccessivo quando avete solo pochi parametri da specificare, ma alcune pipeline si aspettano dozzine di parametri.
In quei casi, utilizzare un file di parametri ci consentirà di fornire valori di parametri al momento dell'esecuzione senza dover digitare righe di comando enormi e senza modificare lo script del workflow.

Rende anche più facile distribuire insiemi di parametri ai collaboratori, o come informazioni di supporto per una pubblicazione, ad esempio.
Questo rende il vostro lavoro più riproducibile da altri.

### Takeaway

Sapete come sfruttare le opzioni di configurazione chiave per gestire gli input del workflow.

### Cosa c'è dopo?

Imparate come gestire dove e come vengono pubblicati gli output del vostro workflow.

---

## 2. Gestire gli output del workflow

??? example "Scenario"

    La vostra pipeline pubblica gli output in una directory hardcoded, ma volete organizzare i risultati per progetto o nome dell'esperimento senza modificare il codice del workflow ogni volta.

Il workflow che abbiamo ereditato utilizza percorsi per le dichiarazioni di output a livello di workflow, il che non è terribilmente flessibile e comporta molta ripetizione.

Vediamo alcuni modi comuni in cui potreste configurarlo per essere più flessibile.

### 2.1. Personalizzare il nome della directory `outputDir`

Ogni versione del workflow che abbiamo eseguito finora ha pubblicato i suoi output in una sottodirectory diversa hardcoded nelle definizioni di output.

Abbiamo cambiato dove si trovava quella sottodirectory nella Parte 1 utilizzando il flag CLI `-output-dir`, ma quella è ancora solo una stringa statica.
Configuriamo invece questo in un file di configurazione, dove possiamo definire percorsi dinamici più complessi.
Potremmo creare un parametro completamente nuovo per questo, ma utilizziamo il parametro `batch` poiché è proprio lì.

#### 2.1.1. Impostare un valore per `outputDir` nel file di configurazione

Il percorso che Nextflow utilizza per pubblicare gli output è controllato dall'opzione `outputDir`.
Per cambiare il percorso per tutti gli output, potete impostare un valore per questa opzione nel file di configurazione `nextflow.config`.

Aggiungete il seguente codice al file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Parametri della pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Impostazioni di output
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Parametri della pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Questo sostituirà il percorso predefinito integrato, `results/`, con `results_config/` più il valore del parametro `batch` come sottodirectory.

Ricordate che potete anche impostare questa opzione dalla riga di comando utilizzando il parametro `-output-dir` nel vostro comando (`-o` in breve), ma in quel caso non potreste utilizzare il valore del parametro `batch`.
L'utilizzo del flag CLI sovrascriverà `outputDir` nella configurazione se è impostato.

#### 2.1.2. Rimuovere la parte ripetuta del percorso hardcoded

Abbiamo ancora una sottodirectory hardcoded nelle opzioni di output, quindi eliminiamola ora.

Apportate le seguenti modifiche al codice nel file del workflow:

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

Avremmo anche potuto semplicemente aggiungere `${params.batch}` a ciascun percorso invece di modificare il valore predefinito di `outputDir`, ma questo è più conciso.

#### 2.1.3. Eseguire la pipeline

Verifichiamo che funzioni correttamente, impostando il nome del batch su `outdir` dalla riga di comando.

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

??? abstract "Contenuto della directory"

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

Potete combinare questo approccio con definizioni di percorsi personalizzati per costruire qualsiasi gerarchia di directory desideriate.

### 2.2. Organizzare gli output per processo

Un modo popolare per organizzare ulteriormente gli output è farlo per processo, _cioè_ creare sottodirectory per ogni processo eseguito nella pipeline.

#### 2.2.1. Sostituire i percorsi di output con un riferimento ai nomi dei processi

Tutto ciò che dovete fare è fare riferimento al nome del processo come `<processo>.name` nella dichiarazione del percorso di output.

Apportate le seguenti modifiche nel file del workflow:

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

#### 2.2.2. Eseguire la pipeline

Verifichiamo che funzioni correttamente, impostando il nome del batch su `pnames` dalla riga di comando.

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

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results_config/pnames/`, e sono raggruppati per processo.

??? abstract "Contenuto della directory"

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

!!! note

    Notate che qui abbiamo cancellato la distinzione tra `intermediates` e output finali al livello superiore.
    Potete mescolare e abbinare questi approcci e persino includere più variabili, ad esempio impostando il percorso del primo output come `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Impostare la modalità di pubblicazione a livello di workflow

Infine, nello spirito di ridurre la quantità di codice ripetitivo, possiamo sostituire le dichiarazioni `mode` per output con una singola riga nella configurazione.

#### 2.3.1. Aggiungere `workflow.output.mode` al file di configurazione

Aggiungete il seguente codice al file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Impostazioni di output
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Impostazioni di output
    */
    outputDir = "results_config/${params.batch}"
    ```

Proprio come l'opzione `outputDir`, dare a `workflow.output.mode` un valore nel file di configurazione sarebbe sufficiente per sovrascrivere ciò che è impostato nel file del workflow, ma rimuoviamo comunque il codice non necessario.

#### 2.3.2. Rimuovere la modalità di output dal file del workflow

Apportate le seguenti modifiche nel file del workflow:

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

#### 2.3.3. Eseguire la pipeline

Verifichiamo che funzioni correttamente, impostando il nome del batch su `outmode` dalla riga di comando.

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
Sono ancora tutte copie corrette, non symlink.

??? abstract "Contenuto della directory"

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

Il motivo principale per cui potreste ancora voler utilizzare il modo per output di impostare la modalità è se volete mescolare e abbinare all'interno dello stesso workflow, _cioè_ avere alcuni output copiati e alcuni collegati tramite symlink.

Ci sono molte altre opzioni che potete personalizzare in questo modo, ma si spera che questo vi dia un'idea della gamma di opzioni e di come utilizzarle efficacemente per adattarle alle vostre preferenze.

### Takeaway

Sapete come controllare la denominazione e la struttura delle directory in cui vengono pubblicati i vostri output, nonché la modalità di pubblicazione dell'output del workflow.

### Cosa c'è dopo?

Imparate come adattare la configurazione del vostro workflow al vostro ambiente di calcolo, iniziando con la tecnologia di packaging del software.

---

## 3. Selezionare una tecnologia di packaging del software

Finora abbiamo esaminato elementi di configurazione che controllano come entrano gli input e dove escono gli output. Ora è il momento di concentrarci più specificamente sull'adattamento della configurazione del vostro workflow al vostro ambiente di calcolo.

Il primo passo su quel percorso è specificare da dove verranno i pacchetti software che verranno eseguiti in ogni passaggio.
Sono già installati nell'ambiente di calcolo locale?
Dobbiamo recuperare immagini ed eseguirle tramite un sistema di container?
O dobbiamo recuperare pacchetti Conda e costruire un ambiente Conda locale?

Nella primissima parte di questo corso di formazione (Parti 1-4) abbiamo semplicemente utilizzato software installato localmente nel nostro workflow.
Poi nella Parte 5, abbiamo introdotto i container Docker e il file `nextflow.config`, che abbiamo utilizzato per abilitare l'uso dei container Docker.

Ora vediamo come possiamo configurare un'opzione di packaging software alternativa tramite il file `nextflow.config`.

### 3.1. Disabilitare Docker e abilitare Conda nel file di configurazione

??? example "Scenario"

    State spostando la vostra pipeline su un cluster HPC dove Docker non è consentito per motivi di sicurezza.
    Il cluster supporta Singularity e Conda, quindi dovete cambiare la vostra configurazione di conseguenza.

Come notato in precedenza, Nextflow supporta più tecnologie di container tra cui Singularity (che è più ampiamente utilizzato su HPC), così come gestori di pacchetti software come Conda.

Possiamo modificare il nostro file di configurazione per utilizzare Conda invece di Docker.
Per farlo, cambiamo il valore di `docker.enabled` in `false` e aggiungiamo una direttiva che abilita l'uso di Conda:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Questo consentirà a Nextflow di creare e utilizzare ambienti Conda per i processi che hanno pacchetti Conda specificati.
Il che significa che ora dobbiamo aggiungerne uno al nostro processo `cowpy`!

### 3.2. Specificare un pacchetto Conda nella definizione del processo

Abbiamo già recuperato l'URI per un pacchetto Conda contenente lo strumento `cowpy`: `conda-forge::cowpy==1.1.5`

Ora aggiungiamo l'URI alla definizione del processo `cowpy` utilizzando la direttiva `conda`:

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

!!! tip

    Ci sono alcuni modi diversi per ottenere l'URI per un dato pacchetto conda.
    Consigliamo di utilizzare la query di ricerca [Seqera Containers](https://seqera.io/containers/), che vi darà un URI che potete copiare e incollare, anche se non state pianificando di creare un container da esso.

### 3.3. Eseguire il workflow per verificare che possa utilizzare Conda

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

Dietro le quinte, Nextflow ha recuperato i pacchetti Conda e creato l'ambiente, il che normalmente richiede un po' di lavoro; quindi è bello che non dobbiamo fare nulla di tutto ciò noi stessi!

!!! info

    Questo viene eseguito rapidamente perché il pacchetto `cowpy` è piuttosto piccolo, ma se state lavorando con pacchetti di grandi dimensioni, potrebbe richiedere un po' più di tempo del solito la prima volta, e potreste vedere l'output della console rimanere 'bloccato' per un minuto circa prima di completarsi.
    Questo è normale ed è dovuto al lavoro extra che Nextflow fa la prima volta che utilizzate un nuovo pacchetto.

Dal nostro punto di vista, sembra che funzioni esattamente come l'esecuzione con Docker, anche se sul backend la meccanica è un po' diversa.

Questo significa che siamo tutti pronti per eseguire con ambienti Conda se necessario.

??? info "Mescolare e abbinare Docker e Conda"

    Poiché queste direttive sono assegnate per processo, è possibile 'mescolare e abbinare', _cioè_ configurare alcuni dei processi nel vostro workflow per essere eseguiti con Docker e altri con Conda, ad esempio, se l'infrastruttura di calcolo che state utilizzando supporta entrambi.
    In tal caso, abilitereste sia Docker che Conda nel vostro file di configurazione.
    Se entrambi sono disponibili per un dato processo, Nextflow darà priorità ai container.

    E come notato in precedenza, Nextflow supporta più altre tecnologie di packaging software e container, quindi non siete limitati solo a queste due.

### Takeaway

Sapete come configurare quale pacchetto software ogni processo dovrebbe utilizzare e come passare da una tecnologia all'altra.

### Cosa c'è dopo?

Imparate come cambiare la piattaforma di esecuzione utilizzata da Nextflow per fare effettivamente il lavoro.

---

## 4. Selezionare una piattaforma di esecuzione

??? example "Scenario"

    Avete sviluppato e testato la vostra pipeline sul vostro laptop, ma ora dovete eseguirla su migliaia di campioni.
    La vostra istituzione ha un cluster HPC con uno scheduler Slurm che vorreste utilizzare invece.

Fino ad ora, abbiamo eseguito la nostra pipeline con l'executor locale.
Questo esegue ogni attività sulla macchina su cui è in esecuzione Nextflow.
Quando Nextflow inizia, esamina le CPU e la memoria disponibili.
Se le risorse delle attività pronte per l'esecuzione superano le risorse disponibili, Nextflow tratterrà le ultime attività dall'esecuzione fino a quando una o più delle attività precedenti non saranno terminate, liberando le risorse necessarie.

L'executor locale è conveniente ed efficiente, ma è limitato a quella singola macchina. Per carichi di lavoro molto grandi, potreste scoprire che la vostra macchina locale è un collo di bottiglia, sia perché avete una singola attività che richiede più risorse di quelle disponibili, sia perché avete così tante attività che aspettare che una singola macchina le esegua richiederebbe troppo tempo.

Nextflow supporta [molti backend di esecuzione diversi](https://nextflow.io/docs/latest/executor.html), inclusi scheduler HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor e altri) così come backend di esecuzione cloud (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e altro).

### 4.1. Indirizzare un backend diverso

La scelta dell'executor è impostata da una direttiva di processo chiamata `executor`.
Per impostazione predefinita è impostato su `local`, quindi la seguente configurazione è implicita:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Per impostare l'executor per indirizzare un backend diverso, specifichereste semplicemente l'executor che desiderate utilizzando una sintassi simile a quella descritta sopra per le allocazioni di risorse (vedere [Executors](https://nextflow.io/docs/latest/executor.html) per tutte le opzioni).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Non possiamo effettivamente testare questo nell'ambiente di formazione perché non è configurato per connettersi a un HPC.

### 4.2. Gestire la sintassi specifica del backend per i parametri di esecuzione

La maggior parte delle piattaforme di calcolo ad alte prestazioni consente (e talvolta richiede) di specificare determinati parametri come richieste e limitazioni di allocazione delle risorse (ad es. numero di CPU e memoria) e nome della coda di lavoro da utilizzare.

Sfortunatamente, ciascuno di questi sistemi utilizza tecnologie, sintassi e configurazioni diverse per definire come un lavoro dovrebbe essere definito e inviato allo scheduler pertinente.

??? abstract "Esempi"

    Ad esempio, lo stesso lavoro che richiede 8 CPU e 4GB di RAM da eseguire sulla coda "my-science-work" deve essere espresso nei seguenti modi diversi a seconda del backend.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Fortunatamente, Nextflow semplifica tutto questo.
Fornisce una sintassi standardizzata in modo da poter specificare le proprietà pertinenti come `cpus`, `memory` e `queue` una sola volta (vedere [Direttive di processo](https://nextflow.io/docs/latest/reference/process.html#process-directives) per tutte le opzioni disponibili).
Quindi, al momento dell'esecuzione, Nextflow utilizzerà quelle impostazioni per generare gli script appropriati specifici del backend in base all'impostazione dell'executor.

Tratteremo quella sintassi standardizzata nella prossima sezione.

### Takeaway

Ora sapete come cambiare l'executor per utilizzare diversi tipi di infrastruttura di calcolo.

### Cosa c'è dopo?

Imparate come valutare ed esprimere allocazioni e limitazioni di risorse in Nextflow.

---

## 5. Controllare le allocazioni delle risorse di calcolo

??? example "Scenario"

    La vostra pipeline continua a fallire sul cluster perché le attività vengono terminate per aver superato i limiti di memoria.
    O forse vi viene addebitato per risorse che non state utilizzando e volete ottimizzare i costi.

La maggior parte delle piattaforme di calcolo ad alte prestazioni consente (e talvolta richiede) di specificare determinati parametri di allocazione delle risorse come numero di CPU e memoria.

Per impostazione predefinita, Nextflow utilizzerà una singola CPU e 2GB di memoria per ogni processo.
Le direttive di processo corrispondenti sono chiamate `cpus` e `memory`, quindi la seguente configurazione è implicita:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Potete modificare questi valori, sia per tutti i processi che per processi specifici denominati, utilizzando direttive di processo aggiuntive nel vostro file di configurazione.
Nextflow li tradurrà nelle istruzioni appropriate per l'executor scelto.

Ma come fate a sapere quali valori utilizzare?

### 5.1. Eseguire il workflow per generare un report di utilizzo delle risorse

??? example "Scenario"

    Non sapete quanta memoria o CPU richiedono i vostri processi e volete evitare di sprecare risorse o di far terminare i lavori.

Se non sapete in anticipo quanta CPU e memoria è probabile che i vostri processi necessitino, potete fare un po' di profilazione delle risorse, il che significa eseguire il workflow con alcune allocazioni predefinite, registrare quanto ha utilizzato ogni processo e da lì, stimare come regolare le allocazioni di base.

Convenientemente, Nextflow include strumenti integrati per farlo e genererà volentieri un report per voi su richiesta.

Per farlo, aggiungete `-with-report <filename>.html` alla vostra riga di comando.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Il report è un file html, che potete scaricare e aprire nel vostro browser. Potete anche fare clic destro su di esso nell'esploratore di file a sinistra e fare clic su `Show preview` per visualizzarlo nell'ambiente di formazione.

Prendetevi qualche minuto per esaminare il report e vedere se riuscite a identificare alcune opportunità per regolare le risorse.
Assicuratevi di fare clic sulle schede che mostrano i risultati di utilizzo come percentuale di ciò che è stato allocato.

Vedere [Reports](https://nextflow.io/docs/latest/reports.html) per la documentazione su tutte le funzionalità disponibili.

### 5.2. Impostare le allocazioni di risorse per tutti i processi

La profilazione mostra che i processi nel nostro workflow di formazione sono molto leggeri, quindi riduciamo l'allocazione di memoria predefinita a 1GB per processo.

Aggiungete quanto segue al vostro file `nextflow.config`, prima della sezione dei parametri della pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Impostazioni di processo
*/
process {
    memory = 1.GB
}
```

Questo aiuterà a ridurre la quantità di calcolo che consumiamo.

### 5.3. Impostare le allocazioni di risorse per un processo specifico

Allo stesso tempo, fingeremo che il processo `cowpy` richieda più risorse degli altri, solo per poter dimostrare come regolare le allocazioni per un singolo processo.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Impostazioni di processo
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
    * Impostazioni di processo
    */
    process {
        memory = 1.GB
    }
    ```

Con questa configurazione, tutti i processi richiederanno 1GB di memoria e una singola CPU (il valore predefinito implicito), tranne il processo `cowpy`, che richiederà 2GB e 2 CPU.

!!! info

    Se avete una macchina con poche CPU e allocate un numero elevato per processo, potreste vedere le chiamate di processo messe in coda una dietro l'altra.
    Questo perché Nextflow assicura che non richiediamo più CPU di quelle disponibili.

### 5.4. Eseguire il workflow con la configurazione aggiornata

Proviamolo, fornendo un nome file diverso per il report di profilazione in modo da poter confrontare le prestazioni prima e dopo le modifiche alla configurazione.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Probabilmente non noterete alcuna differenza reale poiché questo è un carico di lavoro così piccolo, ma questo è l'approccio che utilizzereste per analizzare le prestazioni e i requisiti di risorse di un workflow reale.

È molto utile quando i vostri processi hanno requisiti di risorse diversi. Vi consente di dimensionare correttamente le allocazioni di risorse che impostate per ogni processo in base a dati effettivi, non a congetture.

!!! tip

    Questo è solo un piccolo assaggio di ciò che potete fare per ottimizzare l'uso delle risorse.
    Nextflow stesso ha una [logica di retry dinamico](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) davvero interessante integrata per riprovare i lavori che falliscono a causa di limitazioni di risorse.
    Inoltre, la Seqera Platform offre strumenti basati sull'IA per ottimizzare automaticamente le vostre allocazioni di risorse.

### 5.5. Aggiungere limiti di risorse

A seconda dell'executor di calcolo e dell'infrastruttura di calcolo che state utilizzando, potrebbero esserci alcuni vincoli su ciò che potete (o dovete) allocare.
Ad esempio, il vostro cluster potrebbe richiedere di rimanere entro determinati limiti.

Potete utilizzare la direttiva `resourceLimits` per impostare le limitazioni pertinenti. La sintassi appare così quando è da sola in un blocco di processo:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow tradurrà questi valori nelle istruzioni appropriate a seconda dell'executor che avete specificato.

Non eseguiremo questo, poiché non abbiamo accesso all'infrastruttura pertinente nell'ambiente di formazione.
Tuttavia, se provaste a eseguire il workflow con allocazioni di risorse che superano questi limiti, quindi cercaste il comando `sbatch` nel file di script `.command.run`, vedreste che le richieste che vengono effettivamente inviate all'executor sono limitate ai valori specificati da `resourceLimits`.

??? info "Configurazioni di riferimento istituzionali"

    Il progetto nf-core ha compilato una [raccolta di file di configurazione](https://nf-co.re/configs/) condivisi da varie istituzioni in tutto il mondo, che coprono un'ampia gamma di executor HPC e cloud.

    Quelle configurazioni condivise sono preziose sia per le persone che lavorano lì e possono quindi semplicemente utilizzare la configurazione della loro istituzione pronta all'uso, sia come modello per le persone che stanno cercando di sviluppare una configurazione per la propria infrastruttura.

### Takeaway

Sapete come generare un report di profilazione per valutare l'utilizzo delle risorse e come modificare le allocazioni di risorse per tutti i processi e/o per singoli processi, nonché impostare limitazioni di risorse per l'esecuzione su HPC.

### Cosa c'è dopo?

Imparate come configurare profili di configurazione preimpostati e passare da uno all'altro al momento dell'esecuzione.

---

## 6. Utilizzare i profili per passare tra configurazioni preimpostate

??? example "Scenario"

    Passate regolarmente dall'esecuzione di pipeline sul vostro laptop per lo sviluppo all'HPC della vostra istituzione per le esecuzioni di produzione.
    Siete stanchi di modificare manualmente le impostazioni di configurazione ogni volta che cambiate ambiente.

Vi abbiamo mostrato diversi modi in cui potete personalizzare la configurazione della vostra pipeline a seconda del progetto su cui state lavorando o dell'ambiente di calcolo che state utilizzando.

Potreste voler passare tra impostazioni alternative a seconda dell'infrastruttura di calcolo che state utilizzando. Ad esempio, potreste voler sviluppare ed eseguire test su piccola scala localmente sul vostro laptop, quindi eseguire carichi di lavoro su larga scala su HPC o cloud.

Nextflow vi consente di configurare un numero qualsiasi di [**profili**](https://nextflow.io/docs/latest/config.html#profiles) che descrivono diverse configurazioni, che potete quindi selezionare al momento dell'esecuzione utilizzando un argomento da riga di comando, piuttosto che dover modificare il file di configurazione stesso.

### 6.1. Creare profili per passare tra sviluppo locale ed esecuzione su HPC

Configuriamo due profili alternativi; uno per eseguire carichi su piccola scala su un computer normale, dove utilizzeremo container Docker, e uno per eseguire su un HPC universitario con uno scheduler Slurm, dove utilizzeremo pacchetti Conda.

#### 6.1.1. Configurare i profili

Aggiungete quanto segue al vostro file `nextflow.config`, dopo la sezione dei parametri della pipeline ma prima delle impostazioni di output:

```groovy title="nextflow.config" linenums="24"
/*
* Profili
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

Vedete che per l'HPC universitario, stiamo anche specificando limitazioni di risorse.

#### 6.1.2. Eseguire il workflow con un profilo

Per specificare un profilo nella nostra riga di comando Nextflow, utilizziamo l'argomento `-profile`.

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

Come potete vedere, questo ci consente di passare tra configurazioni molto comodamente al momento dell'esecuzione.

!!! warning

    Il profilo `univ_hpc` non verrà eseguito correttamente nell'ambiente di formazione poiché non abbiamo accesso a uno scheduler Slurm.

Se in futuro troviamo altri elementi di configurazione che si verificano sempre insieme a questi, possiamo semplicemente aggiungerli al/ai profilo/i corrispondente/i.
Possiamo anche creare profili aggiuntivi se ci sono altri elementi di configurazione che vogliamo raggruppare insieme.

### 6.2. Creare un profilo di parametri di test

??? example "Scenario"

    Volete che altri possano provare rapidamente la vostra pipeline senza raccogliere i propri dati di input.

I profili non sono solo per la configurazione dell'infrastruttura.
Possiamo anche usarli per impostare valori predefiniti per i parametri del workflow, per rendere più facile per gli altri provare il workflow senza dover raccogliere valori di input appropriati da soli.
Potete considerare questa un'alternativa all'utilizzo di un file di parametri.

#### 6.2.1. Configurare il profilo

La sintassi per esprimere valori predefiniti in questo contesto appare così, per un profilo che chiamiamo `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Se aggiungiamo un profilo di test per il nostro workflow, il blocco `profiles` diventa:

```groovy title="nextflow.config" linenums="24"
/*
* Profili
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

Proprio come per i profili di configurazione tecnica, potete configurare più profili diversi specificando parametri con qualsiasi nome arbitrario desideriate.

#### 6.2.2. Eseguire il workflow localmente con il profilo di test

Convenientemente, i profili non si escludono a vicenda, quindi possiamo specificare più profili nella nostra riga di comando utilizzando la seguente sintassi `-profile <profilo1>,<profilo2>` (per qualsiasi numero di profili).

Se combinate profili che impostano valori per gli stessi elementi di configurazione e sono descritti nello stesso file di configurazione, Nextflow risolverà il conflitto utilizzando qualsiasi valore abbia letto per ultimo (_cioè_ qualunque cosa venga dopo nel file).
Se le impostazioni in conflitto sono impostate in diverse fonti di configurazione, si applica l'[ordine di precedenza](https://www.nextflow.io/docs/latest/config.html) predefinito.

Proviamo ad aggiungere il profilo di test al nostro comando precedente:

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

Questo utilizzerà Docker dove possibile e produrrà output sotto `results_config/test`, e questa volta il carattere è il duo comico `dragonandcow`.

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

Questo significa che finché distribuiamo eventuali file di dati di test con il codice del workflow, chiunque può provare rapidamente il workflow senza dover fornire i propri input tramite la riga di comando o un file di parametri.

!!! tip

    Possiamo puntare a URL per file più grandi che sono archiviati esternamente.
    Nextflow li scaricherà automaticamente finché c'è una connessione aperta.

    Per maggiori dettagli, vedere la Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Utilizzare `nextflow config` per vedere la configurazione risolta

Come notato sopra, a volte lo stesso parametro può essere impostato su valori diversi in profili che volete combinare.
E più in generale, ci sono numerosi posti dove gli elementi di configurazione possono essere archiviati, e talvolta le stesse proprietà possono essere impostate su valori diversi in posti diversi.

Nextflow applica un [ordine di precedenza](https://nextflow.io/docs/latest/config.html#configuration-file) impostato per risolvere eventuali conflitti, ma può essere difficile determinarlo da soli.
E anche se nulla è in conflitto, può essere noioso cercare tutti i possibili posti dove le cose potrebbero essere configurate.

Fortunatamente, Nextflow include uno strumento di utilità conveniente chiamato `config` che può automatizzare l'intero processo per voi.

Lo strumento `config` esplorerà tutti i contenuti nella vostra directory di lavoro corrente, raccoglierà eventuali file di configurazione e produrrà la configurazione completamente risolta che Nextflow utilizzerebbe per eseguire il workflow.
Questo vi consente di scoprire quali impostazioni verranno utilizzate senza dover avviare nulla.

#### 6.3.1. Risolvere la configurazione predefinita

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

Questo vi mostra la configurazione di base che ottenete se non specificate nulla di extra nella riga di comando.

#### 6.3.2. Risolvere la configurazione con impostazioni specifiche attivate

Se fornite parametri da riga di comando, ad es. abilitando uno o più profili o caricando un file di parametri, il comando li prenderà inoltre in considerazione.

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

Questo diventa particolarmente utile per progetti complessi che coinvolgono più livelli di configurazione.

### Takeaway

Sapete come utilizzare i profili per selezionare una configurazione preimpostata al momento dell'esecuzione con il minimo fastidio.
Più in generale, sapete come configurare le vostre esecuzioni di workflow per adattarle a diverse piattaforme di calcolo e migliorare la riproducibilità delle vostre analisi.

### Cosa c'è dopo?

Imparate come eseguire pipeline direttamente da repository remoti come GitHub.

---

## 7. Eseguire pipeline da repository remoti

??? example "Scenario"

    Volete eseguire una pipeline ben consolidata come quelle di nf-core senza dover scaricare e gestire il codice da soli.

Finora abbiamo eseguito script di workflow situati nella directory corrente.
In pratica, vorrete spesso eseguire pipeline archiviate in repository remoti, come GitHub.

Nextflow rende questo semplice: potete eseguire qualsiasi pipeline direttamente da un URL di repository Git senza scaricarla manualmente prima.

### 7.1. Eseguire una pipeline da GitHub

La sintassi di base per eseguire una pipeline remota è `nextflow run <repository>`, dove `<repository>` può essere un percorso di repository GitHub come `nextflow-io/hello`, un URL completo o un percorso a GitLab, Bitbucket o altri servizi di hosting Git.

Provate a eseguire la pipeline demo ufficiale "hello" di Nextflow:

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

La prima volta che eseguite una pipeline remota, Nextflow la scarica e la memorizza nella cache localmente.
Le esecuzioni successive utilizzano la versione memorizzata nella cache a meno che non richiediate esplicitamente un aggiornamento.

### 7.2. Specificare una versione per la riproducibilità

Per impostazione predefinita, Nextflow esegue l'ultima versione dal branch predefinito.
Potete specificare una versione particolare (tag), branch o commit utilizzando il flag `-r`:

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

### Takeaway

Sapete come eseguire pipeline direttamente da GitHub e altri repository remoti, e come specificare versioni per la riproducibilità.

### Cosa c'è dopo?

Datevi una grande pacca sulla spalla!
Sapete tutto ciò che dovete sapere per iniziare a eseguire e gestire pipeline Nextflow.

Questo conclude questo corso, ma se siete desiderosi di continuare a imparare, abbiamo due raccomandazioni principali:

- Se volete approfondire lo sviluppo delle vostre pipeline, date un'occhiata a [Hello Nextflow](../hello_nextflow/index.md), un corso per principianti che copre la stessa progressione generale di questo ma entra molto più nel dettaglio sui canali e gli operatori.
- Se desiderate continuare a imparare come eseguire pipeline Nextflow senza approfondire il codice, date un'occhiata alla prima parte di [Hello nf-core](../hello_nf-core/index.md), che introduce gli strumenti per trovare ed eseguire pipeline dal progetto [nf-core](https://nf-co.re/) estremamente popolare.

Buon divertimento!

---

## Quiz

<quiz>
Quando i valori dei parametri sono impostati sia nel file del workflow che in `nextflow.config`, quale ha la precedenza?
- [ ] Il valore del file del workflow
- [x] Il valore del file di configurazione
- [ ] Il primo valore incontrato
- [ ] Causa un errore

Per saperne di più: [1.1. Impostare i valori in `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Qual è la differenza di sintassi tra l'impostazione di un valore predefinito di un parametro in un file di workflow rispetto a un file di configurazione?
- [ ] Usano la stessa sintassi
- [x] Il workflow usa dichiarazione tipizzata (`#!groovy param: Type = value`), la configurazione usa assegnazione (`#!groovy param = value`)
- [ ] La configurazione usa dichiarazione tipizzata, il workflow usa assegnazione
- [ ] Solo i file di configurazione possono impostare valori predefiniti

Per saperne di più: [1.1. Impostare i valori in `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Come si specifica un file di parametri quando si esegue un workflow?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Per saperne di più: [1.3. Utilizzare un file di parametri](#13-use-a-parameter-file)
</quiz>

<quiz>
Cosa controlla l'opzione di configurazione `outputDir`?
- [ ] La posizione della directory di lavoro
- [x] Il percorso base dove vengono pubblicati gli output del workflow
- [ ] La directory per i file di log
- [ ] La posizione dei file dei moduli

Per saperne di più: [2.1. Personalizzare il nome della directory outputDir](#21-customize-the-outputdir-directory-name)
</quiz>

<quiz>
Come si fa riferimento dinamicamente a un nome di processo nella configurazione del percorso di output?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<processo>.name"`
- [x] `#!groovy path { <processo>.name }`
- [ ] `@processName`

Per saperne di più: [2.2. Organizzare gli output per processo](#22-organize-outputs-by-process)
</quiz>

<quiz>
Se sia Docker che Conda sono abilitati e un processo ha entrambe le direttive, quale ha la priorità?
- [x] Docker (container)
- [ ] Conda
- [ ] Il primo definito nel processo
- [ ] Causa un errore

Per saperne di più: [3. Selezionare una tecnologia di packaging del software](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Qual è l'executor predefinito in Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Per saperne di più: [4. Selezionare una piattaforma di esecuzione](#4-select-an-execution-platform)
</quiz>

<quiz>
Quale comando genera un report di utilizzo delle risorse?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Per saperne di più: [5.1. Eseguire il workflow per generare un report di utilizzo delle risorse](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
Come si impostano i requisiti di risorse per un processo specifico denominato `cowpy` nel file di configurazione?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Per saperne di più: [5.3. Impostare le allocazioni di risorse per un processo specifico](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
Cosa fa la direttiva `resourceLimits`?
- [ ] Imposta i requisiti minimi di risorse
- [ ] Alloca risorse ai processi
- [x] Limita le risorse massime che possono essere richieste
- [ ] Monitora l'utilizzo delle risorse in tempo reale

Per saperne di più: [5.5. Aggiungere limiti di risorse](#55-add-resource-limits)
</quiz>

<quiz>
Come si specificano più profili in un singolo comando?
- [ ] `-profile profilo1 -profile profilo2`
- [ ] `-profiles profilo1,profilo2`
- [x] `-profile profilo1,profilo2`
- [ ] `--profile profilo1 --profile profilo2`

Per saperne di più: [6. Utilizzare i profili per passare tra configurazioni preimpostate](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
Quale comando mostra la configurazione completamente risolta che Nextflow utilizzerebbe?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Per saperne di più: [6.3. Utilizzare `nextflow config` per vedere la configurazione risolta](#63-use-nextflow-config-to-see-the-resolved-configuration)
</quiz>

<quiz>
Per cosa possono essere utilizzati i profili? (Selezionare tutte le opzioni applicabili)
- [x] Definire impostazioni specifiche dell'infrastruttura (executor, container)
- [x] Impostare limiti di risorse per diversi ambienti
- [x] Fornire parametri di test per facilitare il test del workflow
- [ ] Definire nuovi processi

Per saperne di più: [6. Utilizzare i profili per passare tra configurazioni preimpostate](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
