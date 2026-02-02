# Parte 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/06_hello_config.md).
///
-->

Questa sezione esplorerà come configurare e gestire la configurazione della vostra pipeline Nextflow in modo che possiate personalizzarne il comportamento, adattarla a diversi ambienti e ottimizzare l'uso delle risorse _senza modificare una singola riga del codice del workflow stesso_.

Ci sono diversi modi per farlo, che possono essere usati in combinazione e vengono interpretati secondo l'ordine di precedenza descritto [qui](https://www.nextflow.io/docs/latest/config.html).

In questa parte del corso, vi mostreremo il meccanismo di file di configurazione più semplice e comune, il file `nextflow.config`, che avete già incontrato nella Parte 5: Hello Containers.

Esamineremo i componenti essenziali della configurazione di Nextflow come le direttive dei processi, gli executor, i profili e i file di parametri.
Imparando a utilizzare efficacemente queste opzioni di configurazione, potete migliorare la flessibilità, la scalabilità e le prestazioni delle vostre pipeline.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato le Parti 1-5 del corso [Hello Nextflow](./index.md) e abbiate una pipeline funzionante completa.

    Se state iniziando il corso da questo punto, dovrete copiare la directory `modules` e il file `nextflow.config` dalle soluzioni:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    Il file `nextflow.config` contiene la riga `docker.enabled = true` che abilita l'uso dei container Docker.

    Se non avete familiarità con la pipeline Hello o avete bisogno di un promemoria, vedete [questa pagina informativa](../info/hello_pipeline.md).

---

## 0. Riscaldamento: Eseguire `hello-config.nf`

Useremo lo script del workflow `hello-config.nf` come punto di partenza.
È equivalente allo script prodotto seguendo la Parte 5 di questo corso di formazione, tranne che abbiamo cambiato le destinazioni dell'output:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Solo per assicurarci che tutto funzioni, eseguite lo script una volta prima di apportare modifiche:

```bash
nextflow run hello-config.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Come in precedenza, troverete i file di output nella directory specificata nel blocco `output` (`results/hello_config/`).

??? abstract "Contenuti della directory"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

L'output ASCII art finale è nella directory `results/hello_config/`, sotto il nome `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contenuti del file"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Se tutto ha funzionato, siete pronti a imparare come configurare le vostre pipeline.

---

## 1. Gestire i parametri di input del workflow

Inizieremo con un aspetto della configurazione che è semplicemente un'estensione di ciò con cui abbiamo lavorato finora: la gestione dei parametri di input.

Attualmente, il nostro workflow è configurato per accettare diversi valori di parametro tramite la riga di comando, con valori predefiniti impostati in un blocco `params` nello script del workflow stesso.
Tuttavia, potreste voler sovrascrivere quei valori predefiniti senza dover specificare parametri sulla riga di comando, o modificare il file di script originale.

Ci sono diversi modi per farlo; vi mostreremo tre modi di base che sono molto comunemente usati.

### 1.1. Spostare i valori predefiniti in `nextflow.config`

Questo è l'approccio più semplice, anche se è forse il meno flessibile dato che il file `nextflow.config` principale non è qualcosa che si vuole modificare per ogni esecuzione.
Ma ha il vantaggio di separare le preoccupazioni della _dichiarazione_ dei parametri nel workflow (che appartiene sicuramente lì) rispetto alla fornitura di _valori predefiniti_, che sono più appropriati in un file di configurazione.

Facciamo questo in due passaggi.

#### 1.1.1. Creare un blocco `params` nel file di configurazione

Effettuate le seguenti modifiche al codice nel file `nextflow.config`:

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

Nota che non abbiamo semplicemente copiato il blocco `params` dal workflow al file di configurazione.
La sintassi è leggermente diversa.
Nel file del workflow, quelle sono dichiarazioni tipizzate.
Nella configurazione, quelle sono assegnazioni di valori.

Tecnicamente, questo è sufficiente per sovrascrivere i valori predefiniti ancora specificati nel file del workflow.
Potreste modificare il personaggio, per esempio, ed eseguire il workflow per verificare che il valore impostato nel file di configurazione sovrascriva quello impostato nel file del workflow.

Ma nello spirito di spostare completamente la configurazione nel file di configurazione, rimuoviamo completamente quei valori dal file del workflow.

#### 1.1.2. Rimuovere i valori dal blocco `params` nel file del workflow

Effettuate le seguenti modifiche al codice nel file del workflow `hello-config.nf`:

=== "Dopo"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
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

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Ora il file del workflow stesso non imposta alcun valore predefinito per questi parametri.

#### 1.1.3. Eseguire la pipeline

Testiamo che funzioni correttamente.

```bash
nextflow run hello-config.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima.

L'output ASCII art finale è nella directory `results/hello_config/`, sotto il nome `cowpy-COLLECTED-batch-output.txt`, come prima.

??? abstract "Contenuti del file"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

### 1.2. Usare un file di configurazione specifico per l'esecuzione

Questo è ottimo, ma a volte potreste voler eseguire alcuni esperimenti temporanei con valori predefiniti diversi senza toccare il file di configurazione principale.
Potete farlo creando un nuovo file `nextflow.config` in una sottodirectory che userete come directory di lavoro per i vostri esperimenti.

#### 1.2.1. Creare la directory di lavoro con una configurazione vuota

Iniziamo creando una nuova directory e spostandoci al suo interno:

```bash
mkdir -p tux-run
cd tux-run
```

Poi, crei un file di configurazione vuoto in quella directory:

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

Nota che il percorso al file di input deve riflettere la struttura della directory.

#### 1.2.3. Eseguire la pipeline

Ora possiamo eseguire la nostra pipeline dall'interno della nostra nuova directory di lavoro.
Assicuratevi di adattare il percorso di conseguenza!

```bash
nextflow run ../hello-config.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Questo creerà un nuovo set di directory sotto `tux-run/` incluse `tux-run/work/` e `tux-run/results/`.

In questa esecuzione, Nextflow combina il `nextflow.config` nella nostra directory corrente con il `nextflow.config` nella directory radice della pipeline, e quindi sovrascrive il personaggio predefinito (turkey) con il personaggio tux.

Il file di output finale dovrebbe contenere il personaggio tux che dice i saluti.

??? abstract "Contenuti del file"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

Ora vediamo un altro modo utile per impostare i valori dei parametri.

### 1.3. Usare un file di parametri

L'approccio della sottodirectory funziona benissimo per sperimentare, ma comporta un po' di configurazione e richiede di adattare i percorsi di conseguenza.
C'è un approccio più semplice per quando volete eseguire la vostra pipeline con un set specifico di valori, o permettere a qualcun altro di farlo con il minimo sforzo.

Nextflow vi permette di specificare parametri tramite un file di parametri in formato YAML o JSON, il che rende molto conveniente gestire e distribuire set alternativi di valori predefiniti, per esempio, così come valori di parametri specifici per l'esecuzione.

#### 1.3.1. Esaminare il file di parametri di esempio

Per dimostrare questo, forniamo un file di parametri di esempio nella directory corrente, chiamato `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
{
  input: "greetings.csv"
  batch: "yaml"
  character: "stegosaurus"
}
```

Questo file di parametri contiene una coppia chiave-valore per ciascuno degli input che vogliamo specificare.
Nota l'uso dei due punti (`:`) invece dei segni di uguale (`=`) se confronta la sintassi con il file di configurazione.
Il file di configurazione è scritto in Groovy, mentre il file di parametri è scritto in YAML.

!!! info "Informazione"

    Forniamo anche una versione JSON del file di parametri come esempio ma non la eseguiremo qui.
    Sentitevi liberi di provare quella da soli.

#### 1.3.2. Eseguire la pipeline

Per eseguire il workflow con questo file di parametri, aggiunga semplicemente `-params-file <filename>` al comando base.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Il file di output finale dovrebbe contenere il personaggio stegosaurus che dice i saluti.

??? abstract "Contenuti del file"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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

Usare un file di parametri può sembrare eccessivo quando si hanno solo pochi parametri da specificare, ma alcune pipeline si aspettano decine di parametri.
In quei casi, usare un file di parametri ci permetterà di fornire valori di parametri a runtime senza dover digitare righe di comando massive e senza modificare lo script del workflow.

Rende anche più facile distribuire set di parametri ai collaboratori, o come informazione di supporto per una pubblicazione, per esempio.
Questo rende il vostro lavoro più riproducibile da altri.

### Takeaway

Sapete come sfruttare le principali opzioni di configurazione per gestire gli input del workflow.

### Cosa c'è dopo?

Imparare come gestire dove e come vengono pubblicati gli output del vostro workflow.

---

## 2. Gestire gli output del workflow

Finora abbiamo hardcodato tutti i percorsi per le dichiarazioni di output a livello di workflow, e come abbiamo notato quando abbiamo iniziato ad aggiungere output multipli, può esserci un po' di ripetizione coinvolta.

Vediamo alcuni modi comuni in cui potreste configurare questo per essere più flessibile.

### 2.1. Personalizzare il nome della directory `outputDir`

Per ogni capitolo di questo corso, abbiamo pubblicato gli output in una sottodirectory diversa hardcodata nelle definizioni di output.

Cambiamo questo per usare un parametro configurabile dall'utente.
Potremmo creare un parametro completamente nuovo per questo, ma useremo il parametro `batch` dato che è già lì.

#### 2.1.1. Impostare un valore per `outputDir` nel file di configurazione

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
    outputDir = "results/${params.batch}"
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

Questo sostituirà il percorso predefinito integrato, `results/`, con `results/` più il valore del parametro `batch` come sottodirectory.
Potreste anche cambiare la parte `results` se lo desideraste.

Per un cambiamento temporaneo, potreste impostare questa opzione dalla riga di comando usando il parametro `-output-dir` nel vostro comando (ma in quel caso non potreste usare il valore del parametro `batch`).

#### 2.1.2. Rimuovere la parte ripetuta del percorso hardcodato

Abbiamo ancora una sottodirectory hardcodata nelle opzioni di output, quindi eliminiamola ora.

Effettuate le seguenti modifiche al codice nel file del workflow:

=== "Dopo"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Avremmo anche potuto semplicemente aggiungere `${params.batch}` a ogni percorso invece di modificare il `outputDir` predefinito, ma questo è più conciso.

#### 2.1.3. Eseguire la pipeline

Testiamo che funzioni correttamente, impostando il nome del batch a `outdir` dalla riga di comando.

```bash
nextflow run hello-config.nf --batch outdir
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results/outdir/`.

??? abstract "Contenuti della directory"

    ```console
    results/outdir/
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

Potete combinare questo approccio con definizioni di percorso personalizzate per costruire qualsiasi gerarchia di directory desideriate.

### 2.2. Organizzare gli output per processo

Un modo popolare per organizzare ulteriormente gli output è farlo per processo, _cioè_ creare sottodirectory per ogni processo eseguito nella pipeline.

#### 2.2.1. Sostituire i percorsi di output con un riferimento ai nomi dei processi

Tutto ciò che dovete fare è riferirvi al nome del processo come `<task>.name` nella dichiarazione del percorso di output.

Effettuate le seguenti modifiche nel file del workflow:

=== "Dopo"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Questo rimuove gli elementi hardcodati rimanenti dalla configurazione del percorso di output.

#### 2.2.2. Eseguire la pipeline

Testiamo che funzioni correttamente, impostando il nome del batch a `pnames` dalla riga di comando.

```bash
nextflow run hello-config.nf --batch pnames
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results/pnames/`, e sono raggruppati per processo.

??? abstract "Contenuti della directory"

    ```console
    results/pnames/
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

Nota che qui abbiamo eliminato la distinzione tra `intermediates` rispetto agli output finali al livello superiore.
Potreste ovviamente mescolare questi approcci, per esempio impostando il percorso del primo output come `intermediates/${sayHello.process}`

### 2.3. Impostare la modalità di pubblicazione a livello di workflow

Infine, nello spirito di ridurre la quantità di codice ripetitivo, possiamo sostituire le dichiarazioni `mode` per-output con una singola riga nella configurazione.

#### 2.3.1. Aggiungere `workflow.output.mode` al file di configurazione

Aggiungete il seguente codice al file `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

Proprio come l'opzione `outputDir`, dare a `workflow.output.mode` un valore nel file di configurazione sarebbe sufficiente per sovrascrivere ciò che è impostato nel file del workflow, ma rimuoviamo comunque il codice non necessario.

#### 2.3.2. Rimuovere la modalità di output dal file del workflow

Effettuate le seguenti modifiche nel file del workflow:

=== "Dopo"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { sayHello.process }
        }
        uppercased {
            path { convertToUpper.process }
        }
        collected {
            path { collectGreetings.process }
        }
        batch_report {
            path { collectGreetings.process }
        }
        cowpy_art {
            path { cowpy.process }
        }
    }
    ```

=== "Prima"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.process }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.process }
            mode 'copy'
        }
        collected {
            path { collectGreetings.process }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.process }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.process }
            mode 'copy'
        }
    }
    ```

Più conciso, vero?

#### 2.3.3. Eseguire la pipeline

Testiamo che funzioni correttamente, impostando il nome del batch a `outmode` dalla riga di comando.

```bash
nextflow run hello-config.nf --batch outmode
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Questo produce ancora lo stesso output di prima, tranne che questa volta troviamo i nostri output sotto `results/outmode/`.
Sono ancora tutte copie appropriate, non symlink.

??? abstract "Contenuti della directory"

    ```console
    results/outmode/
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

Il motivo principale per cui potreste ancora voler usare il modo per-output di impostare la modalità è se volete mescolare all'interno dello stesso workflow, _cioè_ avere alcuni output copiati e alcuni come symlink.

Ci sono molte altre opzioni che potete personalizzare in questo modo, ma speriamo che questo vi dia un senso della gamma di opzioni e di come utilizzarle efficacemente per adattarle alle vostre preferenze.

### Takeaway

Sapete come controllare la denominazione e la struttura delle directory dove vengono pubblicati i vostri output, così come la modalità di pubblicazione dell'output del workflow.

### Cosa c'è dopo?

Imparare come adattare la configurazione del vostro workflow al vostro ambiente di calcolo, partendo dalla tecnologia di packaging del software.

---

## 3. Selezionare una tecnologia di packaging del software

Finora abbiamo esaminato elementi di configurazione che controllano come gli input entrano e dove escono gli output. Ora è il momento di concentrarci più specificamente sull'adattamento della configurazione del vostro workflow al vostro ambiente di calcolo.

Il primo passo su quel percorso è specificare da dove provengono i pacchetti software che verranno eseguiti in ogni passaggio.
Sono già installati nell'ambiente di calcolo locale?
Dobbiamo recuperare immagini ed eseguirle tramite un sistema container?
O dobbiamo recuperare pacchetti Conda e costruire un ambiente Conda locale?

Nella primissima parte di questo corso di formazione (Parti 1-4) abbiamo semplicemente usato software installato localmente nel nostro workflow.
Poi nella Parte 5, abbiamo introdotto i container Docker e il file `nextflow.config`, che abbiamo usato per abilitare l'uso dei container Docker.

Ora vediamo come potete configurare un'opzione alternativa di packaging del software tramite il file `nextflow.config`.

### 3.1. Disabilitare Docker e abilitare Conda nel file di configurazione

Facciamo finta di lavorare su un cluster HPC e l'amministratore non permette l'uso di Docker per motivi di sicurezza.
Fortunatamente per noi, Nextflow supporta diverse altre tecnologie container come Singularity (che è più ampiamente usato su HPC), e gestori di pacchetti software come Conda.

Potete cambiare il vostro file di configurazione per usare Conda invece di Docker.
Per farlo, cambiate il valore di `docker.enabled` a `false`, e aggiungete una direttiva che abilita l'uso di Conda:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Questo permetterà a Nextflow di creare e utilizzare ambienti Conda per i processi che hanno pacchetti Conda specificati.
Il che significa che ora dobbiamo aggiungere uno di quelli al nostro processo `cowpy`!

### 3.2. Specificare un pacchetto Conda nella definizione del processo

Abbiamo già recuperato l'URI per un pacchetto Conda contenente lo strumento `cowpy`: `conda-forge::cowpy==1.1.5`

Ora aggiungiamo l'URI alla definizione del processo `cowpy` usando la direttiva `conda`:

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

    Ci sono diversi modi per ottenere l'URI per un dato pacchetto conda.
    Raccomandiamo di usare la query di ricerca di [Seqera Containers](https://seqera.io/containers/), che vi darà un URI che potete copiare e incollare, anche se non avete intenzione di creare un container da esso.

### 3.3. Eseguire il workflow per verificare che possa usare Conda

Proviamo.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Output del comando"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Questo dovrebbe funzionare senza problemi e produrre gli stessi output di prima sotto `results/conda`.

Dietro le quinte, Nextflow ha recuperato i pacchetti Conda e creato l'ambiente, il che normalmente richiede un po' di lavoro; quindi è bello che non dobbiamo fare nulla di tutto ciò noi stessi!

!!! note "Nota"

    Questo viene eseguito rapidamente perché il pacchetto `cowpy` è abbastanza piccolo, ma se state lavorando con pacchetti grandi, potrebbe richiedere un po' più di tempo del solito la prima volta, e potreste vedere l'output della console rimanere 'bloccato' per un minuto o così prima di completare.
    Questo è normale ed è dovuto al lavoro extra che Nextflow fa la prima volta che usa un nuovo pacchetto.

Dal nostro punto di vista, sembra che funzioni esattamente come l'esecuzione con Docker, anche se sul backend i meccanismi sono un po' diversi.

Questo significa che siamo pronti a eseguire con ambienti Conda se necessario.

??? info "Mescolare Docker e Conda"

    Dato che queste direttive sono assegnate per processo, è possibile 'mescolare', _cioè_ configurare alcuni dei processi nel vostro workflow per essere eseguiti con Docker e altri con Conda, per esempio, se l'infrastruttura di calcolo che state usando supporta entrambi.
    In quel caso, abilitereste sia Docker che Conda nel vostro file di configurazione.
    Se entrambi sono disponibili per un dato processo, Nextflow darà priorità ai container.

    E come notato in precedenza, Nextflow supporta diverse altre tecnologie di packaging software e container, quindi non è limitato a quelle due.

### Takeaway

Sapete come configurare quale pacchetto software ogni processo dovrebbe usare, e come passare da una tecnologia all'altra.

### Cosa c'è dopo?

Imparare come cambiare la piattaforma di esecuzione usata da Nextflow per effettivamente fare il lavoro.

---

## 4. Selezionare una piattaforma di esecuzione

Finora, abbiamo eseguito la nostra pipeline con l'executor locale.
Questo esegue ogni task sulla macchina su cui è in esecuzione Nextflow.
Quando Nextflow inizia, guarda le CPU e la memoria disponibili.
Se le risorse dei task pronti per l'esecuzione superano le risorse disponibili, Nextflow tratterrà gli ultimi task dall'esecuzione fino a quando uno o più dei task precedenti non sono terminati, liberando le risorse necessarie.

L'executor locale è conveniente ed efficiente, ma è limitato a quella singola macchina. Per carichi di lavoro molto grandi, potreste scoprire che la vostra macchina locale è un collo di bottiglia, sia perché avete un singolo task che richiede più risorse di quelle disponibili, sia perché avete così tanti task che aspettare che una singola macchina li esegua richiederebbe troppo tempo.

Nextflow supporta [molti backend di esecuzione diversi](https://www.nextflow.io/docs/latest/executor.html), inclusi scheduler HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor e altri) così come backend di esecuzione cloud (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e altro).

### 4.1. Puntare a un backend diverso

La scelta dell'executor è impostata da una direttiva di processo chiamata `executor`.
Per impostazione predefinita è impostato su `local`, quindi la seguente configurazione è implicita:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Per impostare l'executor per puntare a un backend diverso, specificate semplicemente l'executor che volete usando una sintassi simile a quella descritta sopra per le allocazioni di risorse (vedete la [documentazione](https://www.nextflow.io/docs/latest/executor.html) per tutte le opzioni).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Avviso"

    Non possiamo effettivamente testare questo nell'ambiente di formazione perché non è configurato per connettersi a un HPC.

### 4.2. Gestire la sintassi specifica del backend per i parametri di esecuzione

La maggior parte delle piattaforme di calcolo ad alte prestazioni permette (e a volte richiede) di specificare certi parametri come le richieste di allocazione delle risorse e le limitazioni (per es. numero di CPU e memoria) e nome della coda di lavoro da usare.

Sfortunatamente, ciascuno di questi sistemi usa tecnologie, sintassi e configurazioni diverse per definire come un lavoro dovrebbe essere definito e inviato allo scheduler rilevante.

??? abstract "Esempi"

    Per esempio, lo stesso lavoro che richiede 8 CPU e 4GB di RAM per essere eseguito sulla coda "my-science-work" deve essere espresso in modi diversi a seconda del backend.

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
Fornisce una sintassi standardizzata in modo che possa specificare le proprietà rilevanti come `cpus`, `memory` e `queue` (veda la documentazione per altre proprietà) una sola volta.
Poi, a runtime, Nextflow userà quelle impostazioni per generare gli script appropriati specifici del backend basandosi sull'impostazione dell'executor.

Tratteremo quella sintassi standardizzata nella prossima sezione.

### Takeaway

Ora sapete come cambiare l'executor per usare diversi tipi di infrastruttura di calcolo.

### Cosa c'è dopo?

Imparare come valutare ed esprimere allocazioni e limitazioni di risorse in Nextflow.

---

## 5. Controllare le allocazioni delle risorse di calcolo

La maggior parte delle piattaforme di calcolo ad alte prestazioni permette (e a volte richiede) di specificare certi parametri di allocazione delle risorse come numero di CPU e memoria.

Per impostazione predefinita, Nextflow userà una singola CPU e 2GB di memoria per ogni processo.
Le corrispondenti direttive di processo si chiamano `cpus` e `memory`, quindi la seguente configurazione è implicita:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Potete modificare questi valori, sia per tutti i processi che per processi specifici nominati, usando direttive di processo aggiuntive nel vostro file di configurazione.
Nextflow le tradurrà nelle istruzioni appropriate per l'executor scelto.

Ma come sapete quali valori usare?

### 5.1. Eseguire il workflow per generare un report di utilizzo delle risorse

Se non sapete in anticipo quanta CPU e memoria i vostri processi probabilmente necessitano, potete fare un po' di profilazione delle risorse, il che significa che eseguite il workflow con alcune allocazioni predefinite, registrate quanto ogni processo ha usato, e da lì, stimate come regolare le allocazioni base.

Convenientemente, Nextflow include strumenti integrati per farlo, e genererà felicemente un report per voi su richiesta.

Per farlo, aggiungete `-with-report <filename>.html` alla vostra riga di comando.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Il report è un file html, che potete scaricare e aprire nel vostro browser. Potete anche fare clic destro su di esso nell'esploratore file a sinistra e cliccare su `Show preview` per visualizzarlo nell'ambiente di formazione.

Prendetevi qualche minuto per esaminare il report e vedere se riuscite a identificare alcune opportunità per regolare le risorse.
Assicuratevi di cliccare sulle schede che mostrano i risultati di utilizzo come percentuale di ciò che è stato allocato.
C'è della [documentazione](https://www.nextflow.io/docs/latest/reports.html) che descrive tutte le funzionalità disponibili.

### 5.2. Impostare allocazioni di risorse per tutti i processi

La profilazione mostra che i processi nel nostro workflow di formazione sono molto leggeri, quindi riduciamo l'allocazione di memoria predefinita a 1GB per processo.

Aggiungete il seguente al vostro file `nextflow.config`, prima della sezione dei parametri della pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

Questo aiuterà a ridurre la quantità di calcolo che consumiamo.

### 5.3. Impostare allocazioni di risorse per un processo specifico

Allo stesso tempo, facciamo finta che il processo `cowpy` richieda più risorse degli altri, solo per dimostrare come regolare le allocazioni per un singolo processo.

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

Con questa configurazione, tutti i processi richiederanno 1GB di memoria e una singola CPU (il valore predefinito implicito), tranne il processo `cowpy`, che richiederà 2GB e 2 CPU.

!!! tip "Suggerimento"

    Se avete una macchina con poche CPU e allocate un numero elevato per processo, potreste vedere le chiamate dei processi accodarsi l'una dietro l'altra.
    Questo perché Nextflow si assicura che non richiediamo più CPU di quelle disponibili.

### 5.4. Eseguire il workflow con la configurazione aggiornata

Proviamo, fornendo un nome di file diverso per il report di profilazione in modo da poter confrontare le prestazioni prima e dopo le modifiche alla configurazione.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Probabilmente non noterete alcuna differenza reale dato che questo è un carico di lavoro così piccolo, ma questo è l'approccio che usereste per analizzare le prestazioni e i requisiti di risorse di un workflow del mondo reale.

È molto utile quando i vostri processi hanno requisiti di risorse diversi. Vi permette di dimensionare correttamente le allocazioni di risorse che impostate per ogni processo basandovi su dati reali, non su congetture.

!!! tip "Suggerimento"

    Questo è solo un piccolo assaggio di ciò che potete fare per ottimizzare il vostro uso delle risorse.
    Nextflow stesso ha una [logica di retry dinamica](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) davvero interessante integrata per riprovare i lavori che falliscono a causa di limitazioni di risorse.
    Inoltre, la Piattaforma Seqera offre strumenti guidati dall'AI per ottimizzare le vostre allocazioni di risorse automaticamente.

### 5.5. Aggiungere limiti di risorse

A seconda di quale executor di calcolo e infrastruttura di calcolo state usando, potrebbero esserci alcuni vincoli su ciò che potete (o dovete) allocare.
Per esempio, il vostro cluster potrebbe richiedere di rimanere entro certi limiti.

Potete usare la direttiva `resourceLimits` per impostare le limitazioni rilevanti. La sintassi appare così quando è da sola in un blocco process:

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

Non eseguiremo questo, dato che non abbiamo accesso all'infrastruttura rilevante nell'ambiente di formazione.
Tuttavia, se provaste a eseguire il workflow con allocazioni di risorse che superano questi limiti, poi cercaste il comando `sbatch` nel file di script `.command.run`, vedreste che le richieste che effettivamente vengono inviate all'executor sono limitate ai valori specificati da `resourceLimits`.

??? info "Configurazioni di riferimento istituzionali"

    Il progetto nf-core ha compilato una [collezione di file di configurazione](https://nf-co.re/configs/) condivisi da varie istituzioni in tutto il mondo, che coprono una vasta gamma di executor HPC e cloud.

    Quelle configurazioni condivise sono preziose sia per le persone che lavorano lì e possono quindi semplicemente utilizzare la configurazione della loro istituzione pronta all'uso, sia come modello per le persone che stanno cercando di sviluppare una configurazione per la propria infrastruttura.

### Takeaway

Sapete come generare un report di profilazione per valutare l'utilizzo delle risorse e come modificare le allocazioni di risorse per tutti i processi e/o per singoli processi, così come impostare limitazioni di risorse per l'esecuzione su HPC.

### Cosa c'è dopo?

Imparare come configurare profili di configurazione preimpostati e passare da uno all'altro a runtime.

---

## 6. Usare profili per passare tra configurazioni preimpostate

Vi abbiamo mostrato diversi modi in cui potete personalizzare la configurazione della vostra pipeline a seconda del progetto su cui state lavorando o dell'ambiente di calcolo che state usando.

Potreste voler passare tra impostazioni alternative a seconda di quale infrastruttura di calcolo state usando. Per esempio, potreste voler sviluppare ed eseguire test su piccola scala localmente sul vostro laptop, poi eseguire carichi di lavoro su scala completa su HPC o cloud.

Nextflow vi permette di configurare qualsiasi numero di profili che descrivono diverse configurazioni, che potete poi selezionare a runtime usando un argomento da riga di comando, piuttosto che dover modificare il file di configurazione stesso.

### 6.1. Creare profili per passare tra sviluppo locale ed esecuzione su HPC

Configuriamo due profili alternativi; uno per eseguire carichi su piccola scala su un computer normale, dove useremo container Docker, e uno per eseguire su un HPC universitario con uno scheduler Slurm, dove useremo pacchetti Conda.

#### 6.1.1. Configurare i profili

Aggiungete il seguente al vostro file `nextflow.config`, dopo la sezione dei parametri della pipeline ma prima delle impostazioni di output:

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

Vedete che per l'HPC universitario, stiamo anche specificando limitazioni di risorse.

#### 6.1.2. Eseguire il workflow con un profilo

Per specificare un profilo nella nostra riga di comando Nextflow, usiamo l'argomento `-profile`.

Proviamo a eseguire il workflow con la configurazione `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Come potete vedere, questo ci permette di alternare tra configurazioni molto comodamente a runtime.

!!! warning "Avviso"

    Il profilo `univ_hpc` non funzionerà correttamente nell'ambiente di formazione dato che non abbiamo accesso a uno scheduler Slurm.

Se in futuro troviamo altri elementi di configurazione che co-occorrono sempre con questi, possiamo semplicemente aggiungerli al/ai profilo/i corrispondente/i.
Possiamo anche creare profili aggiuntivi se ci sono altri elementi di configurazione che vogliamo raggruppare insieme.

### 6.2. Creare un profilo di parametri di test

I profili non sono solo per la configurazione dell'infrastruttura.
Possiamo anche usarli per impostare valori predefiniti per i parametri del workflow, per rendere più facile per altri provare il workflow senza dover raccogliere valori di input appropriati da soli.
Potete considerare questo un'alternativa all'uso di un file di parametri.

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
        params.greeting = 'greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Proprio come per i profili di configurazione tecnica, potete configurare più profili diversi specificando parametri sotto qualsiasi nome arbitrario che desiderate.

#### 6.2.2. Eseguire il workflow localmente con il profilo di test

Convenientemente, i profili non sono mutuamente esclusivi, quindi possiamo specificare più profili nella nostra riga di comando usando la seguente sintassi `-profile <profile1>,<profile2>` (per qualsiasi numero di profili).

Se combinate profili che impostano valori per gli stessi elementi di configurazione e sono descritti nello stesso file di configurazione, Nextflow risolverà il conflitto usando qualsiasi valore abbia letto per ultimo (_cioè_ qualsiasi cosa venga dopo nel file).
Se le impostazioni in conflitto sono impostate in diverse fonti di configurazione, si applica l'[ordine di precedenza](https://www.nextflow.io/docs/latest/config.html) predefinito.

Proviamo ad aggiungere il profilo di test al nostro comando precedente:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Questo userà Docker dove possibile e produrrà output sotto `results/test`, e questa volta il personaggio è il duo comico `dragonandcow`.

??? abstract "Contenuti del file"

    ```console title="results/test/"
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

Questo significa che finché distribuiamo qualsiasi file di dati di test con il codice del workflow, chiunque può provare rapidamente il workflow senza dover fornire i propri input tramite la riga di comando o un file di parametri.

!!! tip "Suggerimento"

    Possiamo puntare a URL per file più grandi che sono memorizzati esternamente.
    Nextflow li scaricherà automaticamente finché c'è una connessione aperta.

    Per maggiori dettagli, vedete la Side Quest [Lavorare con i file](../side_quests/working_with_files.md)

### 6.3. Usare `nextflow config` per vedere la configurazione risolta

Come notato sopra, a volte lo stesso parametro può essere impostato a valori diversi in profili che volete combinare.
E più in generale, ci sono numerosi posti dove elementi di configurazione possono essere memorizzati, e a volte le stesse proprietà possono essere impostate a valori diversi in posti diversi.

Nextflow applica un [ordine di precedenza](https://www.nextflow.io/docs/latest/config.html) stabilito per risolvere qualsiasi conflitto, ma può essere complicato da determinare da soli.
E anche se nulla è in conflitto, può essere tedioso cercare tutti i possibili posti dove le cose potrebbero essere configurate:

Fortunatamente, Nextflow include uno strumento utility conveniente chiamato `config` che può automatizzare l'intero processo per voi.

Lo strumento `config` esplorerà tutti i contenuti nella vostra directory di lavoro corrente, raccoglierà qualsiasi file di configurazione, e produrrà la configurazione completamente risolta che Nextflow userebbe per eseguire il workflow.
Questo vi permette di scoprire quali impostazioni verranno usate senza dover lanciare nulla.

#### 6.3.1. Risolvere la configurazione predefinita

Eseguite questo comando per risolvere la configurazione che verrebbe applicata per impostazione predefinita.

```bash
nextflow config
```

??? success "Output del comando"

    ```groovy
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

    params {
      input = 'greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }
    ```

Questo vi mostra la configurazione base che ottenete se non specificate nulla di extra nella riga di comando.

#### 6.3.2. Risolvere la configurazione con impostazioni specifiche attivate

Se fornite parametri da riga di comando, es. abilitando uno o più profili o caricando un file di parametri, il comando prenderà in considerazione anche quelli.

```bash
nextflow config -profile my_laptop,test
```

??? success "Output del comando"

    ```groovy
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

    params {
      input = 'greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }
    ```

Questo diventa particolarmente utile per progetti complessi che coinvolgono più livelli di configurazione.

### Takeaway

Sapete come usare i profili per selezionare una configurazione preimpostata a runtime con il minimo sforzo.
Più in generale, sapete come configurare le esecuzioni del vostro workflow per adattarsi a diverse piattaforme di calcolo e migliorare la riproducibilità delle vostre analisi.

### Cosa c'è dopo?

Festeggiate e datevi una bella pacca sulla spalla! Avete completato il vostro primo corso di sviluppo Nextflow.

Passate al [riepilogo finale del corso](./next_steps.md) per rivedere ciò che avete imparato e scoprire cosa viene dopo.

---

## Quiz

<quiz>
Qual è il nome del file di configurazione che Nextflow carica automaticamente?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
Cosa ha la precedenza quando lo stesso parametro è impostato sia nel file di configurazione che nella riga di comando?
- [ ] Il valore del file di configurazione
- [x] Il valore della riga di comando
- [ ] Il primo valore incontrato
- [ ] Nessuno dei due; causa un errore

Per approfondire: [1.1. Spostare i valori predefiniti in `nextflow.config`](#11-spostare-i-valori-predefiniti-in-nextflowconfig)
</quiz>

<quiz>
Si possono avere sia Docker che Conda abilitati nella stessa configurazione?
- [x] Sì, Nextflow potete usare entrambi a seconda delle direttive del processo
- [ ] No, solo uno può essere abilitato alla volta
- [ ] Sì, ma solo nei profili
- [ ] No, sono mutuamente esclusivi
</quiz>

<quiz>
Se sia Docker che Conda sono abilitati e un processo ha entrambe le direttive, quale ha la priorità?
- [x] Docker (container)
- [ ] Conda
- [ ] Il primo definito
- [ ] Causa un errore

Per approfondire: [3. Selezionare una tecnologia di packaging del software](#3-selezionare-una-tecnologia-di-packaging-del-software)
</quiz>

<quiz>
Qual è l'allocazione di memoria predefinita per i processi Nextflow?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Nessun limite
</quiz>

<quiz>
Come si impostano i requisiti di risorse per un processo specifico nel file di configurazione?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Per approfondire: [5.3. Impostare allocazioni di risorse per un processo specifico](#53-impostare-allocazioni-di-risorse-per-un-processo-specifico)
</quiz>

<quiz>
Quale opzione da riga di comando genera un report di utilizzo delle risorse?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Per approfondire: [5.1. Eseguire il workflow per generare un report di utilizzo delle risorse](#51-eseguire-il-workflow-per-generare-un-report-di-utilizzo-delle-risorse)
</quiz>

<quiz>
Cosa fa la direttiva `resourceLimits`?
- [ ] Imposta i requisiti minimi di risorse
- [ ] Alloca risorse ai processi
- [x] Limita le risorse massime che possono essere richieste
- [ ] Monitora l'uso delle risorse

Per approfondire: [5.5. Aggiungere limiti di risorse](#55-aggiungere-limiti-di-risorse)
</quiz>

<quiz>
Qual è l'executor predefinito in Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Per approfondire: [4. Selezionare una piattaforma di esecuzione](#4-selezionare-una-piattaforma-di-esecuzione)
</quiz>

<quiz>
Come si specifica un file di parametri quando si esegue Nextflow?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Per approfondire: [1.3. Usare un file di parametri](#13-usare-un-file-di-parametri)
</quiz>

<quiz>
Per cosa possono essere usati i profili? (Selezioni tutte le risposte applicabili)
- [x] Definire impostazioni specifiche dell'infrastruttura
- [x] Impostare limiti di risorse per ambienti diversi
- [x] Fornire parametri di test
- [ ] Definire nuovi processi

Per approfondire: [6. Usare profili per passare tra configurazioni preimpostate](#6-usare-profili-per-passare-tra-configurazioni-preimpostate)
</quiz>

<quiz>
Come si specificano più profili in un singolo comando?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Per approfondire: [6. Usare profili per passare tra configurazioni preimpostate](#6-usare-profili-per-passare-tra-configurazioni-preimpostate)
</quiz>
