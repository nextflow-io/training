# Parte 2: Eseguire pipeline reali

Nella Parte 1 di questo corso (Eseguire Operazioni di Base), abbiamo iniziato con un workflow di esempio che aveva solo funzionalità minime per mantenere bassa la complessità del codice.
Ad esempio, `1-hello.nf` utilizzava un parametro da riga di comando (`--input`) per fornire un singolo valore alla volta.

Tuttavia, la maggior parte delle pipeline del mondo reale utilizza funzionalità più sofisticate per consentire l'elaborazione efficiente di grandi quantità di dati su larga scala e applicare più passaggi di elaborazione concatenati da logiche talvolta complesse.

In questa parte della formazione, dimostriamo le caratteristiche chiave delle pipeline del mondo reale provando versioni espanse della pipeline originale Hello World.

## 1. Elaborare dati di input da un file

In una pipeline del mondo reale, tipicamente vogliamo elaborare più punti dati (o serie di dati) contenuti in uno o più file di input.
E ovunque possibile, vogliamo eseguire l'elaborazione di dati indipendenti in parallelo, per ridurre il tempo di attesa per l'analisi.

Per dimostrare come Nextflow fa questo, abbiamo preparato un file CSV chiamato `greetings.csv` che contiene diversi saluti di input, imitando il tipo di dati colonnari che potreste voler elaborare in un'analisi dati reale.
Notate che i numeri non sono significativi, sono lì solo a scopo illustrativo.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Abbiamo anche scritto una versione migliorata del workflow originale, ora chiamata `2a-inputs.nf`, che leggerà il file CSV, estrarrà i saluti e scriverà ciascuno di essi in un file separato.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Eseguiamo prima il workflow, e poi daremo un'occhiata al codice Nextflow rilevante.

### 1.1. Eseguire il workflow

Eseguite il seguente comando nel vostro terminale.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Entusiasticamente, questo sembra indicare che sono state effettuate '3 di 3' chiamate per il processo, il che è incoraggiante, dato che c'erano tre righe di dati nel CSV che abbiamo fornito come input.
Questo suggerisce che il processo `sayHello()` è stato chiamato tre volte, una volta su ogni riga di input.

### 1.2. Trovare gli output pubblicati nella directory `results`

Diamo un'occhiata alla directory 'results' per vedere se il nostro workflow sta ancora scrivendo una copia dei nostri output lì.

??? abstract "Directory contents"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Sì! Vediamo una nuova directory chiamata `2a-inputs` con tre file di output con nomi diversi, abbastanza convenientemente.

Potete aprire ciascuno di essi per verificare che contengano la stringa di saluto appropriata.

??? abstract "File contents"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Questo conferma che ogni saluto nel file di input è stato elaborato appropriatamente.

### 1.3. Trovare gli output originali e i registri

Potreste aver notato che l'output della console sopra si riferiva a una sola directory di attività.
Significa che tutte e tre le chiamate a `sayHello()` sono state eseguite all'interno di quella singola directory di attività?

#### 1.3.1. Esaminare la directory di attività indicata nel terminale

Diamo un'occhiata all'interno di quella directory di attività `8e/0eb066`.

??? abstract "Directory contents"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Troviamo solo l'output corrispondente a uno dei saluti (così come i file accessori se abilitiamo la visualizzazione dei file nascosti).

Quindi cosa sta succedendo qui?

Per impostazione predefinita, il sistema di registrazione ANSI scrive le informazioni di stato per tutte le chiamate allo stesso processo sulla stessa riga.
Di conseguenza, ci ha mostrato solo uno dei tre percorsi delle directory di attività (`8e/0eb066`) nell'output della console.
Ce ne sono altri due che non sono elencati lì.

#### 1.3.2. Far mostrare più dettagli al terminale

Possiamo modificare il comportamento di registrazione per vedere l'elenco completo delle chiamate ai processi aggiungendo `-ansi-log false` al comando come segue:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Output del comando"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Questa volta vediamo tutte e tre le esecuzioni del processo e le loro sottodirectory di lavoro associate elencate nell'output.
Disabilitare la registrazione ANSI ha anche impedito a Nextflow di usare i colori nell'output del terminale.

Notate che il modo in cui viene riportato lo stato è un po' diverso tra le due modalità di registrazione.
Nella modalità condensata, Nextflow riporta se le chiamate sono state completate con successo o meno.
In questa modalità espansa, riporta solo che sono state inviate.

Questo conferma che il processo `sayHello()` viene chiamato tre volte e viene creata una directory di attività separata per ciascuna.

Se guardiamo all'interno di ciascuna delle directory di attività elencate lì, possiamo verificare che ciascuna corrisponde a uno dei saluti.

??? abstract "Directory contents"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Questo conferma che ogni chiamata al processo viene eseguita in isolamento da tutte le altre.
Questo ha molti vantaggi, incluso evitare collisioni se il processo produce file intermedi con nomi non univoci.

!!! tip

    Per un workflow complesso, o un gran numero di input, avere l'elenco completo nell'output del terminale potrebbe diventare un po' opprimente, quindi le persone normalmente non usano `-ansi-log false` nell'uso di routine.

### 1.4. Esaminare il codice del workflow

Quindi questa versione del workflow è in grado di leggere un file CSV di input, elaborare gli input separatamente e nominare gli output in modo univoco.

Diamo un'occhiata a cosa rende possibile questo nel codice del workflow.

??? full-code "File di codice completo"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Ancora una volta, non è necessario memorizzare la sintassi del codice, ma è bene imparare a riconoscere i componenti chiave del workflow che forniscono funzionalità importanti.

#### 1.4.1. Caricare i dati di input dal CSV

Questa è la parte più interessante: come siamo passati dal prendere un singolo valore dalla riga di comando, al prendere un file CSV, analizzarlo ed elaborare i singoli saluti che contiene?

In Nextflow, lo facciamo con un [**canale**](https://nextflow.io/docs/latest/channel.html): un costrutto di coda progettato per gestire gli input in modo efficiente e trasportarli da un passaggio all'altro in workflow multi-step, fornendo al contempo parallelismo integrato e molti vantaggi aggiuntivi.

Analizziamolo.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

Questo codice crea un canale chiamato `greeting_ch` che legge il file CSV, lo analizza ed estrae la prima colonna da ogni riga.
Il risultato è un canale contenente `Hello`, `Bonjour` e `Holà`.

??? tip "Come funziona questo?"

    Ecco cosa significa quella riga in italiano semplice:

    - `channel.fromPath` è una **fabbrica di canali** che crea un canale da percorsi di file
    - `(params.input)` specifica che il percorso del file è fornito da `--input` sulla riga di comando

    In altre parole, quella riga dice a Nextflow: prendi il percorso del file fornito con `--input` e preparati a trattare il suo contenuto come dati di input.

    Poi le due righe successive applicano **operatori** che effettuano l'analisi effettiva del file e il caricamento dei dati nella struttura dati appropriata:

    - `.splitCsv()` dice a Nextflow di analizzare il file CSV in un array che rappresenta righe e colonne
    - `.map { line -> line[0] }` dice a Nextflow di prendere solo l'elemento nella prima colonna da ogni riga

    Quindi in pratica, partendo dal seguente file CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Abbiamo trasformato questo in un array che appare così:

    ```txt title="Array contents"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    E poi abbiamo preso il primo elemento da ciascuna delle tre righe e li abbiamo caricati in un canale Nextflow che ora contiene: `Hello`, `Bonjour` e `Holà`.

    Se volete comprendere i canali e gli operatori in profondità, incluso come scriverli voi stessi, consultate [Hello Nextflow Parte 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Chiamare il processo su ogni saluto

Successivamente, nell'ultima riga del blocco `main:` del workflow, forniamo il canale `greeting_ch` caricato come input al processo `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

Questo dice a Nextflow di eseguire il processo individualmente su ogni elemento nel canale, _cioè_ su ogni saluto.
E poiché Nextflow è intelligente in questo modo, eseguirà queste chiamate ai processi in parallelo se possibile, a seconda dell'infrastruttura di calcolo disponibile.

È così che potete ottenere un'elaborazione efficiente e scalabile di molti dati (molti campioni, o punti dati, qualunque sia la vostra unità di ricerca) con relativamente pochissimo codice.

#### 1.4.3. Come vengono nominati gli output

Infine, vale la pena dare una rapida occhiata al codice del processo per vedere come otteniamo che i file di output siano nominati in modo univoco.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

Vedete che, rispetto alla versione di questo processo in `1-hello.nf`, la dichiarazione di output e il bit rilevante del comando sono cambiati per includere il valore del saluto nel nome del file di output.

Questo è un modo per garantire che i nomi dei file di output non entrino in collisione quando vengono pubblicati nella directory dei risultati comune.

E questo è l'unico cambiamento che abbiamo dovuto fare all'interno della dichiarazione del processo!

### Takeaway

Comprendete a un livello di base come i canali e gli operatori ci consentono di elaborare più input in modo efficiente.

### Cosa c'è dopo?

Scoprite come vengono costruiti i workflow multi-step e come operano.

---

## 2. Eseguire workflow multi-step

La maggior parte dei workflow del mondo reale coinvolge più di un passaggio.
Costruiamo su ciò che abbiamo appena imparato sui canali e vediamo come Nextflow utilizza canali e operatori per connettere i processi insieme in un workflow multi-step.

A tal fine, vi forniamo un workflow di esempio che concatena tre passaggi separati e dimostra quanto segue:

1. Far fluire i dati da un processo al successivo
2. Raccogliere gli output da più chiamate ai processi in una singola chiamata al processo

Nello specifico, abbiamo realizzato una versione espansa del workflow chiamata `2b-multistep.nf` che prende ogni saluto di input, lo converte in maiuscolo, quindi raccoglie tutti i saluti in maiuscolo in un singolo file di output.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Come in precedenza, eseguiremo prima il workflow e poi guarderemo il codice per vedere cosa c'è di nuovo.

### 2.1. Eseguire il workflow

Eseguite il seguente comando nel vostro terminale:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Output del comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Vedete che come promesso, sono stati eseguiti più passaggi come parte del workflow; i primi due (`sayHello` e `convertToUpper`) sono stati presumibilmente eseguiti su ogni singolo saluto, e il terzo (`collectGreetings`) sarà stato eseguito solo una volta, sugli output di tutte e tre le chiamate a `convertToUpper`.

### 2.2. Trovare gli output

Verifichiamo che sia effettivamente quello che è successo dando un'occhiata alla directory `results`.

??? abstract "Directory contents"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Come potete vedere, abbiamo una nuova directory chiamata `2b-multistep`, e contiene molti più file di prima.
Alcuni dei file sono stati raggruppati in una sottodirectory chiamata `intermediates`, mentre due file si trovano al livello superiore.

Questi due sono i risultati finali del workflow multi-step.
Prendetevi un minuto per guardare i nomi dei file e controllare i loro contenuti per confermare che siano ciò che vi aspettate.

??? abstract "File contents"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

Il primo contiene i nostri tre saluti, in maiuscolo e raccolti di nuovo in un singolo file come promesso.
Il secondo è un file di report che riassume alcune informazioni sull'esecuzione.

### 2.3. Esaminare il codice

Diamo un'occhiata al codice e identifichiamo i pattern chiave per i workflow multi-step.

??? full-code "File di codice completo"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

C'è molto in corso lì, ma la differenza più ovvia rispetto alla versione precedente del workflow è che ora ci sono più definizioni di processi e, di conseguenza, diverse chiamate ai processi nel blocco workflow.

Diamo un'occhiata più da vicino e vediamo se riusciamo a identificare i pezzi più interessanti.

#### 2.3.1. Visualizzare la struttura del workflow

Se state usando VSCode con l'estensione Nextflow, potete ottenere un diagramma utile di come i processi sono connessi cliccando sul piccolo link `DAG preview` visualizzato appena sopra il blocco workflow in qualsiasi script Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Questo vi dà una bella panoramica di come i processi sono connessi e cosa producono.

Vedete che oltre al processo originale `sayHello`, ora abbiamo anche `convertToUpper` e `collectGreetings`, che corrispondono ai nomi dei processi che abbiamo visto nell'output della console.
Le due nuove definizioni di processi sono strutturate allo stesso modo del processo `sayHello`, tranne che `collectGreetings` prende un parametro di input aggiuntivo chiamato `batch` e produce due output.

Non entreremo nei dettagli del codice per ciascuno, ma se siete curiosi, potete consultare i dettagli nella [Parte 2 di Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Per ora, approfondiamo come i processi sono connessi l'uno all'altro.

#### 2.3.2. Come i processi sono connessi

La cosa davvero interessante da guardare qui è come le chiamate ai processi sono concatenate nel blocco `main:` del workflow.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Potete vedere che la prima chiamata al processo, `sayHello(greeting_ch)`, è invariata.
Poi la chiamata al processo successivo, a `convertToUpper`, si riferisce all'output di `sayHello` come `sayHello.out`.

Il pattern è semplice: `processName.out` si riferisce al canale di output di un processo, che può essere passato direttamente al processo successivo.
È così che trasportiamo i dati da un passaggio all'altro in Nextflow.

#### 2.3.3. Un processo può prendere più input

La terza chiamata al processo, a `collectGreetings`, è un po' diversa.

```groovy title="2b-multistep.nf" linenums="77"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Vedete che questa chiamata riceve due input, `convertToUpper.out.collect()` e `params.batch`.
Ignorando il bit `.collect()` per ora, possiamo generalizzare questo come `collectGreetings(input1, input2)`.

Questo corrisponde alle due dichiarazioni di input nel modulo del processo:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Quando Nextflow analizza questo, assegnerà il primo input nella chiamata a `path input_files`, e il secondo a `val batch_name`.

Quindi ora sapete che un processo può prendere più input, e come appare la chiamata nel blocco workflow.

Ora diamo un'occhiata più da vicino a quel primo input, `convertToUpper.out.collect()`.

#### 2.3.4. Cosa fa `collect()` nella chiamata a `collectGreetings`

Per passare l'output di `sayHello` a `convertToUpper`, ci siamo semplicemente riferiti al canale di output di `sayHello` come `sayHello.out`. Ma per il passaggio successivo, stiamo vedendo un riferimento a `convertToUpper.out.collect()`.

Cos'è questo bit `.collect()` e cosa fa?

È un operatore, ovviamente. Proprio come gli operatori `splitCsv` e `map` che abbiamo incontrato prima.
Questa volta l'operatore si chiama `collect`, ed è applicato al canale di output prodotto da `convertToUpper`.

L'operatore `collect` viene utilizzato per raccogliere gli output da più chiamate allo stesso processo e impacchettarli in un singolo elemento del canale.

Nel contesto di questo workflow, sta prendendo i tre saluti in maiuscolo nel canale `convertToUpper.out` (che sono tre elementi separati del canale, e normalmente sarebbero gestiti in chiamate separate dal processo successivo) e impacchettandoli in un singolo elemento.
È così che riportiamo tutti i saluti nello stesso file.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

Al contrario, se non applicassimo `collect()` all'output di `convertToUpper()` prima di alimentarlo a `collectGreetings()`, Nextflow eseguirebbe semplicemente `collectGreetings()` indipendentemente su ogni saluto, il che non raggiungerebbe il nostro obiettivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Ci sono molti altri [operatori](https://nextflow.io/docs/latest/reference/operator.html) disponibili per applicare trasformazioni ai contenuti dei canali tra le chiamate ai processi.

Questo dà agli sviluppatori di pipeline molta flessibilità per personalizzare la logica di flusso della loro pipeline.
Lo svantaggio è che a volte può rendere più difficile decifrare cosa sta facendo la pipeline.

#### 2.3.5. Un parametro di input può avere un valore predefinito

Potreste aver notato che `collectGreetings` prende un secondo input, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Questo passa un parametro CLI chiamato `--batch` al workflow.
Tuttavia, quando abbiamo lanciato il workflow prima, non abbiamo specificato un parametro `--batch`.

Cosa sta succedendo lì?
Date un'occhiata al blocco `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

C'è un valore predefinito configurato nel workflow, quindi non dobbiamo fornirlo.
Ma se ne forniamo uno sulla riga di comando, il valore che specifichiamo verrà utilizzato al posto del predefinito.

Provatelo:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Output del comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Dovreste vedere nuovi output finali nominati con il vostro nome batch personalizzato.

??? abstract "Directory contents"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Questo è un aspetto della configurazione degli input, che tratteremo più in dettaglio nella Parte 3, ma per ora la cosa importante è sapere che i parametri di input possono avere valori predefiniti.

#### 2.3.6. Un processo può produrre più output

Nella definizione del processo `collectGreetings`, vediamo le seguenti dichiarazioni di output:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Che vengono poi riferiti dal nome dato con `emit:` nel blocco `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Questo rende facile passare poi output specifici individualmente ad altri processi nel workflow, in combinazione con vari operatori.

#### 2.3.7. Gli output pubblicati possono essere organizzati

Nel blocco `output`, abbiamo usato percorsi personalizzati per raggruppare i risultati intermedi al fine di rendere più facile individuare solo gli output finali del workflow.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Ci sono modi più sofisticati per organizzare gli output pubblicati; ne toccheremo alcuni nella parte sulla configurazione.

!!! tip "Volete saperne di più sulla costruzione di workflow?"

    Per una copertura dettagliata della costruzione di workflow multi-step, consultate [Hello Nextflow Parte 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Takeaway

Comprendete a un livello di base come i workflow multi-step sono costruiti usando canali e operatori e come operano.
Avete anche visto che i processi possono prendere più input e produrre più output, e che questi possono essere pubblicati in modo strutturato.

### Cosa c'è dopo?

Imparate come le pipeline Nextflow possono essere modularizzate per promuovere il riutilizzo del codice e la manutenibilità.

---

## 3. Eseguire pipeline modularizzate

Finora, tutti i workflow che abbiamo esaminato consistevano in un singolo file di workflow contenente tutto il codice rilevante.

Tuttavia, le pipeline del mondo reale tipicamente beneficiano dall'essere _modularizzate_, il che significa che il codice è suddiviso in file diversi.
Questo può rendere il loro sviluppo e manutenzione più efficienti e sostenibili.

Qui dimostreremo la forma più comune di modularità del codice in Nextflow, che è l'uso di **moduli**.

In Nextflow, un [**modulo**](https://nextflow.io/docs/latest/module.html) è una singola definizione di processo che è incapsulata da sola in un file di codice autonomo.
Per usare un modulo in un workflow, basta aggiungere un'istruzione di importazione di una sola riga al file di codice del workflow; poi potete integrare il processo nel workflow nello stesso modo in cui fareste normalmente.
Questo rende possibile riutilizzare le definizioni dei processi in più workflow senza produrre più copie del codice.

Fino ad ora abbiamo eseguito workflow che avevano tutti i loro processi inclusi in un file di codice monolitico.
Ora vedremo come appare quando i processi sono memorizzati in moduli individuali.

Abbiamo ovviamente ancora una volta preparato un workflow adatto a scopo dimostrativo, chiamato `2c-modules.nf`, insieme a un set di moduli situati nella directory `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Directory contents"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Vedete che ci sono quattro file Nextflow, ciascuno nominato dopo uno dei processi.
Potete ignorare il file `cowpy.nf` per ora; ci arriveremo più tardi.

### 3.1. Esaminare il codice

Questa volta guarderemo prima il codice.
Iniziate aprendo il file di workflow `2c-modules.nf`.

??? full-code "File di codice completo"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Vedete che la logica del workflow è esattamente la stessa della versione precedente del workflow.
Tuttavia, il codice del processo è scomparso dal file del workflow, e invece ci sono istruzioni `include` che puntano a file separati sotto `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Aprite uno di quei file e troverete il codice per il processo corrispondente.

??? full-code "File di codice completo"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

Come potete vedere, il codice del processo non è cambiato; è stato semplicemente copiato in un file di modulo individuale invece di essere nel file del workflow principale.
Lo stesso vale per gli altri due processi.

Quindi vediamo come appare eseguire questa nuova versione.

### 3.2. Eseguire il workflow

Eseguite questo comando nel vostro terminale, con il flag `-resume`:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Noterete che le esecuzioni dei processi sono state tutte memorizzate nella cache con successo, il che significa che Nextflow ha riconosciuto di aver già fatto il lavoro richiesto, anche se il codice è stato suddiviso e il file del workflow principale è stato rinominato.

Niente di tutto ciò importa a Nextflow; ciò che importa è lo script di lavoro che viene generato una volta che tutto il codice è stato assemblato e valutato.

!!! tip

    È anche possibile incapsulare una sezione di un workflow come 'subworkflow' che può essere importato in una pipeline più grande, ma questo è al di fuori dell'ambito di questo corso.

    Potete saperne di più sullo sviluppo di workflow componibili nella Side Quest su [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Takeaway

Sapete come i processi possono essere memorizzati in moduli autonomi per promuovere il riutilizzo del codice e migliorare la manutenibilità.

### Cosa c'è dopo?

Imparate a usare i container per gestire le dipendenze software.

---

## 4. Usare software containerizzato

Finora i workflow che abbiamo usato come esempi avevano solo bisogno di eseguire operazioni di elaborazione del testo molto basilari usando strumenti UNIX disponibili nel nostro ambiente.

Tuttavia, le pipeline del mondo reale tipicamente richiedono strumenti e pacchetti specializzati che non sono inclusi per impostazione predefinita nella maggior parte degli ambienti.
Di solito, dovreste installare questi strumenti, gestire le loro dipendenze e risolvere eventuali conflitti.

Tutto questo è molto noioso e fastidioso.
Un modo molto migliore per affrontare questo problema è usare i **container**.

Un **container** è un'unità di software leggera, autonoma ed eseguibile creata da un'**immagine** container che include tutto il necessario per eseguire un'applicazione, inclusi codice, librerie di sistema e impostazioni.

!!! Tip

    Insegniamo questo usando la tecnologia [Docker](https://www.docker.com/get-started/), ma Nextflow supporta anche diverse altre tecnologie di container.
    Potete saperne di più sul supporto di Nextflow per i container [qui](https://nextflow.io/docs/latest/container.html).

### 4.1. Usare un container direttamente

Prima, proviamo a interagire direttamente con un container.
Questo aiuterà a consolidare la vostra comprensione di cosa sono i container prima di iniziare a usarli in Nextflow.

#### 4.1.1. Scaricare l'immagine del container

Per usare un container, di solito scaricate o "pull" un'immagine container da un registro di container, e poi eseguite l'immagine del container per creare un'istanza del container.

La sintassi generale è la seguente:

```bash title="Syntax"
docker pull '<container>'
```

- `docker pull` è l'istruzione al sistema di container per scaricare un'immagine container da un repository.
- `'<container>'` è l'indirizzo URI dell'immagine del container.

Come esempio, scarichiamo un'immagine container che contiene [cowpy](https://github.com/jeffbuttars/cowpy), un'implementazione python di uno strumento chiamato `cowsay` che genera arte ASCII per visualizzare input di testo arbitrari in modo divertente.

Ci sono vari repository dove potete trovare container pubblicati.
Abbiamo usato il servizio [Seqera Containers](https://seqera.io/containers/) per generare questa immagine container Docker dal pacchetto Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Eseguite il comando pull completo:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Output del comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Questo dice al sistema di scaricare l'immagine specificata.
Una volta completato il download, avete una copia locale dell'immagine del container.

#### 4.1.2. Avviare il container

I container possono essere eseguiti come comando una tantum, ma potete anche usarli in modo interattivo, il che vi dà un prompt della shell all'interno del container e vi permette di giocare con il comando.

La sintassi generale è la seguente:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` è l'istruzione al sistema di container per avviare un'istanza del container da un'immagine container ed eseguire un comando in essa.
- `--rm` dice al sistema di arrestare l'istanza del container dopo che il comando è stato completato.

Completamente assemblato, il comando di esecuzione del container appare così:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Eseguite quel comando, e dovreste vedere il vostro prompt cambiare in qualcosa come `(base) root@b645838b3314:/tmp#`, che indica che ora siete all'interno del container.

Potete verificare questo eseguendo `ls` per elencare i contenuti della directory:

```bash
ls /
```

??? success "Output del comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Vedete che il filesystem all'interno del container è diverso dal filesystem sul vostro sistema host.

!!! Tip

    Quando eseguite un container, è isolato dal sistema host per impostazione predefinita.
    Questo significa che il container non può accedere a nessun file sul sistema host a meno che non lo permettiate esplicitamente specificando che volete montare un volume come parte del comando `docker run` usando la seguente sintassi:

    ```bash title="Syntax"
    -v <outside_path>:<inside_path>
    ```

    Questo stabilisce effettivamente un tunnel attraverso il muro del container che potete usare per accedere a quella parte del vostro filesystem.

    Questo è trattato più in dettaglio nella [Parte 5 di Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Eseguire lo strumento `cowpy`

Dall'interno del container, potete eseguire il comando `cowpy` direttamente.

```bash
cowpy "Hello Containers"
```

??? success "Output del comando"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Questo produce arte ASCII del personaggio mucca predefinito (o 'cowacter') con una bolla di dialogo contenente il testo che abbiamo specificato.

Ora che avete testato l'uso di base, potete provare a dargli alcuni parametri.
Ad esempio, la documentazione dello strumento dice che possiamo impostare il personaggio con `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Output del comando"

    ```console
    __________________
    < Hello Containers >
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

Questa volta l'output dell'arte ASCII mostra il pinguino Linux, Tux, perché abbiamo specificato il parametro `-c tux`.

Poiché siete all'interno del container, potete eseguire il comando cowpy tutte le volte che volete, variando i parametri di input, senza dovervi preoccupare di installare librerie sul vostro sistema stesso.

??? tip "Altri personaggi disponibili"

    Usate il flag '-c' per scegliere un personaggio diverso, inclusi:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Sentitevi liberi di giocare con questo.
Quando avete finito, uscite dal container usando il comando `exit`:

```bash
exit
```

Vi ritroverete nella vostra shell normale.

### 4.2. Usare un container in un workflow

Quando eseguiamo una pipeline, vogliamo essere in grado di dire a Nextflow quale container usare ad ogni passaggio, e soprattutto, vogliamo che gestisca tutto quel lavoro che abbiamo appena fatto: scaricare il container, avviarlo, eseguire il comando e smontare il container quando ha finito.

Buone notizie: è esattamente quello che Nextflow farà per noi.
Dobbiamo solo specificare un container per ogni processo.

Per dimostrare come funziona questo, abbiamo realizzato un'altra versione del nostro workflow che esegue `cowpy` sul file di saluti raccolti prodotto nel terzo passaggio.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Questo dovrebbe produrre un file contenente l'arte ASCII con i tre saluti nella bolla di dialogo.

#### 4.2.1. Esaminare il codice

Il workflow è molto simile al precedente, più il passaggio extra per eseguire `cowpy`.

??? full-code "File di codice completo"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Vedete che questo workflow importa un processo `cowpy` da un file di modulo, e lo chiama sull'output della chiamata a `collectGreetings()`, più un parametro di input chiamato `params.character`.

```groovy title="2d-container.nf" linenums="31"
// generate ASCII art of the greetings with cowpy
cowpy(collectGreetings.out.outfile, params.character)
```

Il processo `cowpy`, che avvolge il comando cowpy per generare arte ASCII, è definito nel modulo `cowpy.nf`.

??? full-code "File di codice completo"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

Il processo `cowpy` richiede due input: il percorso a un file di input contenente il testo da mettere nella bolla di dialogo (`input_file`), e un valore per la variabile character.

Importante, include anche la riga `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, che punta all'URI del container che abbiamo usato prima.

#### 4.2.2. Verificare che Docker sia abilitato nella configurazione

Anticiperemo leggermente la Parte 3 di questo corso di formazione introducendo il file di configurazione `nextflow.config`, che è uno dei modi principali che Nextflow offre per configurare l'esecuzione del workflow.
Quando un file chiamato `nextflow.config` è presente nella directory corrente, Nextflow lo caricherà automaticamente e applicherà qualsiasi configurazione contenga.

A tal fine, abbiamo incluso un file `nextflow.config` con una singola riga di codice che abilita Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Questa configurazione dice a Nextflow di usare Docker per qualsiasi processo che specifichi un container compatibile.

!!! tip

    È tecnicamente possibile abilitare l'esecuzione Docker dalla riga di comando, su base per-esecuzione, usando il parametro `-with-docker <container>`.
    Tuttavia, questo ci permette solo di specificare un container per l'intero workflow, mentre l'approccio che vi abbiamo appena mostrato ci permette di specificare un container diverso per processo.
    Quest'ultimo è molto migliore per modularità, manutenzione del codice e riproducibilità.

#### 4.2.3. Eseguire il workflow

Solo per ricapitolare, questo è ciò che stiamo per eseguire:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Pensate che funzionerà?

Eseguiamo il workflow con il flag `-resume`, e specifichiamo che vogliamo che il personaggio sia il tacchino.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

I primi tre passaggi sono stati memorizzati nella cache poiché li abbiamo già eseguiti prima, ma il processo `cowpy` è nuovo quindi viene effettivamente eseguito.

Potete trovare l'output del passaggio `cowpy` nella directory `results`.

??? abstract "File contents"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
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

Vedete che il personaggio sta dicendo tutti i saluti, poiché è stato eseguito sul file di saluti raccolti in maiuscolo.

Più importante, siamo stati in grado di eseguire questo come parte della nostra pipeline senza dover fare un'installazione corretta di cowpy e tutte le sue dipendenze.
E ora possiamo condividere la pipeline con i collaboratori e farla eseguire sulla loro infrastruttura senza che debbano installare nulla neanche loro, a parte Docker o una delle sue alternative (come Singularity/Apptainer) come menzionato sopra.

#### 4.2.4. Ispezionare come Nextflow ha lanciato l'attività containerizzata

Come coda finale a questa sezione, diamo un'occhiata alla sottodirectory di lavoro per una delle chiamate al processo `cowpy` per ottenere un po' più di informazioni su come Nextflow lavora con i container sotto il cofano.

Controllate l'output dal vostro comando `nextflow run` per trovare il percorso alla sottodirectory di lavoro per il processo `cowpy`.
Guardando ciò che abbiamo ottenuto per l'esecuzione mostrata sopra, la riga di log della console per il processo `cowpy` inizia con `[7f/caf718]`.
Questo corrisponde al seguente percorso di directory troncato: `work/7f/caf718`.

In quella directory, troverete il file `.command.run` che contiene tutti i comandi che Nextflow ha eseguito per vostro conto nel corso dell'esecuzione della pipeline.

??? abstract "File contents"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Se cercate `nxf_launch` in questo file, dovreste vedere qualcosa del genere:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Questo comando di lancio mostra che Nextflow sta usando un comando `docker run` molto simile per lanciare la chiamata al processo come abbiamo fatto quando l'abbiamo eseguito manualmente.
Monta anche la sottodirectory di lavoro corrispondente nel container, imposta la directory di lavoro all'interno del container di conseguenza ed esegue il nostro script bash templato nel file `.command.sh`.

Questo conferma che tutto il duro lavoro che abbiamo dovuto fare manualmente nella sezione precedente è ora fatto per noi da Nextflow!

### Takeaway

Comprendete quale ruolo giocano i container nella gestione delle versioni degli strumenti software e nel garantire la riproducibilità.

Più in generale, avete una comprensione di base di quali sono i componenti principali delle pipeline Nextflow del mondo reale e come sono organizzati.
Conoscete i fondamenti di come Nextflow può elaborare più input in modo efficiente, eseguire workflow composti da più passaggi connessi insieme, sfruttare componenti di codice modulari e utilizzare container per una maggiore riproducibilità e portabilità.

### Cosa c'è dopo?

Prendetevi un'altra pausa! È stata una grande quantità di informazioni su come funzionano le pipeline Nextflow.

Nell'ultima sezione di questa formazione, approfondiremo l'argomento della configurazione.
Imparerete come configurare l'esecuzione della vostra pipeline per adattarla alla vostra infrastruttura e gestire la configurazione di input e parametri.

---

## Quiz

<quiz>
Perché Nextflow crea una directory di attività separata per ogni chiamata al processo?
- [ ] Per migliorare la velocità di esecuzione
- [ ] Per ridurre l'uso della memoria
- [x] Per isolare le esecuzioni ed evitare collisioni tra gli output
- [ ] Per abilitare la compressione parallela dei file

Approfondisci: [1.3. Trovare gli output originali e i registri](#13-find-the-original-outputs-and-logs)
</quiz>

<quiz>
Cosa fa l'opzione `-ansi-log false` quando si esegue un workflow?
- [ ] Disabilita tutto l'output della console
- [x] Rimuove il colore dall'output
- [x] Mostra tutti i percorsi delle directory di attività invece di condensarli su una riga
- [ ] Abilita la modalità di debug verbosa

Approfondisci: [1.3.2. Far mostrare più dettagli al terminale](#132-make-the-terminal-show-more-details)

Potete anche usare una delle seguenti variabili d'ambiente se preferite questo stile:

```bash
export NXF_ANSI_LOG=0
# oppure
export NO_COLOR=1
```

</quiz>

<quiz>
Nel codice `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, cosa fa `#!groovy .map { line -> line[0] }`?
- [ ] Filtra le righe vuote
- [ ] Ordina le righe alfabeticamente
- [x] Estrae la prima colonna da ogni riga CSV
- [ ] Conta il numero di righe

Approfondisci: [1.4.1. Caricare i dati di input dal CSV](#141-loading-the-input-data-from-the-csv)
</quiz>

<quiz>
Perché è importante includere il valore di input nei nomi dei file di output (ad es., `#!groovy "${greeting}-output.txt"`)?
- [ ] Per migliorare la velocità di elaborazione
- [ ] Per abilitare la funzionalità di resume
- [x] Per impedire che i file di output si sovrascrivano a vicenda quando si elaborano più input
- [ ] Per rendere i file più facili da comprimere

Approfondisci: [1.4.3. Come vengono nominati gli output](#143-how-the-outputs-are-named)
</quiz>

<quiz>
Qual è lo scopo dell'istruzione `include` in un workflow modularizzato?
- [ ] Copiare il codice del processo nel file del workflow
- [x] Importare una definizione di processo da un file di modulo esterno
- [ ] Includere impostazioni di configurazione
- [ ] Aggiungere commenti di documentazione

Approfondisci: [3. Eseguire pipeline modularizzate](#3-running-modularized-pipelines)
</quiz>

<quiz>
Quando modularizzate un workflow e lo eseguite con `-resume`, cosa succede?
- [ ] La cache è disabilitata per i processi modulari
- [ ] Tutte le attività devono essere rieseguite
- [x] La cache funziona normalmente in base agli script di lavoro generati
- [ ] Solo il file del workflow principale viene memorizzato nella cache

Approfondisci: [3.2. Eseguire il workflow](#32-run-the-workflow)
</quiz>

<quiz>
Cosa specifica la direttiva `container` in una definizione di processo?
- [ ] La directory di lavoro per il processo
- [ ] L'allocazione massima di memoria
- [x] L'URI dell'immagine del container da usare per eseguire il processo
- [ ] Il formato del file di output

Approfondisci: [4.2. Usare un container in un workflow](#42-use-a-container-in-a-workflow)
</quiz>

<quiz>
Nel file `.command.run`, cosa contiene la funzione `nxf_launch`?
- [ ] Le informazioni sulla versione di Nextflow
- [ ] I parametri del workflow
- [x] Il comando `docker run` con i mount dei volumi e le impostazioni del container
- [ ] Le dichiarazioni di input del processo

Approfondisci: [4.2.4. Ispezionare come Nextflow ha lanciato l'attività containerizzata](#424-inspect-how-nextflow-launched-the-containerized-task)
</quiz>

<quiz>
Cosa gestisce automaticamente Nextflow quando esegue un processo containerizzato? (Selezionate tutte le opzioni applicabili)
- [x] Scaricare l'immagine del container se necessario
- [x] Montare la directory di lavoro nel container
- [x] Eseguire lo script del processo all'interno del container
- [x] Pulire l'istanza del container dopo l'esecuzione

Approfondisci: [4. Usare software containerizzato](#4-using-containerized-software)
</quiz>
