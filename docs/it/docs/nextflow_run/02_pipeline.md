# Parte 2: Eseguire pipeline reali

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella Parte 1 di questo corso (Operazioni base), abbiamo iniziato con un workflow di esempio che aveva solo funzionalità minime per mantenere bassa la complessità del codice.
Per esempio, `1-hello.nf` usava un parametro da riga di comando (`--input`) per fornire un singolo valore alla volta.

Tuttavia, la maggior parte delle pipeline del mondo reale usa funzionalità più sofisticate per consentire l'elaborazione efficiente di grandi quantità di dati su scala, e applicare step di elaborazione multipli concatenati insieme da logiche a volte complesse.

In questa parte della formazione, dimostriamo le funzionalità chiave delle pipeline del mondo reale provando versioni espanse della pipeline Hello World originale.

## 1. Elaborare dati di input da un file

In una pipeline del mondo reale, tipicamente vogliamo elaborare punti dati multipli (o serie di dati) contenuti in uno o più file di input.
E dove possibile, vogliamo eseguire l'elaborazione di dati indipendenti in parallelo, per ridurre il tempo di attesa per l'analisi.

Per dimostrare come Nextflow fa questo, abbiamo preparato un file CSV chiamato `greetings.csv` che contiene diversi saluti di input, simulando il tipo di dati in colonne che potresti voler elaborare in un'analisi dati reale.
Nota che i numeri non sono significativi, sono lì solo a scopo illustrativo.

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

### 1.1. Esegui il workflow

Esegui il seguente comando nel tuo terminale.

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

Eccitante, questo sembra indicare che sono state fatte '3 of 3' chiamate per il process, il che è incoraggiante, dato che c'erano tre righe di dati nel CSV che abbiamo fornito come input.
Questo suggerisce che il process `sayHello()` è stato chiamato tre volte, una volta per ogni riga di input.

### 1.2. Trova gli output pubblicati nella directory `results`

Guardiamo la directory 'results' per vedere se il nostro workflow sta ancora scrivendo una copia dei nostri output lì.

??? abstract "Contenuti della directory"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Sì! Vediamo una nuova directory chiamata `2a-inputs` con tre file di output con nomi diversi, convenientemente.

Puoi aprire ciascuno di essi per verificare che contengano la stringa di saluto appropriata.

??? abstract "Contenuto dei file"

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

### 1.3. Trova gli output originali e i log

Potresti aver notato che l'output della console sopra si riferiva a una sola directory di attività.
Significa che tutte e tre le chiamate a `sayHello()` sono state eseguite all'interno di quella singola directory di attività?

#### 1.3.1. Esamina la directory di attività data nel terminale

Diamo un'occhiata dentro quella directory di attività `8e/0eb066`.

??? abstract "Contenuti della directory"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Troviamo solo l'output corrispondente a uno dei saluti (così come i file accessori se abilitiamo la visualizzazione dei file nascosti).

Quindi cosa sta succedendo?

Per impostazione predefinita, il sistema di logging ANSI scrive le informazioni di stato per tutte le chiamate allo stesso process sulla stessa riga.
Di conseguenza, ci ha mostrato solo uno dei tre percorsi delle directory di attività (`8e/0eb066`) nell'output della console.
Ce ne sono altri due che non sono elencati lì.

#### 1.3.2. Fai mostrare al terminale più dettagli

Possiamo modificare il comportamento del logging per vedere la lista completa delle chiamate ai process aggiungendo `-ansi-log false` al comando come segue:

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

Questa volta vediamo tutti e tre i process eseguiti e le loro sottodirectory di lavoro associate elencate nell'output.
Disabilitare il logging ANSI ha anche impedito a Nextflow di usare i colori nell'output del terminale.

Nota che il modo in cui viene riportato lo stato è un po' diverso tra le due modalità di logging.
Nella modalità condensata, Nextflow riporta se le chiamate sono state completate con successo o meno.
In questa modalità espansa, riporta solo che sono state sottomesse.

Questo conferma che il process `sayHello()` viene chiamato tre volte, e una directory di attività separata viene creata per ciascuna.

Se guardiamo dentro ciascuna delle directory di attività elencate lì, possiamo verificare che ognuna corrisponde a uno dei saluti.

??? abstract "Contenuti della directory"

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

Questo conferma che ogni chiamata al process viene eseguita in isolamento da tutte le altre.
Questo ha molti vantaggi, incluso evitare collisioni se il process produce file intermedi con nomi non univoci.

!!! tip "Suggerimento"

    Per un workflow complesso, o un grande numero di input, avere la lista completa stampata nel terminale potrebbe diventare un po' opprimente, quindi le persone normalmente non usano `-ansi-log false` nell'utilizzo di routine.

### 1.4. Esamina il codice del workflow

Quindi questa versione del workflow è capace di leggere un file CSV di input, elaborare gli input separatamente, e nominare gli output in modo univoco.

Diamo un'occhiata a cosa rende questo possibile nel codice del workflow.

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
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
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

Ancora una volta, non devi memorizzare la sintassi del codice, ma è bene imparare a riconoscere i componenti chiave del workflow che forniscono funzionalità importanti.

#### 1.4.1. Caricare i dati di input dal CSV

Questa è la parte più interessante: come siamo passati dal prendere un singolo valore dalla riga di comando, al prendere un file CSV, analizzarlo ed elaborare i singoli saluti che contiene?

In Nextflow, lo facciamo con un [**canale**](https://nextflow.io/docs/latest/channel.html): un costrutto progettato per gestire gli input efficientemente e trasferirli da uno step all'altro in workflow multi-step, fornendo al contempo parallelismo integrato e molti altri benefici.

Analizziamolo.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // crea un canale per gli input da un file CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emette un saluto
    sayHello(greeting_ch)
```

Questo codice crea un canale chiamato `greeting_ch` che legge il file CSV, lo analizza, e estrae la prima colonna da ogni riga.
Il risultato è un canale contenente `Hello`, `Bonjour`, e `Holà`.

??? tip "Come funziona?"

    Ecco cosa significa quella riga in italiano semplice:

    - `channel.fromPath` è una **fabbrica di canali** che crea un canale da percorsi di file
    - `(params.input)` specifica che il percorso del file è fornito da `--input` sulla riga di comando

    In altre parole, quella riga dice a Nextflow: prendi il percorso del file dato con `--input` e preparati a trattare i suoi contenuti come dati di input.

    Poi le due righe successive applicano **operatori** che fanno l'effettivo parsing del file e il caricamento dei dati nella struttura dati appropriata:

    - `.splitCsv()` dice a Nextflow di analizzare il file CSV in un array che rappresenta righe e colonne
    - `.map { line -> line[0] }` dice a Nextflow di prendere solo l'elemento nella prima colonna da ogni riga

    Quindi in pratica, partendo dal seguente file CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Abbiamo trasformato quello in un array che appare così:

    ```txt title="Contenuti dell'array"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    E poi abbiamo preso il primo elemento da ciascuna delle tre righe e li abbiamo caricati in un canale Nextflow che ora contiene: `Hello`, `Bonjour`, e `Holà`.

    Se vuoi capire i canali e gli operatori in profondità, incluso come scriverli tu stesso, vedi [Hello Nextflow Parte 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Chiama il process su ogni saluto

Successivamente, nell'ultima riga del blocco `main:` del workflow, forniamo il canale `greeting_ch` caricato come input al process `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // crea un canale per gli input da un file CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emette un saluto
    sayHello(greeting_ch)
```

Questo dice a Nextflow di eseguire il process individualmente su ogni elemento nel canale, _cioè_ su ogni saluto.
E poiché Nextflow è intelligente così, eseguirà queste chiamate al process in parallelo se possibile, a seconda dell'infrastruttura di calcolo disponibile.

È così che puoi ottenere un'elaborazione efficiente e scalabile di molti dati (molti campioni, o punti dati, qualunque sia la tua unità di ricerca) con relativamente pochissimo codice.

#### 1.4.3. Come vengono nominati gli output

Infine, vale la pena dare una rapida occhiata al codice del process per vedere come facciamo a nominare i file di output in modo univoco.

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

Vedi che, rispetto alla versione di questo process in `1-hello.nf`, la dichiarazione dell'output e la parte rilevante del comando sono cambiate per includere il valore del saluto nel nome del file di output.

Questo è un modo per assicurarsi che i nomi dei file di output non collideranno quando vengono pubblicati nella directory dei risultati comune.

E questa è l'unica modifica che abbiamo dovuto fare dentro la dichiarazione del process!

### Takeaway

Capisci a livello base come i canali e gli operatori ci permettono di elaborare input multipli efficientemente.

### Cosa c'è dopo?

Scopri come sono costruiti i workflow multi-step e come operano.

---

## 2. Eseguire workflow multi-step

La maggior parte dei workflow del mondo reale coinvolge più di uno step.
Costruiamo su ciò che abbiamo appena imparato sui canali, e guardiamo come Nextflow usa canali e operatori per connettere i process insieme in un workflow multi-step.

A tal fine, ti forniamo un workflow di esempio che concatena tre step separati e dimostra quanto segue:

1. Far fluire i dati da un process al successivo
2. Raccogliere output da chiamate multiple al process in una singola chiamata al process

Specificamente, abbiamo fatto una versione espansa del workflow chiamata `2b-multistep.nf` che prende ogni saluto di input, lo converte in maiuscolo, poi raccoglie tutti i saluti maiuscoli in un singolo file di output.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Come in precedenza, eseguiremo prima il workflow poi guarderemo il codice per vedere cosa c'è di nuovo.

### 2.1. Esegui il workflow

Esegui il seguente comando nel tuo terminale:

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

Vedi che come promesso, step multipli sono stati eseguiti come parte del workflow; i primi due (`sayHello` e `convertToUpper`) sono stati presumibilmente eseguiti su ogni singolo saluto, e il terzo (`collectGreetings`) sarà stato eseguito solo una volta, sugli output di tutte e tre le chiamate di `convertToUpper`.

### 2.2. Trova gli output

Verifichiamo che sia effettivamente ciò che è successo dando un'occhiata nella directory `results`.

??? abstract "Contenuti della directory"

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

Come puoi vedere, abbiamo una nuova directory chiamata `2b-multistep`, e contiene molti più file di prima.
Alcuni dei file sono stati raggruppati in una sottodirectory chiamata `intermediates`, mentre due file si trovano al livello superiore.

Questi due sono i risultati finali del workflow multi-step.
Prenditi un minuto per guardare i nomi dei file e controllare i loro contenuti per confermare che sono ciò che ti aspetti.

??? abstract "Contenuto dei file"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

Il primo contiene i nostri tre saluti, convertiti in maiuscolo e raccolti di nuovo in un singolo file come promesso.
Il secondo è un file di report che riassume alcune informazioni sull'esecuzione.

### 2.3. Esamina il codice

Guardiamo il codice e identifichiamo i pattern chiave per i workflow multi-step.

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
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
        sayHello(greeting_ch)
        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)
        // raccoglie tutti i saluti in un file
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

C'è molto da vedere lì dentro, ma la differenza più ovvia rispetto alla versione precedente del workflow è che ora ci sono definizioni di process multiple, e corrispondentemente, diverse chiamate ai process nel blocco workflow.

Diamo uno sguardo più da vicino e vediamo se possiamo identificare i pezzi più interessanti.

#### 2.3.1. Visualizzare la struttura del workflow

Se stai usando VSCode con l'estensione Nextflow, puoi ottenere un diagramma utile di come i process sono connessi cliccando sul piccolo link `DAG preview` visualizzato appena sopra il blocco workflow in qualsiasi script Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Questo ti dà una bella panoramica di come i process sono connessi e cosa producono.

Vedi che oltre al process originale `sayHello`, ora abbiamo anche `convertToUpper` e `collectGreetings`, che corrispondono ai nomi dei process che abbiamo visto nell'output della console.
Le due nuove definizioni di process sono strutturate nello stesso modo del process `sayHello`, tranne che `collectGreetings` prende un parametro di input aggiuntivo chiamato `batch` e produce due output.

Non entreremo nel codice per ciascuno in dettaglio, ma se sei curioso, puoi cercare i dettagli nella [Parte 2 di Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Per ora, approfondiamo come i process sono connessi l'uno all'altro.

#### 2.3.2. Come sono connessi i process

La cosa veramente interessante da guardare qui è come le chiamate ai process sono concatenate insieme nel blocco `main:` del workflow.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // crea un canale per gli input da un file CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emette un saluto
    sayHello(greeting_ch)
    // converte il saluto in maiuscolo
    convertToUpper(sayHello.out)
    // raccoglie tutti i saluti in un file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Puoi vedere che la prima chiamata al process, `sayHello(greeting_ch)`, è invariata.
Poi la chiamata successiva al process, a `convertToUpper`, si riferisce all'output di `sayHello` come `sayHello.out`.

Il pattern è semplice: `processName.out` si riferisce al canale di output di un process, che può essere passato direttamente al process successivo.
È così che trasferiamo i dati da uno step al successivo in Nextflow.

#### 2.3.3. Un process può prendere input multipli

La terza chiamata al process, a `collectGreetings`, è un po' diversa.

```groovy title="2b-multistep.nf" linenums="77"
    // raccoglie tutti i saluti in un file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Vedi che a questa chiamata vengono dati due input, `convertToUpper.out.collect()` e `params.batch`.
Ignorando la parte `.collect()` per ora, possiamo generalizzare questo come `collectGreetings(input1, input2)`.

Questo corrisponde alle due dichiarazioni di input nel modulo del process:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Quando Nextflow analizza questo, assegnerà il primo input nella chiamata a `path input_files`, e il secondo a `val batch_name`.

Quindi ora sai che un process può prendere input multipli, e come appare la chiamata nel blocco workflow.

Ora diamo uno sguardo più da vicino a quel primo input, `convertToUpper.out.collect()`.

#### 2.3.4. Cosa fa `collect()` nella chiamata `collectGreetings`

Per passare l'output di `sayHello` a `convertToUpper`, abbiamo semplicemente fatto riferimento al canale di output di `sayHello` come `sayHello.out`. Ma per lo step successivo, stiamo vedendo un riferimento a `convertToUpper.out.collect()`.

Cos'è questa parte `collect()` e cosa fa?

È un operatore, naturalmente. Proprio come gli operatori `splitCsv` e `map` che abbiamo incontrato prima.
Questa volta l'operatore si chiama `collect`, ed è applicato al canale di output prodotto da `convertToUpper`.

L'operatore `collect` è usato per raccogliere gli output da chiamate multiple allo stesso process e impacchettarli in un singolo elemento del canale.

Nel contesto di questo workflow, sta prendendo i tre saluti maiuscoli nel canale `convertToUpper.out` (che sono tre elementi separati del canale, e normalmente verrebbero gestiti in chiamate separate dal process successivo) e li impacchetta in un singolo elemento.
È così che otteniamo tutti i saluti di nuovo nello stesso file.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

Al contrario, se non applicassimo `collect()` all'output di `convertToUpper()` prima di passarlo a `collectGreetings()`, Nextflow eseguirebbe semplicemente `collectGreetings()` indipendentemente su ogni saluto, il che non raggiungerebbe il nostro obiettivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Ci sono molti altri [operatori](https://www.nextflow.io/docs/latest/reference/operator.html#operator-page) disponibili per applicare trasformazioni ai contenuti dei canali tra le chiamate ai process.

Questo dà agli sviluppatori di pipeline molta flessibilità per personalizzare la logica di flusso della loro pipeline.
Il lato negativo è che a volte può rendere più difficile decifrare cosa sta facendo la pipeline.

#### 2.3.5. Un parametro di input può avere un valore predefinito

Potresti aver notato che `collectGreetings` prende un secondo input, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // raccoglie tutti i saluti in un file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Questo passa un parametro CLI chiamato `--batch` al workflow.
Tuttavia, quando abbiamo lanciato il workflow prima, non abbiamo specificato un parametro `--batch`.

Cosa sta succedendo?
Dai un'occhiata al blocco `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

C'è un valore predefinito configurato nel workflow, quindi non dobbiamo fornirlo.
Ma se ne forniamo uno sulla riga di comando, il valore che specifichiamo verrà usato al posto di quello predefinito.

Provalo:

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

Dovresti vedere nuovi output finali nominati con il tuo nome batch personalizzato.

??? abstract "Contenuti della directory"

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

Questo è un aspetto della configurazione degli input, che copriremo più in dettaglio nella Parte 3, ma per ora la cosa importante è sapere che i parametri di input possono avere valori predefiniti.

#### 2.3.6. Un process può produrre output multipli

Nella definizione del process `collectGreetings`, vediamo le seguenti dichiarazioni di output:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

A cui poi si fa riferimento con il nome dato con `emit:` nel blocco `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Questo rende facile passare output specifici individualmente ad altri process nel workflow, in combinazione con vari operatori.

#### 2.3.7. Gli output pubblicati possono essere organizzati

Nel blocco `output`, abbiamo usato percorsi personalizzati per raggruppare i risultati intermedi in modo da rendere più facile individuare solo gli output finali del workflow.

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

!!! tip "Vuoi saperne di più sulla costruzione di workflow?"

    Per una copertura dettagliata sulla costruzione di workflow multi-step, vedi [Hello Nextflow Parte 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Takeaway

Capisci a livello base come i workflow multi-step sono costruiti usando canali e operatori e come operano.
Hai anche visto che i process possono prendere input multipli e produrre output multipli, e che questi possono essere pubblicati in modo strutturato.

### Cosa c'è dopo?

Impara come le pipeline Nextflow possono essere modularizzate per promuovere il riuso del codice e la manutenibilità.

---

## 3. Eseguire pipeline modularizzate

Finora, tutti i workflow che abbiamo guardato consistevano in un singolo file di workflow contenente tutto il codice rilevante.

Tuttavia, le pipeline del mondo reale tipicamente beneficiano dall'essere _modularizzate_, il che significa che il codice è diviso in file diversi.
Questo può rendere il loro sviluppo e manutenzione più efficienti e sostenibili.

Qui dimostreremo la forma più comune di modularità del codice in Nextflow, che è l'uso dei **moduli**.

In Nextflow, un [**modulo**](https://nextflow.io/docs/latest/module.html) è una singola definizione di process che è incapsulata da sola in un file di codice standalone.
Per usare un modulo in un workflow, aggiungi semplicemente una dichiarazione import di una riga al tuo file di codice del workflow; poi puoi integrare il process nel workflow nello stesso modo in cui lo faresti normalmente.
Questo rende possibile riutilizzare le definizioni dei process in workflow multipli senza produrre copie multiple del codice.

Fino ad ora abbiamo eseguito workflow che avevano tutti i loro process inclusi in un file di codice monolitico.
Ora vedremo come appare quando i process sono memorizzati in moduli individuali.

Abbiamo ovviamente preparato ancora una volta un workflow adatto per scopi dimostrativi, chiamato `2c-modules.nf`, insieme a un set di moduli situati nella directory `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Contenuti della directory"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Vedi che ci sono quattro file Nextflow, ciascuno nominato dopo uno dei process.
Puoi ignorare il file `cowpy.nf` per ora; ci arriveremo dopo.

### 3.1. Esamina il codice

Questa volta guarderemo prima il codice.
Inizia aprendo il file di workflow `2c-modules.nf`.

??? full-code "File di codice completo"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Include i moduli
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
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
        sayHello(greeting_ch)
        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)
        // raccoglie tutti i saluti in un file
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

Vedi che la logica del workflow è esattamente la stessa della versione precedente del workflow.
Tuttavia, il codice del process non c'è più nel file del workflow, e invece ci sono dichiarazioni `include` che puntano a file separati sotto `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Include i moduli
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Apri uno di quei file e troverai il codice per il process corrispondente.

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

Come puoi vedere, il codice del process non è cambiato; è stato semplicemente copiato in un file modulo individuale invece di essere nel file del workflow principale.
Lo stesso vale per gli altri due process.

Quindi vediamo come appare eseguire questa nuova versione.

### 3.2. Esegui il workflow

Esegui questo comando nel tuo terminale, con il flag `-resume`:

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

Noterai che le esecuzioni dei process sono state tutte salvate nella cache con successo, il che significa che Nextflow ha riconosciuto che ha già fatto il lavoro richiesto, anche se il codice è stato diviso e il file del workflow principale è stato rinominato.

Niente di tutto ciò importa a Nextflow; ciò che importa è lo script del job che viene generato una volta che tutto il codice è stato messo insieme e valutato.

!!! tip "Suggerimento"

    È anche possibile incapsulare una sezione di un workflow come un 'subworkflow' che può essere importato in una pipeline più grande, ma questo è al di fuori dello scopo di questo corso.

    Puoi saperne di più sullo sviluppo di workflow componibili nella Side Quest su [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Takeaway

Sai come i process possono essere memorizzati in moduli standalone per promuovere il riuso del codice e migliorare la manutenibilità.

### Cosa c'è dopo?

Impara a usare i container per gestire le dipendenze software.

---

## 4. Usare software containerizzato

Finora i workflow che abbiamo usato come esempi avevano solo bisogno di eseguire operazioni di elaborazione di testo molto basilari usando strumenti UNIX disponibili nel nostro ambiente.

Tuttavia, le pipeline del mondo reale tipicamente richiedono strumenti e pacchetti specializzati che non sono inclusi di default nella maggior parte degli ambienti.
Di solito, dovresti installare questi strumenti, gestire le loro dipendenze, e risolvere eventuali conflitti.

Tutto questo è molto noioso e fastidioso.
Un modo molto migliore per affrontare questo problema è usare i **container**.

Un **container** è un'unità software leggera, standalone ed eseguibile creata da un'**immagine** container che include tutto il necessario per eseguire un'applicazione incluso codice, librerie di sistema e impostazioni.

!!! Tip "Suggerimento"

    Insegniamo questo usando la tecnologia [Docker](https://www.docker.com/get-started/), ma Nextflow supporta diverse altre tecnologie container altrettanto bene.
    Puoi saperne di più sul supporto di Nextflow per i container [qui](https://www.nextflow.io/docs/latest/container.html#).

### 4.1. Usa un container direttamente

Prima, proviamo a interagire direttamente con un container.
Questo aiuterà a consolidare la tua comprensione di cosa sono i container prima di iniziare a usarli in Nextflow.

#### 4.1.1. Scarica l'immagine container

Per usare un container, di solito scarichi o "pulli" un'immagine container da un registry di container, e poi esegui l'immagine container per creare un'istanza del container.

La sintassi generale è la seguente:

```bash title="Sintassi"
docker pull '<container>'
```

- `docker pull` è l'istruzione al sistema container per scaricare un'immagine container da un repository.
- `'<container>'` è l'indirizzo URI dell'immagine container.

Come esempio, scarichiamo un'immagine container che contiene [cowpy](https://github.com/jeffbuttars/cowpy), un'implementazione python di uno strumento chiamato `cowsay` che genera arte ASCII per visualizzare input di testo arbitrari in modo divertente.

Ci sono vari repository dove puoi trovare container pubblicati.
Abbiamo usato il servizio [Seqera Containers](https://seqera.io/containers/) per generare questa immagine container Docker dal pacchetto Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Esegui il comando pull completo:

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
Una volta completato il download, hai una copia locale dell'immagine container.

#### 4.1.2. Avvia il container

I container possono essere eseguiti come comando singolo, ma puoi anche usarli interattivamente, il che ti dà un prompt shell all'interno del container e ti permette di giocare con il comando.

La sintassi generale è la seguente:

```bash title="Sintassi"
docker run --rm '<container>' [comando dello strumento]
```

- `docker run --rm '<container>'` è l'istruzione al sistema container per avviare un'istanza del container da un'immagine container ed eseguire un comando al suo interno.
- `--rm` dice al sistema di spegnere l'istanza del container dopo che il comando è stato completato.

Completamente assemblato, il comando di esecuzione del container appare così:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Esegui quel comando, e dovresti vedere il tuo prompt cambiare in qualcosa come `(base) root@b645838b3314:/tmp#`, che indica che sei ora dentro il container.

Puoi verificare questo eseguendo `ls` per elencare i contenuti della directory:

```bash
ls /
```

??? success "Output del comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Vedi che il filesystem dentro il container è diverso dal filesystem sul tuo sistema host.

!!! Tip "Suggerimento"

    Quando esegui un container, è isolato dal sistema host di default.
    Questo significa che il container non può accedere a nessun file sul sistema host a meno che tu non gli permetta esplicitamente di farlo specificando che vuoi montare un volume come parte del comando `docker run` usando la seguente sintassi:

    ```bash title="Sintassi"
    -v <percorso_esterno>:<percorso_interno>
    ```

    Questo effettivamente stabilisce un tunnel attraverso la parete del container che puoi usare per accedere a quella parte del tuo filesystem.

    Questo è coperto più in dettaglio nella [Parte 5 di Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Esegui lo strumento `cowpy`

Dall'interno del container, puoi eseguire il comando `cowpy` direttamente.

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

Questo produce arte ASCII della mucca predefinita (o 'cowacter') con un fumetto contenente il testo che abbiamo specificato.

Ora che hai testato l'utilizzo base, puoi provare a dargli alcuni parametri.
Per esempio, la documentazione dello strumento dice che possiamo impostare il personaggio con `-c`.

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

Dato che sei dentro il container, puoi eseguire il comando cowpy quante volte vuoi, variando i parametri di input, senza doverti preoccupare di installare librerie sul tuo sistema stesso.

??? tip "Altri personaggi disponibili"

    Usa il flag '-c' per scegliere un personaggio diverso, inclusi:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Sentiti libero di giocare con questo.
Quando hai finito, esci dal container usando il comando `exit`:

```bash
exit
```

Ti ritroverai di nuovo nella tua shell normale.

### 4.2. Usa un container in un workflow

Quando eseguiamo una pipeline, vogliamo essere in grado di dire a Nextflow quale container usare a ogni step, e importante, vogliamo che gestisca tutto quel lavoro che abbiamo appena fatto: scaricare il container, avviarlo, eseguire il comando e chiudere il container quando ha finito.

Buone notizie: è esattamente ciò che Nextflow farà per noi.
Dobbiamo solo specificare un container per ogni process.

Per dimostrare come funziona, abbiamo fatto un'altra versione del nostro workflow che esegue `cowpy` sul file di saluti raccolti prodotto nel terzo step.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Questo dovrebbe produrre un file contenente l'arte ASCII con i tre saluti nel fumetto.

#### 4.2.1. Esamina il codice

Il workflow è molto simile a quello precedente, più lo step extra per eseguire `cowpy`.

??? full-code "File di codice completo"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Include i moduli
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
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
        sayHello(greeting_ch)
        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // genera arte ASCII dei saluti con cowpy
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

Vedi che questo workflow importa un process `cowpy` da un file modulo, e lo chiama sull'output della chiamata `collectGreetings()`, più un parametro di input chiamato `params.character`.

```groovy title="2d-container.nf" linenums="31"
// genera arte ASCII dei saluti con cowpy
cowpy(collectGreetings.out.outfile, params.character)
```

Il process `cowpy`, che incapsula il comando cowpy per generare arte ASCII, è definito nel modulo `cowpy.nf`.

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

Il process `cowpy` richiede due input: il percorso di un file di input contenente il testo da mettere nel fumetto (`input_file`), e un valore per la variabile character.

Importante, include anche la riga `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, che punta all'URI del container che abbiamo usato prima.

#### 4.2.2. Verifica che Docker sia abilitato nella configurazione

Anticiperemo leggermente la Parte 3 di questo corso di formazione introducendo il file di configurazione `nextflow.config`, che è uno dei modi principali che Nextflow offre per configurare l'esecuzione del workflow.
Quando un file chiamato `nextflow.config` è presente nella directory corrente, Nextflow lo caricherà automaticamente e applicherà qualsiasi configurazione contenga.

A tal fine, abbiamo incluso un file `nextflow.config` con una singola riga di codice che abilita Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Questa configurazione dice a Nextflow di usare Docker per qualsiasi process che specifica un container compatibile.

!!! tip "Suggerimento"

    È tecnicamente possibile abilitare l'esecuzione Docker dalla riga di comando, su base per-esecuzione, usando il parametro `-with-docker <container>`.
    Tuttavia, questo ci permette solo di specificare un container per l'intero workflow, mentre l'approccio che ti abbiamo appena mostrato ci permette di specificare un container diverso per process.
    Quest'ultimo è molto meglio per la modularità, la manutenzione del codice e la riproducibilità.

#### 4.2.3. Esegui il workflow

Solo per riepilogare, questo è ciò che stiamo per eseguire:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Pensi che funzionerà?

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

I primi tre step sono stati salvati nella cache dato che li abbiamo già eseguiti prima, ma il process `cowpy` è nuovo quindi viene effettivamente eseguito.

Puoi trovare l'output dello step `cowpy` nella directory `results`.

??? abstract "Contenuto del file"

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

Vedi che il personaggio sta dicendo tutti i saluti, dato che è stato eseguito sul file di saluti maiuscoli raccolti.

Ancora più importante, siamo stati in grado di eseguire questo come parte della nostra pipeline senza dover fare un'installazione vera e propria di cowpy e tutte le sue dipendenze.
E ora possiamo condividere la pipeline con i collaboratori e farla eseguire sulla loro infrastruttura senza che loro debbano installare nulla, a parte Docker o una delle sue alternative (come Singularity/Apptainer) come menzionato sopra.

#### 4.2.4. Ispeziona come Nextflow ha lanciato l'attività containerizzata

Come coda finale a questa sezione, diamo un'occhiata alla sottodirectory di lavoro per una delle chiamate al process `cowpy` per ottenere un po' più di insight su come Nextflow lavora con i container sotto il cofano.

Controlla l'output dal tuo comando `nextflow run` per trovare il percorso alla sottodirectory di lavoro per il process `cowpy`.
Guardando ciò che abbiamo ottenuto per l'esecuzione mostrata sopra, la riga del log della console per il process `cowpy` inizia con `[7f/caf718]`.
Questo corrisponde al seguente percorso di directory troncato: `work/7f/caf718`.

In quella directory, troverai il file `.command.run` che contiene tutti i comandi che Nextflow ha eseguito per te nel corso dell'esecuzione della pipeline.

??? abstract "Contenuto del file"

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
    ...
    ```

Se cerchi `nxf_launch` in questo file, dovresti vedere qualcosa come questo:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Questo comando di lancio mostra che Nextflow sta usando un comando `docker run` molto simile per lanciare la chiamata al process come abbiamo fatto quando l'abbiamo eseguito manualmente.
Monta anche la corrispondente sottodirectory di lavoro nel container, imposta la directory di lavoro all'interno del container di conseguenza, e esegue il nostro script bash templato nel file `.command.sh`.

Questo conferma che tutto il duro lavoro che abbiamo dovuto fare manualmente nella sezione precedente viene ora fatto per noi da Nextflow!

### Takeaway

Capisci quale ruolo giocano i container nella gestione delle versioni degli strumenti software e nell'assicurare la riproducibilità.

Più in generale, hai una comprensione di base di quali sono i componenti principali delle pipeline Nextflow del mondo reale e come sono organizzati.
Conosci i fondamenti di come Nextflow può elaborare input multipli efficientemente, eseguire workflow composti da step multipli connessi insieme, sfruttare componenti di codice modulari, e utilizzare container per una maggiore riproducibilità e portabilità.

### Cosa c'è dopo?

Prenditi un'altra pausa! Era un grande mucchio di informazioni su come funzionano le pipeline Nextflow.

Nell'ultima sezione di questa formazione, approfondiremo il tema della configurazione.
Imparerai come configurare l'esecuzione della tua pipeline per adattarla alla tua infrastruttura e come gestire la configurazione di input e parametri.

---

## Quiz

<quiz>
Perché Nextflow crea una directory di attività separata per ogni chiamata al process?
- [ ] Per migliorare la velocità di esecuzione
- [ ] Per ridurre l'utilizzo della memoria
- [x] Per isolare le esecuzioni ed evitare collisioni tra gli output
- [ ] Per abilitare la compressione dei file in parallelo

Approfondisci: [1.3. Trova gli output originali e i log](#13-trova-gli-output-originali-e-i-log)
</quiz>

<quiz>
Cosa fa l'opzione `-ansi-log false` quando si esegue un workflow?
- [ ] Disabilita tutto l'output della console
- [x] Rimuove i colori dall'output
- [x] Mostra tutti i percorsi delle directory di attività invece di condensarli su una riga
- [ ] Abilita la modalità di debug verbose

Approfondisci: [1.3.2. Fai mostrare al terminale più dettagli](#132-fai-mostrare-al-terminale-piu-dettagli)

Puoi anche usare una delle seguenti variabili d'ambiente se preferisci questo stile:

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

Approfondisci: [1.4.1. Caricare i dati di input dal CSV](#141-caricare-i-dati-di-input-dal-csv)
</quiz>

<quiz>
Perché è importante includere il valore di input nei nomi dei file di output (es. `#!groovy "${greeting}-output.txt"`)?
- [ ] Per migliorare la velocità di elaborazione
- [ ] Per abilitare la funzionalità resume
- [x] Per evitare che i file di output si sovrascrivano a vicenda quando si elaborano input multipli
- [ ] Per rendere i file più facili da comprimere

Approfondisci: [1.4.3. Come vengono nominati gli output](#143-come-vengono-nominati-gli-output)
</quiz>

<quiz>
Qual è lo scopo della dichiarazione `include` in un workflow modularizzato?
- [ ] Copiare il codice del process nel file del workflow
- [x] Importare una definizione di process da un file modulo esterno
- [ ] Includere impostazioni di configurazione
- [ ] Aggiungere commenti di documentazione

Approfondisci: [3. Eseguire pipeline modularizzate](#3-eseguire-pipeline-modularizzate)
</quiz>

<quiz>
Quando modularizzi un workflow e lo esegui con `-resume`, cosa succede?
- [ ] Il caching è disabilitato per i process modulari
- [ ] Tutte le attività devono essere rieseguite
- [x] Il caching funziona normalmente in base agli script di job generati
- [ ] Solo il file del workflow principale viene salvato nella cache

Approfondisci: [3.2. Esegui il workflow](#32-esegui-il-workflow)
</quiz>

<quiz>
Cosa specifica la direttiva `container` in una definizione di process?
- [ ] La directory di lavoro per il process
- [ ] L'allocazione massima di memoria
- [x] L'URI dell'immagine container da usare per eseguire il process
- [ ] Il formato del file di output

Approfondisci: [4.2. Usa un container in un workflow](#42-usa-un-container-in-un-workflow)
</quiz>

<quiz>
Nel file `.command.run`, cosa contiene la funzione `nxf_launch`?
- [ ] Le informazioni sulla versione di Nextflow
- [ ] I parametri del workflow
- [x] Il comando `docker run` con i mount dei volumi e le impostazioni del container
- [ ] Le dichiarazioni di input del process

Approfondisci: [4.2.4. Ispeziona come Nextflow ha lanciato l'attività containerizzata](#424-ispeziona-come-nextflow-ha-lanciato-lattivita-containerizzata)
</quiz>

<quiz>
Cosa gestisce automaticamente Nextflow quando esegue un process containerizzato? (Seleziona tutte le risposte corrette)
- [x] Scaricare l'immagine container se necessario
- [x] Montare la directory di lavoro nel container
- [x] Eseguire lo script del process all'interno del container
- [x] Pulire l'istanza del container dopo l'esecuzione

Approfondisci: [4. Usare software containerizzato](#4-usare-software-containerizzato)
</quiz>
