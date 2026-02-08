# Parte 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=it" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guardate [l'intera playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/02_hello_channels.md).
///

Nella Parte 1 di questo corso (Hello World), vi abbiamo mostrato come fornire un input variabile a un processo passandolo direttamente nella chiamata al processo: `sayHello(params.input)`.
Quello era un approccio volutamente semplificato.
In pratica, quell'approccio ha limitazioni importanti; in particolare, funziona solo per casi molto semplici in cui vogliamo eseguire il processo una sola volta, su un singolo valore.
Nella maggior parte dei casi d'uso realistici dei flussi di lavoro, vogliamo elaborare più valori (dati sperimentali per più campioni, per esempio), quindi abbiamo bisogno di un modo più sofisticato per gestire gli input.

Ecco a cosa servono i [**canali**](https://nextflow.io/docs/latest/channel.html) di Nextflow.
I canali sono code progettate per gestire gli input in modo efficiente e trasferirli da un passaggio all'altro nei flussi di lavoro multi-step, fornendo al contempo parallelismo integrato e molti altri vantaggi.

In questa parte del corso, imparerete come utilizzare un canale per gestire più input da diverse fonti.
Imparerete anche a utilizzare gli [**operatori**](https://nextflow.io/docs/latest/reference/operator.html) per trasformare i contenuti dei canali secondo necessità.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato la Parte 1 del corso [Hello Nextflow](./index.md), ma se avete familiarità con i concetti base trattati in quella sezione, potete iniziare da qui senza fare nulla di speciale.

---

## 0. Riscaldamento: Eseguire `hello-channels.nf`

Useremo lo script del flusso di lavoro `hello-channels.nf` come punto di partenza.
È equivalente allo script prodotto seguendo la Parte 1 di questo corso di formazione, tranne che abbiamo cambiato la destinazione dell'output:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Solo per assicurarci che tutto funzioni, eseguite lo script una volta prima di apportare modifiche:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Come in precedenza, troverete il file di output chiamato `output.txt` nella directory `results/hello_channels` (come specificato nel blocco `output` dello script del flusso di lavoro, mostrato sopra).

??? abstract "Directory contents"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contenuti del file"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Se tutto ha funzionato, siete pronti a imparare i canali.

---

## 1. Fornire input variabili tramite un canale esplicitamente

Creeremo un **canale** per passare l'input variabile al processo `sayHello()` invece di affidarci alla gestione implicita, che ha alcune limitazioni.

### 1.1. Creare un canale di input

Esistono diverse [**fabbriche di canali**](https://nextflow.io/docs/latest/reference/channel.html) che possiamo usare per configurare un canale.
Per mantenere le cose semplici per ora, useremo la fabbrica di canali più basilare, chiamata [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of), che creerà un canale contenente un singolo valore.
Funzionalmente sarà simile a come lo avevamo configurato prima, ma invece di far creare a Nextflow un canale implicitamente, ora lo stiamo facendo esplicitamente.

Questa è la riga di codice che useremo:

```console title="Sintassi"
greeting_ch = channel.of('Hello Channels!')
```

Questo crea un canale chiamato `greeting_ch` usando la fabbrica di canali `channel.of()`, che configura un semplice canale di coda, e carica la stringa `'Hello Channels!'` da usare come valore di saluto.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Nota"

    Stiamo temporaneamente tornando alle stringhe hardcoded invece di usare un parametro CLI per motivi di leggibilità. Torneremo a usare i parametri CLI una volta che avremo spiegato cosa succede a livello di canale.

Nel blocco workflow, aggiungete il codice della fabbrica di canali:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // crea un canale per gli input
        greeting_ch = channel.of('Hello Channels!')
        // emette un saluto
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // emette un saluto
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Questo non è ancora funzionale poiché non abbiamo ancora cambiato l'input nella chiamata al processo.

### 1.2. Aggiungere il canale come input alla chiamata del processo

Ora dobbiamo effettivamente collegare il nostro canale appena creato alla chiamata del processo `sayHello()`, sostituendo il parametro CLI che fornivamo direttamente prima.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canale per gli input
        greeting_ch = channel.of('Hello Channels!')
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canale per gli input
        greeting_ch = channel.of('Hello Channels!')
        // emette un saluto
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Questo dice a Nextflow di eseguire il processo `sayHello` sui contenuti del canale `greeting_ch`.

Ora il nostro flusso di lavoro è correttamente funzionale; è l'equivalente esplicito di scrivere `sayHello('Hello Channels!')`.

### 1.3. Eseguire il flusso di lavoro

Eseguiamolo!

```bash
nextflow run hello-channels.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Se avete effettuato entrambe le modifiche correttamente, dovreste ottenere un'esecuzione riuscita.
Potete controllare la directory dei risultati per verificare che l'output sia ancora lo stesso di prima.

??? abstract "Contenuti del file"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Quindi abbiamo aumentato la flessibilità del nostro flusso di lavoro ottenendo lo stesso risultato finale.
Questo può sembrare scrivere più codice senza un beneficio tangibile, ma il valore diventerà chiaro non appena inizieremo a gestire più input.

Come anteprima, vediamo un'altra cosa prima di procedere: un piccolo ma conveniente vantaggio dell'uso di un canale esplicito per gestire l'input dei dati.

### 1.4. Usare `view()` per ispezionare i contenuti del canale

I canali di Nextflow sono costruiti in modo da permetterci di operare sui loro contenuti usando operatori, che tratteremo in dettaglio più avanti in questo capitolo.

Per ora, vi mostreremo semplicemente come usare un operatore super semplice chiamato [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) per ispezionare i contenuti di un canale.
Potete pensare a `view()` come uno strumento di debug, come un'istruzione `print()` in Python, o il suo equivalente in altri linguaggi.

Aggiungete questa piccola riga al blocco workflow:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canale per gli input
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canale per gli input
        greeting_ch = channel.of('Hello Channels!')
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

La quantità esatta di spazi non importa purché sia un multiplo di 4; stiamo solo cercando di allineare l'inizio dell'istruzione `.view()` alla parte `.of()` della costruzione del canale.

Ora eseguite di nuovo il flusso di lavoro:

```bash
nextflow run hello-channels.nf
```

??? success "Output del comando"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Come potete vedere, questo mostra i contenuti del canale sulla console.
Qui abbiamo solo un elemento, ma quando inizieremo a caricare più valori nel canale nella prossima sezione, vedrete che questo è impostato per mostrare un elemento per riga.

### Takeaway

Sapete come usare una fabbrica di canali di base per fornire un input a un processo.

### Cosa c'è dopo?

Imparate come usare i canali per far iterare il flusso di lavoro su più valori di input.

---

## 2. Modificare il flusso di lavoro per eseguire su più valori di input

I flussi di lavoro tipicamente vengono eseguiti su lotti di input che devono essere elaborati in blocco, quindi vogliamo aggiornare il flusso di lavoro per accettare più valori di input.

### 2.1. Caricare più saluti nel canale di input

Convenientemente, la fabbrica di canali `channel.of()` che abbiamo usato è perfettamente in grado di accettare più di un valore, quindi non abbiamo bisogno di modificarla affatto.
Possiamo semplicemente caricare più valori nel canale.

Usiamo `'Hello'`, `'Bonjour'` e `'Holà'`.

#### 2.1.1. Aggiungere più saluti

Prima del blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // crea un canale per gli input
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // crea un canale per gli input
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

La documentazione ci dice che questo dovrebbe funzionare. Può davvero essere così semplice?

#### 2.1.2. Eseguire il comando e guardare l'output del log

Proviamo.

```bash
nextflow run hello-channels.nf
```

??? success "Output del comando"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Sembra certamente che abbia funzionato correttamente.
Il monitor di esecuzione mostra che sono state effettuate `3 of 3` chiamate per il processo `sayHello`, e vediamo i tre saluti elencati dall'istruzione `view()`, uno per riga come promesso.

Tuttavia, c'è ancora un solo output nella directory dei risultati:

??? abstract "Directory contents"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contenuti del file"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Dovreste vedere uno dei tre saluti lì dentro, ma quello che avete ottenuto potrebbe essere diverso da quello mostrato qui.
Riuscite a pensare al motivo?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Nel diagramma, il canale è rappresentato in verde, e l'ordine degli elementi è rappresentato come biglie in un tubo: il primo caricato è a destra, poi il secondo al centro, poi il terzo è a sinistra._

Guardando indietro al monitor di esecuzione, ci ha dato solo un percorso di sottodirectory (`f4/c9962c`).
Diamo un'occhiata lì dentro.

??? abstract "Directory contents"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contenuti del file"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Questo non è nemmeno lo stesso saluto che abbiamo ottenuto nella directory dei risultati! Cosa sta succedendo?

A questo punto, dobbiamo dirvi che per impostazione predefinita, il sistema di logging ANSI scrive il logging di più chiamate allo stesso processo sulla stessa riga.
Quindi lo stato di tutte e tre le chiamate al processo sayHello() finiscono nello stesso punto.

Fortunatamente, possiamo disabilitare quel comportamento per vedere l'elenco completo delle chiamate ai processi.

#### 2.1.3. Eseguire di nuovo il comando con l'opzione `-ansi-log false`

Per espandere il logging e visualizzare una riga per chiamata di processo, aggiungete `-ansi-log false` al comando.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Output del comando"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Questa volta vediamo tutte e tre le esecuzioni dei processi e le loro sottodirectory di lavoro associate elencate nell'output.

Molto meglio, almeno per un flusso di lavoro semplice.
Per un flusso di lavoro complesso, o un grande numero di input, avere l'elenco completo in output sul terminale diventerebbe un po' opprimente.
Ecco perché `-ansi-log false` non è il comportamento predefinito.

!!! tip "Suggerimento"

    Il modo in cui viene riportato lo stato è un po' diverso tra le due modalità di logging.
    Nella modalità condensata, Nextflow riporta se le chiamate sono state completate con successo o meno.
    In questa modalità espansa, riporta solo che sono state inviate.

Comunque, ora che abbiamo le sottodirectory di ogni chiamata di processo, possiamo cercare i loro log e output.

??? abstract "Directory contents"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contenuti del file"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Questo mostra che tutti e tre i processi sono stati eseguiti con successo (evviva).

Detto questo, abbiamo ancora il problema che c'è un solo file di output nella directory dei risultati.

Potete ricordare che abbiamo hardcodato il nome del file di output per il processo `sayHello`, quindi tutte e tre le chiamate hanno prodotto un file chiamato `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Finché i file di output rimangono nelle sottodirectory di lavoro, isolati dagli altri processi, va bene.
Ma quando vengono pubblicati nella stessa directory dei risultati, quello che viene copiato per primo viene sovrascritto dal successivo, e così via.

### 2.2. Assicurarsi che i nomi dei file di output siano unici

Possiamo continuare a pubblicare tutti gli output nella stessa directory dei risultati, ma dobbiamo assicurarci che abbiano nomi unici.
Nello specifico, dobbiamo modificare il primo processo per generare un nome di file dinamicamente in modo che i nomi dei file finali siano unici.

Quindi come rendiamo i nomi dei file unici?
Un modo comune per farlo è usare qualche pezzo unico di metadati dagli input (ricevuti dal canale di input) come parte del nome del file di output.
Qui, per comodità, useremo semplicemente il saluto stesso poiché è solo una stringa corta, e lo anteporremo al nome base del file di output.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Costruire un nome di file di output dinamico

Nel blocco process, effettuate le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

Assicuratevi di sostituire `output.txt` sia nella definizione dell'output che nel blocco del comando `script:`.

!!! tip "Suggerimento"

    Nella definizione dell'output, DEVE usare le virgolette doppie attorno all'espressione del nome del file (NON le virgolette singole), altrimenti fallirà.

Questo dovrebbe produrre un nome di file di output unico ogni volta che il processo viene chiamato, in modo che possa essere distinto dagli output di altre chiamate allo stesso processo nella directory di output.

#### 2.2.2. Eseguire il flusso di lavoro

Eseguiamolo. Nota che siamo tornati a eseguire con le impostazioni predefinite del log ANSI.

```bash
nextflow run hello-channels.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Tornando alla vista riassuntiva, l'output è di nuovo riassunto su una riga.
Date un'occhiata alla directory `results` per vedere se tutti i saluti di output sono presenti.

??? abstract "Directory contents"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Sì! E ognuno ha i contenuti attesi.

??? abstract "Contenuti del file"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Successo! Ora possiamo aggiungere quanti saluti vogliamo senza preoccuparci che i file di output vengano sovrascritti.

!!! tip "Suggerimento"

    In pratica, nominare i file basandosi sui dati di input stessi è quasi sempre poco pratico.
    Il modo migliore per generare nomi di file dinamici è passare metadati a un processo insieme ai file di input.
    I metadati sono tipicamente forniti tramite un 'sample sheet' o equivalenti.
    Imparerete come farlo più avanti nella vostra formazione su Nextflow (vedete la [side quest sui metadati](../side_quests/metadata.md)).

### Takeaway

Sapete come fornire più elementi di input attraverso un canale.

### Cosa c'è dopo?

Imparate a usare un operatore per trasformare i contenuti di un canale.

---

## 3. Fornire più input tramite un array

Vi abbiamo appena mostrato come gestire più elementi di input che erano hardcodati direttamente nella fabbrica di canali.
E se volessimo fornire questi input multipli in modo diverso?

Per esempio, immaginate di configurare una variabile di input contenente un array di elementi come questo:

`greetings_array = ['Hello','Bonjour','Holà']`

Possiamo caricarla nel nostro canale di output e aspettarci che funzioni?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Scopriamolo.

### 3.1. Fornire un array di valori come input al canale

Il buon senso suggerisce che dovremmo essere in grado di passare semplicemente un array di valori invece di un singolo valore.
Proviamo; dovremo configurare la variabile di input e caricarla nella fabbrica di canali.

#### 3.1.1. Configurare la variabile di input

Prendiamo la variabile `greetings_array` che abbiamo appena immaginato e rendiamola realtà aggiungendola al blocco workflow:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canale per gli input
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Questo non è ancora funzionale, abbiamo solo aggiunto una dichiarazione per l'array.

#### 3.1.2. Impostare l'array di saluti come input alla fabbrica di canali

Ora sostituiremo i valori `'Hello','Bonjour','Holà'` attualmente hardcodati nella fabbrica di canali con il `greetings_array` che abbiamo appena creato.

Nel blocco workflow, effettuate la seguente modifica:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Questo dovrebbe essere funzionale ora.

#### 3.1.3. Eseguire il flusso di lavoro

Proviamo a eseguirlo:

```bash
nextflow run hello-channels.nf
```

??? failure "Output del comando"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh no! C'è un errore!

Guardate l'output di `view()` e i messaggi di errore.

Sembra che Nextflow abbia provato a eseguire una singola chiamata di processo, usando `[Hello, Bonjour, Holà]` come valore stringa, invece di usare le tre stringhe nell'array come valori separati.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Quindi è il 'confezionamento' che sta causando il problema.
Come facciamo a far sì che Nextflow spacchetti l'array e carichi le singole stringhe nel canale?

### 3.2. Usare un operatore per trasformare i contenuti del canale

Qui entrano in gioco gli [**operatori**](https://www.nextflow.io/docs/latest/reference/operator.html).
Avete già usato l'operatore `.view()`, che semplicemente guarda cosa c'è dentro.
Ora vedremo operatori che ci permettono di agire sui contenuti di un canale.

Se scorrete la [lista degli operatori](https://www.nextflow.io/docs/latest/reference/operator.html) nella documentazione di Nextflow, troverete [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten), che fa esattamente ciò di cui abbiamo bisogno: spacchetta i contenuti di un array e li emette come elementi individuali.

#### 3.2.1. Aggiungere l'operatore `flatten()`

Per applicare l'operatore `flatten()` al nostro canale di input, lo aggiungiamo alla dichiarazione della fabbrica di canali.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Qui abbiamo aggiunto l'operatore sulla riga successiva per leggibilità, ma potete aggiungere operatori sulla stessa riga della fabbrica di canali se preferite, così:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Affinare le istruzioni `view()`

Potremmo eseguirlo subito per testare se funziona, ma già che ci siamo, affineremo come ispezioniamo i contenuti del canale.

Vogliamo poter confrontare come appaiono i contenuti prima e dopo l'applicazione dell'operatore `flatten()`, quindi ne aggiungeremo un secondo, E aggiungeremo un po' di codice per etichettarli più chiaramente nell'output.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Prima di flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Dopo flatten: $greeting" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Vedete che abbiamo aggiunto una seconda istruzione `.view`, e per ciascuna di esse, abbiamo sostituito le parentesi vuote (`()`) con parentesi graffe contenenti del codice, come `{ greeting -> "Prima di flatten: $greeting" }`.

Queste si chiamano _closure_. Il codice che contengono verrà eseguito per ogni elemento nel canale.
Definiamo una variabile temporanea per il valore interno, qui chiamata `greeting` (ma potrebbe essere qualsiasi nome arbitrario), che viene usata solo nell'ambito di quella closure.

In questo esempio, `$greeting` rappresenta ogni singolo elemento caricato nel canale.
Questo risulterà in un output della console ordinatamente etichettato.

!!! info "Informazione"

    In alcune pipeline potrete vedere una variabile speciale chiamata `$it` usata all'interno delle closure degli operatori.
    Questa è una variabile _implicita_ che permette un accesso abbreviato alla variabile interna,
    senza doverla definire con un `->`.

    Preferiamo essere espliciti per aiutare la chiarezza del codice, quindi la sintassi `$it` è sconsigliata e sarà gradualmente eliminata dal linguaggio Nextflow.

#### 3.2.3. Eseguire il flusso di lavoro

Finalmente, potete provare a eseguire di nuovo il flusso di lavoro!

```bash
nextflow run hello-channels.nf
```

??? success "Output del comando"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Prima di flatten: [Hello, Bonjour, Holà]
    Dopo flatten: Hello
    Dopo flatten: Bonjour
    Dopo flatten: Holà
    ```

Questa volta funziona E ci dà la comprensione aggiuntiva di come appaiono i contenuti del canale prima e dopo l'esecuzione dell'operatore `flatten()`.

- Vedete che otteniamo una singola istruzione `Prima di flatten:` perché a quel punto il canale contiene un elemento, l'array originale.
  Poi otteniamo tre istruzioni separate `Dopo flatten:`, una per ogni saluto, che ora sono elementi individuali nel canale.

Importante, questo significa che ogni elemento può ora essere elaborato separatamente dal flusso di lavoro.

!!! tip "Suggerimento"

    È tecnicamente possibile ottenere gli stessi risultati usando una fabbrica di canali diversa, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), che include un passo di mapping implicito nella sua operazione.
    Qui abbiamo scelto di non usarla per dimostrare l'uso di un operatore su un caso d'uso semplice.

### Takeaway

Sapete come usare un operatore come `flatten()` per trasformare i contenuti di un canale, e come usare l'operatore `view()` per ispezionare i contenuti del canale prima e dopo l'applicazione di un operatore.

### Cosa c'è dopo?

Imparate come far prendere al flusso di lavoro un file come fonte di valori di input.

---

## 4. Leggere i valori di input da un file CSV

Realisticamente, raramente o mai partiremo da un array di valori.
Molto probabilmente, avremo uno o più file contenenti i dati che devono essere elaborati, in qualche tipo di formato strutturato.

Abbiamo preparato un file CSV chiamato `greetings.csv` che contiene diversi saluti di input, imitando il tipo di dati colonnari che potreste voler elaborare in un'analisi dati reale, memorizzato sotto `data/`.
(I numeri non sono significativi, sono lì solo a scopo illustrativo.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Il nostro prossimo compito è adattare il nostro flusso di lavoro per leggere i valori da questo file.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Vediamo come possiamo realizzarlo.

### 4.1. Modificare lo script per aspettarsi un file CSV come fonte di saluti

Per iniziare, dovremo apportare due modifiche chiave allo script:

- Cambiare il parametro di input per puntare al file CSV
- Cambiare la fabbrica di canali con una progettata per gestire un file

#### 4.1.1. Cambiare il parametro di input per puntare al file CSV

Ricordate il parametro `params.input` che abbiamo configurato nella Parte 1?
Lo aggiorneremo per puntare al file CSV contenente i nostri saluti.

Effettuate la seguente modifica alla dichiarazione del parametro:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    input: String = 'Holà mundo!'
    ```

Questo presuppone che il file sia nella stessa posizione del codice del flusso di lavoro.
Imparerete come gestire altre posizioni dei dati più avanti nel vostro percorso con Nextflow.

#### 4.1.2. Passare a una fabbrica di canali progettata per gestire un file

Poiché ora vogliamo usare un file invece di semplici stringhe come input, non possiamo usare la fabbrica di canali `channel.of()` di prima.
Dobbiamo passare a usare una nuova fabbrica di canali, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath), che ha alcune funzionalità integrate per gestire i percorsi dei file.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Prima di flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Dopo flatten: $greeting" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // dichiara un array di saluti di input
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canale per gli input
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Prima di flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Dopo flatten: $greeting" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Noterete che abbiamo riportato l'input del canale a `param.input`, e cancellato la dichiarazione `greetings_array` poiché non ne avremo più bisogno.
Abbiamo anche commentato `flatten()` e la seconda istruzione `view()`.

#### 4.1.3. Eseguire il flusso di lavoro

Proviamo a eseguire il flusso di lavoro con la nuova fabbrica di canali e il file di input.

```bash
nextflow run hello-channels.nf
```

??? failure "Output del comando"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Prima di flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh no, non funziona. Date un'occhiata all'inizio dell'output della console e al messaggio di errore.
La parte `Command executed:` è particolarmente utile qui.

Questo potrebbe sembrare un po' familiare.
Sembra che Nextflow abbia provato a eseguire una singola chiamata di processo usando il percorso del file stesso come valore stringa.
Quindi ha risolto correttamente il percorso del file, ma non ha effettivamente analizzato i suoi contenuti, che è quello che volevamo.

Come facciamo a far sì che Nextflow apra il file e carichi i suoi contenuti nel canale?

Sembra che abbiamo bisogno di un altro [operatore](https://nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Usare l'operatore `splitCsv()` per analizzare il file

Guardando di nuovo la lista degli operatori, troviamo [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), che è progettato per analizzare e dividere testo formattato CSV.

#### 4.2.1. Applicare `splitCsv()` al canale

Per applicare l'operatore, lo aggiungiamo alla riga della fabbrica di canali come fatto in precedenza.

Nel blocco workflow, effettuate la seguente modifica al codice per sostituire `flatten()` con `splitcsv()` (senza commento):

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Prima di splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Dopo splitCsv: $csv" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Prima di flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Dopo flatten: $greeting" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Come potete vedere, abbiamo anche aggiornato le istruzioni `view()` prima/dopo.
Tecnicamente avremmo potuto usare lo stesso nome di variabile (`greeting`) ma l'abbiamo aggiornato a qualcosa di più appropriato (`csv`) per rendere il codice più leggibile agli altri.

#### 4.2.2. Eseguire di nuovo il flusso di lavoro

Proviamo a eseguire il flusso di lavoro con la logica di parsing CSV aggiunta.

```bash
nextflow run hello-channels.nf
```

??? failure "Output del comando"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Prima di splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Dopo splitCsv: [Hello, English, 123]
    Dopo splitCsv: [Bonjour, French, 456]
    Dopo splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Interessante, fallisce anche questo, ma con un errore diverso.
Questa volta Nextflow ha analizzato i contenuti del file (evviva!) ma ha caricato ogni riga come un array, e ogni array è un elemento nel canale.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Dobbiamo dirgli di prendere solo la prima colonna di ogni riga.
Quindi come spacchettamo questo?

Abbiamo usato precedentemente `flatten()` per spacchettare i contenuti di un canale, ma non funzionerebbe qui perché flatten spacchetta _tutto_ (potete provarlo se volete vedere di persona).

Invece, useremo un altro operatore chiamato `map()` che è davvero utile e compare spesso nelle pipeline Nextflow.

### 4.3. Usare l'operatore `map()` per estrarre i saluti

L'operatore [`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) è un piccolo strumento molto utile che ci permette di fare ogni tipo di mappatura sui contenuti di un canale.

In questo caso, lo useremo per estrarre quell'unico elemento che vogliamo da ogni riga nel nostro file di dati.
Ecco come appare la sintassi:

```groovy title="Sintassi"
.map { row -> row[0] }
```

Questo significa 'per ogni riga nel canale, prendi l'elemento 0 (primo) che contiene'.

Quindi applichiamolo al nostro parsing CSV.

#### 4.3.1. Applicare `map()` al canale

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Prima di splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Dopo splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "Dopo map: $csv" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Prima di splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Dopo splitCsv: $csv" }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Vedete che abbiamo aggiunto un'altra chiamata `view()` per confermare che l'operatore fa quello che ci aspettiamo.

#### 4.3.2. Eseguire il flusso di lavoro

Eseguiamolo ancora una volta:

```bash
nextflow run hello-channels.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Prima di splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Dopo splitCsv: [Hello, English, 123]
    Dopo splitCsv: [Bonjour, French, 456]
    Dopo splitCsv: [Holà, Spanish, 789]
    Dopo map: Hello
    Dopo map: Bonjour
    Dopo map: Holà
    ```

Questa volta dovrebbe funzionare senza errori.

Guardando l'output delle istruzioni `view()`, vedete il seguente:

- Una singola istruzione `Prima di splitCsv:`: a quel punto il canale contiene un elemento, il percorso del file originale.
- Tre istruzioni separate `Dopo splitCsv:`: una per ogni saluto, ma ciascuna è contenuta in un array che corrisponde a quella riga nel file.
- Tre istruzioni separate `Dopo map:`: una per ogni saluto, che ora sono elementi individuali nel canale.

_Nota che le righe potrebbero apparire in un ordine diverso nel vostro output._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

Potete anche guardare i file di output per verificare che ogni saluto sia stato correttamente estratto ed elaborato attraverso il flusso di lavoro.

Abbiamo ottenuto lo stesso risultato di prima, ma ora abbiamo molta più flessibilità per aggiungere più elementi al canale di saluti che vogliamo elaborare modificando un file di input, senza modificare alcun codice.
Imparerete approcci più sofisticati per gestire input complessi in una formazione successiva.

### Takeaway

Sapete come usare il costruttore di canali `.fromPath()` e gli operatori `splitCsv()` e `map()` per leggere un file di valori di input e gestirli appropriatamente.

Più in generale, avete una comprensione di base di come Nextflow usa i **canali** per gestire gli input ai processi e gli **operatori** per trasformare i loro contenuti.
Avete anche visto come i canali gestiscono l'esecuzione parallela implicitamente.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### Cosa c'è dopo?

Prendetevi una bella pausa, avete lavorato duramente in questa sezione!

Quando siete pronti, passate alla [**Parte 3: Hello Workflow**](./03_hello_workflow.md) per imparare come aggiungere più passaggi e connetterli insieme in un flusso di lavoro vero e proprio.

---

## Quiz

<quiz>
Cos'è un canale in Nextflow?
- [ ] Una specifica di percorso file
- [ ] Una definizione di processo
- [x] Una struttura simile a una coda per passare dati tra processi
- [ ] Un'impostazione di configurazione

Per approfondire: [1.1. Creare un canale di input](#11-creare-un-canale-di-input)
</quiz>

<quiz>
Cosa produrrà questo codice?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (una singola lista)
- [x] Ogni elemento su una riga separata: `Hello`, `Bonjour`, `Hola`
- [ ] Niente (i canali non stampano per impostazione predefinita)
- [ ] Un errore (sintassi non valida)

Per approfondire: [1.1. Creare un canale di input](#11-creare-un-canale-di-input)
</quiz>

<quiz>
Quando un canale contiene più valori, come gestisce Nextflow l'esecuzione del processo?
- [ ] Il processo viene eseguito una volta con tutti i valori
- [x] Il processo viene eseguito una volta per ogni valore nel canale
- [ ] Il processo viene eseguito solo con il primo valore
- [ ] Il processo viene eseguito solo con l'ultimo valore

Per approfondire: [2. Modificare il flusso di lavoro per eseguire su più valori di input](#2-modificare-il-flusso-di-lavoro-per-eseguire-su-piu-valori-di-input)
</quiz>

<quiz>
Cosa fa l'operatore `flatten()`?
- [ ] Combina più canali in uno
- [ ] Ordina gli elementi del canale
- [x] Spacchetta gli array in elementi individuali
- [ ] Rimuove gli elementi duplicati

Per approfondire: [3.2.1. Aggiungere l'operatore `flatten()`](#321-aggiungere-loperatore-flatten)
</quiz>

<quiz>
Qual è lo scopo dell'operatore `view()`?
- [ ] Filtrare i contenuti del canale
- [ ] Trasformare gli elementi del canale
- [x] Ispezionare e fare debug dei contenuti del canale
- [ ] Salvare i contenuti del canale in un file

Per approfondire: [1.4. Usare `view()` per ispezionare i contenuti del canale](#14-usare-view-per-ispezionare-i-contenuti-del-canale)
</quiz>

<quiz>
Cosa fa `splitCsv()`?
- [ ] Crea un file CSV dai contenuti del canale
- [ ] Divide una stringa per virgole
- [x] Analizza un file CSV in array che rappresentano ogni riga
- [ ] Unisce più file CSV

Per approfondire: [4.2. Usare l'operatore `splitCsv()` per analizzare il file](#42-usare-loperatore-splitcsv-per-analizzare-il-file)
</quiz>

<quiz>
Qual è lo scopo dell'operatore `map()`?
- [ ] Filtrare elementi da un canale
- [ ] Combinare più canali
- [x] Trasformare ogni elemento in un canale
- [ ] Contare gli elementi in un canale

Per approfondire: [4.3. Usare l'operatore `map()` per estrarre i saluti](#43-usare-loperatore-map-per-estrarre-i-saluti)
</quiz>

<quiz>
Perché è importante usare nomi di file di output dinamici quando si elaborano più input?
- [ ] Per migliorare le prestazioni
- [ ] Per ridurre lo spazio su disco
- [x] Per evitare che i file di output si sovrascrivano a vicenda
- [ ] Per abilitare la funzionalità di resume

Per approfondire: [2.2. Assicurarsi che i nomi dei file di output siano unici](#22-assicurarsi-che-i-nomi-dei-file-di-output-siano-unici)
</quiz>
