# Parte 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/01_hello_world.md).
///
-->

In questa prima parte del corso di formazione Hello Nextflow, introduciamo l'argomento con un esempio Hello World molto semplice e indipendente dal dominio, che costruiremo progressivamente per dimostrare l'uso della logica e dei componenti fondamentali di Nextflow.

??? info "Cos'è un esempio Hello World?"

    Un "Hello World!" è un esempio minimalista pensato per dimostrare la sintassi e la struttura di base di un linguaggio di programmazione o framework software.
    L'esempio consiste tipicamente nella stampa della frase "Hello, World!" su un dispositivo di output, come la console o il terminale, oppure nella sua scrittura su un file.

---

## 0. Riscaldamento: Eseguire un esempio Hello World direttamente

Dimostriamolo con un semplice comando che eseguiamo direttamente nel terminale, per mostrare cosa fa prima di incapsularlo in Nextflow.

!!! tip "Suggerimento"

    Ricordate che ora dovreste trovarvi nella directory `hello-nextflow/` come descritto nella pagina [Per iniziare](00_orientation.md).

### 0.1. Far dire hello al terminale

Eseguite il seguente comando nel terminale.

```bash
echo 'Hello World!'
```

??? success "Output del comando"

    ```console
    Hello World!
    ```

Questo produce il testo 'Hello World' direttamente nel terminale.

### 0.2. Scrivere l'output su un file

L'esecuzione di pipeline coinvolge principalmente la lettura di dati da file e la scrittura di risultati su altri file, quindi modifichiamo il comando per scrivere l'output testuale su un file per rendere l'esempio un po' più pertinente.

```bash
echo 'Hello World!' > output.txt
```

??? success "Output del comando"

    ```console

    ```

Questo non produce alcun output nel terminale.

### 0.3. Trovare l'output

Il testo 'Hello World' dovrebbe ora trovarsi nel file di output che abbiamo specificato, chiamato `output.txt`.
Potete aprirlo nell'esplora file o dalla riga di comando utilizzando l'utility `cat`, ad esempio.

??? abstract "Contenuti del file"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Questo è ciò che cercheremo di replicare con il nostro primo workflow Nextflow.

### Takeaway

Ora sapete come eseguire un semplice comando nel terminale che produce del testo, e opzionalmente, come farlo scrivere l'output su un file.

### Cosa c'è dopo?

Scoprite come apparirebbe scritto come workflow Nextflow.

---

## 1. Esaminare lo script ed eseguirlo

Vi forniamo uno script di workflow completamente funzionale, anche se minimalista, chiamato `hello-world.nf` che fa la stessa cosa di prima (scrive 'Hello World!') ma con Nextflow.

Per iniziare, apriamo lo script del workflow così potete farvi un'idea di come è strutturato.
Poi lo eseguiremo e cercheremo i suoi output.

### 1.1. Esaminare il codice

Troverete lo script `hello-world.nf` nella vostra directory corrente, che dovrebbe essere `hello-nextflow`. Apritelo nel riquadro dell'editor.

??? full-code "File di codice completo"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usa echo per stampare 'Hello World!' su un file
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // emette un saluto
        sayHello()
    }
    ```

Uno script di workflow Nextflow include tipicamente una o più definizioni di **process** e il **workflow** stesso, più alcuni blocchi opzionali (non presenti qui) che introdurremo più avanti.

Ogni **process** descrive quale/i operazione/i il corrispondente step nella pipeline dovrebbe eseguire, mentre il **workflow** descrive la logica del flusso di dati che connette i vari step.

Esamineremo prima più da vicino il blocco **process**, poi guarderemo il blocco **workflow**.

#### 1.1.1. La definizione del `process`

Il primo blocco di codice descrive un **process**.

La definizione del process inizia con la parola chiave `process`, seguita dal nome del process e infine dal corpo del process delimitato da parentesi graffe.
Il corpo del process deve contenere un blocco script che specifica il comando da eseguire, che può essere qualsiasi cosa si possa eseguire in un terminale a riga di comando.

```groovy title="hello-world.nf" linenums="3"
/*
* Usa echo per stampare 'Hello World!' su un file
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Qui abbiamo un **process** chiamato `sayHello` che scrive il suo **output** su un file chiamato `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Questa è una definizione di process molto minimale che contiene solo una definizione di `output` e lo `script` da eseguire.

La definizione di `output` include il qualificatore `path`, che dice a Nextflow che questo dovrebbe essere gestito come un percorso (include sia percorsi di directory che file).
Un altro qualificatore comune è `val`.

È importante notare che la definizione di output non _determina_ quale output verrà creato.
Semplicemente _dichiara_ qual è l'output atteso, così che Nextflow possa cercarlo una volta completata l'esecuzione.
Questo è necessario per verificare che il comando sia stato eseguito con successo e per passare l'output ai processi a valle se necessario. L'output prodotto che non corrisponde a quanto dichiarato nel blocco output non verrà passato ai processi a valle.

!!! warning "Avviso"

    Questo esempio è fragile perché abbiamo codificato in modo fisso il nome del file di output in due posti separati (lo script e i blocchi output).
    Se modifichiamo uno ma non l'altro, lo script si romperà.
    Più avanti imparerete modi per usare le variabili per mitigare questo problema.

In una pipeline reale, un process contiene di solito blocchi aggiuntivi come direttive e input, che introdurremo tra poco.

#### 1.1.2. La definizione del `workflow`

Il secondo blocco di codice descrive il **workflow** stesso.
La definizione del workflow inizia con la parola chiave `workflow`, seguita da un nome opzionale, poi dal corpo del workflow delimitato da parentesi graffe.

Qui abbiamo un **workflow** che consiste in un blocco `main:` (che dice 'questo è il corpo principale del workflow') contenente una chiamata al process `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emette un saluto
    sayHello()
}
```

Questa è una definizione di **workflow** molto minimale.
In una pipeline reale, il workflow contiene tipicamente multiple chiamate a **processi** connessi da **channel**, e i processi si aspettano uno o più **input** variabili.

Imparerete come aggiungere input variabili più avanti in questo modulo di formazione; e imparerete come aggiungere più processi e connetterli tramite channel nella Parte 3 di questo corso.

!!! tip "Suggerimento"

    Tecnicamente la riga `main:` non è richiesta per workflow semplici come questo, quindi potreste incontrare workflow che non ce l'hanno.
    Ma ne avremo bisogno per sfruttare gli output a livello di workflow, quindi tanto vale includerla dall'inizio.

### 1.2. Eseguire il workflow

Guardare il codice non è divertente quanto eseguirlo, quindi proviamolo in pratica.

#### 1.2.1. Avviare il workflow e monitorare l'esecuzione

Nel terminale, eseguite il seguente comando:

```bash
nextflow run hello-world.nf
```

??? success "Output del comando"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Se l'output della console assomiglia a questo, congratulazioni, avete appena eseguito il vostro primo workflow Nextflow!

L'output più importante qui è l'ultima riga, che è evidenziata nell'output sopra:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Questo ci dice che il process `sayHello` è stato eseguito con successo una volta (`1 of 1 ✔`).

È importante notare che questa riga indica anche dove trovare l'output della chiamata al process `sayHello`.
Guardiamolo ora.

#### 1.2.2. Trovare l'output e i log nella directory `work`

Quando eseguite Nextflow per la prima volta in una data directory, crea una directory chiamata `work` dove scriverà tutti i file (e tutti i symlink) generati nel corso dell'esecuzione.

All'interno della directory `work`, Nextflow organizza output e log per ogni chiamata di process.
Per ogni chiamata di process, Nextflow crea una sottodirectory nidificata, denominata con un hash per renderla unica, dove preparerà tutti gli input necessari (usando symlink per impostazione predefinita), scriverà file di supporto e scriverà log e tutti gli output del process.

Il percorso di quella sottodirectory è mostrato in forma troncata tra parentesi quadre nell'output della console.
Guardando cosa abbiamo ottenuto per l'esecuzione mostrata sopra, la riga di log della console per il process sayHello inizia con `[65/7be2fa]`. Questo corrisponde al seguente percorso di directory: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Diamo un'occhiata a cosa c'è dentro.

??? abstract "Contenuti della directory"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Non vedete la stessa cosa?"

    I nomi esatti delle sottodirectory saranno diversi sul vostro sistema.

    Se naviga nei contenuti della sottodirectory dell'attività nell'esplora file di VSCode, vedrà tutti i file immediatamente.
    Tuttavia, i file di log sono impostati per essere invisibili nel terminale, quindi se vuole usare `ls` o `tree` per visualizzarli, dovrà impostare l'opzione pertinente per visualizzare i file invisibili.

    ```bash
    tree -a work
    ```

La prima cosa che vorrete guardare è l'output effettivo del workflow, cioè il file `output.txt` prodotto dal process `sayHello`.
Apritelo e troverete il saluto `Hello World!`, che era lo scopo del nostro workflow minimalista.

??? abstract "Contenuti del file"

    ```console title="output.txt"
    Hello World!
    ```

Ha funzionato!

Certo, potrebbe sembrare molto codice di wrapper per un risultato così piccolo, ma il valore di tutto quel codice di wrapper diventerà più ovvio una volta che inizieremo a leggere file di input e a concatenare più step insieme.

Detto questo, diamo anche un'occhiata agli altri file in quella directory. Sono file di supporto e di log prodotti da Nextflow come parte dell'esecuzione dell'attività.

- **`.command.begin`**: Metadati relativi all'inizio dell'esecuzione della chiamata del process
- **`.command.err`**: Messaggi di errore (`stderr`) emessi dalla chiamata del process
- **`.command.log`**: Output di log completo emesso dalla chiamata del process
- **`.command.out`**: Output regolare (`stdout`) della chiamata del process
- **`.command.run`**: Script completo eseguito da Nextflow per eseguire la chiamata del process
- **`.command.sh`**: Il comando che è stato effettivamente eseguito dalla chiamata del process
- **`.exitcode`**: Il codice di uscita risultante dal comando

Il file `.command.sh` è particolarmente utile perché indica il comando principale che Nextflow ha eseguito, non includendo tutta la contabilità e la configurazione dell'attività/ambiente.

??? abstract "Contenuti del file"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Questo corrisponde a quanto abbiamo eseguito manualmente prima.

In questo caso è molto semplice perché il comando del process era codificato in modo fisso, ma più avanti nel corso vedrà comandi di process che coinvolgono qualche interpolazione di variabili.
Questo rende particolarmente prezioso poter vedere esattamente come Nextflow ha interpretato il codice e quale comando è stato prodotto quando si sta risolvendo un problema in un'esecuzione fallita.

### 1.3. Eseguire di nuovo il workflow

Provate a rieseguire il workflow alcune volte, poi guardate le directory delle attività sotto `work/`.

??? abstract "Contenuti della directory"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Vede che è stata creata una nuova sottodirectory con un set completo di file di output e log per ogni esecuzione.
Questo mostra che eseguire lo stesso workflow più volte non sovrascriverà i risultati delle esecuzioni precedenti.

### Takeaway

Sapete come decifrare un semplice script Nextflow, eseguirlo e trovare l'output e i file di log rilevanti nella directory work.

### Cosa c'è dopo?

Imparate come pubblicare gli output del workflow in una posizione più conveniente.

---

## 2. Pubblicare gli output

Come ha appena appreso, l'output prodotto dalla nostra pipeline è sepolto in una directory di lavoro a diversi livelli di profondità.
Questo è fatto di proposito; Nextflow ha il controllo di questa directory e non dovremmo interagire con essa.
Tuttavia, questo rende scomodo recuperare gli output che ci interessano.

Fortunatamente, Nextflow fornisce un modo per pubblicare gli output in una directory designata usando le [definizioni di output a livello di workflow](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Uso di base

Questo coinvolgerà due nuovi pezzi di codice:

1. Un blocco `publish:` all'interno del corpo del `workflow`, che dichiara gli output dei processi.
2. Un blocco `output` nello script che specifica le opzioni di output come modalità e posizione.

#### 2.1.1. Dichiarare l'output del process `sayHello`

Dobbiamo aggiungere un blocco `publish:` al corpo del workflow (stesso tipo di elemento di codice del blocco `main:`) e listare l'output del process `sayHello()`.

Nel file dello script del workflow `hello-world.nf`, aggiungete le seguenti righe di codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emette un saluto
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emette un saluto
        sayHello()
    }
    ```

Vedete che possiamo riferirci all'output del process semplicemente facendo `sayHello().out`, e assegnargli un nome arbitrario, `first_output`.

#### 2.1.2. Aggiungere un blocco `output:` allo script

Ora dobbiamo solo aggiungere il blocco `output:` dove verrà specificato il percorso della directory di output. Noti che questo nuovo blocco si trova **fuori** e **sotto** il blocco `workflow` all'interno dello script.

Nel file dello script del workflow `hello-world.nf`, aggiungete le seguenti righe di codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emette un saluto
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emette un saluto
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Possiamo usare questo per assegnare percorsi specifici a qualsiasi output di process dichiarato nel blocco `workflow`.
Più avanti imparerete modi per generare strutture di directory di output sofisticate, ma per ora stiamo semplicemente codificando in modo fisso un percorso minimo per semplicità.

#### 2.1.3. Eseguire il workflow

Ora eseguite lo script del workflow modificato:

```bash
nextflow run hello-world.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

L'output del terminale dovrebbe sembrare familiare. Esternamente, nulla è cambiato.

Tuttavia, controllate il vostro esplora file: questa volta, Nextflow ha creato una nuova directory chiamata `results/`.

??? abstract "Contenuti della directory"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

All'interno della directory `results`, troviamo un link simbolico al file `output.txt` prodotto nella directory work dal comando che abbiamo appena eseguito.

Questo ci permette di recuperare facilmente i file di output senza dover scavare nella sottodirectory work.

### 2.2. Impostare una posizione personalizzata

Avere una posizione predefinita è ottimo, ma potreste voler personalizzare dove vengono salvati i risultati e come sono organizzati.

Per esempio, potreste voler organizzare i vostri output in sottodirectory.
Il modo più semplice per farlo è assegnare percorsi di output specifici per ogni output.

#### 2.2.1. Modificare il percorso di output

Ancora una volta, modificare il comportamento di pubblicazione per un output specifico è davvero semplice.
Per impostare una posizione personalizzata, basta modificare il `path` di conseguenza:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Poiché questo è impostato a livello del singolo output, potete specificare posizioni e sottodirectory diverse per soddisfare le vostre esigenze.

#### 2.2.2. Eseguire di nuovo il workflow

Proviamo.

```bash
nextflow run hello-world.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Questa volta il risultato viene scritto nella sottodirectory specificata.

??? abstract "Contenuti della directory"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Vedete che il risultato dell'esecuzione precedente è ancora lì.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Potete usare quanti livelli di nidificazione desiderate.
È anche possibile usare il nome del process o altre variabili per nominare le directory usate per organizzare i risultati, ed è possibile cambiare il nome predefinito della directory di output di primo livello (che è controllata dalla variabile speciale `outputDir`).
Copriremo queste opzioni in formazioni successive.

### 2.3. Impostare la modalità di pubblicazione su copy

Per impostazione predefinita, gli output vengono pubblicati come link simbolici dalla directory `work`.
Questo significa che c'è un solo file nel filesystem.

Questo è ottimo quando si hanno a che fare con file molto grandi, per i quali non si vogliono memorizzare copie multiple.
Tuttavia, se eliminate la directory work in qualsiasi momento (copriremo le operazioni di pulizia a breve), perderete l'accesso al file.
Quindi dovete avere un piano per salvare copie di tutti i file importanti in un posto sicuro.

Un'opzione facile è cambiare la modalità di pubblicazione su copy per gli output che vi interessano.

#### 2.3.1. Aggiungere la direttiva mode

Questa parte è davvero semplice.
Basta aggiungere `mode 'copy'` alla definizione di output a livello di workflow pertinente:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Questo imposta la modalità di pubblicazione per quello specifico output.

#### 2.3.2. Eseguire di nuovo il workflow

Proviamo.

```bash
nextflow run hello-world.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Questa volta, se guarda i risultati, il file è una copia vera e propria invece di un semplice symlink.

??? abstract "Contenuti della directory"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Poiché anche questo è impostato a livello del singolo output, vi permette di impostare la modalità di pubblicazione in modo granulare.
Questo sarà particolarmente utile più avanti quando passeremo a pipeline multi-step, dove potrebbe voler copiare solo gli output finali e lasciare gli output intermedi come symlink, per esempio.

Come notato prima, ci sono altre opzioni più sofisticate per controllare come vengono pubblicati gli output.
Vi mostreremo come usarle a tempo debito nel vostro percorso con Nextflow.

### 2.4. Nota sulle direttive `publishDir` a livello di process

Fino a molto recentemente, il modo stabilito per pubblicare gli output era farlo a livello di ogni singolo process usando una direttiva `publishDir`.

Per ottenere quanto abbiamo appena fatto per gli output del process `sayHello`, avremmo invece aggiunto la seguente riga alla definizione del process:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Troverete ancora questo pattern di codice ovunque nelle pipeline Nextflow più vecchie e nei moduli di process, quindi è importante esserne consapevoli.
Tuttavia, non raccomandiamo di usarlo in nessun nuovo lavoro poiché sarà eventualmente vietato nelle versioni future del linguaggio Nextflow.

### Takeaway

Sapete come pubblicare gli output del workflow in una posizione più conveniente.

### Cosa c'è dopo?

Imparate a fornire un input variabile tramite un parametro da riga di comando e utilizzare efficacemente i valori predefiniti.

---

## 3. Usare un input variabile passato dalla riga di comando

Nel suo stato attuale, il nostro workflow usa un saluto codificato in modo fisso nel comando del process.
Vogliamo aggiungere un po' di flessibilità usando una variabile di input, così da poter cambiare più facilmente il saluto a runtime.

Questo richiede di apportare tre serie di modifiche al nostro script:

1. Modificare il process per aspettarsi un input variabile
2. Impostare un parametro da riga di comando per catturare l'input dell'utente
3. Passare l'input al process nel corpo del workflow

Facciamo queste modifiche una alla volta.

### 3.1. Modificare il process `sayHello` per aspettarsi un input variabile

Dobbiamo modificare la definizione del process per (1) accettare una variabile di input e (2) usare quella variabile nella riga di comando.

#### 3.1.1. Aggiungere un blocco input alla definizione del process

Prima, adattiamo la definizione del process per accettare un input chiamato `greeting`.

Nel blocco del process, fate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

La variabile `greeting` è prefissata da `val` per dire a Nextflow che è un valore (non un percorso).

#### 3.1.2. Modificare il comando del process per usare la variabile di input

Ora scambiamo il valore originale codificato in modo fisso con il valore della variabile di input che ci aspettiamo di ricevere.

Nel blocco del process, fate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

Il simbolo `$` e le parentesi graffe (`{ }`) dicono a Nextflow che questo è un nome di variabile che deve essere sostituito con il valore di input effettivo (=interpolato).

!!! tip "Suggerimento"

    Le parentesi graffe (`{ }`) erano tecnicamente opzionali nelle versioni precedenti di Nextflow, quindi potreste vedere workflow più vecchi dove questo è scritto come `echo '$greeting' > output.txt`.

Ora che il process `sayHello()` è pronto ad accettare un input variabile, abbiamo bisogno di un modo per fornire un valore di input alla chiamata del process a livello di workflow.

### 3.2. Impostare un parametro da riga di comando per catturare l'input dell'utente

Potremmo semplicemente codificare in modo fisso un input direttamente facendo la chiamata al process `sayHello('Hello World!')`.
Tuttavia, quando facciamo lavoro reale con il nostro workflow, vorremo essere in grado di controllare i suoi input dalla riga di comando.

Buone notizie: Nextflow ha un sistema di parametri di workflow integrato chiamato `params`, che rende facile dichiarare e usare parametri CLI.

La sintassi generale è dichiarare `params.<nome_parametro>` per dire a Nextflow di aspettarsi un parametro `--<nome_parametro>` sulla riga di comando.

Qui vogliamo creare un parametro chiamato `--input`, quindi dobbiamo dichiarare `params.input` da qualche parte nel workflow.
In linea di principio possiamo scriverlo ovunque; ma poiché vorremo darlo alla chiamata del process `sayHello()`, possiamo inserirlo direttamente lì scrivendo `sayHello(params.input)`.

Nel blocco del workflow, fate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emette un saluto
    sayHello(params.input)
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emette un saluto
    sayHello()
    ```

Questo dice a Nextflow di eseguire il process `sayHello` sul valore fornito tramite il parametro `--input`.

In effetti, abbiamo compiuto i passaggi (2) e (3) delineati all'inizio della sezione in un colpo solo.

### 3.3. Eseguire il comando del workflow

Eseguiamolo!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Se avete fatto tutte queste modifiche correttamente, dovreste ottenere un'altra esecuzione riuscita.

Assicuratevi di aprire il file di output per verificare che ora abbiate la nuova versione del saluto.

??? abstract "Contenuti del file"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà!

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Noti come la nuova esecuzione ha sovrascritto il file di output pubblicato nella directory `results`.
Tuttavia, i risultati delle esecuzioni precedenti sono ancora conservati nelle directory delle attività sotto `work`.

!!! tip "Suggerimento"

    Potete facilmente distinguere i parametri a livello di Nextflow dai parametri a livello di pipeline.

    - I parametri che si applicano a una pipeline prendono sempre un doppio trattino (`--`).
    - I parametri che modificano un'impostazione di Nextflow, _ad es._ la funzione `-resume` che abbiamo usato prima, prendono un singolo trattino (`-`).

### 3.4. Usare valori predefiniti per i parametri da riga di comando

Ok, quello era conveniente, ma in molti casi ha senso fornire un valore predefinito per un dato parametro così da non doverlo specificare per ogni esecuzione.

#### 3.4.1. Impostare un valore predefinito per il parametro CLI

Diamo al parametro `input` un valore predefinito dichiarandolo prima della definizione del workflow.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String = 'Holà mundo!'
}
```

Come vedete, possiamo specificare il tipo di input che il workflow si aspetta (Nextflow 25.10.2 e successivi).
La sintassi è `nome: Tipo = valore_predefinito`.
I tipi supportati includono `String`, `Integer`, `Float`, `Boolean` e `Path`.

!!! info "Nota"

    Nei workflow più vecchi, potreste vedere quel intero blocco `params` scritto semplicemente come `input = 'Holà mundo!'`.

Man mano che aggiungete più parametri alla vostra pipeline, dovreste aggiungerli tutti a questo blocco, che debbano o meno avere un valore predefinito.
Questo renderà facile trovare tutti i parametri configurabili a colpo d'occhio.

#### 3.4.2. Eseguire di nuovo il workflow senza specificare il parametro

Ora che avete un valore predefinito impostato, potete eseguire di nuovo il workflow senza dover specificare un valore nella riga di comando.

```bash
nextflow run hello-world.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

L'output sarà nella stessa posizione di prima, ma i contenuti dovrebbero essere aggiornati con il nuovo testo.

??? abstract "Contenuti del file"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow ha usato il valore predefinito del parametro greeting per creare l'output.

#### 3.4.3. Sovrascrivere il valore predefinito

Se fornite il parametro sulla riga di comando, il valore CLI sovrascriverà il valore predefinito.

Provate:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Ancora una volta, dovreste trovare l'output aggiornato corrispondente nella vostra directory results.

??? abstract "Contenuti del file"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Nota"

    In Nextflow, ci sono più posti dove può specificare valori per i parametri.
    Se lo stesso parametro è impostato su valori diversi in più posti, Nexflow determinerà quale valore usare in base all'ordine di precedenza descritto [qui](https://www.nextflow.io/docs/latest/config.html).

    Copriremo questo in maggior dettaglio nella Parte 6 (Configuration).

### Takeaway

Sapete come usare un semplice input variabile fornito a runtime tramite un parametro da riga di comando, così come impostare, usare e sovrascrivere valori predefiniti.

### Cosa c'è dopo?

Imparate a gestire le esecuzioni in modo più conveniente.

---

## 4. Gestire le esecuzioni del workflow

Sapere come avviare workflow e recuperare output è ottimo, ma scoprirete rapidamente che ci sono alcuni altri aspetti della gestione del workflow che vi renderanno la vita più facile, specialmente se state sviluppando i vostri workflow.

Qui vi mostriamo come usare la funzione `resume` per quando dovete ri-avviare lo stesso workflow, come ispezionare il log delle esecuzioni passate con `nextflow log`, e come eliminare le directory work più vecchie con `nextflow clean`.

<!-- Any other cool options we should include? Added log -->

### 4.1. Ri-avviare un workflow con `-resume`

A volte vorrà ri-eseguire una pipeline che ha già avviato in precedenza senza rifare alcuno step che è già stato completato con successo.

Nextflow ha un'opzione chiamata `-resume` che vi permette di fare questo.
Specificamente, in questa modalità, tutti i processi che sono già stati eseguiti con lo stesso identico codice, impostazioni e input verranno saltati.
Questo significa che Nextflow eseguirà solo i processi che ha aggiunto o modificato dalla precedente esecuzione, o a cui sta fornendo nuove impostazioni o input.

Ci sono due vantaggi chiave nel fare questo:

- Se state sviluppando la vostra pipeline, potete iterare più rapidamente poiché dovete solo eseguire il/i process su cui state lavorando attivamente per testare le vostre modifiche.
- Se state eseguendo una pipeline in produzione e qualcosa va storto, in molti casi potete risolvere il problema e ri-avviare la pipeline, e riprenderà dal punto di fallimento, il che può farvi risparmiare molto tempo e calcolo.

Per usarlo, aggiungete semplicemente `-resume` al vostro comando ed eseguite:

```bash
nextflow run hello-world.nf -resume
```

??? success "Output del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

L'output della console dovrebbe sembrare familiare, ma c'è una cosa leggermente diversa rispetto a prima.

Cercate la parte `cached:` che è stata aggiunta nella riga di stato del process (riga 5), che significa che Nextflow ha riconosciuto che ha già fatto questo lavoro e ha semplicemente riutilizzato il risultato dall'esecuzione precedente riuscita.

Potete anche vedere che l'hash della sottodirectory work è lo stesso dell'esecuzione precedente.
Nextflow vi sta letteralmente indicando l'esecuzione precedente e dicendo "Ho già fatto quello lì."

!!! tip "Suggerimento"

    Quando ri-esegue una pipeline con `resume`, Nextflow non sovrascrive alcun file pubblicato fuori dalla directory work da esecuzioni che sono state eseguite con successo in precedenza.

### 4.2. Ispezionare il log delle esecuzioni passate

Che stiate sviluppando una nuova pipeline o eseguendo pipeline in produzione, a un certo punto avrete probabilmente bisogno di cercare informazioni sulle esecuzioni passate.
Ecco come fare.

Ogni volta che avvia un workflow nextflow, una riga viene scritta in un file di log chiamato `history`, sotto una directory nascosta chiamata `.nextflow` nella directory di lavoro corrente.

??? abstract "Contenuti del file"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Questo file fornisce il timestamp, il nome dell'esecuzione, lo stato, l'ID di revisione, l'ID di sessione e la riga di comando completa per ogni esecuzione Nextflow che è stata avviata dalla directory di lavoro corrente.

Un modo più conveniente per accedere a queste informazioni è usare il comando `nextflow log`.

```bash
nextflow log
```

??? success "Output del comando"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Questo produrrà i contenuti del file di log nel terminale, arricchiti con una riga di intestazione.

Noterete che l'ID di sessione cambia ogni volta che eseguite un nuovo comando `nextflow run`, ECCETTO se state usando l'opzione `-resume`.
In quel caso, l'ID di sessione rimane lo stesso.

Nextflow usa l'ID di sessione per raggruppare le informazioni di cache dell'esecuzione nella directory `cache`, anch'essa situata sotto `.nextflow`.

### 4.3. Eliminare le directory work più vecchie

Durante il processo di sviluppo, tipicamente eseguirete la vostra bozza di pipeline un gran numero di volte, il che può portare a un accumulo di molti file in molte sottodirectory.

Fortunatamente Nextflow include un utile sottocomando `clean` che potete eliminare automaticamente le sottodirectory work per le esecuzioni passate che non Le interessano più.

#### 4.3.1. Determinare i criteri di eliminazione

Ci sono multiple [opzioni](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per determinare cosa eliminare.

Qui vi mostriamo un esempio che elimina tutte le sottodirectory dalle esecuzioni precedenti a una data esecuzione, specificata usando il suo nome di esecuzione.

Cercate l'esecuzione riuscita più recente dove non avete usato `-resume`; nel nostro caso il nome dell'esecuzione era `golden_cantor`.

Il nome dell'esecuzione è la stringa in due parti generata dalla macchina mostrata tra parentesi quadre nell'output della console `Launching (...)`.
Potete anche usare il log di Nextflow per cercare un'esecuzione in base al suo timestamp e/o riga di comando.

#### 4.3.2. Fare una prova a secco

Prima usiamo il flag di prova a secco `-n` per verificare cosa verrà eliminato dato il comando:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Output del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Il vostro output avrà nomi di directory di attività diversi e potrebbe avere un numero diverso di righe, ma dovrebbe sembrare simile all'esempio.

Se non vedete alcuna riga di output, o non avete fornito un nome di esecuzione valido o non ci sono esecuzioni passate da eliminare. Assicuratevi di cambiare `golden_cantor` nel comando di esempio con qualunque sia il nome dell'esecuzione più recente corrispondente nel vostro log.

#### 4.3.3. Procedere con l'eliminazione

Se l'output sembra come previsto e volete procedere con l'eliminazione, ri-eseguite il comando con il flag `-f` invece di `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Output del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

L'output dovrebbe essere simile a prima, ma ora dicendo 'Removed' invece di 'Would remove'.
Noti che questo non rimuove le sottodirectory a due caratteri (come `a3/` sopra) ma ne svuota il contenuto.

!!! Warning "Avviso"

    Eliminare le sottodirectory work dalle esecuzioni passate le rimuove dalla cache di Nextflow e cancella tutti gli output che erano memorizzati in quelle directory.
    Questo significa che interrompe la capacità di Nextflow di riprendere l'esecuzione senza ri-eseguire i processi corrispondenti.

    È vostra responsabilità salvare tutti gli output che vi interessano o su cui contate! Questo è il motivo principale per cui preferiamo usare la modalità `copy` piuttosto che la modalità `symlink` per la direttiva `publish`.

### Takeaway

Sapete come pubblicare output in una directory specifica, ri-avviare una pipeline senza ripetere step che erano già stati eseguiti in modo identico, e usare il comando `nextflow clean` per pulire le directory work vecchie.

Più in generale, sapete come interpretare un semplice workflow Nextflow, gestire la sua esecuzione e recuperare gli output.

### Cosa c'è dopo?

Prendetevi una piccola pausa, ve la siete meritata!

Quando siete pronti, passate alla [**Parte 2: Hello Channels**](./02_hello_channels.md) per imparare come usare i channel per alimentare gli input nel vostro workflow, il che vi permetterà di sfruttare il parallelismo del flusso dati integrato di Nextflow e altre potenti funzionalità.

---

## Quiz

<quiz>
Quali sono i componenti minimi richiesti di un process Nextflow?
- [ ] Solo blocchi input e output
- [x] Blocchi output e script
- [ ] Blocchi input, output e script
- [ ] Solo un blocco script

Approfondisci: [1.1.1. La definizione del process](#111-la-definizione-del-process)
</quiz>

<quiz>
Qual è lo scopo del blocco output in un process?
- [ ] Stampare i risultati sulla console
- [ ] Salvare i file nella directory work
- [x] Dichiarare gli output attesi dal process
- [ ] Definire variabili d'ambiente

Approfondisci: [1.1.1. La definizione del process](#111-la-definizione-del-process)
</quiz>

<quiz>
Quale comando viene usato per eseguire un workflow Nextflow?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Guardando la directory work di un'attività, quale file contiene il comando effettivo che è stato eseguito?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Approfondisci: [1.2.2. Trovare l'output e i log nella directory `work`](#122-trovare-loutput-e-i-log-nella-directory-work)
</quiz>

<quiz>
Cosa fa il flag `-resume`?
- [ ] Riavvia il workflow dall'inizio
- [ ] Mette in pausa il workflow
- [x] Salta i processi che sono già stati completati con successo
- [ ] Crea un backup del workflow

Approfondisci: [4.1. Ri-avviare un workflow con `-resume`](#41-ri-avviare-un-workflow-con--resume)
</quiz>

<quiz>
Qual è la modalità predefinita per pubblicare gli output del workflow?
- [ ] Copiare i file nella directory di output
- [x] Creare link simbolici nella directory di output
- [ ] Spostare i file nella directory di output
- [ ] Comprimere i file nella directory di output

Approfondisci: [2.3. Impostare la modalità di pubblicazione su copy](#23-impostare-la-modalità-di-pubblicazione-su-copy)
</quiz>

<quiz>
Come si passa un valore di parametro a un workflow Nextflow dalla riga di comando?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Approfondisci: [3.2. Impostare un parametro da riga di comando per catturare l'input dell'utente](#32-impostare-un-parametro-da-riga-di-comando-per-catturare-linput-dellutente)
</quiz>

<quiz>
Come si fa riferimento a una variabile all'interno di un blocco script Nextflow?
- [ ] Usare la sintassi `%variable%`
- [x] Usare la sintassi `#!groovy ${variable}`
- [ ] Usare la sintassi `{{variable}}`
- [ ] Usare la sintassi `[variable]`
</quiz>
