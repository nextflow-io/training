# Parte 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guardate [l'intera playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/01_hello_world.md).
///

In questa prima parte del corso di formazione Hello Nextflow, ci addentriamo nell'argomento con un esempio Hello World molto semplice e indipendente dal dominio, che svilupperemo progressivamente per dimostrare l'uso della logica e dei componenti fondamentali di Nextflow.

??? info "Cos'è un esempio Hello World?"

    Un "Hello World!" è un esempio minimalista pensato per dimostrare la sintassi e la struttura di base di un linguaggio di programmazione o framework software.
    L'esempio consiste tipicamente nello stampare la frase "Hello, World!" su un dispositivo di output, come la console o il terminale, oppure nello scriverla in un file.

---

## 0. Riscaldamento: Eseguire un esempio Hello World direttamente

Dimostriamo questo con un semplice comando che eseguiamo direttamente nel terminale, per mostrare cosa fa prima di inserirlo in Nextflow.

!!! tip

    Ricordate che ora dovreste trovarvi all'interno della directory `hello-nextflow/` come descritto nella pagina [Getting Started](00_orientation.md).

### 0.1. Fare in modo che il terminale dica hello

Eseguite il seguente comando nel vostro terminale.

```bash
echo 'Hello World!'
```

??? success "Output del comando"

    ```console
    Hello World!
    ```

Questo stampa il testo 'Hello World' direttamente nel terminale.

### 0.2. Scrivere l'output in un file

L'esecuzione di pipeline consiste principalmente nel leggere dati da file e scrivere risultati in altri file, quindi modifichiamo il comando per scrivere l'output di testo in un file per rendere l'esempio un po' più rilevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Output del comando"

    ```console

    ```

Questo non produce alcun output nel terminale.

### 0.3. Trovare l'output

Il testo 'Hello World' dovrebbe ora trovarsi nel file di output che abbiamo specificato, chiamato `output.txt`.
Potete aprirlo nell'esploratore di file o dalla riga di comando usando l'utility `cat`, per esempio.

??? abstract "Contenuto del file"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Questo è ciò che cercheremo di replicare con il nostro primo flusso di lavoro Nextflow.

### Takeaway

Ora sapete come eseguire un semplice comando nel terminale che produce del testo in output e, opzionalmente, come farlo scrivere l'output in un file.

### Cosa c'è dopo?

Scoprite come apparirebbe scritto come flusso di lavoro Nextflow.

---

## 1. Esaminare lo script ed eseguirlo

Vi forniamo uno script di flusso di lavoro completamente funzionale, anche se minimalista, chiamato `hello-world.nf` che fa la stessa cosa di prima (scrivere 'Hello World!') ma con Nextflow.

Per iniziare, apriamo lo script del flusso di lavoro in modo che possiate farvi un'idea di come è strutturato.
Poi lo eseguiremo e cercheremo i suoi output.

### 1.1. Esaminare il codice

Troverete lo script `hello-world.nf` nella vostra directory corrente, che dovrebbe essere `hello-nextflow`. Apritelo nel pannello dell'editor.

??? full-code "File di codice completo"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usa echo per stampare 'Hello World!' in un file
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
        // emit a greeting
        sayHello()
    }
    ```

Uno script di flusso di lavoro Nextflow include tipicamente una o più definizioni di [**process**](https://nextflow.io/docs/latest/process.html) e il [**workflow**](https://nextflow.io/docs/latest/workflow.html) stesso, più alcuni blocchi opzionali (non presenti qui) che introdurremo più avanti.

Ogni **process** descrive quali operazioni deve compiere il passo corrispondente nella pipeline, mentre il **workflow** descrive la logica di flusso dati che collega i vari passi.

Esamineremo prima più da vicino il blocco **process**, poi guarderemo il blocco **workflow**.

#### 1.1.1. La definizione del `process`

Il primo blocco di codice descrive un **process**.

La definizione del processo inizia con la parola chiave `process`, seguita dal nome del processo e infine dal corpo del processo delimitato da parentesi graffe.
Il corpo del processo deve contenere un blocco script che specifica il comando da eseguire, che può essere qualsiasi cosa possiate eseguire in un terminale a riga di comando.

```groovy title="hello-world.nf" linenums="3"
/*
* Usa echo per stampare 'Hello World!' in un file
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

Qui abbiamo un **process** chiamato `sayHello` che scrive il suo **output** in un file chiamato `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Questa è una definizione di processo molto minimale che contiene solo una definizione di `output` e lo `script` da eseguire.

La definizione di `output` include il qualificatore `path`, che dice a Nextflow che questo dovrebbe essere gestito come un percorso (include sia percorsi di directory che file).
Un altro qualificatore comune è `val`.

È importante notare che la definizione di output non _determina_ quale output verrà creato.
Semplicemente _dichiara_ qual è l'output atteso, in modo che Nextflow possa cercarlo una volta completata l'esecuzione.
Questo è necessario per verificare che il comando sia stato eseguito con successo e per passare l'output ai processi a valle se necessario. L'output prodotto che non corrisponde a quanto dichiarato nel blocco output non verrà passato ai processi a valle.

!!! warning

    Questo esempio è fragile perché abbiamo codificato il nome del file di output in due posti separati (i blocchi script e output).
    Se cambiamo uno ma non l'altro, lo script si romperà.
    Più avanti imparerete modi per usare variabili per mitigare questo problema.

In una pipeline reale, un processo contiene solitamente blocchi aggiuntivi come direttive e input, che introdurremo tra poco.

#### 1.1.2. La definizione del `workflow`

Il secondo blocco di codice descrive il **workflow** stesso.
La definizione del workflow inizia con la parola chiave `workflow`, seguita da un nome opzionale, poi dal corpo del workflow delimitato da parentesi graffe.

Qui abbiamo un **workflow** che consiste in un blocco `main:` (che dice 'questo è il corpo principale del flusso di lavoro') contenente una chiamata al processo `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emit a greeting
    sayHello()
}
```

Questa è una definizione di **workflow** molto minimale.
In una pipeline reale, il workflow contiene tipicamente chiamate multiple a **processi** collegati da **canali**, e i processi si aspettano uno o più **input** variabili.

Imparerete come aggiungere input variabili più avanti in questo modulo di formazione; e imparerete come aggiungere più processi e collegarli tramite canali nella Parte 3 di questo corso.

!!! tip

    Tecnicamente la riga `main:` non è richiesta per flussi di lavoro semplici come questo, quindi potreste incontrare flussi di lavoro che non ce l'hanno.
    Ma ne avremo bisogno per sfruttare gli output a livello di workflow, quindi tanto vale includerla dall'inizio.

### 1.2. Eseguire il flusso di lavoro

Guardare il codice non è così divertente come eseguirlo, quindi proviamolo in pratica.

#### 1.2.1. Lanciare il flusso di lavoro e monitorare l'esecuzione

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

Se l'output della vostra console assomiglia a quello, allora congratulazioni, avete appena eseguito il vostro primo flusso di lavoro Nextflow!

L'output più importante qui è l'ultima riga, che è evidenziata nell'output sopra:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Questo ci dice che il processo `sayHello` è stato eseguito con successo una volta (`1 of 1 ✔`).

È importante notare che questa riga vi dice anche dove trovare l'output della chiamata al processo `sayHello`.
Guardiamolo ora.

#### 1.2.2. Trovare l'output e i registri nella directory `work`

Quando eseguite Nextflow per la prima volta in una determinata directory, crea una directory chiamata `work` dove scriverà tutti i file (e qualsiasi symlink) generati nel corso dell'esecuzione.

All'interno della directory `work`, Nextflow organizza output e registri per chiamata di processo.
Per ogni chiamata di processo, Nextflow crea una sottodirectory nidificata, denominata con un hash per renderla unica, dove preparerà tutti gli input necessari (usando symlink per impostazione predefinita), scriverà file di supporto e scriverà registri e qualsiasi output del processo.

Il percorso a quella sottodirectory è mostrato in forma troncata tra parentesi quadre nell'output della console.
Guardando ciò che abbiamo ottenuto per l'esecuzione mostrata sopra, la riga di log della console per il processo sayHello inizia con `[65/7be2fa]`. Questo corrisponde al seguente percorso di directory: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Diamo un'occhiata a cosa c'è lì dentro.

??? abstract "Directory contents"

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

    Se navigate i contenuti della sottodirectory dell'attività nell'esploratore file di VSCode, vedrete tutti i file immediatamente.
    Tuttavia, i file di log sono impostati per essere invisibili nel terminale, quindi se volete usare `ls` o `tree` per visualizzarli, dovrete impostare l'opzione rilevante per visualizzare i file invisibili.

    ```bash
    tree -a work
    ```

La prima cosa che volete guardare è l'output effettivo del flusso di lavoro, cioè il file `output.txt` prodotto dal processo `sayHello`.
Apritelo e troverete il saluto `Hello World!`, che era lo scopo del nostro flusso di lavoro minimalista.

??? abstract "Contenuto del file"

    ```console title="output.txt"
    Hello World!
    ```

Ha funzionato!

Certo, potrebbe sembrare molto codice wrapper per un risultato così piccolo, ma il valore di tutto quel codice wrapper diventerà più ovvio una volta che inizieremo a leggere file di input e a concatenare più passi.

Detto questo, guardiamo anche gli altri file in quella directory. Questi sono file di supporto e di log prodotti da Nextflow come parte dell'esecuzione dell'attività.

- **`.command.begin`**: Metadati relativi all'inizio dell'esecuzione della chiamata al processo
- **`.command.err`**: Messaggi di errore (`stderr`) emessi dalla chiamata al processo
- **`.command.log`**: Output di log completo emesso dalla chiamata al processo
- **`.command.out`**: Output regolare (`stdout`) della chiamata al processo
- **`.command.run`**: Script completo eseguito da Nextflow per eseguire la chiamata al processo
- **`.command.sh`**: Il comando che è stato effettivamente eseguito dalla chiamata al processo
- **`.exitcode`**: Il codice di uscita risultante dal comando

Il file `.command.sh` è particolarmente utile perché vi dice il comando principale che Nextflow ha eseguito, senza includere tutta la contabilità e la configurazione dell'attività/ambiente.

??? abstract "Contenuto del file"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Questo corrisponde a ciò che abbiamo eseguito prima manualmente.

In questo caso è molto semplice perché il comando del processo era codificato, ma più avanti nel corso vedrete comandi di processo che coinvolgono un po' di interpolazione di variabili.
Questo rende particolarmente prezioso poter vedere esattamente come Nextflow ha interpretato il codice e quale comando è stato prodotto quando state risolvendo problemi di un'esecuzione fallita.

### 1.3. Eseguire nuovamente il flusso di lavoro

Provate a rieseguire il flusso di lavoro alcune volte, poi guardate le directory delle attività sotto `work/`.

??? abstract "Directory contents"

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

Vedete che è stata creata una nuova sottodirectory con un set completo di file di output e log per ogni esecuzione.
Questo vi mostra che eseguire lo stesso flusso di lavoro più volte non sovrascriverà i risultati delle esecuzioni precedenti.

### Takeaway

Sapete come decifrare un semplice script Nextflow, eseguirlo e trovare l'output e i file di log rilevanti nella directory work.

### Cosa c'è dopo?

Imparate come pubblicare gli output del flusso di lavoro in una posizione più comoda.

---

## 2. Pubblicare gli output

Come avete appena imparato, l'output prodotto dalla nostra pipeline è sepolto in una directory di lavoro a diversi livelli di profondità.
Questo è fatto apposta; Nextflow ha il controllo di questa directory e non dovremmo interagire con essa.
Tuttavia, questo rende scomodo recuperare gli output che ci interessano.

Fortunatamente, Nextflow fornisce un modo per pubblicare gli output in una directory designata usando le [definizioni di output del workflow](https://nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Uso di base

Questo coinvolgerà due nuovi pezzi di codice:

1. Un blocco `publish:` all'interno del corpo del `workflow`, che dichiara gli output del processo.
2. Un blocco `output` nello script che specifica opzioni di output come modalità e posizione.

#### 2.1.1. Dichiarare l'output del processo `sayHello`

Dobbiamo aggiungere un blocco `publish:` al corpo del workflow (stesso tipo di elemento di codice del blocco `main:`) e elencare l'output del processo `sayHello()`.

Nel file dello script del flusso di lavoro `hello-world.nf`, aggiungete le seguenti righe di codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emit a greeting
        sayHello()
    }
    ```

Vedete che possiamo riferirci all'output del processo semplicemente facendo `sayHello().out`, e assegnargli un nome arbitrario, `first_output`.

#### 2.1.2. Aggiungere un blocco `output:` allo script

Ora dobbiamo solo aggiungere il blocco `output:` dove verrà specificato il percorso della directory di output. Notate che questo nuovo blocco si trova **fuori** e **sotto** il blocco `workflow` all'interno dello script.

Nel file dello script del flusso di lavoro `hello-world.nf`, aggiungete le seguenti righe di codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emit a greeting
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
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Possiamo usare questo per assegnare percorsi specifici a qualsiasi output di processo dichiarato nel blocco `workflow`.
Più avanti imparerete modi per generare strutture di directory di output sofisticate, ma per ora stiamo solo codificando un percorso minimale per semplicità.

#### 2.1.3. Eseguire il flusso di lavoro

Ora eseguite lo script del flusso di lavoro modificato:

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

Tuttavia, controllate il vostro esploratore file: questa volta, Nextflow ha creato una nuova directory chiamata `results/`.

??? abstract "Directory contents"

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

All'interno della directory `results`, troviamo un link simbolico all'`output.txt` prodotto nella directory work dal comando che abbiamo appena eseguito.

Questo ci permette di recuperare facilmente i file di output senza dover scavare nella sottodirectory work.

### 2.2. Impostare una posizione personalizzata

Avere una posizione predefinita è ottimo, ma potreste voler personalizzare dove vengono salvati i risultati e come sono organizzati.

Per esempio, potreste voler organizzare i vostri output in sottodirectory.
Il modo più semplice per farlo è assegnare un percorso di output specifico per output.

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

#### 2.2.2. Eseguire nuovamente il flusso di lavoro

Proviamolo.

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

??? abstract "Directory contents"

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

Potete usare tutti i livelli di nidificazione che volete.
È anche possibile usare il nome del processo o altre variabili per nominare le directory usate per organizzare i risultati, ed è possibile cambiare il nome predefinito della directory di output di livello superiore (che è controllato dal flag CLI `-o` o dalla variabile di configurazione `outputDir`).
Copriremo queste opzioni più avanti nella formazione.

### 2.3. Impostare la modalità di pubblicazione su copy

Per impostazione predefinita, gli output sono pubblicati come link simbolici dalla directory `work`.
Questo significa che c'è solo un singolo file sul filesystem.

Questo è ottimo quando si ha a che fare con file molto grandi, per i quali non si vogliono memorizzare copie multiple.
Tuttavia, se eliminate la directory work in qualsiasi momento (copriremo le operazioni di pulizia a breve), perderete l'accesso al file.
Quindi dovete avere un piano per salvare copie di qualsiasi file importante in un luogo sicuro.

Un'opzione facile è cambiare la modalità di pubblicazione in copy per gli output che vi interessano.

#### 2.3.1. Aggiungere la direttiva mode

Questo bit è davvero semplice.
Basta aggiungere `mode 'copy'` alla definizione di output a livello di workflow rilevante:

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

Questo imposta la modalità di pubblicazione per quell'output specifico.

#### 2.3.2. Eseguire nuovamente il flusso di lavoro

Proviamolo.

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

Questa volta, se guardate i risultati, il file è una copia vera e propria invece di un semplice symlink.

??? abstract "Directory contents"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Poiché anche questo è impostato a livello del singolo output, vi permette di impostare la modalità di pubblicazione in modo granulare.
Questo tornerà particolarmente utile più avanti quando passeremo a pipeline multi-step, dove potreste voler copiare solo gli output finali e lasciare gli output intermedi come symlink, per esempio.

Come notato in precedenza, ci sono altre opzioni più sofisticate per controllare come vengono pubblicati gli output.
Vi mostreremo come usarle a tempo debito nel vostro percorso con Nextflow.

### 2.4. Nota sulle direttive `publishDir` a livello di processo

Fino a molto recentemente, il modo consolidato per pubblicare gli output era farlo a livello di ogni singolo processo usando una direttiva `publishDir`.

Per ottenere ciò che abbiamo appena fatto per gli output del processo `sayHello`, avremmo invece aggiunto la seguente riga alla definizione del processo:

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

Troverete ancora questo pattern di codice ovunque in pipeline Nextflow e moduli di processo più vecchi, quindi è importante esserne consapevoli.
Tuttavia, non raccomandiamo di usarlo in alcun nuovo lavoro poiché sarà eventualmente vietato nelle versioni future del linguaggio Nextflow.

### Takeaway

Sapete come pubblicare gli output del flusso di lavoro in una posizione più comoda.

### Cosa c'è dopo?

Imparate a fornire un input variabile tramite un parametro da riga di comando e utilizzare i valori predefiniti in modo efficace.

---

## 3. Usare un input variabile passato dalla riga di comando

Nel suo stato attuale, il nostro flusso di lavoro usa un saluto codificato nel comando del processo.
Vogliamo aggiungere un po' di flessibilità usando una variabile di input, in modo da poter cambiare più facilmente il saluto a runtime.

Questo richiede di apportare tre serie di modifiche al nostro script:

1. Cambiare il processo per aspettarsi un input variabile
2. Impostare un parametro da riga di comando per catturare l'input dell'utente
3. Passare l'input al processo nel corpo del workflow

Facciamo queste modifiche una alla volta.

### 3.1. Cambiare il processo `sayHello` per aspettarsi un input variabile

Dobbiamo modificare la definizione del processo per (1) accettare una variabile di input e (2) usare quella variabile nella riga di comando.

#### 3.1.1. Aggiungere un blocco input alla definizione del processo

Prima, adattiamo la definizione del processo per accettare un input chiamato `greeting`.

Nel blocco process, fate la seguente modifica al codice:

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

La variabile `greeting` è preceduta da `val` per dire a Nextflow che è un valore (non un percorso).

#### 3.1.2. Modificare il comando del processo per usare la variabile di input

Ora sostituiamo il valore originale codificato con il valore della variabile di input che ci aspettiamo di ricevere.

Nel blocco process, fate la seguente modifica al codice:

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

!!! tip

    Le parentesi graffe (`{ }`) erano tecnicamente opzionali nelle versioni precedenti di Nextflow, quindi potreste vedere flussi di lavoro più vecchi dove questo è scritto come `echo '$greeting' > output.txt`.

Ora che il processo `sayHello()` è pronto ad accettare un input variabile, abbiamo bisogno di un modo per fornire un valore di input alla chiamata del processo a livello di workflow.

### 3.2. Impostare un parametro da riga di comando per catturare l'input dell'utente

Potremmo semplicemente codificare un input direttamente facendo la chiamata al processo `sayHello('Hello World!')`.
Tuttavia, quando facciamo lavoro reale con il nostro flusso di lavoro, vorremo essere in grado di controllare i suoi input dalla riga di comando, quindi possiamo fare qualcosa del genere:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Fortunatamente, Nextflow ha un sistema di parametri del flusso di lavoro integrato chiamato [`params`](https://nextflow.io/docs/latest/config.html#params) che rende facile dichiarare e usare parametri CLI.

La sintassi generale è dichiarare `params.<nome_parametro>` per dire a Nextflow di aspettarsi un parametro `--<nome_parametro>` sulla riga di comando.

Qui, vogliamo creare un parametro chiamato `--input`, quindi dobbiamo dichiarare `params.input` da qualche parte nel flusso di lavoro.
In linea di principio possiamo scriverlo ovunque; ma poiché vogliamo darlo alla chiamata del processo `sayHello()`, possiamo inserirlo lì direttamente scrivendo `sayHello(params.input)`.

Nel blocco workflow, fate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emit a greeting
    sayHello(params.input)
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emit a greeting
    sayHello()
    ```

Questo dice a Nextflow di eseguire il processo `sayHello` sul valore fornito tramite il parametro `--input`.

In effetti, abbiamo compiuto i passi (2) e (3) delineati all'inizio della sezione in un colpo solo.

### 3.3. Eseguire il comando del flusso di lavoro

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

Assicuratevi di aprire il file di output per verificare che ora avete la nuova versione del saluto.

??? abstract "Contenuto del file"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà!

Notate come la nuova esecuzione ha sovrascritto il file di output pubblicato nella directory `results`.
Tuttavia, i risultati delle esecuzioni precedenti sono ancora preservati nelle directory delle attività sotto `work`.

!!! tip

    Potete facilmente distinguere i parametri a livello di Nextflow dai parametri a livello di pipeline.

    - I parametri che si applicano a una pipeline prendono sempre un doppio trattino (`--`).
    - I parametri che modificano un'impostazione di Nextflow, _ad es._ la funzionalità `-resume` che abbiamo usato prima, prendono un singolo trattino (`-`).

### 3.4. Usare valori predefiniti per i parametri da riga di comando

Ok, questo era comodo, ma in molti casi ha senso fornire un valore predefinito per un dato parametro in modo da non doverlo specificare per ogni esecuzione.

#### 3.4.1. Impostare un valore predefinito per il parametro CLI

Diamo al parametro `input` un valore predefinito dichiarandolo prima della definizione del workflow.

```groovy title="hello-world.nf" linenums="20"
/*
 * Parametri della pipeline
 */
params {
    input: String = 'Holà mundo!'
}
```

Come vedete, possiamo specificare il tipo di input che il flusso di lavoro si aspetta (Nextflow 25.10.2 e successivi).
La sintassi è `nome: Tipo = valore_predefinito`.
I tipi supportati includono `String`, `Integer`, `Float`, `Boolean` e `Path`.

!!! info

    Nei flussi di lavoro più vecchi, potreste vedere quel blocco `params` intero scritto come solo `input = 'Holà mundo!'`.

Man mano che aggiungete più parametri alla vostra pipeline, dovreste aggiungerli tutti a questo blocco, che abbiate o meno bisogno di dare loro un valore predefinito.
Questo renderà facile trovare tutti i parametri configurabili a colpo d'occhio.

#### 3.4.2. Eseguire nuovamente il flusso di lavoro senza specificare il parametro

Ora che avete un valore predefinito impostato, potete eseguire nuovamente il flusso di lavoro senza dover specificare un valore nella riga di comando.

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

L'output sarà nello stesso posto di prima, ma i contenuti dovrebbero essere aggiornati con il nuovo testo.

??? abstract "Contenuto del file"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow ha usato il valore predefinito del parametro greeting per creare l'output.

#### 3.4.3. Sovrascrivere il valore predefinito

Se fornite il parametro sulla riga di comando, il valore CLI sovrascriverà il valore predefinito.

Provatelo:

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

??? abstract "Contenuto del file"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note

    In Nextflow, ci sono più posti dove potete specificare valori per i parametri.
    Se lo stesso parametro è impostato su valori diversi in più posti, Nextflow determinerà quale valore usare in base all'ordine di precedenza descritto [qui](https://www.nextflow.io/docs/latest/config.html).

    Copriremo questo in maggior dettaglio nella Parte 6 (Configurazione).

### Takeaway

Sapete come usare un semplice input variabile fornito a runtime tramite un parametro da riga di comando, così come impostare, usare e sovrascrivere valori predefiniti.

### Cosa c'è dopo?

Imparate come gestire le esecuzioni in modo più comodo.

---

## 4. Gestire le esecuzioni del flusso di lavoro

Sapere come lanciare flussi di lavoro e recuperare output è ottimo, ma scoprirete rapidamente che ci sono alcuni altri aspetti della gestione del flusso di lavoro che renderanno la vostra vita più facile, specialmente se state sviluppando i vostri flussi di lavoro.

Qui vi mostriamo come usare la funzionalità [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) per quando dovete rilanciare lo stesso flusso di lavoro, come ispezionare il log delle esecuzioni passate con [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), e come eliminare le directory work più vecchie con [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

### 4.1. Rilanciare un flusso di lavoro con `-resume`

A volte vorrete rieseguire una pipeline che avete già lanciato in precedenza senza rifare alcun passo che è già stato completato con successo.

Nextflow ha un'opzione chiamata [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) che vi permette di farlo.
Nello specifico, in questa modalità, qualsiasi processo che è già stato eseguito con esattamente lo stesso codice, impostazioni e input verrà saltato.
Questo significa che Nextflow eseguirà solo i processi che avete aggiunto o modificato dall'ultima esecuzione, o ai quali state fornendo nuove impostazioni o input.

Ci sono due vantaggi chiave nel farlo:

- Se siete nel mezzo dello sviluppo della vostra pipeline, potete iterare più rapidamente poiché dovete solo eseguire il/i processo/i su cui state lavorando attivamente per testare le vostre modifiche.
- Se state eseguendo una pipeline in produzione e qualcosa va storto, in molti casi potete correggere il problema e rilanciare la pipeline, e riprenderà l'esecuzione dal punto di fallimento, il che può farvi risparmiare molto tempo e calcolo.

Per usarlo, basta aggiungere `-resume` al vostro comando ed eseguirlo:

```bash
nextflow run hello-world.nf -resume
```

??? success "Output del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

L'output della console dovrebbe sembrare familiare, ma c'è una cosa che è un po' diversa rispetto a prima.

Cercate il bit `cached:` che è stato aggiunto nella riga di stato del processo (riga 5), che significa che Nextflow ha riconosciuto di aver già fatto questo lavoro e ha semplicemente riutilizzato il risultato dell'esecuzione precedente riuscita.

Potete anche vedere che l'hash della sottodirectory work è lo stesso dell'esecuzione precedente.
Nextflow vi sta letteralmente indicando l'esecuzione precedente e dicendo "L'ho già fatto laggiù."

!!! tip

    Quando rieseguite una pipeline con `resume`, Nextflow non sovrascrive alcun file pubblicato fuori dalla directory work da qualsiasi esecuzione che è stata eseguita con successo in precedenza.

### 4.2. Ispezionare il log delle esecuzioni passate

Che stiate sviluppando una nuova pipeline o eseguendo pipeline in produzione, a un certo punto probabilmente dovrete cercare informazioni sulle esecuzioni passate.
Ecco come farlo.

Ogni volta che lanciate un flusso di lavoro nextflow, viene scritta una riga in un file di log chiamato `history`, sotto una directory nascosta chiamata `.nextflow` nella directory di lavoro corrente.

??? abstract "Contenuto del file"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Questo file vi dà il timestamp, il nome dell'esecuzione, lo stato, l'ID di revisione, l'ID di sessione e la riga di comando completa per ogni esecuzione Nextflow che è stata lanciata dall'interno della directory di lavoro corrente.

Un modo più comodo per accedere a queste informazioni è usare il comando `nextflow log`.

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

Questo produrrà i contenuti del file di log nel terminale, aumentati con una riga di intestazione.

Noterete che l'ID di sessione cambia ogni volta che eseguite un nuovo comando `nextflow run`, TRANNE se state usando l'opzione `-resume`.
In quel caso, l'ID di sessione rimane lo stesso.

Nextflow usa l'ID di sessione per raggruppare le informazioni di caching dell'esecuzione sotto la directory `cache`, anch'essa situata sotto `.nextflow`.

### 4.3. Eliminare le directory work più vecchie

Durante il processo di sviluppo, eseguirete tipicamente la vostra bozza di pipeline un gran numero di volte, il che può portare a un accumulo di molti file attraverso molte sottodirectory.

Fortunatamente Nextflow include un utile sottocomando `clean` che può eliminare automaticamente le sottodirectory work per esecuzioni passate che non vi interessano più.

#### 4.3.1. Determinare i criteri di eliminazione

Ci sono multiple [opzioni](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per determinare cosa eliminare.

Qui vi mostriamo un esempio che elimina tutte le sottodirectory da esecuzioni prima di una data esecuzione, specificata usando il suo nome di esecuzione.

Cercate l'esecuzione riuscita più recente dove non avete usato `-resume`; nel nostro caso il nome dell'esecuzione era `golden_cantor`.

Il nome dell'esecuzione è la stringa in due parti generata dalla macchina mostrata tra parentesi quadre nella riga di output della console `Launching (...)`.
Potete anche usare il log di Nextflow per cercare un'esecuzione in base al suo timestamp e/o riga di comando.

#### 4.3.2. Fare un'esecuzione di prova

Prima usiamo il flag di esecuzione di prova `-n` per controllare cosa verrà eliminato dato il comando:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Output del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Il vostro output avrà nomi di directory delle attività diversi e potrebbe avere un numero diverso di righe, ma dovrebbe sembrare simile all'esempio.

Se non vedete alcuna riga in output, o non avete fornito un nome di esecuzione valido o non ci sono esecuzioni passate da eliminare. Assicuratevi di cambiare `golden_cantor` nel comando di esempio con qualunque sia il nome dell'esecuzione più recente corrispondente nel vostro log.

#### 4.3.3. Procedere con l'eliminazione

Se l'output sembra come previsto e volete procedere con l'eliminazione, rieseguite il comando con il flag `-f` invece di `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Output del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

L'output dovrebbe essere simile a prima, ma ora dicendo 'Removed' invece di 'Would remove'.
Notate che questo non rimuove le sottodirectory di due caratteri (come `a3/` sopra) ma ne svuota i contenuti.

!!! Warning

    Eliminare le sottodirectory work da esecuzioni passate le rimuove dalla cache di Nextflow ed elimina qualsiasi output che era memorizzato in quelle directory.
    Questo significa che rompe la capacità di Nextflow di riprendere l'esecuzione senza rieseguire i processi corrispondenti.

    Siete responsabili di salvare qualsiasi output che vi interessa o su cui pianificate di fare affidamento! Questa è la ragione principale per cui preferiamo usare la modalità `copy` piuttosto che la modalità `symlink` per la direttiva `publish`.

### Takeaway

Sapete come pubblicare output in una directory specifica, rilanciare una pipeline senza ripetere passi che sono già stati eseguiti in modo identico, e usare il comando `nextflow clean` per pulire le vecchie directory work.

Più in generale, sapete come interpretare un semplice flusso di lavoro Nextflow, gestire la sua esecuzione e recuperare gli output.

### Cosa c'è dopo?

Prendetevi una piccola pausa, ve la siete meritata!

Quando siete pronti, passate alla [**Parte 2: Hello Channels**](./02_hello_channels.md) per imparare come usare i canali per alimentare input nel vostro flusso di lavoro, il che vi permetterà di sfruttare il parallelismo di flusso dati integrato di Nextflow e altre funzionalità potenti.

---

## Quiz

<quiz>
Quali sono i componenti minimi richiesti di un processo Nextflow?
- [ ] Solo blocchi input e output
- [x] Blocchi output e script
- [ ] Blocchi input, output e script
- [ ] Solo un blocco script

Per saperne di più: [1.1.1. La definizione del processo](#111-the-process-definition)
</quiz>

<quiz>
Qual è lo scopo del blocco output in un processo?
- [ ] Stampare risultati sulla console
- [ ] Salvare file nella directory work
- [x] Dichiarare gli output attesi dal processo
- [ ] Definire variabili d'ambiente

Per saperne di più: [1.1.1. La definizione del processo](#111-the-process-definition)
</quiz>

<quiz>
Quale comando viene usato per eseguire un flusso di lavoro Nextflow?
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

Per saperne di più: [1.2.2. Trovare l'output e i registri nella directory `work`](#122-find-the-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Cosa fa il flag `-resume`?
- [ ] Riavvia il flusso di lavoro dall'inizio
- [ ] Mette in pausa il flusso di lavoro
- [x] Salta i processi che sono già stati completati con successo
- [ ] Crea un backup del flusso di lavoro

Per saperne di più: [4.1. Rilanciare un flusso di lavoro con `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Qual è la modalità predefinita per pubblicare gli output del flusso di lavoro?
- [ ] Copiare i file nella directory di output
- [x] Creare link simbolici nella directory di output
- [ ] Spostare i file nella directory di output
- [ ] Comprimere i file nella directory di output

Per saperne di più: [2.3. Impostare la modalità di pubblicazione su copy](#23-set-the-publish-mode-to-copy)
</quiz>

<quiz>
Come si passa un valore di parametro a un flusso di lavoro Nextflow dalla riga di comando?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Per saperne di più: [3.2. Impostare un parametro da riga di comando per catturare l'input dell'utente](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
Come si fa riferimento a una variabile all'interno di un blocco script Nextflow?
- [ ] Usare la sintassi `%variable%`
- [x] Usare la sintassi `#!groovy ${variable}`
- [ ] Usare la sintassi `{{variable}}`
- [ ] Usare la sintassi `[variable]`
</quiz>
