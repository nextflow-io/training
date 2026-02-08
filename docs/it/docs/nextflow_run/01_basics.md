# Parte 1: Operazioni base

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa prima parte del corso di formazione Nextflow Run, ci introduciamo all'argomento con un esempio Hello World molto basilare e indipendente dal dominio, che useremo per dimostrare le operazioni essenziali e indicare i corrispondenti componenti del codice Nextflow.

??? info "Cos'è un esempio Hello World?"

    Un "Hello World!" è un esempio minimalista pensato per dimostrare la sintassi base e la struttura di un linguaggio di programmazione o framework software.
    L'esempio tipicamente consiste nello stampare la frase "Hello, World!" sul dispositivo di output, come la console o il terminale, o nello scriverla su un file.

---

## 1. Esegui un Hello World direttamente

Dimostriamo questo concetto con un semplice comando che eseguiamo direttamente nel terminale, per mostrare cosa fa prima di incapsularlo in Nextflow.

!!! tip "Suggerimento"

    Ricorda che dovresti ora trovarti all'interno della directory `nextflow-run/` come descritto nella pagina [Per iniziare](00_orientation.md).

### 1.1. Fai dire hello al terminale

Esegui il seguente comando nel tuo terminale.

```bash
echo 'Hello World!'
```

??? success "Output del comando"

    ```console
    Hello World!
    ```

Questo stampa il testo 'Hello World' proprio lì nel terminale.

### 1.2. Scrivi l'output su un file

Eseguire pipeline comporta principalmente la lettura di dati da file e la scrittura di risultati su altri file, quindi modifichiamo il comando per scrivere l'output di testo su un file per rendere l'esempio un po' più rilevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Output del comando"

    ```console

    ```

Questo non stampa nulla nel terminale.

### 1.3. Trova l'output

Il testo 'Hello World' dovrebbe ora essere nel file di output che abbiamo specificato, chiamato `output.txt`.
Puoi aprirlo nell'esploratore file o dalla riga di comando usando l'utility `cat`, per esempio.

??? abstract "Contenuto del file"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Questo è ciò che cercheremo di replicare con il nostro primo workflow Nextflow.

### Riepilogo

Ora sai come eseguire un semplice comando nel terminale che produce del testo, e opzionalmente, come fargli scrivere l'output su un file.

### Cosa c'è dopo?

Scopri cosa serve per eseguire un workflow Nextflow che ottiene lo stesso risultato.

---

## 2. Esegui il workflow

Ti forniamo uno script di workflow chiamato `1-hello.nf` che prende un saluto in input tramite un argomento da riga di comando chiamato `--input` e produce un file di testo contenente quel saluto.

Non guarderemo ancora il codice; prima vediamo come appare eseguirlo.

### 2.1. Avvia il workflow e monitora l'esecuzione

Nel terminale, esegui il seguente comando:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Output del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Se l'output della tua console appare più o meno così, congratulazioni, hai appena eseguito il tuo primo workflow Nextflow!

L'output più importante qui è l'ultima riga, che è evidenziata nell'output sopra:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Questo ci dice che il process `sayHello` è stato eseguito con successo una volta (`1 of 1 ✔`).

Ottimo, ma potresti chiederti: dov'è l'output?

### 2.2. Trova il file di output nella directory `results`

Questo workflow è configurato per pubblicare il suo output in una directory dei risultati.
Se guardi la tua directory corrente, vedrai che quando hai eseguito il workflow, Nextflow ha creato una nuova directory chiamata `results`, oltre a una sottodirectory chiamata `1-hello` sotto di essa, contenente un file chiamato `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Apri il file; il contenuto dovrebbe corrispondere alla stringa che hai specificato nella riga di comando.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Ottimo, il nostro workflow ha fatto ciò che doveva fare!

### 2.3. Salva i risultati in una directory diversa

Per impostazione predefinita, Nextflow salva gli output della pipeline in una directory chiamata `results` nel tuo percorso corrente.
Per cambiare dove i tuoi file vengono pubblicati, usa il flag CLI `-output-dir` (o `-o` in breve)

!!! danger "Attenzione"

    Nota che `--input` ha due trattini e `-output-dir` ne ha uno!
    Questo perché `--input` è un _parametro_ della pipeline e `-output-dir` è un flag CLI principale di Nextflow.
    Approfondiremo questo più avanti.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Dovresti vedere che i tuoi output sono ora pubblicati in una directory chiamata `hello_results` invece di `results`:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

I file all'interno di questa directory sono esattamente gli stessi di prima, è solo la directory di primo livello ad essere diversa.
Tuttavia, tieni presente in entrambi i casi che il risultato 'pubblicato' è una copia (o in alcuni casi un link simbolico) dell'output effettivo prodotto da Nextflow quando ha eseguito il workflow.

Quindi ora daremo un'occhiata sotto il cofano per vedere dove Nextflow ha effettivamente eseguito il lavoro.

!!! Warning "Avviso"

    Non tutti i workflow saranno configurati per pubblicare output in una directory dei risultati, e/o i nomi delle directory e la struttura potrebbero essere diversi.
    Un po' più avanti in questa sezione, ti mostreremo come scoprire dove questo comportamento viene specificato.

### 2.4. Trova l'output originale e i log nella directory `work/`

Quando esegui un workflow, Nextflow crea una distinta 'directory di attività' per ogni singola invocazione di ogni process nel workflow (=ogni step nella pipeline).
Per ciascuna, metterà in staging gli input necessari, eseguirà le istruzioni rilevanti e scriverà output e file di log all'interno di quella singola directory, che viene nominata automaticamente usando un hash per renderla unica.

Tutte queste directory di attività risiederanno sotto una directory chiamata `work` nella tua directory corrente (dove stai eseguendo il comando).

Potrebbe sembrare confuso, quindi vediamo come appare in pratica.

Tornando all'output della console per il workflow che abbiamo eseguito prima, avevamo questa riga:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Vedi come la riga inizia con `[a3/1e1535]`?
Questa è una forma troncata del percorso della directory di attività per quella chiamata di process, e ti dice dove trovare l'output della chiamata al process `sayHello` all'interno del percorso della directory `work/`.

Puoi trovare il percorso completo digitando il seguente comando (sostituendo `a3/1e1535` con quello che vedi nel tuo terminale) e premendo il tasto tab per l'autocompletamento del percorso o aggiungendo un asterisco:

```bash
ls work/a3/1e1535*
```

Questo dovrebbe produrre il percorso completo della directory: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Diamo un'occhiata a cosa c'è dentro.

??? abstract "Contenuti della directory"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Non vedi la stessa cosa?"

    I nomi esatti delle sottodirectory saranno diversi sul tuo sistema.

    Se esplori i contenuti della sottodirectory di attività nell'esploratore file di VSCode, vedrai tutti i file subito.
    Tuttavia, i file di log sono impostati per essere invisibili nel terminale, quindi se vuoi usare `ls` o `tree` per visualizzarli, dovrai impostare l'opzione rilevante per mostrare i file invisibili.

    ```bash
    tree -a work
    ```

Ci sono due set di directory in `work/`, dalle due diverse esecuzioni della pipeline che abbiamo fatto.
Ogni esecuzione di attività ottiene la propria directory isolata in cui lavorare.
In questo caso la pipeline ha fatto la stessa cosa entrambe le volte, quindi i contenuti di ciascuna directory di attività sono identici.

Dovresti riconoscere immediatamente il file `output.txt`, che è infatti l'output originale del process `sayHello` che è stato pubblicato nella directory `results`.
Se lo apri, troverai di nuovo il saluto `Hello World!`.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

E tutti quegli altri file?

Questi sono i file helper e di log che Nextflow ha scritto come parte dell'esecuzione dell'attività:

- **`.command.begin`**: File sentinella creato non appena l'attività viene lanciata.
- **`.command.err`**: Messaggi di errore (`stderr`) emessi dalla chiamata del process
- **`.command.log`**: Output di log completo emesso dalla chiamata del process
- **`.command.out`**: Output regolare (`stdout`) della chiamata del process
- **`.command.run`**: Script completo eseguito da Nextflow per eseguire la chiamata del process
- **`.command.sh`**: Il comando che è stato effettivamente eseguito dalla chiamata del process
- **`.exitcode`**: Il codice di uscita risultante dal comando

Il file `.command.sh` è particolarmente utile perché ti mostra il comando principale che Nextflow ha eseguito, non includendo tutta la gestione e il setup dell'attività/ambiente.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Quindi questo conferma che il workflow ha composto lo stesso comando che abbiamo eseguito direttamente sulla riga di comando prima.

Quando qualcosa va storto e devi risolvere il problema, può essere utile guardare lo script `command.sh` per controllare esattamente quale comando Nextflow ha composto in base alle istruzioni del workflow, all'interpolazione delle variabili e così via.

### 2.5. Riesegui il workflow con saluti diversi

Prova a rieseguire il workflow alcune volte con valori diversi per l'argomento `--input`, poi guarda le directory di attività.

??? abstract "Contenuti della directory"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Vedi che una nuova sottodirectory con un set completo di file di output e log è stata creata per ogni esecuzione.

Al contrario, se guardi la directory `results`, c'è ancora solo un set di risultati, e il contenuto del file di output corrisponde a qualunque cosa tu abbia eseguito per ultima.

??? abstract "Contenuti della directory"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Questo ti mostra che i risultati pubblicati verranno sovrascritti dalle esecuzioni successive, mentre le directory di attività sotto `work/` vengono preservate.

### Riepilogo

Sai come eseguire un semplice script Nextflow, monitorare la sua esecuzione e trovare i suoi output.

### Cosa c'è dopo?

Impara a leggere uno script Nextflow base e identificare come i suoi componenti si relazionano alla sua funzionalità.

---

## 3. Esamina lo script starter del workflow Hello World

Ciò che abbiamo fatto finora era trattare lo script del workflow come una scatola nera.
Ora che abbiamo visto cosa fa, apriamo la scatola e guardiamo dentro.

Il nostro obiettivo qui non è memorizzare la sintassi del codice Nextflow, ma formare un'intuizione di base su quali sono i componenti principali e come sono organizzati.

### 3.1. Esamina la struttura generale del codice

Troverai lo script `1-hello.nf` nella tua directory corrente, che dovrebbe essere `nextflow-run`. Aprilo nel pannello dell'editor.

??? full-code "File di codice completo"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
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

Uno script di workflow Nextflow tipicamente include una o più definizioni di **process**, il **workflow** stesso, e alcuni blocchi opzionali come **params** e **output**.

Ogni **process** descrive quale/i operazione/i lo step corrispondente nella pipeline dovrebbe compiere, mentre il **workflow** descrive la logica del dataflow che connette i vari step.

Diamo prima uno sguardo più da vicino al blocco **process**, poi guarderemo il blocco **workflow**.

### 3.2. La definizione del `process`

Il primo blocco di codice descrive un [**process**](https://nextflow.io/docs/latest/process.html).
La definizione del process inizia con la parola chiave `process`, seguita dal nome del process e infine dal corpo del process delimitato da parentesi graffe.
Il corpo del process deve contenere un blocco script che specifica il comando da eseguire, che può essere qualsiasi cosa tu possa eseguire in un terminale da riga di comando.

```groovy title="1-hello.nf" linenums="3"
/*
* Use echo to print a greeting to a file
*/
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

Qui abbiamo un **process** chiamato `sayHello` che prende una variabile di **input** chiamata `greeting` e scrive il suo **output** in un file chiamato `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Questa è una definizione di process molto minimale che contiene solo una definizione `input`, una definizione `output` e lo `script` da eseguire.

La definizione `input` include il qualificatore `val`, che dice a Nextflow di aspettarsi un valore di qualche tipo (può essere una stringa, un numero, qualsiasi cosa).

La definizione `output` include il qualificatore `path`, che dice a Nextflow che questo dovrebbe essere gestito come un percorso (include sia percorsi di directory che file).

### 3.3. La definizione del `workflow`

Il secondo blocco di codice descrive il [**workflow**](https://nextflow.io/docs/latest/workflow.html) stesso.
La definizione del workflow inizia con la parola chiave `workflow`, seguita da un nome opzionale, poi il corpo del workflow delimitato da parentesi graffe.

Qui abbiamo un **workflow** che consiste in un blocco `main:` e un blocco `publish:`.
Il blocco `main:` è il corpo principale del workflow e il blocco `publish:` elenca gli output che dovrebbero essere pubblicati nella directory `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emette un saluto
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

In questo caso il blocco `main:` contiene una chiamata al process `sayHello` e gli fornisce un input chiamato `params.input` da usare come saluto.

Come discuteremo più in dettaglio tra un momento, `params.input` contiene il valore che abbiamo dato al parametro `--input` nella nostra riga di comando.

Il blocco `publish:` elenca l'output della chiamata al process `sayHello()`, a cui si riferisce come `sayHello.out` e dà il nome `first_output` (questo può essere qualsiasi cosa l'autore del workflow voglia).

Questa è una definizione di **workflow** molto minimale.
In una pipeline del mondo reale, il workflow tipicamente contiene chiamate multiple a **processes** connessi da **channels**, e potrebbero esserci valori predefiniti impostati per gli input variabili.

Approfondiremo questo nella Parte 2 del corso.
Per ora, diamo uno sguardo più da vicino a come il nostro workflow gestisce input e output.

### 3.4. Il sistema `params` per i parametri da riga di comando

Il `params.input` che forniamo alla chiamata del process `sayHello()` è un utile pezzo di codice Nextflow e vale la pena spendere un minuto in più.

Come menzionato sopra, è così che passiamo il valore del parametro da riga di comando `--input` alla chiamata del process `sayHello()`.
In realtà, semplicemente dichiarare `params.someParameterName` è sufficiente per dare al workflow un parametro chiamato `--someParameterName` dalla riga di comando.

Qui abbiamo formalizzato quella dichiarazione di parametro impostando un blocco `params` che specifica il tipo di input che il workflow si aspetta (Nextflow 25.10.2 e successivi).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

I tipi supportati includono `String`, `Integer`, `Float`, `Boolean` e `Path`.
Per saperne di più, consulta [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) nella documentazione di riferimento di Nextflow.

!!! tip "Suggerimento"

    Ricorda che i parametri del _workflow_ dichiarati usando il sistema `params` prendono sempre due trattini nella riga di comando (`--`).
    Questo li distingue dai flag CLI _a livello di Nextflow_, che prendono solo un trattino (`-`).

### 3.5. La direttiva `publish`

All'altra estremità del workflow, abbiamo già dato un'occhiata al blocco `publish:`.
Questa è una metà del sistema di gestione dell'output; l'altra metà è il blocco `output` situato sotto.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Questo specifica che l'output `first_output` elencato nel blocco `publish:` dovrebbe essere copiato in una sottodirectory chiamata `1-hello` sotto la directory di output predefinita `results`.

La riga `mode 'copy'` sovrascrive il comportamento predefinito del sistema, che è di creare un link simbolico (o symlink) al file originale nella directory `work/` invece di una copia vera e propria.

Ci sono più opzioni di quelle mostrate qui per controllare il comportamento di pubblicazione; ne copriremo alcune più avanti.
Vedrai anche che quando un workflow genera output multipli, ognuno viene elencato in questo modo nel blocco `output`.

Per saperne di più, consulta [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) nella documentazione di riferimento di Nextflow.

??? info "Sintassi precedente per pubblicare output usando `publishDir`"

    Fino a molto recentemente, il modo stabilito per pubblicare gli output era farlo a livello di ogni singolo process usando una direttiva `publishDir`.

    Troverai ancora questo pattern di codice ovunque nelle pipeline Nextflow più vecchie e nei moduli di process, quindi è importante esserne consapevoli.

    Invece di avere un blocco `publish:` nel workflow e un blocco `output` al livello superiore, vedresti una riga `publishDir` nella definizione del process `sayHello`:

    ```groovy title="Esempio di sintassi" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Tuttavia, non raccomandiamo di usarlo in nessun nuovo lavoro poiché verrà eventualmente vietato nelle future versioni del linguaggio Nextflow.

### Riepilogo

Ora sai come è strutturato un semplice workflow Nextflow, e come i componenti base si relazionano alla sua funzionalità.

### Cosa c'è dopo?

Impara a gestire le esecuzioni del tuo workflow in modo conveniente.

---

## 4. Gestisci le esecuzioni del workflow

Sapere come lanciare workflow e recuperare output è ottimo, ma scoprirai presto che ci sono alcuni altri aspetti della gestione dei workflow che ti renderanno la vita più facile.

Qui ti mostriamo come sfruttare la funzionalità `resume` per quando devi rilanciare lo stesso workflow, come ispezionare i log di esecuzione con `nextflow log`, e come eliminare le vecchie directory di lavoro con `nextflow clean`.

### 4.1. Rilancia un workflow con `-resume`

A volte, vorrai rieseguire una pipeline che hai già lanciato in precedenza senza rifare alcun lavoro che è già stato completato con successo.

Nextflow ha un'opzione chiamata `-resume` che ti permette di farlo.
Specificamente, in questa modalità, tutti i process che sono già stati eseguiti con lo stesso identico codice, impostazioni e input verranno saltati.
Questo significa che Nextflow eseguirà solo i process che hai aggiunto o modificato dall'ultima esecuzione, o a cui stai fornendo nuove impostazioni o input.

Ci sono due vantaggi chiave nel farlo:

- Se sei nel mezzo dello sviluppo di una pipeline, puoi iterare più rapidamente poiché devi eseguire solo il/i process su cui stai attivamente lavorando per testare le tue modifiche.
- Se stai eseguendo una pipeline in produzione e qualcosa va storto, in molti casi puoi risolvere il problema e rilanciare la pipeline, e riprenderà l'esecuzione dal punto di fallimento, il che può farti risparmiare molto tempo e calcolo.

Per usarlo, aggiungi semplicemente `-resume` al tuo comando ed eseguilo:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Output del comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

L'output della console dovrebbe sembrare familiare, ma c'è una cosa che è un po' diversa rispetto a prima.

Cerca la parte `cached:` che è stata aggiunta nella riga di stato del process (riga 5), il che significa che Nextflow ha riconosciuto che ha già fatto questo lavoro e ha semplicemente riutilizzato il risultato dell'esecuzione precedente riuscita.

Puoi anche vedere che l'hash della sottodirectory di lavoro è lo stesso dell'esecuzione precedente.
Nextflow ti sta letteralmente indicando l'esecuzione precedente e dicendo "L'ho già fatto laggiù."

!!! tip "Suggerimento"

    Quando riesegui una pipeline con `resume`, Nextflow non sovrascrive alcun file pubblicato al di fuori della directory di lavoro da esecuzioni che sono state eseguite con successo in precedenza.

    Per saperne di più, consulta [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) nella documentazione di riferimento di Nextflow.

### 4.2. Ispeziona il log delle esecuzioni passate

Ogni volta che lanci un workflow nextflow, una riga viene scritta in un file di log chiamato `history`, sotto una directory nascosta chiamata `.nextflow` nella directory di lavoro corrente.

??? abstract "Contenuto del file"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Questo file ti dà il timestamp, il nome dell'esecuzione, lo stato, l'ID di revisione, l'ID di sessione e la riga di comando completa per ogni esecuzione Nextflow che è stata lanciata dalla directory di lavoro corrente.

Un modo più conveniente per accedere a queste informazioni è usare il comando [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

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

Questo stamperà i contenuti del file di log nel terminale, arricchiti con una riga di intestazione.

Noterai che l'ID di sessione cambia ogni volta che esegui un nuovo comando `nextflow run`, TRANNE se stai usando l'opzione `-resume`.
In quel caso, l'ID di sessione rimane lo stesso.

Nextflow usa l'ID di sessione per raggruppare le informazioni di caching dell'esecuzione sotto la directory `cache`, anch'essa situata sotto `.nextflow`.

### 4.3. Elimina le vecchie directory di lavoro

Se esegui molte pipeline, potresti finire per accumulare moltissimi file in molte sottodirectory.
Poiché le sottodirectory sono nominate casualmente, è difficile capire dai loro nomi quali sono esecuzioni più vecchie vs. più recenti.

Fortunatamente Nextflow include un utile comando chiamato [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) che può automaticamente eliminare le sottodirectory di lavoro per le esecuzioni passate che non ti interessano più.

#### 4.3.1. Determina i criteri di eliminazione

Ci sono molteplici opzioni per determinare cosa eliminare, che puoi esplorare nella documentazione linkata sopra.
Qui ti mostriamo un esempio che elimina tutte le sottodirectory dalle esecuzioni prima di una data esecuzione, specificata usando il suo nome di esecuzione.

Cerca l'esecuzione riuscita più recente dove non hai usato `-resume`; nel nostro caso il nome dell'esecuzione era `backstabbing_swartz`.

Il nome dell'esecuzione è la stringa in due parti generata automaticamente mostrata tra parentesi quadre nell'output della console nella riga `Launching (...)`.
Puoi anche usare il log di Nextflow per cercare un'esecuzione in base al suo timestamp e/o riga di comando.

#### 4.3.2. Fai un'esecuzione di prova

Prima usiamo il flag dry run `-n` per controllare cosa verrà eliminato dato il comando:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Output del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Il tuo output avrà nomi di directory di attività diversi e potrebbe avere un numero diverso di righe, ma dovrebbe apparire simile all'esempio.

Se non vedi righe nell'output, o non hai fornito un nome di esecuzione valido o non ci sono esecuzioni passate da eliminare. Assicurati di cambiare `backstabbing_swartz` nel comando di esempio con qualsiasi sia il nome dell'esecuzione più recente corrispondente nel tuo log.

#### 4.3.3. Procedi con l'eliminazione

Se l'output sembra come previsto e vuoi procedere con l'eliminazione, riesegui il comando con il flag `-f` invece di `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Output del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

L'output dovrebbe essere simile a prima, ma ora dicendo 'Removed' invece di 'Would remove'.
Nota che questo non rimuove le sottodirectory di due caratteri (come `eb/` sopra) ma ne svuota i contenuti.

!!! Warning "Avviso"

    Eliminare le sottodirectory di lavoro dalle esecuzioni passate le rimuove dalla cache di Nextflow ed elimina tutti gli output che erano memorizzati in quelle directory.
    Questo significa che rompe la capacità di Nextflow di riprendere l'esecuzione senza rieseguire i process corrispondenti.

    Sei responsabile di salvare tutti gli output a cui tieni! Questo è il motivo principale per cui preferiamo usare la modalità `copy` piuttosto che la modalità `symlink` per la direttiva `publish`.

### Riepilogo

Sai come rilanciare una pipeline senza ripetere step che sono già stati eseguiti in modo identico, ispezionare il log di esecuzione, e usare il comando `nextflow clean` per pulire le vecchie directory di lavoro.

### Cosa c'è dopo?

Prenditi una piccola pausa! Hai appena assorbito i blocchi costitutivi della sintassi Nextflow e le istruzioni di utilizzo base.

Nella prossima sezione di questa formazione, guarderemo quattro versioni successivamente più realistiche della pipeline Hello World che dimostreranno come Nextflow ti permette di elaborare input multipli efficientemente, eseguire workflow composti da step multipli connessi insieme, sfruttare componenti di codice modulari, e utilizzare container per una maggiore riproducibilità e portabilità.

---

## Quiz

<quiz>
Nella riga di output della console `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, cosa rappresenta `[a3/7be2fa]`?
- [ ] Il numero di versione del process
- [ ] Un identificatore univoco dell'esecuzione
- [x] Il percorso troncato alla directory di lavoro dell'attività
- [ ] Il checksum del file di output

Approfondisci: [2.4. Trova l'output originale e i log nella directory `work/`](#24-trova-loutput-originale-e-i-log-nella-directory-work)
</quiz>

<quiz>
Qual è lo scopo del file `.command.sh` in una directory di attività?
- [ ] Memorizza le impostazioni di configurazione dell'attività
- [x] Mostra il comando effettivo che è stato eseguito dal process
- [ ] Contiene messaggi di errore dalle attività fallite
- [ ] Elenca i file di input messi in staging per l'attività

Approfondisci: [2.4. Trova l'output originale e i log nella directory `work/`](#24-trova-loutput-originale-e-i-log-nella-directory-work)
</quiz>

<quiz>
Cosa succede ai risultati pubblicati quando riesegui un workflow senza `-resume`?
- [ ] Vengono preservati in directory separate con timestamp
- [x] Vengono sovrascritti dalla nuova esecuzione
- [ ] Nextflow previene la sovrascrittura e fallisce
- [ ] Vengono automaticamente backed up

Approfondisci: [2.5. Riesegui il workflow con saluti diversi](#25-riesegui-il-workflow-con-saluti-diversi)
</quiz>

<quiz>
Cosa indica questo output della console?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] L'attività è fallita ed è stata saltata
- [ ] L'attività è in attesa in una coda
- [x] Nextflow ha riutilizzato i risultati di una precedente esecuzione identica
- [ ] L'attività è stata cancellata manualmente

Approfondisci: [4.1. Rilancia un workflow con `-resume`](#41-rilancia-un-workflow-con--resume)
</quiz>

<quiz>
Dove memorizza Nextflow la cronologia delle esecuzioni che il comando `nextflow log` visualizza?
- [ ] Nella directory results
- [ ] Nella directory work
- [x] Nel file `.nextflow/history`
- [ ] In `nextflow.config`

Approfondisci: [4.2. Ispeziona il log delle esecuzioni passate](#42-ispeziona-il-log-delle-esecuzioni-passate)
</quiz>

<quiz>
Qual è lo scopo del blocco `params` in un file di workflow?
- [ ] Definire i requisiti di risorse dei process
- [ ] Configurare l'executor
- [x] Dichiarare e tipizzare i parametri di input del workflow
- [ ] Specificare le opzioni di pubblicazione dell'output

Approfondisci: [3.4. Il sistema params per i parametri da riga di comando](#34-il-sistema-params-per-i-parametri-da-riga-di-comando)
</quiz>

<quiz>
Nel blocco `output` del workflow, cosa fa `mode 'copy'`?
- [ ] Crea un backup della directory di lavoro
- [x] Fa una copia completa dei file invece di link simbolici
- [ ] Copia lo script del workflow nei risultati
- [ ] Abilita la copia incrementale dei file

Approfondisci: [3.5. La direttiva publish](#35-la-direttiva-publish)
</quiz>

<quiz>
Qual è il flag raccomandato da usare con il comando `nextflow clean` prima di eliminare effettivamente i file?
- [x] `-n` (dry run) per vedere in anteprima cosa verrebbe eliminato
- [ ] `-v` (verbose) per vedere output dettagliato
- [ ] `-a` (all) per selezionare tutte le directory
- [ ] `-q` (quiet) per sopprimere gli avvisi

Approfondisci: [4.3. Elimina le vecchie directory di lavoro](#43-elimina-le-vecchie-directory-di-lavoro)
</quiz>
