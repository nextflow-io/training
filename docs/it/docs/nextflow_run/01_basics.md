# Parte 1: Eseguire operazioni di base

In questa prima parte del corso di formazione Nextflow Run, ci addentriamo nell'argomento con un esempio Hello World molto semplice e indipendente dal dominio, che useremo per dimostrare le operazioni essenziali e indicare i componenti corrispondenti del codice Nextflow.

??? info "Cos'è un esempio Hello World?"

    Un "Hello World!" è un esempio minimalista pensato per dimostrare la sintassi di base e la struttura di un linguaggio di programmazione o framework software.
    L'esempio consiste tipicamente nella stampa della frase "Hello, World!" sul dispositivo di output, come la console o il terminale, oppure nella sua scrittura in un file.

---

## 1. Eseguire un Hello World direttamente

Dimostriamo questo concetto con un semplice comando che eseguiamo direttamente nel terminale, per mostrare cosa fa prima di inserirlo in Nextflow.

!!! tip

    Ricordate che ora dovreste trovarvi all'interno della directory `nextflow-run/` come descritto nella pagina [Primi Passi](00_orientation.md).

### 1.1. Far dire hello al terminale

Eseguite il seguente comando nel vostro terminale.

```bash
echo 'Hello World!'
```

??? success "Output del comando"

    ```console
    Hello World!
    ```

Questo produce il testo 'Hello World' direttamente nel terminale.

### 1.2. Scrivere l'output in un file

L'esecuzione di pipeline consiste principalmente nella lettura di dati da file e nella scrittura di risultati in altri file, quindi modifichiamo il comando per scrivere l'output di testo in un file per rendere l'esempio un po' più rilevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Output del comando"

    ```console

    ```

Questo non produce alcun output nel terminale.

### 1.3. Trovare l'output

Il testo 'Hello World' dovrebbe ora trovarsi nel file di output che abbiamo specificato, chiamato `output.txt`.
Potete aprirlo nell'esploratore di file o dalla riga di comando usando l'utility `cat`, per esempio.

??? abstract "Contenuto del file"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Questo è ciò che cercheremo di replicare con il nostro primissimo flusso di lavoro Nextflow.

### Takeaway

Ora sapete come eseguire un semplice comando nel terminale che produce del testo, e opzionalmente, come farlo scrivere l'output in un file.

### Cosa c'è dopo?

Scoprite cosa serve per eseguire un flusso di lavoro Nextflow che ottiene lo stesso risultato.

---

## 2. Eseguire il flusso di lavoro

Vi forniamo uno script di flusso di lavoro chiamato `1-hello.nf` che prende un saluto di input tramite un argomento da riga di comando chiamato `--input` e produce un file di testo contenente quel saluto.

Non esamineremo ancora il codice; prima vediamo come si presenta la sua esecuzione.

### 2.1. Lanciare il flusso di lavoro e monitorare l'esecuzione

Nel terminale, eseguite il seguente comando.

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

Se l'output della vostra console assomiglia a quello, allora congratulazioni, avete appena eseguito il vostro primo flusso di lavoro Nextflow!

L'output più importante qui è l'ultima riga, che è evidenziata nell'output sopra:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Questo ci dice che il processo `sayHello` è stato eseguito con successo una volta (`1 of 1 ✔`).

È fantastico, ma potreste chiedervi: dove si trova l'output?

### 2.2. Trovare il file di output nella directory `results`

Questo flusso di lavoro è configurato per pubblicare il suo output in una directory results.
Se guardate la vostra directory corrente, vedrete che quando avete eseguito il flusso di lavoro, Nextflow ha creato una nuova directory chiamata `results`, così come una sottodirectory chiamata `1-hello` sotto di essa, contenente un file chiamato `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Aprite il file; il contenuto dovrebbe corrispondere alla stringa che avete specificato sulla riga di comando.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Fantastico, il nostro flusso di lavoro ha fatto quello che doveva fare!

### 2.3. Salvare i risultati in una directory diversa

Per impostazione predefinita, Nextflow salverà gli output della pipeline in una directory chiamata `results` nel vostro percorso corrente.
Per cambiare dove vengono pubblicati i vostri file, usate il flag CLI `-output-dir` (o `-o` in forma breve)

!!! danger

    Notate che `--input` ha due trattini e `-output-dir` ne ha uno!
    Questo perché `--input` è un _parametro_ della pipeline e `-output-dir` è un flag CLI core di Nextflow.
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

Dovreste vedere che i vostri output sono ora pubblicati in una directory chiamata `hello_results` invece di `results`:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

I file all'interno di questa directory sono esattamente gli stessi di prima, è solo la directory di primo livello che è diversa.
Tuttavia, siate consapevoli che in entrambi i casi il risultato 'pubblicato' è una copia (o in alcuni casi un link simbolico) dell'output effettivo prodotto da Nextflow quando ha eseguito il flusso di lavoro.

Quindi ora, daremo un'occhiata sotto il cofano per vedere dove Nextflow ha effettivamente eseguito il lavoro.

!!! Warning

    Non tutti i flussi di lavoro saranno configurati per pubblicare gli output in una directory results, e/o i nomi delle directory e la struttura potrebbero essere diversi.
    Un po' più avanti in questa sezione, vi mostreremo come scoprire dove è specificato questo comportamento.

### 2.4. Trovare l'output originale e i registri nella directory `work/`

Quando eseguite un flusso di lavoro, Nextflow crea una 'directory di attività' distinta per ogni singola invocazione di ciascun processo nel flusso di lavoro (=ogni passo nella pipeline).
Per ciascuna, preparerà gli input necessari, eseguirà le istruzioni rilevanti e scriverà output e file di log all'interno di quella directory, che viene nominata automaticamente usando un hash per renderla unica.

Tutte queste directory di attività vivranno sotto una directory chiamata `work` all'interno della vostra directory corrente (dove state eseguendo il comando).

Questo può sembrare confuso, quindi vediamo come appare in pratica.

Tornando all'output della console per il flusso di lavoro che abbiamo eseguito prima, avevamo questa riga:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Vedete come la riga inizia con `[a3/1e1535]`?
Questa è una forma troncata del percorso della directory di attività per quella chiamata di processo, e vi dice dove trovare l'output della chiamata del processo `sayHello` all'interno del percorso della directory `work/`.

Potete trovare il percorso completo digitando il seguente comando (sostituendo `a3/1e1535` con quello che vedete nel vostro terminale) e premendo il tasto tab per completare automaticamente il percorso o aggiungendo un asterisco:

```bash
ls work/a3/1e1535*
```

Questo dovrebbe produrre il percorso completo della directory: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Diamo un'occhiata a cosa c'è dentro.

??? abstract "Contenuto della directory"

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

??? question "Non vedete la stessa cosa?"

    I nomi esatti delle sottodirectory saranno diversi sul vostro sistema.

    Se navigate il contenuto della sottodirectory di attività nell'esploratore di file di VSCode, vedrete tutti i file immediatamente.
    Tuttavia, i file di log sono impostati per essere invisibili nel terminale, quindi se volete usare `ls` o `tree` per visualizzarli, dovrete impostare l'opzione rilevante per visualizzare i file invisibili.

    ```bash
    tree -a work
    ```

Ci sono due set di directory in `work/`, dalle due diverse esecuzioni della pipeline che abbiamo fatto.
Ogni esecuzione di attività ottiene la propria directory isolata in cui lavorare.
In questo caso la pipeline ha fatto la stessa cosa entrambe le volte, quindi il contenuto di ciascuna directory di attività è identico

Dovreste riconoscere immediatamente il file `output.txt`, che è infatti l'output originale del processo `sayHello` che è stato pubblicato nella directory `results`.
Se lo aprite, troverete di nuovo il saluto `Hello World!`.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

E tutti quegli altri file?

Questi sono i file di supporto e di log che Nextflow ha scritto come parte dell'esecuzione dell'attività:

- **`.command.begin`**: File sentinella creato non appena l'attività viene lanciata.
- **`.command.err`**: Messaggi di errore (`stderr`) emessi dalla chiamata del processo
- **`.command.log`**: Output completo del log emesso dalla chiamata del processo
- **`.command.out`**: Output regolare (`stdout`) della chiamata del processo
- **`.command.run`**: Script completo eseguito da Nextflow per eseguire la chiamata del processo
- **`.command.sh`**: Il comando che è stato effettivamente eseguito dalla chiamata del processo
- **`.exitcode`**: Il codice di uscita risultante dal comando

Il file `.command.sh` è particolarmente utile perché mostra il comando principale che Nextflow ha eseguito, senza includere tutta la contabilità e la configurazione dell'attività/ambiente.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Quindi questo conferma che il flusso di lavoro ha composto lo stesso comando che abbiamo eseguito direttamente sulla riga di comando in precedenza.

Quando qualcosa va storto e dovete risolvere i problemi su cosa è successo, può essere utile guardare lo script `command.sh` per verificare esattamente quale comando Nextflow ha composto in base alle istruzioni del flusso di lavoro, all'interpolazione delle variabili e così via.

### 2.5. Rieseguire il flusso di lavoro con saluti diversi

Provate a rieseguire il flusso di lavoro alcune volte con valori diversi per l'argomento `--input`, quindi guardate le directory di attività.

??? abstract "Contenuto della directory"

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

Vedete che è stata creata una nuova sottodirectory con un set completo di file di output e log per ogni esecuzione.

Al contrario, se guardate la directory `results`, c'è ancora solo un set di risultati, e il contenuto del file di output corrisponde a qualunque cosa abbiate eseguito per ultima.

??? abstract "Contenuto della directory"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Questo vi mostra che i risultati pubblicati verranno sovrascritti dalle esecuzioni successive, mentre le directory di attività sotto `work/` vengono preservate.

### Takeaway

Sapete come eseguire un semplice script Nextflow, monitorare la sua esecuzione e trovare i suoi output.

### Cosa c'è dopo?

Imparate come leggere uno script Nextflow di base e identificare come i suoi componenti si relazionano alla sua funzionalità.

---

## 3. Esaminare lo script iniziale del flusso di lavoro Hello World

Quello che abbiamo fatto lì era fondamentalmente trattare lo script del flusso di lavoro come una scatola nera.
Ora che abbiamo visto cosa fa, apriamo la scatola e guardiamo dentro.

Il nostro obiettivo qui non è memorizzare la sintassi del codice Nextflow, ma formare un'intuizione di base su quali sono i componenti principali e come sono organizzati.

### 3.1. Esaminare la struttura complessiva del codice

Troverete lo script `1-hello.nf` nella vostra directory corrente, che dovrebbe essere `nextflow-run`. Apritelo nel pannello dell'editor.

??? full-code "File di codice completo"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usa echo per stampare 'Hello World!' in un file
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
    * Parametri della pipeline
    */
    params {
        input: String
    }

    workflow {

        main:
        // emit a greeting
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

Uno script di flusso di lavoro Nextflow include tipicamente una o più definizioni di **process**, il **workflow** stesso, e alcuni blocchi opzionali come **params** e **output**.

Ogni **process** descrive quali operazioni il passo corrispondente nella pipeline dovrebbe compiere, mentre il **workflow** descrive la logica di flusso dati che connette i vari passi.

Diamo un'occhiata più da vicino prima al blocco **process**, poi esamineremo il blocco **workflow**.

### 3.2. La definizione del `process`

Il primo blocco di codice descrive un [**process**](https://nextflow.io/docs/latest/process.html).
La definizione del processo inizia con la parola chiave `process`, seguita dal nome del processo e infine dal corpo del processo delimitato da parentesi graffe.
Il corpo del processo deve contenere un blocco script che specifica il comando da eseguire, che può essere qualsiasi cosa sareste in grado di eseguire in un terminale a riga di comando.

```groovy title="1-hello.nf" linenums="3"
/*
* Usa echo per stampare un saluto in un file
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

Questa è una definizione di processo molto minimale che contiene solo una definizione `input`, una definizione `output` e lo `script` da eseguire.

La definizione `input` include il qualificatore `val`, che dice a Nextflow di aspettarsi un valore di qualche tipo (può essere una stringa, un numero, qualsiasi cosa).

La definizione `output` include il qualificatore `path`, che dice a Nextflow che questo dovrebbe essere gestito come un percorso (include sia percorsi di directory che file).

### 3.3. La definizione del `workflow`

Il secondo blocco di codice descrive il [**workflow**](https://nextflow.io/docs/latest/workflow.html) stesso.
La definizione del flusso di lavoro inizia con la parola chiave `workflow`, seguita da un nome opzionale, poi dal corpo del flusso di lavoro delimitato da parentesi graffe.

Qui abbiamo un **workflow** che consiste in un blocco `main:` e un blocco `publish:`.
Il blocco `main:` è il corpo principale del flusso di lavoro e il blocco `publish:` elenca gli output che dovrebbero essere pubblicati nella directory `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

In questo caso il blocco `main:` contiene una chiamata al processo `sayHello` e gli fornisce un input chiamato `params.input` da usare come saluto.

Come discuteremo più in dettaglio tra un momento, `params.input` contiene il valore che abbiamo dato al parametro `--input` nella nostra riga di comando.

Il blocco `publish:` elenca l'output della chiamata del processo `sayHello()`, che viene indicato come `sayHello.out` e gli viene dato il nome `first_output` (questo può essere qualsiasi cosa l'autore del flusso di lavoro voglia).

Questa è una definizione di **workflow** molto minimale.
In una pipeline del mondo reale, il flusso di lavoro contiene tipicamente più chiamate a **process** connessi da **channel**, e potrebbero esserci valori predefiniti impostati per gli input variabili.

Approfondiremo questo nella Parte 2 del corso.
Per ora, diamo un'occhiata più da vicino a come il nostro flusso di lavoro gestisce input e output.

### 3.4. Il sistema `params` dei parametri da riga di comando

Il `params.input` che forniamo alla chiamata del processo `sayHello()` è un pezzo interessante di codice Nextflow e vale la pena spendere un minuto in più su di esso.

Come menzionato sopra, è così che passiamo il valore del parametro da riga di comando `--input` alla chiamata del processo `sayHello()`.
Infatti, semplicemente dichiarare `params.someParameterName` è sufficiente per dare al flusso di lavoro un parametro chiamato `--someParameterName` dalla riga di comando.

Qui abbiamo formalizzato quella dichiarazione di parametro impostando un blocco `params` che specifica il tipo di input che il flusso di lavoro si aspetta (Nextflow 25.10.2 e successivi).

```groovy title="1-hello.nf" linenums="20"
/*
 * Parametri della pipeline
 */
params {
    input: String
}
```

I tipi supportati includono `String`, `Integer`, `Float`, `Boolean` e `Path`.
Per saperne di più, consultate [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) nella documentazione di riferimento di Nextflow.

!!! tip

    Ricordate che i parametri del _flusso di lavoro_ dichiarati usando il sistema `params` prendono sempre due trattini sulla riga di comando (`--`).
    Questo li distingue dai flag CLI a _livello Nextflow_, che prendono solo un trattino (`-`).

### 3.5. La direttiva `publish`

Dall'altra parte del flusso di lavoro, abbiamo già dato un'occhiata al blocco `publish:`.
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

La riga `mode 'copy'` sovrascrive il comportamento predefinito del sistema, che è fare un link simbolico (o symlink) al file originale nella directory `work/` invece di una copia vera e propria.

Ci sono più opzioni di quelle mostrate qui per controllare il comportamento di pubblicazione; ne tratteremo alcune più avanti.
Vedrete anche che quando un flusso di lavoro genera più output, ciascuno viene elencato in questo modo nel blocco `output`.

Per saperne di più, consultate [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) nella documentazione di riferimento di Nextflow.

??? info "Sintassi più vecchia per pubblicare output usando `publishDir`"

    Fino a molto recentemente, il modo consolidato per pubblicare output era farlo a livello di ogni singolo processo usando una direttiva `publishDir`.

    Troverete ancora questo pattern di codice ovunque nelle pipeline Nextflow più vecchie e nei moduli di processo, quindi è importante esserne consapevoli.

    Invece di avere un blocco `publish:` nel flusso di lavoro e un blocco `output` al livello superiore, vedreste una riga `publishDir` nella definizione del processo `sayHello`:

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
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

    Tuttavia, non raccomandiamo di usare questo in nessun nuovo lavoro poiché sarà eventualmente vietato nelle versioni future del linguaggio Nextflow.

### Takeaway

Ora sapete come è strutturato un semplice flusso di lavoro Nextflow, e come i componenti di base si relazionano alla sua funzionalità.

### Cosa c'è dopo?

Imparate a gestire comodamente le vostre esecuzioni di flussi di lavoro.

---

## 4. Gestire le esecuzioni dei flussi di lavoro

Sapere come lanciare flussi di lavoro e recuperare output è fantastico, ma scoprirete rapidamente che ci sono alcuni altri aspetti della gestione dei flussi di lavoro che vi renderanno la vita più facile.

Qui vi mostriamo come sfruttare la funzionalità `resume` per quando dovete rilanciare lo stesso flusso di lavoro, come ispezionare i registri di esecuzione con `nextflow log`, e come eliminare le directory di lavoro più vecchie con `nextflow clean`.

### 4.1. Rilanciare un flusso di lavoro con `-resume`

A volte, vorrete rieseguire una pipeline che avete già lanciato in precedenza senza rifare alcun lavoro che è già stato completato con successo.

Nextflow ha un'opzione chiamata `-resume` che vi permette di farlo.
Specificamente, in questa modalità, qualsiasi processo che è già stato eseguito con esattamente lo stesso codice, impostazioni e input verrà saltato.
Questo significa che Nextflow eseguirà solo i processi che avete aggiunto o modificato dall'ultima esecuzione, o ai quali state fornendo nuove impostazioni o input.

Ci sono due vantaggi chiave nel fare questo:

- Se state sviluppando una pipeline, potete iterare più rapidamente poiché dovete eseguire solo i processi su cui state lavorando attivamente per testare le vostre modifiche.
- Se state eseguendo una pipeline in produzione e qualcosa va storto, in molti casi potete correggere il problema e rilanciare la pipeline, e riprenderà l'esecuzione dal punto di fallimento, il che può farvi risparmiare molto tempo e calcolo.

Per usarlo, aggiungete semplicemente `-resume` al vostro comando ed eseguitelo:

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

Cercate la parte `cached:` che è stata aggiunta nella riga di stato del processo (riga 5), che significa che Nextflow ha riconosciuto di aver già fatto questo lavoro e ha semplicemente riutilizzato il risultato dall'esecuzione precedente riuscita.

Potete anche vedere che l'hash della sottodirectory di lavoro è lo stesso dell'esecuzione precedente.
Nextflow vi sta letteralmente indicando l'esecuzione precedente e dicendo "L'ho già fatto laggiù."

!!! tip

    Quando rieseguite una pipeline con `resume`, Nextflow non sovrascrive alcun file pubblicato al di fuori della directory di lavoro da qualsiasi esecuzione che è stata eseguita con successo in precedenza.

    Per saperne di più, consultate [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) nella documentazione di riferimento di Nextflow.

### 4.2. Ispezionare il log delle esecuzioni passate

Ogni volta che lanciate un flusso di lavoro nextflow, viene scritta una riga in un file di log chiamato `history`, sotto una directory nascosta chiamata `.nextflow` nella directory di lavoro corrente.

??? abstract "Contenuto del file"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Questo file vi fornisce il timestamp, il nome dell'esecuzione, lo stato, l'ID di revisione, l'ID di sessione e la riga di comando completa per ogni esecuzione Nextflow che è stata lanciata dall'interno della directory di lavoro corrente.

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

Questo produrrà il contenuto del file di log nel terminale, arricchito con una riga di intestazione.

Noterete che l'ID di sessione cambia ogni volta che eseguite un nuovo comando `nextflow run`, TRANNE se state usando l'opzione `-resume`.
In quel caso, l'ID di sessione rimane lo stesso.

Nextflow usa l'ID di sessione per raggruppare le informazioni di caching dell'esecuzione sotto la directory `cache`, anch'essa situata sotto `.nextflow`.

### 4.3. Eliminare le directory di lavoro più vecchie

Se eseguite molte pipeline, potreste finire per accumulare moltissimi file attraverso molte sottodirectory.
Poiché le sottodirectory sono nominate casualmente, è difficile dire dai loro nomi quali sono esecuzioni più vecchie rispetto a quelle più recenti.

Fortunatamente Nextflow include un comando utile chiamato [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) che può eliminare automaticamente le sottodirectory di lavoro per esecuzioni passate che non vi interessano più.

#### 4.3.1. Determinare i criteri di eliminazione

Ci sono più opzioni per determinare cosa eliminare, che potete esplorare nella documentazione linkata sopra.
Qui vi mostriamo un esempio che elimina tutte le sottodirectory da esecuzioni precedenti a una data esecuzione, specificata usando il suo nome di esecuzione.

Cercate l'esecuzione riuscita più recente in cui non avete usato `-resume`; nel nostro caso il nome dell'esecuzione era `backstabbing_swartz`.

Il nome dell'esecuzione è la stringa generata automaticamente in due parti mostrata tra parentesi quadre nella riga di output della console `Launching (...)`.
Potete anche usare il log Nextflow per cercare un'esecuzione in base al suo timestamp e/o alla riga di comando.

#### 4.3.2. Fare una prova a secco

Prima usiamo il flag di prova a secco `-n` per verificare cosa verrà eliminato dato il comando:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Output del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Il vostro output avrà nomi di directory di attività diversi e potrebbe avere un numero diverso di righe, ma dovrebbe sembrare simile all'esempio.

Se non vedete alcuna riga di output, o non avete fornito un nome di esecuzione valido o non ci sono esecuzioni passate da eliminare. Assicuratevi di cambiare `backstabbing_swartz` nel comando di esempio con qualunque sia il nome dell'esecuzione più recente corrispondente nel vostro log.

#### 4.3.3. Procedere con l'eliminazione

Se l'output sembra come previsto e volete procedere con l'eliminazione, rieseguite il comando con il flag `-f` invece di `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Output del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

L'output dovrebbe essere simile a prima, ma ora dice 'Removed' invece di 'Would remove'.
Notate che questo non rimuove le sottodirectory di due caratteri (come `eb/` sopra) ma ne svuota il contenuto.

!!! Warning

    L'eliminazione delle sottodirectory di lavoro da esecuzioni passate le rimuove dalla cache di Nextflow ed elimina qualsiasi output che era memorizzato in quelle directory.
    Ciò significa che interrompe la capacità di Nextflow di riprendere l'esecuzione senza rieseguire i processi corrispondenti.

    Siete responsabili di salvare qualsiasi output che vi interessa! Questa è la ragione principale per cui preferiamo usare la modalità `copy` piuttosto che la modalità `symlink` per la direttiva `publish`.

### Takeaway

Sapete come rilanciare una pipeline senza ripetere passi che sono già stati eseguiti in modo identico, ispezionare il log di esecuzione e usare il comando `nextflow clean` per pulire le vecchie directory di lavoro.

### Cosa c'è dopo?

Prendetevi una piccola pausa! Avete appena assorbito i blocchi costitutivi della sintassi Nextflow e le istruzioni di utilizzo di base.

Nella prossima sezione di questa formazione, esamineremo quattro versioni successivamente più realistiche della pipeline Hello World che dimostreranno come Nextflow vi permette di elaborare più input in modo efficiente, eseguire flussi di lavoro composti da più passi connessi insieme, sfruttare componenti di codice modulari e utilizzare container per maggiore riproducibilità e portabilità.

---

## Quiz

<quiz>
Nella riga di output della console `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, cosa rappresenta `[a3/7be2fa]`?
- [ ] Il numero di versione del processo
- [ ] Un identificatore univoco dell'esecuzione
- [x] Il percorso troncato alla directory di lavoro dell'attività
- [ ] Il checksum del file di output

Per saperne di più: [2.3. Trovare l'output originale e i registri nella directory `work/`](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Qual è lo scopo del file `.command.sh` in una directory di attività?
- [ ] Memorizza le impostazioni di configurazione dell'attività
- [x] Mostra il comando effettivo che è stato eseguito dal processo
- [ ] Contiene messaggi di errore da attività fallite
- [ ] Elenca i file di input preparati per l'attività

Per saperne di più: [2.3. Trovare l'output originale e i registri nella directory `work/`](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Cosa succede ai risultati pubblicati quando rieseguite un flusso di lavoro senza `-resume`?
- [ ] Vengono preservati in directory separate con timestamp
- [x] Vengono sovrascritti dalla nuova esecuzione
- [ ] Nextflow impedisce la sovrascrittura e fallisce
- [ ] Vengono automaticamente salvati come backup

Per saperne di più: [2.4. Rieseguire il flusso di lavoro con saluti diversi](#24-re-run-the-workflow-with-different-greetings)
</quiz>

<quiz>
Cosa indica questo output della console?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] L'attività è fallita ed è stata saltata
- [ ] L'attività è in attesa in una coda
- [x] Nextflow ha riutilizzato i risultati da un'esecuzione identica precedente
- [ ] L'attività è stata annullata manualmente

Per saperne di più: [4.1. Rilanciare un flusso di lavoro con `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Dove memorizza Nextflow la cronologia di esecuzione che il comando `nextflow log` visualizza?
- [ ] Nella directory results
- [ ] Nella directory work
- [x] Nel file `.nextflow/history`
- [ ] In `nextflow.config`

Per saperne di più: [4.2. Ispezionare il log delle esecuzioni passate](#42-inspect-the-log-of-past-executions)
</quiz>

<quiz>
Qual è lo scopo del blocco `params` in un file di flusso di lavoro?
- [ ] Definire i requisiti di risorse del processo
- [ ] Configurare l'executor
- [x] Dichiarare e tipizzare i parametri di input del flusso di lavoro
- [ ] Specificare le opzioni di pubblicazione dell'output

Per saperne di più: [3.4. Il sistema params dei parametri da riga di comando](#34-the-params-system-of-command-line-parameters)
</quiz>

<quiz>
Nel blocco `output` del flusso di lavoro, cosa fa `mode 'copy'`?
- [ ] Crea un backup della directory di lavoro
- [x] Crea una copia completa dei file invece di link simbolici
- [ ] Copia lo script del flusso di lavoro nei results
- [ ] Abilita la copia incrementale dei file

Per saperne di più: [3.5. La direttiva publish](#35-the-publish-directive)
</quiz>

<quiz>
Qual è il flag raccomandato da usare con il comando `nextflow clean` prima di eliminare effettivamente i file?
- [x] `-n` (prova a secco) per visualizzare in anteprima cosa verrebbe eliminato
- [ ] `-v` (verbose) per vedere l'output dettagliato
- [ ] `-a` (all) per selezionare tutte le directory
- [ ] `-q` (quiet) per sopprimere gli avvisi

Per saperne di più: [4.3. Eliminare le directory di lavoro più vecchie](#43-delete-older-work-directories)
</quiz>
