# Parte 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [tutta la playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale Youtube Nextflow.

:green_book: La trascrizione di questo video è disponibile [qui](./transcripts/01_hello_world.md).
///

In questa prima parte del corso di formazione Hello Nextflow, ci addentriamo nell'argomento con un esempio di Hello World molto elementare e indipendente dal dominio, che svilupperemo progressivamente per dimostrare l'uso della logica e dei componenti fondamentali di Nextflow.

!!! note

      Un “Hello World!” è un esempio minimalista che ha lo scopo di dimostrare la sintassi e la struttura di base di un linguaggio di programmazione o di un framework software. L'esempio consiste tipicamente nello stampare la frase “Hello, World!” su un dispositivo di output, come la console o il terminale, o nello scriverla su un file.

---

## 0. Riscaldamento: Eseguire direttamente Hello World

Possiamo farlo con un semplice comando eseguito direttamente nel terminale, per mostrare cosa fa prima di adentrarci in Nextflow.

!!! tip

    Ricordate che ora dovreste trovarvi all'interno della cartella `hello-nextflow/`, come descritto nella sezione Orientamento.

### 0.1. Facciamo dire al terminale "Hello"

```bash
echo 'Hello World!'
```

Questo testo viene inviato al terminale con la scritta 'Hello World'.

```console title="Output"
Hello World!
```

### 0.2. Ora fate in modo che scriva il testo output in un file

```bash
echo 'Hello World!' > output.txt
```

Questo non produce alcun output nel terminale.

```console title="Output"

```

### 0.3. Mostra il contenuto del file

```bash
cat output.txt
```

Adesso, il testo 'Hello World' si trova nel file output che abbiamo specificato.

```console title="output.txt" linenums="1"
Hello World!
```

!!! tip

     Nell'ambiente di addestramento, è possibile trovare il file di output nell'esploratore di file e visualizzarne il contenuto facendo clic su di esso. In alternativa, si può usare il comando `code` per aprire il file per visualizzarlo

    ```bash
    code output.txt
    ```

### Takeaway

Ora sapete come eseguire un semplice comando nel terminale che produce un testo e, facoltativamente, come far scrivere l'output in un file.

### Cosa c'è dopo?

Scoprire come potrebbe essere scritto un flusso di lavoro Nextflow.

---

## 1. Esaminare lo script di avvio del flusso di lavoro Hello World

Come accennato nella guida, vi forniamo uno script di flusso di lavoro completamente funzionale, anche se minimalista, chiamato `hello-world.nf` che fa la stessa cosa di prima (scrivere “Hello World!”) ma con Nextflow.

Per iniziare, apriremo prima lo script del flusso di lavoro, in modo da avere un'idea di come è strutturato.

### 1.1. Esaminare la struttura complessiva del codice

Apriamo lo script `hello-world.nf` nel pannello dell'editor.

!!! note

    Il file si trova nella cartella `hello-nextflow`, che dovrebbe essere la cartella di lavoro corrente.
    È possibile fare clic sul file nell'esploratore di file, oppure digitare `ls` nel terminale e fare Cmd+clic (MacOS) o Ctrl+clic (PC) sul file per aprirlo.

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

    // emette un saluto
    sayHello()
}
```

Come si può vedere, uno script Nextflow comprende due tipi principali di componenti fondamentali: uno o più **process** e il **workflow** stesso.
Ogni **process** descrive le operazioni che il passo corrispondente della pipeline deve compiere, mentre il **workflow** descrive la logica del flusso di dati che collega i vari passi.

Vediamo prima il blocco **process** e poi il blocco **workflow**.

### 1.2. Definzione di `process`

Il primo blocco di codice descrive un **process**.
La definizione del process inizia con la parola chiave `process`, seguita dal nome del processo e infine dal corpo del processo delimitato da parentesi graffe.
Il corpo del processo deve contenere un blocco di script che specifica il comando da eseguire, che può essere qualsiasi cosa si possa eseguire in un terminale a riga di comando.
Qui abbiamo un **process** chiamato `sayHello` che produce un **output** in un file chiamato `output.txt`.

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

Si tratta di una definizione di processo molto minimale, che contiene solo una definizione di `output' e lo `script' da eseguire

La definizione `output` include il qualificatore `path`, che indica a Nextflow che questo deve essere gestito come un percorso (include sia percorsi di directory che di file).
Un altro qualificatore comune è `val`.

!!! note

    La definzione dell'output non _determina_ quale output verrà creato.
    Semplicemente _dichiara_ qual è l'output atteso, in modo che Nextflow possa cercarlo al termine dell'esecuzione.
    Questo è necessario per verificare che il comando sia stato eseguito correttamente e per passare l'output ai processi a valle, se necessario. L'output prodotto che non corrisponde a quanto dichiarato nel blocco di output non verrà passato ai processi a valle.

!!! warning

    Questo esempio è fragile perché abbiamo codificato il nome del file di output in due punti separati (lo script e i blocchi di output).
    Se ne cambiamo uno ma non l'altro, lo script si rompe.
    Più avanti, si imparerà a usare le variabili per evitare questo problema.

In una pipeline reale, un processo di solito contiene blocchi aggiuntivi come le direttive e gli input, che introdurremo tra poco.

### 1.3. La definizione di `workflow`

Il secondo blocco di codice descrive il **workflow** stesso.
La definizione di workflow comincia con la parola chiabe `workflow`, seguita da un nome opzionale, e dal corpo del workflow delimitato da parentesi graffe.

Qui abbiamo un **workflow** che consiste in una chiamata al processo `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    // emette un saluto
    sayHello()
}
```

Questa è definizione minimale del **workflow**.
In una pipeline reale, the workflow contiene tipicamente chiamate multiple a **processes** connesse da **channels**, e i processi prevedono una o più variabili in **input(s)**.

Si apprenderà come aggiungere ingressi variabili più avanti in questo modulo di formazione; e si apprenderà come aggiungere altri processi e collegarli tramite canali nella Parte 3 di questo corso.

### Takeaway

Adesso sai com'è strutturato un semplice worflow di Nextflow.

### Cosa c'è dopo?

Impara come lanciare un workflow e a monitorare l'esecuzione e a trovare i vostri output.

---

## 2. Esecuzione del workflow

Guarda il codice non è così divertente come eseguirlo, quindi proviamao a metterlo in pratica.

### 2.1. Lancio del workflow e monitorare l'esecuzione

Nel terminale, esegui i seguenti comandi:

```bash
nextflow run hello-world.nf
```

L'output della console dovrebbe apparire simile a questo:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Congratulazione, hai appena eseguito il tuo primo workflow i Nextflow.

L'output più importante è l'ultima riga (riga 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Questo ci dice che il processo `sayHello' è stato eseguito con successo una volta (`1 di 1 ✔`).

È importante che questa riga indichi anche dove trovare l'output della chiamata al processo `sayHello`.
Vediamolo ora.

### 2.2. Trovare l'output e i registri nella directory `work`

Quando si esegue Nextflow per la prima volta in una determinata directory, viene creata una directory chiamata `work` in cui verranno scritti tutti i file (e gli eventuali symlinks) generati nel corso dell'esecuzione.

All'interno della directory `work`, Nextflow organizza gli output e i registri per ogni chiamata di processo.
Per ogni chiamata di processo, Nextflow crea una sottodirectory nidificata, denominata con un hash per renderla unica, in cui inserisce tutti gli input necessari (usando i symlinks per impostazione predefinita), scrive i file di aiuto e scrive i log e gli output del processo.

Il percorso della subdirectory viene mostrato in forma tronca tra parentesi quadre nell'output della console.
Osservando i risultati dell'esecuzione mostrata sopra, la riga di log della console per il processo sayHello inizia con `[a3/7be2fa]`. Ciò corrisponde al seguente percorso di directory: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Diamo un occhiata a cosa c'è dentro.

!!! tip

    Se si sfoglia il contenuto della subdirectory task nel file explorer di VSCode, si vedranno subito tutti i file.
    Tuttavia, i file di log sono impostati per essere invisibili nel terminale, quindi se si vuole usare `ls` o `tree` per visualizzarli, è necessario impostare l'opzione corrispondente per la visualizzazione dei file invisibili.

     ```bash
    tree -a work
    ```

Si dovrebbe vedere qualcosa di simile, anche se i nomi esatti delle subdirectory saranno diversi sul vostro sistema:

```console title="Directory contents"
work
└── a3
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

Questi sono i file di aiuto e di registro:

- **`.command.begin`**: Metadati relativi all'inizio dell'esecuzione della chiamata di processo.
- **`.command.err`**: Messaggi di errore (`stderr`) emessi dalla chiamata al processo.
- **`.command.log`**: Output di log completo emesso dalla chiamata al processo
- **`.command.out`**: Output regolare (`stdout`) emesso dalla chiamata al processo
- **`.command.run`**: Script completo eseguito da Nextflow per eseguire la chiamata al processo
- **`.command.sh`**: Il comando che è stato effettivamente eseguito dalla chiamata di processo
- **`.exitcode`**: Il codice di uscita risultante dal comando

Il file `.command.sh` è particolarmente utile perché dice quale comando Nextflow ha effettivamente eseguito.
In questo caso è molto semplice, ma più avanti nel corso si vedranno comandi che comportano un'interpolazione di variabili.
Quando si ha a che fare con questi comandi, è necessario essere in grado di controllare esattamente cosa è stato eseguito, soprattutto quando si risolve un problema.

L'output effettivo del processo `sayHello` è `output.txt`.
Aprendola, troverete il saluto `Hello World!`, che era il risultato atteso del nostro flusso di lavoro minimalista.

```console title="output.txt" linenums="1"
Hello World!
```

### Takeaway

Sapete come decifrare un semplice script di Nextflow, eseguirlo e trovare l'output e i relativi file di log nella directory di lavoro.

### Cosa c'è dopo?

Imparate a gestire comodamente le esecuzioni dei vostri workflow.

---

## 3. Gestire le esecuzioni dei workflow

Sapere come lanciare workflow e recuperarne i risultati è ottimo, ma scoprirete subito che ci sono altri aspetti della gestione dei workflow che vi renderanno la vita più facile, soprattutto se ne state sviluppando uno di personale.

Qui mostriamo come usare la direttiva `publishDir` per memorizzare in una cartella di output tutti i risultati principali dell'esecuzione della pipeline, la funzione `resume` per quando è necessario rilanciare lo stesso workflow e come cancellare le vecchie directory di lavoro con `nextflow clean`.

### 3.1. Pubblicare i risultati

Come si è appena appreso, l'output prodotto dalla nostra pipeline è sepolto in una directory di lavoro a diversi livelli di profondità.
Questo è stato fatto di proposito; Nextflow ha il controllo di questa directory e noi non dobbiamo interagire con essa.

Tuttavia, questo rende scomodo recuperare gli output che ci interessano.

Fortunatamente, Nextflow offre un modo per gestire più comodamente questo aspetto, chiamato direttiva `publishDir`, che agisce a livello di processo.
Questa direttiva indica a Nextflow di pubblicare gli output del processo in una directory di output designata. Per impostazione predefinita, gli output sono pubblicati come collegamenti simbolici dalla directory `work`.
Permette di recuperare il file di output desiderato senza dover scavare nella directory di lavoro.

#### 3.1.1. Aggiungere una direttiva `publishDir` al processo `sayHello`

Nel file di script del flusso di lavoro `hello-world.nf`, apportare la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3"
    process sayHello {

        publishDir 'results', mode: 'copy'

        output:
            path 'output.txt'
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
            path 'output.txt'
    ```

#### 3.1.2. Eseguire nuovamente il workflow

Eseguite ora lo script del workflow modificato:

```bash
nextflow run hello-world.nf
```

L'output del log dovrebbe essere molto familiare

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

executor >  local (1)
[62/49a1f8] sayHello | 1 of 1 ✔
```

Questa volta, Nextflow ha creato una nuova directory chiamata `results/`.
Il nostro file `output.txt` si trova in questa directory.
Se si controlla il contenuto, dovrebbe corrispondere all'output della sottodirectory di lavoro.
In questo modo si pubblicano i file dei risultati al di fuori delle directory di lavoro.

Quando si ha a che fare con file molto grandi che non devono essere conservati a lungo, si può preferire impostare la direttiva `publishDir` per creare un collegamento simbolico al file invece di copiarlo.
Tuttavia, se si elimina la directory di lavoro come parte di un'operazione di pulizia, si perderà l'accesso al file, quindi assicuratevi sempre di avere copie effettive di tutto ciò che vi interessa prima di eliminare qualcosa.

!!! note

    È stata proposta una nuova opzione di sintassi documentata [qui](https://www.nextflow.io/docs/latest/workflow.html#publishing-outputs) per rendere possibile la dichiarazione e la pubblicazione di output a livello di flusso di lavoro.
    Questo renderà l'uso di `publishDir`, a livello del processo, ridondante per le pipeline completate.
    Tuttavia, ci aspettiamo che `publishDir` rimanga molto utile durante lo sviluppo della pipeline.

### 3.2. Rilanciare un workflow con `-resume`

A volte si desidera rieseguire una pipeline già avviata in precedenza, senza ripetere le fasi già completate con successo.

Nextflow ha un'opzione chiamata `resume` che permette di farlo.
In particolare, in questa modalità vengono saltati tutti i processi già eseguiti con lo stesso codice, le stesse impostazioni e gli stessi input.
Ciò significa che Nextflow eseguirà solo i processi aggiunti o modificati dall'ultima volta che sono stati eseguiti o per i quali sono state fornite nuove impostazioni o input.

Ci sono due vantaggi principali nel fare ciò:

- Se siete nel bel mezzo dello sviluppo della vostra pipeline, potete iterare più rapidamente, poiché dovete eseguire solo i processi su cui state lavorando attivamente per testare le vostre modifiche.
- Se si sta eseguendo una pipeline in produzione e qualcosa va storto, in molti casi è possibile risolvere il problema e rilanciare la pipeline, che riprenderà a funzionare dal punto di guasto, con un notevole risparmio di tempo e di calcolo.

Per usarlo, aggiungete semplicemente `-resume` al vostro comando ed eseguitelo:

```bash
nextflow run hello-world.nf -resume
```

L'output della console dovresse essere simile.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

[62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
```

Cercate il bit `cached:` che è stato aggiunto nella riga di stato del processo (riga 5), il che significa che Nextflow ha riconosciuto di aver già svolto questo lavoro e ha semplicemente riutilizzato il risultato della precedente esecuzione andata a buon fine.

Si può anche notare che l'hash della sottodirectory di lavoro è lo stesso dell'esecuzione precedente.
Nextflow indica letteralmente l'esecuzione precedente e dice: “L'ho già fatto lì”.

!!! note

    Quando si riesegue una pipeline con `resume`, Nextflow non sovrascrive alcun file scritto in una directory `publishDir` da qualsiasi chiamata di processo precedentemente eseguita con successo.

### 3.3. Eliminare vecchie work directories

Durante il processo di sviluppo, di solito si eseguono le pipeline in bozza un gran numero di volte, il che può portare all'accumulo di moltissimi file in molte sottodirectory.
Poiché le sottodirectory sono denominate in modo casuale, è difficile capire dai loro nomi quali sono le esecuzioni più vecchie e quali quelle più recenti.

Nextflow include un utile sottocomando `clean` che può cancellare automaticamente le sottodirectory di lavoro per le esecuzioni passate che non interessano più, con diverse [opzioni](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per controllare cosa verrà cancellato.

Qui viene mostrato un esempio che cancella tutte le sottodirectory delle esecuzioni precedenti a una determinata esecuzione, specificata con il suo nome di esecuzione.
Il nome dell'esecuzione è la stringa in due parti generata dalla macchina e mostrata tra parentesi quadre nella riga di output della console `Launching (...)`.

Per prima cosa usiamo il flag `-n` per verificare cosa verrà cancellato con il comando:

```bash
nextflow clean -before golden_cantor -n
```

L'output dovrebbe risultare il seguente:

```console title="Output"
Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
```

Se non viene visualizzata alcuna riga, significa che non è stato fornito un nome di sessione valido o che non ci sono sessioni precedenti da eliminare.

Se l'output appare come previsto e si vuole procedere con la cancellazione, rieseguire il comando con il flag `-f` invece di `-n`:

```bash
nextflow clean -before golden_cantor -f
```

A questo punto si dovrebbe vedere quanto segue:

```console title="Output"
Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
```

!!! Warning

    L'eliminazione delle sottodirectory di lavoro delle esecuzioni precedenti le rimuove dalla cache di Nextflow e cancella tutti gli output memorizzati in tali directory.
    Ciò significa che interrompe la capacità di Nextflow di riprendere l'esecuzione senza rieseguire i processi corrispondenti.

    L'utente è responsabile del salvataggio di tutti gli output a cui tiene o su cui intende fare affidamento! Se si usa la direttiva `publishDir` a tale scopo, assicurarsi di usare la modalità `copy`, non la modalità `symlink`.

### Takeaway

Sapete come pubblicare gli output in una directory specifica, rilanciare una pipeline senza ripetere i passaggi già eseguiti in modo identico e usare il comando `nextflow clean` per ripulire le vecchie directory di lavoro.

### Cosa c'è dopo?

Imparare a fornire una variabile in ingresso tramite un parametro della riga di comando e a utilizzare efficacemente i valori predefiniti.

---

## 4. Usare una variabile di input passata alla riga di comando

Nel suo stato attuale, il nostro flusso di lavoro utilizza un saluto codificato nel comando di processo.
Vogliamo aggiungere un po' di flessibilità utilizzando una variabile di input, in modo da poter cambiare più facilmente il saluto in fase di esecuzione.

### 4.1. Modificare il flusso di lavoro per accettare e utilizzare un input variabile

Ciò richiede di apportare tre modifiche al nostro script:

1. Indicare al processo di aspettarsi un input variabile aggiungendo un blocco `input:`
2. Modificare il processo per utilizzare l'input
3. Impostare un parametro della riga di comando e fornire il suo valore come input alla chiamata al processo

Apportiamo queste modifiche una alla volta.

#### 4.1.1. Aggiungere un blocco di input alla definizione del processo

Per prima cosa dobbiamo adattare la definizione del processo in modo che accetti un input chiamato `greeting`.

Nel blocco del processo, apportare la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="5 6"
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path 'output.txt'
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        publishDir 'results', mode: 'copy'

        output:
            path 'output.txt'
    ```

La variabile `greeting` è preceduta da `val` per indicare a Nextflow che si tratta di un valore (non di un percorso).

#### 4.1.2. Modificare il comando di processo per utilizzare la variabile di input

Ora scambiamo il valore originale codificato con il valore della variabile di ingresso che ci aspettiamo di ricevere.

Nel blocco del processo, apportare la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="16" hl_lines="3"
    script:
    """
    echo '$greeting' > output.txt
    """
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="16"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

Assicuratevi di aggiungere il simbolo `$` per indicare a Nextflow che si tratta di un nome di variabile che deve essere sostituito con il valore effettivo (=interpolato).

#### 4.1.3. Impostare un parametro CLI e fornirlo come input alla chiamata al processo

Ora dobbiamo impostare un modo per fornire un valore di input alla chiamata di processo `sayHello()`.

Si potrebbe semplicemente codificare direttamente scrivendo `sayHello('Hello World!')`.
Tuttavia, quando si lavora davvero con il flusso di lavoro, spesso si desidera poter controllare gli input dalla riga di comando.

Buone notizie: Nextflow ha un sistema di parametri incorporato nel flusso di lavoro chiamato `params`, che rende facile dichiarare e utilizzare i parametri della CLI. La sintassi generale è dichiarare `params.<nome_parametro>` per dire a Nextflow di aspettarsi un parametro `--<nome_parametro>` sulla riga di comando.

In questo caso, vogliamo creare un parametro chiamato `--greeting`, quindi dobbiamo dichiarare `params.greeting` da qualche parte nel flusso di lavoro.
In linea di principio, possiamo scriverla ovunque; ma poiché vogliamo darla alla chiamata di processo `sayHello()`, possiamo inserirla direttamente scrivendo `sayHello(params.greeting)`.

!!! note

    Il nome del parametro (a livello di flusso di lavoro) non deve necessariamente corrispondere al nome della variabile di input (a livello di processo).
    Usiamo la stessa parola perché è quella che ha senso e mantiene il codice leggibile.

Nel blocco del flusso di lavoro, apportare la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-world.nf" linenums="24" hl_lines="2"
    // emit a greeting
    sayHello(params.greeting)
    ```

=== "Prima"

    ```groovy title="hello-world.nf" linenums="24"
    // emit a greeting
    sayHello()
    ```

Indica a Nextflow di eseguire il processo `sayHello` sul valore fornito attraverso il parametro `--greeting`.

#### 4.1.4. Eseguite nuovamente il comando del flusso di lavoro

Eseguiamolo!

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

Se tutte e tre le modifiche sono state eseguite correttamente, si dovrebbe ottenere un'altra esecuzione di successo:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Avvio di `hello-world.nf` [elated_lavoisier] DSL2 - revisione: 7c031b42ea

executor >  local (1)
[4b/654319] sayHello | 1 of 1 ✔
```

Assicurarsi di aprire il file di output per verificare che sia disponibile la nuova versione del saluto.

```console title="results/output.txt" linenums="1"
Bonjour le monde!
```

Voilà!

!!! tip

    È possibile distinguere facilmente i parametri a livello di Nextflow da quelli a livello di pipeline.

    - I parametri che si applicano a una pipeline hanno sempre un doppio trattino (`--`).
    - I parametri che modificano un'impostazione di Nextflow, ad esempio la funzione `riprendi' usata in precedenza, sono caratterizzati da un singolo trattino (`-`).

### 4.2. Utilizzare i valori predefiniti per i parametri della riga di comando

In molti casi, ha senso fornire un valore predefinito per un determinato parametro, in modo da non doverlo specificare per ogni esecuzione.

#### 4.2.1. Impostare un valore predefinito per il parametro CLI

Diamo al parametro `greeting` un valore predefinito, dichiarandolo prima della definizione del flusso di lavoro.

```groovy title="hello-world.nf" linenums="22"
/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo!'
```

!!! tip

    Se si preferisce, si può inserire la dichiarazione dei parametri all'interno del blocco del flusso di lavoro. Qualunque sia la scelta, si cerchi di raggruppare le cose simili nello stesso posto, in modo da non ritrovarsi con dichiarazioni sparse ovunque.

#### 4.2.2. Eseguire nuovamente il flusso di lavoro senza specificare il parametro

Ora che è stato impostato un valore predefinito, è possibile eseguire nuovamente il flusso di lavoro senza dover specificare un valore nella riga di comando.

```bash
nextflow run hello-world.nf
```

L'output della console dovrebbe essere lo stesso.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

executor >  local (1)
[72/394147] sayHello | 1 of 1 ✔
```

Controllare l'output nella directory dei risultati:

```console title="results/output.txt" linenums="1"
Holà mundo!
```

Nextflow ha utilizzato il valore predefinito per assegnare un nome all'output

#### 4.2.3. Eseguire nuovamente il flusso di lavoro con il parametro per sovrascrivere il valore predefinito

Se si fornisce il parametro sulla riga di comando, il valore della CLI sovrascrive il valore predefinito.

Provate:

```bash
nextflow run hello-world.nf --greeting 'Konnichiwa!'
```

L'output della console dovrebbe essere lo stesso.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

executor >  local (1)
[6f/a12a91] sayHello | 1 of 1 ✔
```

Ora si avrà il nuovo output corrispondente nella cartella dei risultati.

```console title="results/output.txt" linenums="1"
Konnichiwa!
```

!!! note

    In Nextflow sono presenti più punti in cui è possibile specificare i valori dei parametri.
    Se lo stesso parametro è impostato su valori diversi in più punti, Nexflow determinerà quale valore utilizzare in base all'ordine di precedenza descritto [qui](https://www.nextflow.io/docs/latest/config.html).

### Takeaway

Sapete come utilizzare una semplice variabile di input fornita in fase di esecuzione tramite un parametro della riga di comando, nonché come impostare, utilizzare e sovrascrivere i valori predefiniti.

Più in generale, sapete come interpretare un semplice flusso di lavoro Nextflow, gestirne l'esecuzione e recuperare i risultati.

### Cosa c'è dopo?

Prendetevi una piccola pausa, ve la siete meritata!!
Quando siete pronti, passate alla Parte 2 per imparare a usare i canali per alimentare gli input nel vostro flusso di lavoro, il che vi permetterà di sfruttare il parallelismo dataflow integrato di Nextflow e altre potenti funzioni.
