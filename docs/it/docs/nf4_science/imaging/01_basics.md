# Parte 1: Eseguire operazioni di base

In questa prima parte del corso di formazione Nextflow per il Bioimaging, useremo un esempio Hello World molto semplice e indipendente dal dominio per dimostrare le operazioni essenziali e indicare i componenti corrispondenti del codice Nextflow.

## 1. Eseguire il workflow

Vi forniamo uno script di workflow chiamato `hello-world.nf` che riceve un input tramite un argomento da riga di comando chiamato `--greeting` e produce un file di testo contenente quel saluto.
Non esamineremo ancora il codice; prima vediamo come si presenta la sua esecuzione.

### 1.1. Lanciare il workflow e monitorare l'esecuzione

Nel terminale, eseguite il seguente comando:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

L'output della console dovrebbe apparire così:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Congratulazioni, avete appena eseguito il vostro primo workflow Nextflow!

L'output più importante qui è l'ultima riga (riga 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Questo ci dice che il processo `sayHello` è stato eseguito con successo una volta (`1 of 1 ✔`).

È fantastico, ma potreste chiedervi: dove si trova l'output?

### 1.2. Trovare il file di output nella directory `results`

Questo workflow è configurato per pubblicare il suo output in una directory chiamata `results`.
Se guardate la vostra directory corrente, vedrete che quando avete eseguito il workflow, Nextflow ha creato una nuova directory chiamata `results`, che contiene un file chiamato `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Aprite il file; il contenuto dovrebbe corrispondere al saluto che avete specificato dalla riga di comando.

<details>
  <summary>Contenuto del file</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Ottimo, il nostro workflow ha fatto quello che doveva fare!

Tuttavia, tenete presente che il risultato 'pubblicato' è una copia (o in alcuni casi un symlink) dell'output effettivo prodotto da Nextflow quando ha eseguito il workflow.

Quindi ora daremo un'occhiata sotto il cofano per vedere dove Nextflow ha effettivamente eseguito il lavoro.

!!! warning

    Non tutti i workflow saranno configurati per pubblicare gli output in una directory results, e/o il nome della directory potrebbe essere diverso.
    Un po' più avanti in questa sezione, vi mostreremo come scoprire dove è specificato questo comportamento.

### 1.3. Trovare l'output originale e i registri nella directory `work/`

Quando eseguite un workflow, Nextflow crea una 'directory di attività' distinta per ogni singola invocazione di ciascun processo nel workflow (=ogni passaggio nella pipeline).
Per ciascuna, preparerà gli input necessari, eseguirà le istruzioni rilevanti e scriverà output e file di log all'interno di quella directory, che viene nominata automaticamente usando un hash per renderla univoca.

Tutte queste directory di attività si troveranno sotto una directory chiamata `work` all'interno della vostra directory corrente (dove state eseguendo il comando).

Questo può sembrare confuso, quindi vediamo come appare in pratica.

Tornando all'output della console per il workflow che abbiamo eseguito prima, avevamo questa riga:

```console title="Excerpt of command output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Vedete come la riga inizia con `[a3/7be2fa]`?
Questa è una forma troncata del percorso della directory di attività per quella singola chiamata di processo, e vi dice dove trovare l'output della chiamata del processo `sayHello` all'interno del percorso della directory `work/`.

Potete trovare il percorso completo digitando il seguente comando (sostituendo `a3/7be2fa` con quello che vedete nel vostro terminale) e premendo il tasto tab per completare automaticamente il percorso o aggiungendo un asterisco:

```bash
tree work/a3/7be2fa*
```

Questo dovrebbe produrre il percorso completo della directory: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Diamo un'occhiata a cosa c'è dentro.

!!! Tip

    Se navigate il contenuto della sottodirectory di attività nell'esploratore file di VSCode, vedrete tutti i file immediatamente.
    Tuttavia, i file di log sono impostati per essere invisibili nel terminale, quindi se volete usare `ls` o `tree` per visualizzarli, dovrete impostare l'opzione rilevante per visualizzare i file invisibili.

    ```bash
    tree -a work
    ```

I nomi esatti delle sottodirectory saranno diversi sul vostro sistema.

<details>
  <summary>Contenuto della directory</summary>

```console title="work/"
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

</details>

Dovreste riconoscere immediatamente il file `output.txt`, che è infatti l'output originale del processo `sayHello` che è stato pubblicato nella directory `results`.
Se lo aprite, troverete di nuovo il saluto `Hello World!`.

<details>
  <summary>Contenuto del file output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

E tutti quegli altri file?

Questi sono i file di supporto e di log che Nextflow ha scritto come parte dell'esecuzione dell'attività:

- **`.command.begin`**: File sentinella creato non appena l'attività viene lanciata.
- **`.command.err`**: Messaggi di errore (`stderr`) emessi dalla chiamata del processo
- **`.command.log`**: Output di log completo emesso dalla chiamata del processo
- **`.command.out`**: Output regolare (`stdout`) della chiamata del processo
- **`.command.run`**: Script completo eseguito da Nextflow per eseguire la chiamata del processo
- **`.command.sh`**: Il comando che è stato effettivamente eseguito dalla chiamata del processo
- **`.exitcode`**: Il codice di uscita risultante dal comando

Il file `.command.sh` è particolarmente utile perché mostra il comando principale che Nextflow ha eseguito senza includere tutta la contabilità e la configurazione dell'attività/ambiente.

<details>
  <summary>Contenuto del file</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip

    Quando qualcosa va storto e dovete risolvere il problema, può essere utile guardare lo script `command.sh` per verificare esattamente quale comando Nextflow ha composto in base alle istruzioni del workflow, all'interpolazione delle variabili e così via.

### 1.4. Esercizio opzionale: rieseguire con saluti diversi

Provate a rieseguire il workflow alcune volte con valori diversi per l'argomento `--greeting`, quindi guardate sia il contenuto della directory `results/` che le directory di attività.

Osservate come gli output e i registri delle directory di attività isolate vengono preservati, mentre il contenuto della directory `results` viene sovrascritto dall'output delle esecuzioni successive.

### Takeaway

Sapete come eseguire un semplice script Nextflow, monitorare la sua esecuzione e trovare i suoi output.

### Cosa c'è dopo?

Imparate a leggere uno script Nextflow di base e identificare come i suoi componenti si relazionano alla sua funzionalità.

---

## 2. Esaminare lo script iniziale del workflow Hello World

Quello che abbiamo fatto è stato essenzialmente trattare lo script del workflow come una scatola nera.
Ora che abbiamo visto cosa fa, apriamo la scatola e guardiamo dentro.

_L'obiettivo qui non è memorizzare la sintassi del codice Nextflow, ma formare un'intuizione di base di quali sono i componenti principali e come sono organizzati._

### 2.1. Esaminare la struttura complessiva del codice

Apriamo lo script `hello-world.nf` nel pannello dell'editor.

<details>
  <summary>Codice</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Usa echo per stampare un saluto in un file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // emit a greeting
    sayHello(params.greeting)
}
```

</details>

Uno script Nextflow coinvolge due tipi principali di componenti fondamentali: uno o più **process**, e il **workflow** stesso.
Ogni **process** descrive quali operazioni il passaggio corrispondente nella pipeline dovrebbe compiere, mentre il **workflow** descrive la logica del flusso di dati che collega i vari passaggi.

Diamo un'occhiata più da vicino prima al blocco **process**, poi esamineremo il blocco **workflow**.

### 2.2. La definizione del `process`

Il primo blocco di codice descrive un **process**.
La definizione del processo inizia con la parola chiave `process`, seguita dal nome del processo e infine dal corpo del processo delimitato da parentesi graffe.
Il corpo del processo deve contenere un blocco script che specifica il comando da eseguire, che può essere qualsiasi cosa sareste in grado di eseguire in un terminale a riga di comando.

Qui abbiamo un **process** chiamato `sayHello` che riceve una variabile di **input** chiamata `greeting` e scrive il suo **output** in un file chiamato `output.txt`.

<details>
  <summary>Codice</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Usa echo per stampare un saluto in un file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

Questa è una definizione di processo molto minimale che contiene solo una definizione di `input`, una definizione di `output` e lo `script` da eseguire.

La definizione di `input` include il qualificatore `val`, che dice a Nextflow di aspettarsi un valore di qualche tipo (può essere una stringa, un numero, qualsiasi cosa).

La definizione di `output` include il qualificatore `path`, che dice a Nextflow che questo dovrebbe essere gestito come un percorso (include sia percorsi di directory che file).

!!! Tip

    La definizione di output non _determina_ quale output verrà creato.
    Semplicemente _dichiara_ dove trovare i file di output attesi, in modo che Nextflow possa cercarli una volta completata l'esecuzione.

    Questo è necessario per verificare che il comando sia stato eseguito con successo e per passare l'output ai processi a valle se necessario.
    L'output prodotto che non corrisponde a quanto dichiarato nel blocco output non verrà passato ai processi a valle.

In una pipeline del mondo reale, un processo di solito contiene informazioni aggiuntive come le direttive di processo, che introdurremo tra poco.

### 2.3. La definizione del `workflow`

Il secondo blocco di codice descrive il **workflow** stesso.
La definizione del workflow inizia con la parola chiave `workflow`, seguita da un nome opzionale, poi dal corpo del workflow delimitato da parentesi graffe.

Qui abbiamo un **workflow** che consiste in una chiamata al processo `sayHello`, che riceve un input, `params.greeting`, che contiene il valore che abbiamo dato al parametro `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emit a greeting
    sayHello(params.greeting)
}
```

Questa è una definizione di **workflow** molto minimale.
In una pipeline del mondo reale, il workflow tipicamente contiene più chiamate a **process** connessi da **channel**, e potrebbero esserci valori predefiniti impostati per gli input variabili.

Vedremo questo in azione quando eseguiremo nf-core/molkart nella Parte 2 del corso.

### 2.4. Il sistema `params` dei parametri da riga di comando

Il `params.greeting` che forniamo alla chiamata del processo `sayHello()` è un pezzo di codice Nextflow interessante e vale la pena spenderci un minuto in più.

Come menzionato sopra, è così che passiamo il valore del parametro da riga di comando `--greeting` alla chiamata del processo `sayHello()`.
Infatti, semplicemente dichiarando `params.someParameterName` ci permetterà di dare al workflow un parametro chiamato `--someParameterName` dalla riga di comando.

!!! Tip

    Questi parametri del workflow dichiarati usando il sistema `params` prendono sempre due trattini (`--`).
    Questo li distingue dai parametri a livello Nextflow, che prendono solo un trattino (`-`).

### Takeaway

Ora sapete come è strutturato un semplice workflow Nextflow, e come i componenti di base si relazionano alla sua funzionalità.

### Cosa c'è dopo?

Imparate a gestire le vostre esecuzioni di workflow in modo conveniente.

---

## 3. Gestire le esecuzioni del workflow

Sapere come lanciare i workflow e recuperare gli output è fantastico, ma scoprirete rapidamente che ci sono alcuni altri aspetti della gestione del workflow che vi renderanno la vita più facile.

Qui vi mostriamo come sfruttare la funzionalità `resume` per quando dovete rilanciare lo stesso workflow, come ispezionare i registri di esecuzione con `nextflow log`, e come eliminare le directory di lavoro più vecchie con `nextflow clean`.

### 3.1. Rilanciare un workflow con `-resume`

A volte, vorrete rieseguire una pipeline che avete già lanciato in precedenza senza rifare alcun lavoro che è già stato completato con successo.

Nextflow ha un'opzione chiamata `-resume` che vi permette di farlo.
Specificamente, in questa modalità, qualsiasi processo che è già stato eseguito con esattamente lo stesso codice, impostazioni e input verrà saltato.
Questo significa che Nextflow eseguirà solo i processi che avete aggiunto o modificato dall'ultima esecuzione, o ai quali state fornendo nuove impostazioni o input.

Ci sono due vantaggi chiave nel fare questo:

- Se state sviluppando una pipeline, potete iterare più rapidamente poiché dovete eseguire solo i processi su cui state lavorando attivamente per testare le vostre modifiche.
- Se state eseguendo una pipeline in produzione e qualcosa va storto, in molti casi potete risolvere il problema e rilanciare la pipeline, e riprenderà l'esecuzione dal punto di fallimento, il che può farvi risparmiare molto tempo e calcolo.

Per usarlo, aggiungete semplicemente `-resume` al vostro comando ed eseguitelo:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Cercate il bit `cached:` che è stato aggiunto nella riga di stato del processo (riga 5), che significa che Nextflow ha riconosciuto di aver già fatto questo lavoro e ha semplicemente riutilizzato il risultato dell'esecuzione precedente riuscita.

Potete anche vedere che l'hash della sottodirectory di lavoro è lo stesso dell'esecuzione precedente.
Nextflow vi sta letteralmente indicando l'esecuzione precedente e dicendo "L'ho già fatto laggiù."

!!! Tip

    Quando rieseguite una pipeline con `resume`, Nextflow non sovrascrive alcun file scritto in una directory `publishDir` da qualsiasi chiamata di processo che è stata precedentemente eseguita con successo.

### 3.2. Ispezionare il log delle esecuzioni passate

Ogni volta che lanciate un workflow nextflow, viene scritta una riga in un file di log chiamato `history`, sotto una directory nascosta chiamata `.nextflow` nella directory di lavoro corrente.

Un modo più conveniente per accedere a queste informazioni è usare il comando `nextflow log`.

```bash
nextflow log
```

Questo produrrà il contenuto del file di log nel terminale, mostrandovi il timestamp, il nome dell'esecuzione, lo stato e la riga di comando completa per ogni esecuzione Nextflow che è stata lanciata dall'interno della directory di lavoro corrente.

### 3.3. Eliminare le directory di lavoro più vecchie

Durante il processo di sviluppo, tipicamente eseguirete le vostre bozze di pipeline un gran numero di volte, il che può portare a un accumulo di moltissimi file attraverso molte sottodirectory.
Poiché le sottodirectory sono nominate casualmente, è difficile dire dai loro nomi quali sono le esecuzioni più vecchie rispetto a quelle più recenti.

Nextflow include un comodo sottocomando `clean` che può eliminare automaticamente le sottodirectory di lavoro per le esecuzioni passate che non vi interessano più, con diverse [opzioni](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per controllare cosa verrà eliminato.

Potete usare il log Nextflow per cercare un'esecuzione in base al suo timestamp e/o riga di comando, quindi usare `nextflow clean -before <run_name> -f` per eliminare le directory di lavoro dalle esecuzioni precedenti.

!!! Warning

    Eliminare le sottodirectory di lavoro dalle esecuzioni passate le rimuove dalla cache di Nextflow ed elimina qualsiasi output che era memorizzato in quelle directory.
    Questo significa che interrompe la capacità di Nextflow di riprendere l'esecuzione senza rieseguire i processi corrispondenti.

    Siete responsabili di salvare qualsiasi output che vi interessa o su cui pianificate di fare affidamento! Se state usando la direttiva `publishDir` per questo scopo, assicuratevi di usare la modalità `copy`, non la modalità `symlink`.

### Takeaway

Sapete come rilanciare una pipeline senza ripetere i passaggi che sono già stati eseguiti in modo identico, ispezionare il log di esecuzione e usare il comando `nextflow clean` per pulire le vecchie directory di lavoro.

### Cosa c'è dopo?

Ora che comprendete le operazioni di base di Nextflow, siete pronti per eseguire una vera pipeline di bioimaging con nf-core/molkart.
