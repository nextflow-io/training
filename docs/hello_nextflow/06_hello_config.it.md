# Parte 6: Hello Config

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [tutta la playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale Youtube Nextflow.

:green_book: La trascrizione di questo video è disponibile [qui](./transcripts/05_hello_containers.md).
///

Questa sezione esplorerà come impostare e gestire la configurazione della pipeline Nextflow in modo da poterne personalizzare il comportamento, adattarla a diversi ambienti e ottimizzare l'utilizzo delle risorse _senza modificare una sola riga del codice del flusso di lavoro stesso_.

Ci sono diversi modi per farlo; qui useremo il meccanismo di file di configurazione più semplice e comune, il file `nextflow.config`.
Ogni volta che c'è un file chiamato `nextflow.config` nella directory corrente, Nextflow caricherà automaticamente la configurazione da esso.

!!! note

    Tutto ciò che inserisci in `nextflow.config` può essere sovrascritto in fase di esecuzione fornendo le direttive di processo o i parametri e i valori pertinenti sulla riga di comando, oppure importando un altro file di configurazione, secondo l'ordine di precedenza descritto [qui](https://www.nextflow.io/docs/latest/config.html).

In questa parte della formazione, utilizzeremo il file `nextflow.config` per illustrare i componenti essenziali della configurazione di Nextflow, quali direttive di processo, esecutori, profili e file di parametri.

Imparando a utilizzare efficacemente queste opzioni di configurazione, puoi migliorare la flessibilità, la scalabilità e le prestazioni delle tue pipeline.

---

## 0. Warmup: verifica che Docker sia abilitato ed esegui il Hello Config workflow

Innanzitutto, un rapido controllo. C'è un file `nextflow.config` nella directory corrente che contiene la riga `docker.enabled = <setting>`, dove `<setting>` è `true` o `false` a seconda che tu abbia lavorato o meno alla Parte 5 di questo corso nello stesso ambiente.

Se è impostato su `true`, non è necessario fare nulla.

Se è impostato su `false`, impostalo ora su `true`.

```console title="nextflow.config" linenums="1"
docker.enabled = true
```

Una volta fatto questo, verifica che il workflow iniziale funzioni correttamente:

```bash
nextflow run hello-config.nf
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-config.nf` [reverent_heisenberg] DSL2 - revision: 028a841db1

executor >  local (8)
[7f/0da515] sayHello (1)       | 3 of 3 ✔
[f3/42f5a5] convertToUpper (3) | 3 of 3 ✔
[04/fe90e4] collectGreetings   | 1 of 1 ✔
[81/4f5fa9] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

Se tutto funziona, sei pronto per imparare come modificare le proprietà di configurazione di base per adattarle ai requisiti del tuo ambiente di elaborazione.

---

## 1. Determinare quale tecnologia di confezionamento software utilizzare

Il primo passo per adattare la configurazione del workflow al tuo ambiente di calcolo è specificare da dove proverranno i pacchetti software che verranno eseguiti in ogni fase.
Sono già installati nell'ambiente di calcolo locale? Dobbiamo recuperare le immagini ed eseguirle tramite un sistema a container? Oppure dobbiamo recuperare i pacchetti Conda e creare un ambiente Conda locale?

Nella primissima parte di questo corso di formazione (Parti 1-4) abbiamo utilizzato solo il software installato localmente nel nostro workflow.
Nella Parte 5 abbiamo poi introdotto i container Docker e il file `nextflow.config`, che abbiamo utilizzato per abilitare l'uso dei container Docker.

Nel warmup di questa sezione, hai verificato che Docker fosse abilitato nel file `nextflow.config` ed eseguito il workflow, che ha utilizzato un container Docker per eseguire il processo `cowpy()`.

!!! note

    Se questo non ti suona familiare, probabilmente dovresti tornare indietro e completare la Parte 5 prima di continuare.

Vediamo ora come possiamo configurare un'opzione alternativa di software packaging tramite il file `nextflow.config`.

### 1.1. Disabilitare Docker e abilitare Conda nel file di configurazione

Supponiamo di lavorare su un cluster HPC e che l'amministratore non consenta l'uso di Docker per motivi di sicurezza.

Fortunatamente per noi, Nextflow supporta molte altre tecnologie di container, tra cui Singularity (più ampiamente utilizzata su HPC) e gestori di pacchetti software come Conda.

Possiamo modificare il nostro file di configurazione per usare Conda invece di Docker.
Per farlo, cambiamo il valore di `docker.enabled` in `false` e ​​aggiungiamo una direttiva che abilita l'uso di Conda:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"

    docker.enabled = true
    ```

Ciò consentirà a Nextflow di creare e utilizzare ambienti Conda per i processi che hanno pacchetti Conda specificati.
Ciò significa che ora dobbiamo aggiungerne uno al nostro processo `cowpy`!

### 1.2. Specificare un pacchetto Conda nella definizione del processo

Abbiamo già recuperato l'URI per un pacchetto Conda contenente lo strumento `cowpy`: `conda-forge::cowpy==1.1.5`

!!! note

    Esistono diversi modi per ottenere l'URI per un determinato pacchetto conda.
    Consigliamo di usare la ricerca su [Seqera Containers](https://seqera.io/containers/), che ti fornirà un URI che puoi copiare e incollare, anche se non hai intenzione di creare un container da esso.

Ora aggiungiamo l'URI alla definizione del processo `cowpy` utilizzando la direttiva `conda`:

=== "Dopo"

    ```console title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        publishDir 'results', mode: 'copy'
    ```

=== "Prima"

    ```console title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        publishDir 'results', mode: 'copy'
    ```

Per essere chiari, non stiamo _sostituendo_ la direttiva `docker`, stiamo _aggiungendo_ un'opzione alternativa.

### 1.3. Eseguire il workflow per verificare l'uso di Conda

Proviamolo.

```bash
nextflow run hello-config.nf
```

Dovrebbe funzionare senza problemi.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-config.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

executor >  local (8)
[ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
[20/2596a7] convertToUpper (1) | 3 of 3 ✔
[b3/e15de5] collectGreetings   | 1 of 1 ✔
[c5/af5f88] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

Dietro le quinte, Nextflow ha recuperato i pacchetti Conda e creato l'ambiente, il che normalmente richiede un po' di lavoro; è bello quindi che non dobbiamo farlo manualmente!

!!! note

    Questa operazione è veloce perchè il pacchetto `cowpy` è piuttosto piccolo. Tuttavia,se si lavora con pacchetti di grandi dimensioni, il processo potrebbe richiedere più tempo al primo utilizzo, perchè si potrebbe vedere l'output della console rimanere "bloccato" per circa un minuto prima di completarsi.
    Ciò è normale ed è dovuto al lavoro extra che Nextflow esegue la prima volta che si utilizza un nuovo pacchetto.

Dal nostro punto di vista, sembra che funzioni esattamente come con Docker, anche se nel backend i meccanismi sono leggermente diversi.

Ciò significa che siamo pronti per l'esecuzione con gli ambienti Conda, se necessario.

!!!note

    Poiché queste direttive vengono assegnate per processo, è possibile "mescolare e abbinare", ovvero configurare alcuni processi nel workflow in modo che vengano eseguiti con Docker e altri con Conda, ad esempio se l'infrastruttura di elaborazione utilizzata supporta entrambi.
    In tal caso, dovresti abilitare sia Docker che Conda nel tuo file di configurazione.
    Se entrambi sono disponibili per un dato processo, Nextflow darà priorità ai container.

    Come accennato in precedenza, Nextflow supporta molte altre tecnologie di packaging e container software, quindi non sei limitato solo a queste due.

### Conclusione

Sai come configurare quale pacchetto software deve utilizzare ogni processo e come passare da una tecnologia all'altra.

### Prossimi passi

Scopri come modificare l'esecutore utilizzato da Nextflow per eseguire il lavoro.

---

## 2. Assegnare risorse di elaborazione con direttive di processo

La maggior parte delle piattaforme di elaborazione ad alte prestazioni consente (e talvolta richiede) di specificare determinati parametri di allocazione delle risorse, come il numero di CPU e la memoria.

Di default, Nextflow utilizzerà una singola CPU e 2 GB di memoria per ogni processo.
Le direttive di processo corrispondenti sono chiamate `cpus` e `memory`, quindi è implicita la seguente configurazione:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Puoi modificare questi valori, sia per tutti i processi che per specifici processi denominati, utilizzando direttive di processo aggiuntive nel tuo file di configurazione.
Nextflow le tradurrà nelle istruzioni appropriate per l'esecutore scelto.

Ma come fai a sapere quali valori utilizzare?

### 2.1. Eseguire il workflow per generare un report sull'utilizzo delle risorse

Se non si sa in anticipo quanta CPU e memoria saranno probabilmente necessarie ai propri processi, è possibile effettuare una profilazione delle risorse, ovvero eseguire il workflow con alcune allocazioni predefinite, registrare la quantità utilizzata da ciascun processo e, da lì, stimare come modificare le allocazioni di base.

Nextflow include strumenti integrati per fare questo e sarà felice di generare un report per te su richiesta.

Per farlo, aggiungi `-with-report <nomefile>.html` alla riga di comando.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Il report è un file html, che puoi scaricare e aprire nel tuo browser. Puoi anche fare clic destro su di esso nell'esploratore file a sinistra e cliccare su `Mostra anteprima` per visualizzarlo nell'ambiente di training.

Prenditi qualche minuto per esaminare il report e vedere se riesci a identificare alcune opportunità per regolare le risorse.
Assicurati di cliccare sulle schede che mostrano i risultati di utilizzo come percentuale di quanto è stato assegnato.
C'è della [documentazione](https://www.nextflow.io/docs/latest/reports.html) che descrive tutte le funzionalità disponibili.

<!-- TODO: insert images -->

### 2.2. Imposta le allocazioni delle risorse per tutti i processi

La profilazione mostra che i processi nel nostro workflow di formazione sono molto leggeri, quindi riduciamo l'allocazione di memoria predefinita a 1 GB per processo.

Aggiungere quanto segue al file `nextflow.config`:

```groovy title="nextflow.config" linenums="4"
process {
    memory = 1.GB
}
```

### 2.3. Imposta le allocazioni delle risorse per un singolo processo

Allo stesso tempo, faremo finta che il processo `cowpy` richieda più risorse degli altri, così da poter dimostrare come adattare le allocazioni per un singolo processo.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="14"
    process {
        memory = 1.GB
    }
    ```

Con questa configurazione, tutti i processi richiederanno 1 GB di memoria e una singola CPU (impostazione predefinita), tranne il processo `cowpy`, che richiederà 2 GB e 2 CPU.

!!! note

    Se hai una macchina con poche CPU e ne assegni un numero elevato per processo, potresti vedere le chiamate di processo accodarsi una dietro l'altra.
    Questo perché Nextflow assicura che non richiediamo più CPU di quelle disponibili.

### 2.4. Esegui il workflow con la configurazione modificata

Proviamo a fare una prova specificando un nome file diverso per il report di profilazione, in modo da poter confrontare le prestazioni prima e dopo le modifiche alla configurazione.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Probabilmente non noterai alcuna differenza reale, poiché si tratta di un carico di lavoro molto ridotto, ma questo è l'approccio che utilizzeresti per analizzare i requisiti di prestazioni e risorse di un workflow reale.

È molto utile quando i tuoi processi hanno requisiti di risorse diversi. Ti consente di dimensionare correttamente le allocazioni di risorse che imposti per ogni processo in base a dati effettivi, non a supposizioni.

!!!note

    Questo è solo un piccolo assaggio di ciò che puoi fare per ottimizzare l'uso delle risorse.
    Nextflow stesso ha una [logica di retry dinamica](https://training.nextflow.io/basic_training/debugging/#dynamic-resources-allocation) davvero interessante, utilizzata per riprovare i lavori che falliscono a causa di limitazioni di risorse.
    Inoltre, la piattaforma Seqera offre anche strumenti basati sull'intelligenza artificiale per ottimizzare automaticamente le allocazioni delle risorse.

    Parleremo di entrambi gli approcci in una prossima parte di questo corso di formazione.

### 2.5. Aggiungere limiti alle risorse

A seconda dell'esecutore di elaborazione e dell'infrastruttura di elaborazione che stai utilizzando, potrebbero esserci dei vincoli su ciò che puoi (o devi) allocare.
Ad esempio, il tuo cluster potrebbe richiedere di rimanere entro determinati limiti.

Puoi utilizzare la direttiva `resourceLimits` per impostare le limitazioni pertinenti. La sintassi appare così quando è da sola in un blocco di processo:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow tradurrà questi valori nelle istruzioni appropriate a seconda dell'esecutore specificato.

Non lo eseguiremo, poiché non abbiamo accesso all'infrastruttura pertinente nell'ambiente di formazione.
Tuttavia, se provassi a eseguire il workflow con allocazioni di risorse che superano questi limiti, quindi cercassi il comando `sbatch` nel file di script `.command.run`, vedresti che le richieste che vengono effettivamente inviate all'esecutore sono limitate ai valori specificati da `resourceLimits`.

!!!note

    Il progetto nf-core ha compilato una [raccolta di file di configurazione](https://nf-co.re/configs/) condivisa da varie istituzioni in tutto il mondo, che copre un'ampia gamma di esecutori HPC e cloud.

    Tali configurazioni condivise sono preziose sia per le persone che lavorano lì e che possono quindi utilizzare la configurazione della propria istituzione così com'è, sia come modello per coloro che desiderano sviluppare una configurazione per la propria infrastruttura.

### Conclusione

Sai come generare un report di profilazione per valutare l'utilizzo delle risorse e come modificare le allocazioni delle risorse per tutti i processi e/o per singoli processi, nonché come impostare limitazioni delle risorse per l'esecuzione su HPC.

### Prossimi passi

Impara a usare un file di parametri per memorizzare i parametri del workflow.

---

## 3. Utilizzare un file di parametri per memorizzare i parametri del workflow

Finora abbiamo esaminato la configurazione dal punto di vista tecnico dell'infrastruttura di elaborazione.
Ora prendiamo in considerazione un altro aspetto della configurazione del workflow che è molto importante per la riproducibilità: la configurazione dei parametri del workflow.

Attualmente, il nostro workflow è impostato per accettare diversi valori di parametri tramite la riga di comando, con valori predefiniti impostati nello script del workflow stesso.
Questo va bene per un semplice workflow con pochissimi parametri che devono essere impostati per una determinata esecuzione.
Tuttavia, molti workflows del mondo reale avranno molti più parametri che potrebbero essere specifici dell'esecuzione e inserirli tutti nella riga di comando sarebbe noioso e soggetto a errori.

Nextflow ci consente di specificare i parametri tramite un file di parametri in formato JSON, il che rende molto comodo gestire e distribuire set alternativi di valori predefiniti, ad esempio, nonché valori di parametri specifici dell'esecuzione.

Forniamo un file di parametri di esempio nella directory corrente, denominato `test-params.json`:

```json title="test-params.json" linenums="1"
{
  "greeting": "greetings.csv",
  "batch": "Trio",
  "character": "turkey"
}
```

Questo file di parametri contiene una coppia chiave-valore per ciascuno degli input previsti dal nostro workflow.

### 3.1. Eseguire il workflow utilizzando un file di parametri

Per eseguire il workflow con questo file di parametri, è sufficiente aggiungere `-params-file <nomefile>` al comando di base.

```bash
nextflow run hello-config.nf -params-file test-params.json
```

Funziona! E come previsto, produce gli stessi output di prima.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

executor >  local (8)
[f0/35723c] sayHello (2)       | 3 of 3 ✔
[40/3efd1a] convertToUpper (3) | 3 of 3 ✔
[17/e97d32] collectGreetings   | 1 of 1 ✔
[98/c6b57b] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

Questo potrebbe sembrare eccessivo quando hai solo pochi parametri da specificare, ma alcune pipeline si aspettano decine di parametri.
In quei casi, usare un file di parametri ci consentirà di fornire valori di parametri in fase di esecuzione senza dover digitare lunghe righe di comando e senza modificare lo script del workflow.

### Conclusione

Sai come gestire i parametri predefiniti e sovrascriverli in fase di esecuzione utilizzando un file di parametri.

### Prossimi passi

Scopri come utilizzare i profili per passare comodamente da una configurazione alternativa all'altra.

---

## 4. Determinare quale esecutore dovrebbe essere utilizzato per svolgere il lavoro

Finora abbiamo eseguito la nostra pipeline con l'esecutore locale.
Questo esegue ogni task sulla macchina su cui è in esecuzione Nextflow.
Quando Nextflow inizia, esamina le CPU e la memoria disponibili.
Se le risorse dei task pronti per l'esecuzione superano le risorse disponibili, Nextflow tratterrà gli ultimi task dall'esecuzione fino al completamento di uno o più task precedenti, liberando le risorse necessarie.

Per carichi di lavoro molto grandi, potresti scoprire che la tua macchina locale è un collo di bottiglia, sia perché hai un singolo task che richiede più risorse di quelle disponibili, sia perché hai così tanti task che aspettare che un singolo computer li esegua richiederebbe troppo tempo.
L'esecutore locale è comodo ed efficiente, ma è limitato a quel singolo computer.
Nextflow supporta [molti backend di esecuzione diversi](https://www.nextflow.io/docs/latest/executor.html), inclusi gli scheduler HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor e altri) così come i backend di esecuzione cloud come (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e altri).

Ognuno di questi sistemi utilizza tecnologie, sintassi e configurazioni diverse per definire come dovrebbe essere definito un job. Ad esempio, /se non avessimo Nextflow/, un job che richiede 8 CPU e 4 GB di RAM per essere eseguito sulla coda "my-science-work" dovrebbe includere la seguente configurazione su SLURM e inviare il job tramite `sbatch`:

```bash
#SBATCH -o /path/to/my/task/directory/my-task-1.log
#SBATCH --no-requeue
#SBATCH -c 8
#SBATCH --mem 4096M
#SBATCH -p my-science-work
```

Se volessi rendere il workflow disponibile a un collega che utilizza PBS, dovrei ricordarmi di utilizzare un programma di invio diverso, `qsub`, e dovrei modificare i miei script per utilizzare una nuova sintassi per le risorse:

```bash
#PBS -o /path/to/my/task/directory/my-task-1.log
#PBS -j oe
#PBS -q my-science-work
#PBS -l nodes=1:ppn=5
#PBS -l mem=4gb
```

Se volessi usare SGE, la configurazione sarebbe leggermente diversa:

```bash
#$ -o /path/to/my/task/directory/my-task-1.log
#$ -j y
#$ -terse
#$ -notify
#$ -q my-science-work
#$ -l slots=5
#$ -l h_rss=4096M,mem_free=4096M
```

L'esecuzione su un singolo motore di esecuzione cloud richiederebbe di nuovo un nuovo approccio, probabilmente utilizzando un SDK che utilizza le API della piattaforma cloud.

Nextflow semplifica la scrittura di un singolo flusso di lavoro che può essere eseguito su ciascuna di queste diverse infrastrutture e sistemi, senza dover modificare il workflow
L'esecutore è soggetto a una direttiva di processo denominata `executor`.
Per impostazione predefinita è impostato su `local`, quindi è implicita la seguente configurazione:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

### 4.1. Targeting di un backend diverso

Per impostazione predefinita, questo ambiente di formazione non include uno scheduler HPC in esecuzione, ma se si esegue su un sistema con SLURM installato, ad esempio, è possibile far sì che Nextflow converta `cpus`, `memory`, `queue` e altre direttive di processo nella sintassi corretta in fase di esecuzione aggiungendo le seguenti righe al file `nextflow.config`:

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

E... questo è tutto! Come detto prima, questo presuppone che Slurm stesso sia già impostato per te, ma questo è tutto ciò che Nextflow stesso deve sapere.

In pratica stiamo dicendo a Nextflow di generare uno script di invio Slurm e di inviarlo usando un comando `sbatch`.

### Conclusione

Ora sai come modificare l'esecutore per utilizzare diversi tipi di infrastrutture informatiche.

### Prossimi passi?

Scopri come controllare le risorse assegnate per l'esecuzione dei processi.

---

## 5. Utilizza i profili per selezionare configurazioni preimpostate

Potresti voler passare da un'impostazione alternativa all'altra a seconda dell'infrastruttura informatica che stai utilizzando. Ad esempio, potresti voler sviluppare ed eseguire test su piccola scala localmente sul tuo laptop, quindi eseguire carichi di lavoro su larga scala su HPC o cloud.

Nextflow ti consente di impostare profili che descrivono diverse configurazioni, che puoi quindi selezionare in fase di esecuzione utilizzando un argomento della riga di comando, anziché dover modificare il file di configurazione stesso.

### 5.1. Creare profili per passare dallo sviluppo locale all'esecuzione su HPC

Impostiamo due profili alternativi: uno per l'esecuzione di carichi su piccola scala su un computer normale, dove utilizzeremo contenitori Docker, e uno per l'esecuzione su un HPC universitario con uno scheduler Slurm, dove utilizzeremo pacchetti Conda.

Aggiungere quanto segue al file `nextflow.config`:

```groovy title="nextflow.config" linenums="3"
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

Come puoi vedere, per l'HPC universitario stiamo anche specificando le limitazioni delle risorse.

### 5.2. Esegui workflow con un profilo.

Per specificare un profilo nella nostra riga di comando Nextflow, utilizziamo l'argomento `-profile`.

Proviamo a eseguire il workflow con la configurazione `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

Questo produce ancora il seguente output:

```
 N E X T F L O W   ~  version 24.10.0

Launching `hello-config.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

executor >  local (8)
[58/da9437] sayHello (3)       | 3 of 3 ✔
[35/9cbe77] convertToUpper (2) | 3 of 3 ✔
[67/857d05] collectGreetings   | 1 of 1 ✔
[37/7b51b5] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

Come potete vedere, questo ci consente di passare da una configurazione all'altra in modo molto pratico durante l'esecuzione.

!!! warning

    Il profilo `univ_hpc` non verrà eseguito correttamente nell'ambiente di training poiché non abbiamo accesso a uno scheduler Slurm.

Se in futuro dovessimo trovare altri elementi di configurazione che si verificano sempre contemporaneamente a questi, potremmo semplicemente aggiungerli al profilo/i corrispondente/i.
Possiamo anche creare profili aggiuntivi se ci sono altri elementi di configurazione che vogliamo raggruppare.

### 5.3. Crea un profilo di prova

I profili non servono solo per la configurazione dell'infrastruttura.
Possiamo anche usarli per impostare valori predefiniti per i parametri del workflow, per rendere più facile per altri provare il workflow senza dover raccogliere autonomamente i valori di input appropriati.
Questo è inteso come alternativa all'uso di un file di parametri.

La sintassi per esprimere i valori predefiniti è la stessa di quando li scriviamo nel file del workflow stesso, tranne per il fatto che li racchiudiamo in un blocco denominato `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Se aggiungiamo un profilo di prova per il nostro workflow, il blocco `profili` diventa:

```groovy title="nextflow.config" linenums="4"
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
        params.batch = 'test-batch'
        params.character = 'turkey'
    }
}
```

Proprio come per i profili di configurazione tecnica, è possibile impostare più profili diversi specificando i parametri con qualsiasi nome arbitrario si desideri.

### 5.4. Esegui il workflow localmente con il profilo di prova

Per comodità, i profili non si escludono a vicenda, quindi possiamo specificare più profili nella nostra riga di comando utilizzando la seguente sintassi `-profile <profile1>,<profile2>` (per qualsiasi numero di profili).

!!! note

    Se si combinano profili che impostano valori per gli stessi elementi di configurazione e sono descritti nello stesso file di configurazione, Nextflow risolverà il conflitto utilizzando il valore letto per ultimo (ovvero, qualsiasi cosa si trovi dopo nel file).
    Se le impostazioni in conflitto sono impostate in diverse origini di configurazione, si applica l'[ordine di precedenza](https://www.nextflow.io/docs/latest/config.html) predefinito.

Proviamo ad aggiungere il profilo di prova al nostro comando precedente:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

Ciò dovrebbe produrre quanto segue:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-config.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

executor >  local (8)
[58/da9437] sayHello (3)       | 3 of 3 ✔
[35/9cbe77] convertToUpper (2) | 3 of 3 ✔
[67/857d05] collectGreetings   | 1 of 1 ✔
[37/7b51b5] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

<!-- migliorare mostrando e variando gli output per tutti questi forse -->

Ciò significa che finché distribuiamo file di dati di prova con il codice del workflow, chiunque può provare rapidamente il workflow senza dover fornire i propri input tramite la riga di comando o un file di parametri.

!!! note

    Possiamo anche puntare agli URL per file più grandi che sono archiviati esternamente.
    Nextflow li scaricherà automaticamente finché c'è una connessione aperta.

### Conclusione

Sai come usare i profili per selezionare una configurazione preimpostata in fase di esecuzione con il minimo sforzo. Più in generale, sai come configurare le esecuzioni del tuo workflow per adattarle a diverse piattaforme di elaborazione e migliorare la riproducibilità delle tue analisi.

### Prossimi passi

Festeggia e datti una bella pacca sulla spalla! Hai completato il tuo primo corso per sviluppatori Nextflow.

Successivamente, ti chiederemo di compilare un brevissimo sondaggio sulla tua esperienza con questo corso di formazione, dopodiché ti indirizzeremo a una pagina con link ad ulteriori risorse di formazione e link utili.
