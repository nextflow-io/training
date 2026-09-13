# Parte 1: Adattare all'ambiente di calcolo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Nextflow Run](../nextflow_run/index.md), hai configurato gli input, i parametri e gli output di una pipeline.
Questo corso copre l'altra metà del quadro: adattare l'esecuzione di una pipeline a qualsiasi ambiente di calcolo su cui si trovi a girare, senza modificare il codice del flusso di lavoro.

!!! example "Scenario"

    Hai sviluppato e testato la tua pipeline sul tuo laptop usando Docker.
    Ora devi passarla ad altri: un collaboratore ha solo Conda disponibile, e il cluster HPC della tua istituzione si aspetta che i job passino attraverso il proprio scheduler con i propri limiti di risorse.
    Niente di tutto ciò dovrebbe richiedere la riscrittura della pipeline stessa.

Lo stesso codice della pipeline può girare in tutti questi ambienti, perché nulla di tutto ciò è incorporato nel flusso di lavoro.
Il packaging del software, la piattaforma di esecuzione e l'allocazione delle risorse sono tutti controllati tramite configurazione, sovrapposta al codice, ed è questo ciò che copre questo corso: come adattare la stessa pipeline a un nuovo ambiente modificando la configurazione, non il codice.

---

## 1. Selezionare una tecnologia di packaging del software

In [Nextflow Run](../nextflow_run/index.md), hai visto un profilo `conda` già configurato in `nextflow.config` come alternativa a Docker.
Qui costruirai tu stesso lo stesso switch e vedrai cosa serve per rendere un processo effettivamente utilizzabile con Conda.

### 1.1. Disabilitare Docker e abilitare Conda

Imposta `docker.enabled` su `false` e aggiungi una direttiva che abilita Conda.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Questo consente a Nextflow di creare e utilizzare ambienti Conda per qualsiasi processo che abbia un pacchetto Conda specificato.
Il processo `cowpy` non ne ha ancora uno, quindi aggiungiamone uno, interamente dalla configurazione.

### 1.2. Aggiungere un pacchetto Conda tramite configurazione

Una direttiva `conda` può essere impostata nella definizione del processo stesso, nello stesso modo in cui `container` è già presente in `modules/cowpy.nf`, ma non è obbligatorio: `withName` permette di impostarla dalla configurazione, limitandola al solo processo `cowpy`.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Questo non sostituisce la direttiva `container` già presente nel codice della pipeline, ma aggiunge un'alternativa accanto ad essa, senza toccare affatto quel codice.

!!! tip "Suggerimento"

    La ricerca su [Seqera Containers](https://seqera.io/containers/) è un modo comodo per trovare l'URI del pacchetto Conda per un dato strumento, anche se non si ha intenzione di costruire un container a partire da esso.

### 1.3. Eseguire il flusso di lavoro per verificare che possa usare Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/execution-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Questo produce lo stesso output dell'esecuzione con Docker, anche se i meccanismi dietro le quinte sono diversi: Nextflow recupera il pacchetto Conda e costruisce un ambiente a partire da esso, invece di scaricare un'immagine container.

!!! info "Info"

    Costruire un nuovo ambiente Conda può richiedere un po' più di tempo rispetto al download di un container la prima volta, ma il pacchetto usato qui è piccolo quindi dovrebbe essere rapido.

Ora torna a Docker per il resto di questo corso.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Mescolare Docker e Conda"

    Poiché queste impostazioni sono limitate per processo, è possibile mescolarle: alcuni processi usano Docker, altri usano Conda, a seconda di ciò che è disponibile per ogni strumento.
    Se sia una direttiva `container` (nel codice della pipeline) che una direttiva `conda` (qui, dalla configurazione) sono impostate per lo stesso processo e entrambi i sistemi di packaging sono abilitati, Nextflow dà priorità ai container.

### Takeaway

Sai come configurare quale tecnologia di packaging del software un processo dovrebbe usare, e come passare da Docker a Conda.

### Cosa c'è dopo?

Scopri come cambiare la piattaforma di esecuzione che Nextflow usa per eseguire effettivamente le tue attività.

---

## 2. Selezionare una piattaforma di esecuzione

Ogni pipeline che hai eseguito finora ha usato l'executor locale: ogni attività viene eseguita sulla stessa macchina di Nextflow stesso.
Nextflow controlla le CPU e la memoria disponibili, e trattiene le attività finché non si liberano risorse sufficienti.

L'executor locale è comodo, ma non scala oltre una singola macchina.
Nextflow supporta [molti altri backend di esecuzione](https://nextflow.io/docs/latest/executor.html), inclusi gli scheduler HPC (Slurm, LSF, SGE, PBS e altri) e le piattaforme cloud (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e altro ancora).

### 2.1. Puntare a un backend diverso

L'executor è impostato da una direttiva di processo chiamata `executor`.
Per impostazione predefinita è `local`, quindi il seguente è implicito:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Per puntare a un backend diverso, imposta la direttiva sull'executor desiderato.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Avviso"

    L'ambiente di formazione non è connesso a un cluster HPC, quindi questo non è qualcosa che puoi eseguire qui.

### 2.2. La sintassi specifica del backend è astratta

La maggior parte delle piattaforme HPC richiede che le submission dei job specifichino le richieste di risorse, come CPU, memoria e nome della coda, usando la propria sintassi.
La stessa richiesta di 8 CPU e 4 GB di RAM su una coda chiamata `my-science-work` appare completamente diversa a seconda dello scheduler.

??? abstract "Esempi"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow astrae tutto questo: si specificano proprietà standardizzate come `cpus`, `memory` e `queue` una sola volta (vedi [direttive di processo](https://nextflow.io/docs/latest/reference/process.html#process-directives) per l'elenco completo), e Nextflow le traduce negli script specifici del backend appropriato a runtime.

### 2.3. Vedere cosa esegue effettivamente Nextflow

Quella traduzione non è solo una comodità del file di configurazione: è supportata da qualcosa di concreto che puoi ispezionare adesso, anche con l'executor locale.
In [Nextflow Run, sezione 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory), hai guardato all'interno di una directory di attività sotto `work/` e hai trovato `.command.sh`, il comando esatto eseguito da Nextflow.
Quella stessa directory contiene anche un file che non hai ancora esaminato: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Output del comando (estratto)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` è lo script reale che Nextflow consegna per l'esecuzione.
Racchiude `.command.sh` con tutto il necessario per eseguirlo effettivamente: configurazione dell'ambiente, staging di input/output e comunicazione del risultato a Nextflow.
Con l'executor `local`, Nextflow esegue semplicemente questo script sulla stessa macchina.

Questo è esattamente ciò che cambia quando si imposta un `executor` diverso.
Per uno scheduler HPC come Slurm o PBS, Nextflow genera lo stesso tipo di script wrapper, aggiunge l'intestazione specifica dello scheduler che hai visto nella [sezione 2.2](#22-backend-specific-syntax-is-abstracted-away) (tradotta dalle impostazioni di `cpus`, `memory` e `queue`), e consegna il risultato al comando di submission dello scheduler, ad esempio `sbatch` per Slurm.
Da lì, Nextflow interroga lo scheduler per lo stato del job invece di monitorare direttamente un processo locale.
I backend cloud batch funzionano in modo leggermente diverso, poiché sono guidati da chiamate API piuttosto che da un comando di submission, ma la stessa idea di fondo si applica: lo stesso script di attività viene eseguito, cambia solo come viene avviato e tracciato.

### Takeaway

Sai come cambiare l'executor per puntare a infrastrutture di calcolo diverse, che Nextflow astrae la sintassi di submission specifica del backend, e cosa accade effettivamente dietro le quinte quando un'attività viene eseguita su un backend diverso.

### Cosa c'è dopo?

Vai alla [Parte 2](./02_resources_and_retries.md), dove imparerai come profilare e allocare le risorse di calcolo, e come gestire i fallimenti delle attività con i retry.

---

## Riepilogo

In questa parte hai imparato a:

- Passare da una tecnologia di packaging del software all'altra tra Docker e Conda
- Aggiungere una direttiva `conda` alla definizione di un processo
- Cambiare la piattaforma di esecuzione con la direttiva `executor`
- Ispezionare ciò che Nextflow genera ed esegue effettivamente per un'attività, e come questo cambia tra i diversi executor
