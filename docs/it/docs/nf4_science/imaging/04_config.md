# Parte 4: Configurazione

Nelle Parti 1-3, abbiamo imparato come eseguire Nextflow, eseguire una pipeline nf-core e gestire gli input con file di parametri e samplesheet.
Ora esploreremo come configurare le pipeline per diversi ambienti di calcolo utilizzando **file di configurazione** e **profili**.

## Obiettivi di apprendimento

Al termine di questa parte, sarete in grado di:

- Comprendere come Nextflow risolve la configurazione da più fonti
- Utilizzare i profili integrati di nf-core per container e test
- Creare profili personalizzati per diversi ambienti di calcolo
- Personalizzare le richieste di risorse utilizzando le etichette dei processi
- Gestire i limiti di risorse in ambienti vincolati
- Ispezionare la configurazione risolta con `nextflow config`

---

## 1. Comprendere la configurazione di Nextflow

### 1.1. Cos'è un file di configurazione?

Nextflow utilizza file di configurazione per separare la **logica del flusso di lavoro** (cosa fare) dalle **impostazioni di esecuzione** (come e dove farlo).

I file di configurazione controllano:

- Motori di container (Docker, Singularity, Conda)
- Risorse di calcolo (CPU, memoria, tempo)
- Piattaforme di esecuzione (locale, HPC, cloud)
- Parametri della pipeline

### 1.2. Precedenza della configurazione

Nextflow carica la configurazione da più fonti, con le fonti successive che sovrascrivono quelle precedenti:

1. **Configurazione della pipeline**: `nextflow.config` nel repository della pipeline
2. **Configurazione della directory**: `nextflow.config` nella directory di lavoro corrente
3. **Configurazione utente**: `~/.nextflow/config`
4. **Riga di comando**: Parametri e opzioni passati direttamente

Questo approccio stratificato consente di mantenere le impostazioni predefinite nella pipeline, sovrascriverle con impostazioni specifiche dell'utente e apportare modifiche rapide dalla riga di comando.

### 1.3. La nostra configurazione attuale

Diamo un'occhiata alla configurazione che abbiamo utilizzato finora:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Commentiamo o modifichiamo la riga `docker.enabled = true` dalla Parte 2, e scopriamo come possiamo ottenere lo stesso risultato utilizzando un profilo in molkart.

---

## 2. Utilizzare i profili

### 2.1. Cosa sono i profili?

I profili sono insiemi denominati di configurazione che possono essere attivati con il flag `-profile` tramite il comando `nextflow run`.
Rendono facile passare da uno scenario di calcolo all'altro senza modificare i file di configurazione.

Tutte le pipeline nf-core includono una serie di profili predefiniti che possiamo utilizzare.

### 2.2. Ispezionare i profili integrati

Ispezioniamoli nel file `molkart/nextflow.config` associato al codice della pipeline:

```bash
code molkart/nextflow.config
```

Cerchiamo il blocco `profiles`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Profili container comuni:

- `docker`: Utilizza container Docker (più comune per lo sviluppo locale)
- `singularity`: Utilizza Singularity/Apptainer (comune su HPC)
- `conda`: Utilizza ambienti Conda
- `apptainer`: Utilizza container Apptainer

### 2.3. Rieseguire con i profili invece di nextflow.config

Ora che abbiamo disabilitato la configurazione docker nel nostro file `nextflow.config` locale e comprendiamo i profili, rieseguiamo la pipeline utilizzando il flag `-profile`.

Precedentemente nella Parte 3, abbiamo creato un file `params.yaml` con i nostri parametri personalizzati.
Possiamo ora combinarlo con il profilo Docker integrato:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Analizziamo cosa fa ciascun flag:

- `-profile docker`: Attiva il profilo Docker dal `nextflow.config` di molkart, che imposta `docker.enabled = true`
- `-params-file params.yaml`: Carica tutti i parametri della pipeline dal nostro file YAML
- `-resume`: Riutilizza i risultati memorizzati nella cache dalle esecuzioni precedenti

Poiché stiamo usando `-resume`, Nextflow verificherà se qualcosa è cambiato dall'ultima esecuzione.
Se i parametri, gli input e il codice sono gli stessi, tutte le attività verranno recuperate dalla cache e la pipeline si completerà quasi istantaneamente.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Notate che tutti i processi mostrano `cached: 2` o `cached: 1` - nulla è stato rieseguito!

### 2.4. Profili di test

I profili di test forniscono modi rapidi per specificare parametri di input predefiniti e file di dati per verificare che la pipeline funzioni.
Le pipeline nf-core includeranno sempre almeno due profili di test:

- `test`: Dataset piccolo con parametri veloci per test rapidi
- `test_full`: Test più completo con dati più grandi

Diamo un'occhiata più da vicino al profilo `test` in molkart che è incluso utilizzando la direttiva `includeConfig`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Questo significa che ogni volta che eseguiamo la pipeline con `-profile test`, Nextflow caricherà la configurazione da `conf/test.config`.

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Notate che questo profilo contiene gli stessi parametri che abbiamo utilizzato nel nostro file `params.yaml` in precedenza.

Potete attivare più profili separandoli con virgole.
Utilizziamo questo per testare la nostra pipeline senza bisogno del nostro file params:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Questo combina:

- `docker`: Abilita i container Docker
- `test`: Utilizza dataset e parametri di test

I profili vengono applicati da sinistra a destra, quindi i profili successivi sovrascrivono quelli precedenti se impostano gli stessi valori.

### Takeaway

Le pipeline nf-core includono profili integrati per container, test e ambienti speciali.
Potete combinare più profili per costruire la configurazione di cui avete bisogno.

### Cosa c'è dopo?

Imparate a creare i vostri profili personalizzati per diversi ambienti di calcolo.

---

## 3. Creare profili personalizzati

### 3.1. Creare profili per passare dallo sviluppo locale all'esecuzione su HPC

Creiamo profili personalizzati per due scenari:

1. Sviluppo locale con Docker
2. HPC universitario con scheduler Slurm e Singularity

Aggiungete quanto segue al vostro `nextflow.config`:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Ora potete passare facilmente da un ambiente all'altro:

```bash
# Per lo sviluppo locale
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Per HPC (quando disponibile)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note "Nota"

    Non possiamo testare il profilo HPC in questo ambiente di formazione poiché non abbiamo accesso a uno scheduler Slurm.
    Ma questo mostra come lo configurereste per un uso reale.

### 3.2. Utilizzare `nextflow config` per ispezionare la configurazione

Il comando `nextflow config` mostra la configurazione completamente risolta senza eseguire la pipeline.

Visualizzate la configurazione predefinita:

```bash
nextflow config ./molkart
```

Visualizzate la configurazione con un profilo specifico:

```bash
nextflow config -profile local_dev ./molkart
```

Questo è estremamente utile per:

- Debug dei problemi di configurazione
- Comprendere quali valori verranno effettivamente utilizzati
- Verificare come interagiscono più profili

### Takeaway

I profili personalizzati consentono di passare da un ambiente di calcolo all'altro con un singolo flag da riga di comando.
Utilizzate `nextflow config` per ispezionare la configurazione risolta prima dell'esecuzione.

### Cosa c'è dopo?

Imparate a personalizzare le richieste di risorse per singoli processi utilizzando il sistema di etichette dei processi di nf-core.

---

## 4. Personalizzare le richieste di risorse

### 4.1. Comprendere le etichette dei processi nelle pipeline nf-core

Per semplicità, le pipeline nf-core utilizzano [**etichette dei processi**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) per standardizzare l'allocazione delle risorse in tutte le pipeline.
Ogni processo è etichettato con un'etichetta come `process_low`, `process_medium` o `process_high` per descrivere rispettivamente requisiti di risorse di calcolo bassi, medi o alti.
Queste etichette vengono convertite in richieste di risorse specifiche in uno dei file di configurazione situati nella directory `conf/` della pipeline.

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Notate il moltiplicatore `task.attempt` - questo consente ai tentativi successivi di attività di richiedere più risorse, se la pipeline è impostata con `process.maxRetries > 1`.

### 4.2. Sovrascrivere le risorse per processi specifici

Per un controllo più preciso, indirizzate singoli processi per nome:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Se proviamo a eseguire questa pipeline con la sovrascrittura sopra, il processo `CELLPOSE` richiederà 16 CPU e 32 GB di memoria invece del valore predefinito definito dalla sua etichetta.
Questo causerà il fallimento della pipeline nel nostro ambiente attuale poiché non abbiamo così tanta RAM disponibile.
Impareremo come prevenire questi tipi di fallimenti nella prossima sezione.

!!! Tip "Suggerimento"

    Per trovare i nomi dei processi, guardate l'output di esecuzione della pipeline o controllate `.nextflow.log`.
    I nomi dei processi seguono il pattern `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Takeaway

Le pipeline nf-core utilizzano etichette dei processi per standardizzare l'allocazione delle risorse.
Potete sovrascrivere le risorse per etichetta (influenza più processi) o per nome (influenza un processo specifico).

### Cosa c'è dopo?

Imparate a gestire i limiti di risorse in ambienti vincolati come GitHub Codespaces.

---

## 5. Gestire le risorse in ambienti vincolati

### 5.1. Il problema dei limiti di risorse

Se provassimo a eseguire molkart con un processo che richiede 16 CPU e 32 GB di memoria (come mostrato nella sezione 4.2), fallirebbe nel nostro ambiente attuale perché non abbiamo così tante risorse disponibili.
In un ambiente cluster con nodi più grandi, tali richieste verrebbero inviate allo scheduler.

In ambienti vincolati come GitHub Codespaces, senza limiti, Nextflow rifiuterebbe di eseguire processi che superano le risorse disponibili.

### 5.2. Impostare i limiti di risorse

La direttiva `resourceLimits` limita le richieste di risorse a valori specificati:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Questo dice a Nextflow: "Se un processo richiede più di 2 CPU o 7 GB di memoria, limitalo a questi valori invece."

### 5.3. Aggiungere limiti di risorse ai profili personalizzati

Aggiornate i vostri profili personalizzati per includere limiti appropriati:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning "Attenzione"

    Impostare limiti di risorse troppo bassi può causare il fallimento dei processi o un'esecuzione lenta.
    La pipeline potrebbe dover utilizzare algoritmi meno intensivi in termini di memoria o elaborare i dati in blocchi più piccoli.

### Takeaway

Utilizzate `resourceLimits` per eseguire pipeline in ambienti con risorse limitate limitando le richieste di risorse dei processi.
Profili diversi possono avere limiti diversi appropriati per il loro ambiente.

### Cosa c'è dopo?

Avete completato la formazione principale su Nextflow per Bioimaging!

---

## Conclusione

Ora comprendete come configurare le pipeline Nextflow per diversi ambienti di calcolo.

Competenze chiave che avete appreso:

- **Precedenza della configurazione**: Come Nextflow risolve le impostazioni da più fonti
- **Profili nf-core**: Utilizzo dei profili integrati per container, test e utilità
- **Profili personalizzati**: Creazione dei vostri profili per diversi ambienti
- **Etichette dei processi**: Comprensione e sovrascrittura delle richieste di risorse per etichetta
- **Limiti di risorse**: Gestione di ambienti vincolati con `resourceLimits`
- **Ispezione della configurazione**: Utilizzo di `nextflow config` per debug e verifica delle impostazioni

Queste competenze di configurazione sono trasferibili a qualsiasi pipeline Nextflow e vi aiuteranno a eseguire flussi di lavoro in modo efficiente su macchine locali, cluster HPC e piattaforme cloud.

### Cosa c'è dopo?

Congratulazioni per aver completato il corso Nextflow per Bioimaging!

Prossimi passi:

- Compilate il questionario del corso per fornire feedback
- Date un'occhiata a [Hello Nextflow](../hello_nextflow/index.md) per saperne di più sullo sviluppo di flussi di lavoro
- Esplorate [Hello nf-core](../hello_nf-core/index.md) per approfondire gli strumenti nf-core
- Sfogliate altri corsi nelle [raccolte di formazione](../training_collections/index.md)
