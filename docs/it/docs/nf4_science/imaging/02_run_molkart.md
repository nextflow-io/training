# Parte 2: Eseguire nf-core/molkart

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella Parte 1, abbiamo eseguito un semplice workflow Hello World per comprendere le basi dell'esecuzione di Nextflow.
Ora eseguiremo una pipeline di bioimaging reale: **nf-core/molkart**.

Questa pipeline elabora dati di trascrittomica spaziale Molecular Cartography da Resolve Bioscience.
Tuttavia, i pattern di Nextflow che apprenderà qui si applicano a qualsiasi pipeline nf-core o workflow di produzione.

## 1. Comprendere le pipeline nf-core

Prima di eseguire la pipeline, comprendiamo cos'è nf-core e perché è importante per l'esecuzione dei workflow.

### 1.1. Cos'è nf-core?

[nf-core](https://nf-co.re/) è una raccolta guidata dalla comunità di pipeline Nextflow di alta qualità.
Tutte le pipeline nf-core seguono la stessa struttura e convenzioni, il che significa che una volta imparato a eseguirne una, potete eseguirle tutte.

Caratteristiche chiave delle pipeline nf-core:

- **Struttura standardizzata**: Tutte le pipeline hanno nomi di parametri e pattern di utilizzo coerenti
- **Dati di test integrati**: Ogni pipeline include profili di test per una validazione rapida
- **Documentazione completa**: Istruzioni dettagliate sull'uso e descrizioni dei parametri
- **Controllo qualità**: Report QC automatizzati utilizzando MultiQC
- **Supporto container**: Container predefiniti per la riproducibilità

!!! tip "Desidera saperne di più su nf-core?"

    Per un'introduzione approfondita allo sviluppo di pipeline nf-core, consulti il corso di formazione [Hello nf-core](../../hello_nf-core/index.md).
    Copre come creare e personalizzare pipeline nf-core da zero.

### 1.2. La pipeline molkart

![nf-core/molkart pipeline](img/molkart.png)

La pipeline [nf-core/molkart](https://nf-co.re/molkart) elabora dati di imaging di trascrittomica spaziale attraverso diverse fasi:

1. **Preprocessing delle immagini**: Riempimento del pattern a griglia e miglioramento opzionale del contrasto
2. **Segmentazione cellulare**: Opzioni di algoritmi multipli (Cellpose, Mesmer, ilastik, Stardist)
3. **Assegnazione degli spot**: Assegnare spot di trascritti a cellule segmentate
4. **Controllo qualità**: Generare report QC completi

Gli output chiave sono:

- Tabelle di conteggio cellula-per-trascritto
- Maschere di segmentazione
- Report di controllo qualità MultiQC

---

## 2. Eseguire molkart con dati di test

Prima di iniziare, cloniamo il repository molkart localmente in modo da poter ispezionare il suo codice:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Questo crea una directory `molkart/` contenente il codice sorgente completo della pipeline.

!!! note "Perché stiamo clonando localmente?"

    Tipicamente, si eseguirebbero le pipeline nf-core direttamente da GitHub usando `nextflow run nf-core/molkart -r 1.2.0`.
    Nextflow scarica automaticamente la versione richiesta della pipeline in `$HOME/.nextflow/assets/nf-core/molkart` e la esegue da lì.
    Tuttavia, per questa formazione, stiamo clonando la pipeline in una directory locale diversa in modo da poter ispezionare più facilmente il codice.

### 2.1. Comprendere i requisiti dei container

Prima di eseguire la pipeline completa, comprendiamo perché i container sono essenziali per le pipeline nf-core.

Proviamo a eseguire la pipeline utilizzando il dataset di test e i parametri dalla configurazione di test di molkart:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

Analizziamo questi parametri:

- `--input`: Percorso al samplesheet contenente i metadati del campione
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: Parametri per il riempimento del pattern a griglia
- `--clahe_pyramid_tile`: Dimensione del kernel per il miglioramento del contrasto
- `--segmentation_method`: Quale/i algoritmo/i utilizzare per la segmentazione cellulare
- `--outdir`: Dove salvare i risultati

!!! Warning "Questo comando fallirà - è intenzionale!"

    Stiamo deliberatamente eseguendo questo senza container per dimostrare perché sono necessari.

Dopo alcuni istanti, vedrà un errore come questo:

??? failure "Output del comando"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**Cosa sta succedendo?**

L'errore `command not found` (stato di uscita 127) significa che Nextflow ha tentato di eseguire `duplicate_finder.py` ma non è riuscito a trovarlo sul vostro sistema.
Questo perché:

1. La pipeline si aspetta che software bioinformatico specializzato sia installato
2. Questi strumenti (come `duplicate_finder.py`, `apply_clahe.dask.py`, ecc.) non fanno parte delle distribuzioni Linux standard
3. Senza container, Nextflow tenta di eseguire i comandi direttamente sulla vostra macchina locale

**Da dove dovrebbero provenire questi strumenti?**

Ispezioniamo uno dei moduli process per vedere come dichiara i suoi requisiti software.

Aprite il modulo di preprocessing CLAHE:

```bash
code molkart/modules/local/clahe/main.nf
```

Guardi alla riga 5 - vedrà:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Questa riga dice a Nextflow: "Per eseguire questo process, utilizzare l'immagine Docker `ghcr.io/schapirolabor/molkart-local:v0.0.4`, che contiene tutto il software richiesto."

Ogni process dichiara quale immagine container fornisce i suoi strumenti richiesti.
Tuttavia, Nextflow utilizza questi container solo se glielo indica!

**La soluzione: Abilitare Docker nella configurazione**

### 2.2. Configurare Docker e lanciare la pipeline

Per abilitare Docker, dobbiamo cambiare `docker.enabled` da `false` a `true` nel file `nextflow.config`.

Aprite il file di configurazione:

```bash
code nextflow.config
```

Cambi `docker.enabled = false` in `docker.enabled = true`:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Ora eseguite nuovamente la pipeline con lo stesso comando:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

Questa volta, Nextflow:

1. Leggerà l'impostazione `docker.enabled = true` dalla configurazione
2. Scaricherà le immagini Docker richieste (solo la prima volta)
3. Eseguirà ogni process all'interno del suo container specificato
4. Eseguirà con successo perché tutti gli strumenti sono disponibili all'interno dei container

!!! Tip "Perché i container sono importanti"

    La maggior parte delle pipeline nf-core **richiede** la containerizzazione (Docker, Singularity, Podman, ecc.) perché:

    - Utilizzano software bioinformatico specializzato non disponibile negli ambienti standard
    - I container garantiscono la riproducibilità - le stesse identiche versioni del software vengono eseguite ovunque
    - Non è necessario installare manualmente dozzine di strumenti e le loro dipendenze

    Per maggiori dettagli sui container in Nextflow, consulti [Hello Containers](../../hello_nextflow/05_hello_containers.md) dalla formazione Hello Nextflow.

### 2.3. Monitorare l'esecuzione

Durante l'esecuzione della pipeline, vedrà un output simile a questo:

??? success "Output del comando"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

Noti come questo output sia più dettagliato rispetto al nostro esempio Hello World grazie alle convenzioni nf-core che la pipeline segue:

- La pipeline mostra la sua versione e il logo
- I parametri di configurazione vengono visualizzati
- Più process vengono eseguiti in parallelo (indicato da più righe di process)
- I nomi dei process includono il percorso completo del modulo (es., `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Comprendere l'esecuzione dei process

La riga executor `executor > local (22)` Le dice:

- **executor**: Quale ambiente di calcolo viene utilizzato (`local` = la vostra macchina)
- **(22)**: Numero totale di attività lanciate

Ogni riga di process mostra:

- **Hash** (`[1a/2b3c4d]`): Identificatore della directory di lavoro (come prima)
- **Nome del process**: Percorso completo del modulo e nome del process
- **Identificatore input**: Nome del campione tra parentesi
- **Progresso**: Percentuale completata e conteggio (es., `1 of 1 ✔`)

### Takeaway

Sa come lanciare una pipeline nf-core con dati di test e interpretare il suo output di esecuzione.

### Prossimi passi

Impari dove trovare i risultati e come interpretarli.

---

## 3. Trovare ed esaminare gli output

Quando la pipeline viene completata con successo, vedrà un messaggio di completamento e un riepilogo dell'esecuzione.

### 3.1. Localizzare la directory dei risultati

Per impostazione predefinita, le pipeline nf-core scrivono gli output in una directory specificata dal parametro `outdir`, che abbiamo impostato su `results/`.

Elenchi i contenuti:

```bash
tree results/
```

Dovrebbe vedere diverse subdirectory:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Ogni subdirectory contiene output da una fase specifica della pipeline:

- **mindagap/**: Immagini riempite a griglia dal passaggio di preprocessing MindaGap
- **clahe/**: Immagini con contrasto migliorato dal preprocessing CLAHE
- **stack/**: Stack di immagini multi-canale creati per la segmentazione
- **segmentation/**: Risultati di segmentazione da diversi algoritmi (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Tabelle di conteggio cellula-per-trascritto
- **anndata/**: Oggetti AnnData contenenti matrici cellula-per-trascritto e coordinate spaziali
- **molkartqc/**: Metriche di controllo qualità per l'assegnazione degli spot
- **multiqc/**: Report completo di controllo qualità
- **pipeline_info/**: Report di esecuzione e log

### 3.2. Esaminare il report MultiQC

Il report MultiQC è un file HTML completo che aggrega le metriche di qualità da tutti i passaggi della pipeline.

Aprite il report nel file browser e poi cliccate sul pulsante "Show Preview" per vederlo renderizzato direttamente in VS Code.

Il report include:

- Statistiche generali per tutti i campioni
- Metriche di preprocessing
- Metriche di qualità della segmentazione
- Numero di cellule e spot rilevati

!!! Tip

    I report MultiQC sono tipicamente inclusi in tutte le pipeline nf-core.
    Forniscono sempre una panoramica ad alto livello dell'esecuzione della pipeline e della qualità dei dati.

### 3.3. Esaminare le tabelle cellula-per-trascritto

L'output scientifico più importante è la tabella di conteggio cellula-per-trascritto.
Questa Le dice quanti trascritti di ciascun tipo sono stati rilevati in ogni cellula.

Navighi alla directory spot2cell:

```bash
ls results/spot2cell/
```

Troverà file come:

- `cellxgene_mem_only_cellpose.csv`: Tabella cellula-per-trascritto usando la segmentazione Cellpose
- `cellxgene_mem_only_mesmer.csv`: Tabella cellula-per-trascritto usando la segmentazione Mesmer
- `cellxgene_mem_only_stardist.csv`: Tabella cellula-per-trascritto usando la segmentazione Stardist

Abbiamo eseguito solo 1 campione in questo dataset di test, ma in un esperimento reale avremmo queste tabelle per ogni campione.
Noti come Nextflow è in grado di elaborare più metodi di segmentazione in parallelo, rendendo facile confrontare i risultati.

### 3.4. Visualizzare i report di esecuzione

Nextflow genera automaticamente diversi report di esecuzione.

Verifichi la directory pipeline_info:

```bash
ls results/pipeline_info/
```

File chiave:

- **execution_report.html**: Timeline e visualizzazione dell'utilizzo delle risorse
- **execution_timeline.html**: Diagramma di Gantt dell'esecuzione dei process
- **execution_trace.txt**: Metriche dettagliate dell'esecuzione delle attività
- **pipeline_dag.html**: Grafico aciclico diretto che mostra la struttura del workflow

Aprite il report di esecuzione per vedere l'utilizzo delle risorse:

```bash
code results/pipeline_info/execution_report.html
```

Questo mostra:

- Quanto tempo ha impiegato ogni process
- Utilizzo di CPU e memoria
- Quali attività sono state memorizzate nella cache vs. eseguite

!!! Tip

    Questi report sono incredibilmente utili per ottimizzare l'allocazione delle risorse e risolvere problemi di prestazioni.

### Takeaway

Sa come localizzare gli output della pipeline, esaminare i report di controllo qualità e accedere alle metriche di esecuzione.

### Prossimi passi

Impari sulla directory di lavoro e come Nextflow gestisce i file intermedi.

---

## 4. Esplorare la directory di lavoro

Proprio come nel nostro esempio Hello World, tutto il lavoro effettivo avviene nella directory `work/`.

### 4.1. Comprendere la struttura della directory di lavoro

La directory di lavoro contiene una subdirectory per ogni attività che è stata eseguita.
Per questa pipeline con 12 attività, ci saranno 12 subdirectory di lavoro.

Elenchi la directory di lavoro:

```bash
ls -d work/*/*/ | head -5
```

Questo mostra le prime 5 directory di attività.

### 4.2. Ispezionare una directory di attività

Scelga uno degli hash dei process di segmentazione dall'output della console (es., `[3m/4n5o6p]`) e guardi all'interno:

```bash
ls -la work/3m/4n5o6p*/
```

Vedrà:

- **File .command.\***: Script di esecuzione Nextflow e log (come prima)
- **File di input staged**: Symlink ai file di input effettivi
- **File di output**: Maschere di segmentazione, risultati intermedi, ecc.

La differenza chiave rispetto a Hello World:

- Le pipeline reali preparano grandi file di input (immagini, dati di riferimento)
- I file di output possono essere piuttosto grandi (maschere di segmentazione, immagini elaborate)
- Più file di input e output per attività

!!! Tip

    Se un process fallisce, potete navigare alla vostra directory di lavoro, esaminare `.command.err` per i messaggi di errore e persino rieseguire `.command.sh` manualmente per il debug del problema.

### 4.3. Pulizia della directory di lavoro

La directory di lavoro può diventare piuttosto grande nel corso di più esecuzioni della pipeline.
Come abbiamo appreso nella Parte 1, potete utilizzare `nextflow clean` per rimuovere le directory di lavoro dalle vecchie esecuzioni.

Tuttavia, per le pipeline nf-core con grandi file intermedi, è particolarmente importante pulire regolarmente.

### Takeaway

Comprende come le pipeline nf-core organizzano le loro directory di lavoro e come ispezionare singole attività per il debugging.

### Prossimi passi

Impari sulla cache di Nextflow e come riprendere esecuzioni della pipeline fallite.

---

## 5. Riprendere un'esecuzione della pipeline

Una delle caratteristiche più potenti di Nextflow è la capacità di riprendere una pipeline dal punto di fallimento.

### 5.1. Il meccanismo di cache

Quando esegue una pipeline con `-resume`, Nextflow:

1. Controlla la cache per ogni attività
2. Se input, codice e parametri sono identici, riutilizza il risultato memorizzato nella cache
3. Riesegue solo le attività che sono cambiate o fallite

Questo è essenziale per pipeline di lunga durata dove i fallimenti potrebbero verificarsi in ritardo nell'esecuzione.

### 5.2. Provare resume con molkart

Esegua nuovamente lo stesso comando, ma aggiunga `-resume`:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

Dovrebbe vedere un output come: <!-- TODO: full output -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

Noti `cached: 2` o `cached: 1` per ogni process - nulla è stato rieseguito!

### 5.3. Quando resume è utile

Resume è particolarmente prezioso quando:

- Una pipeline fallisce a causa di limiti di risorse (memoria esaurita, limite di tempo superato)
- È necessario modificare i process a valle senza rieseguire i passaggi a monte
- La vostra connessione di rete si interrompe durante il download dei dati
- Desidera aggiungere output aggiuntivi senza rifare il calcolo

!!! Warning

    Resume funziona solo se non ha modificato i dati di input, il codice della pipeline o i parametri.
    Se modifica uno di questi, Nextflow rieseguirà correttamente le attività interessate.

### Takeaway

Sa come utilizzare `-resume` per rieseguire efficientemente le pipeline senza ripetere attività riuscite.

### Prossimi passi

Ora che potete eseguire nf-core/molkart con dati di test, siete pronti per imparare come configurarla per i vostri dataset.
