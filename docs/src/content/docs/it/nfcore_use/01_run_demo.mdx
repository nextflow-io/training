# Parte 1: Eseguire una pipeline demo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa prima parte del corso Use nf-core, vi mostriamo come trovare una pipeline nf-core e provarla utilizzando il suo profilo di test integrato.

Utilizzeremo una pipeline chiamata nf-core/demo, mantenuta dal progetto nf-core come parte del suo inventario di pipeline per scopi dimostrativi e di formazione.

Assicuratevi che la vostra directory di lavoro sia impostata su `nfcore-use/` come indicato nella pagina [Getting started](./00_orientation.md).

---

## 1. Trovare e recuperare la pipeline nf-core/demo

Iniziamo individuando la pipeline nf-core/demo sul sito web del progetto all'indirizzo [nf-co.re](https://nf-co.re), che centralizza tutte le informazioni come: documentazione generale e articoli di supporto, documentazione per ciascuna pipeline, post del blog, annunci di eventi e molto altro.

### 1.1. Trovare la pipeline sul sito web

Nel vostro browser web, andate su [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) e digitate `demo` nella barra di ricerca.

![risultati della ricerca](./img/search-results.png)

Cliccate sul nome della pipeline, `demo`, per accedere alla pagina di documentazione della pipeline.

Ogni pipeline rilasciata ha una pagina dedicata che include le seguenti sezioni di documentazione:

- **Introduction:** Un'introduzione e una panoramica della pipeline
- **Usage:** Descrizioni di come eseguire la pipeline
- **Parameters:** Parametri della pipeline raggruppati con descrizioni
- **Output:** Descrizioni ed esempi dei file di output attesi
- **Results:** File di output di esempio generati dal dataset di test completo
- **Releases & Statistics:** Cronologia delle versioni della pipeline e statistiche

Ogni volta che prendete in considerazione l'adozione di una nuova pipeline, dovreste leggere attentamente la documentazione prima per capire cosa fa e come deve essere configurata, prima di tentare di eseguirla.

Date un'occhiata ora e cercate di scoprire:

- Quali strumenti eseguirà la pipeline (Controllate la scheda: `Introduction`)
- Quali input e parametri la pipeline accetta o richiede (Controllate la scheda: `Parameters`)
- Quali sono gli output prodotti dalla pipeline (Controllate la scheda: `Output`)

#### 1.1.1. Panoramica della pipeline

La scheda `Introduction` fornisce una panoramica della pipeline, inclusa una rappresentazione visiva (chiamata subway map) e un elenco degli strumenti eseguiti come parte della pipeline.

![subway map della pipeline](./img/nf-core-demo-subway-cropped.png)

1. Controllo qualità delle letture ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Trimming degli adattatori e della qualità ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Presentazione del controllo qualità per le letture grezze ([MULTIQC](http://multiqc.info/))
4. Generazione di un messaggio di testo spiritoso da una mucca ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Esempio di riga di comando

La documentazione fornisce anche un file di input di esempio (discusso più avanti) e un esempio di riga di comando.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Noterete che il comando di esempio NON specifica un file del flusso di lavoro, ma solo il riferimento al repository della pipeline, `nf-core/demo`.

Quando invocato in questo modo, Nextflow presuppone che il codice sia organizzato in un certo modo.
Recuperiamo il codice così da poter esaminare questa struttura.

### 1.2. Recuperare il codice della pipeline

Una volta determinato che la pipeline sembra adatta ai nostri scopi, proviamola.
Fortunatamente Nextflow rende facile recuperare pipeline da repository correttamente formattati senza dover scaricare nulla manualmente.

#### 1.2.1. Usare `nextflow pull`

Torniamo al terminale ed eseguiamo il seguente comando:

```bash
nextflow pull nf-core/demo
```

??? success "Output del comando"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow esegue un `pull` del codice della pipeline, ovvero scarica l'intero repository sul vostro disco locale.

Per essere chiari, potete farlo con qualsiasi pipeline Nextflow correttamente configurata su GitHub, non solo con le pipeline nf-core.
Tuttavia nf-core è la più grande raccolta open-source di pipeline Nextflow.

#### 1.2.2. Usare `nextflow list`

Potete chiedere a Nextflow di fornirvi un elenco delle pipeline che avete recuperato in questo modo:

```bash
nextflow list
```

??? success "Output del comando"

    ```console
    nf-core/demo
    ```

Potete provare a fare il pull di alcune altre pipeline per vedere come vengono elencate quando ne avete più di una.

#### 1.2.3. Trovare dove è stata scaricata la pipeline

Noterete che i file non si trovano nella vostra directory di lavoro corrente.
Per impostazione predefinita, Nextflow salva le pipeline scaricate tramite pull in `$NXF_HOME/assets`.

Per trovare dove si trova una pipeline specifica, chiedete direttamente a Nextflow:

```bash
nextflow info nf-core/demo
```

??? success "Output del comando"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Info"

    Il percorso completo potrebbe essere diverso sul vostro sistema se non state utilizzando il nostro ambiente di formazione.

Nextflow mantiene il codice sorgente scaricato intenzionalmente "fuori dalla vista", secondo il principio che queste pipeline dovrebbero essere utilizzate più come librerie che come codice con cui interagire direttamente.

Internamente, Nextflow memorizza ogni pipeline scaricata tramite pull come repository git in `$NXF_HOME/assets/.repos/`, e fa il checkout del codice per ogni revisione in una sottodirectory `clones/<commit>/`.
Poiché `.repos` è una directory nascosta, un semplice `tree -L 2 $NXF_HOME/assets/` apparirà vuoto.

#### 1.2.4. Creare un symlink per accedere facilmente al codice sorgente

Non esamineremo il codice in dettaglio, ma diamo una rapida occhiata per farci un'idea di come appare l'organizzazione generale.

Per rendere più facile sfogliare il codice sorgente della pipeline, create un collegamento simbolico che punta alla copia estratta della pipeline:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Questo crea un collegamento rapido così potete esplorare il codice con `tree -L 2 pipelines/nf-core/demo` o aprire i file direttamente.

#### 1.2.5. Panoramica dell'organizzazione del codice

Potete usare `tree` oppure il file explorer per trovare e aprire la directory `nf-core/demo`.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Contenuto della directory"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows

    7 directories, 12 files
    ```

Come potete vedere, c'è molto in quella directory, ma la maggior parte non è qualcosa di cui dovete preoccuparvi.

In breve, notiamo che al livello superiore si trovano un file README con informazioni di riepilogo, oltre a file accessori che riassumono le informazioni del progetto come licenza, linee guida per i contributi, citazioni e codice di condotta.
La documentazione dettagliata della pipeline si trova nella directory `docs`.
Tutti questi contenuti vengono utilizzati per generare le pagine web sul sito nf-core in modo programmatico, quindi sono sempre aggiornati con il codice.

Per il resto, possiamo distinguere tre gruppi funzionali di file di codice:

1. Componenti del codice della pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configurazione della pipeline
3. Parametri / input della pipeline, validazione

Non esamineremo i componenti del codice della pipeline in questa parte del corso, ma toccheremo elementi di configurazione e validazione che probabilmente saranno rilevanti per voi come utenti finali delle pipeline nf-core.

!!! tip "Suggerimento"

    Potete anche sfogliare il codice sorgente di qualsiasi pipeline nf-core su GitHub, ad esempio [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Ogni pipeline nf-core segue la stessa struttura di directory, quindi una volta che conoscete la struttura, potete trovare file di configurazione, moduli e flussi di lavoro per qualsiasi pipeline nello stesso modo.

Per ora, passiamo all'esecuzione della pipeline!

### Takeaway

Ora sapete come trovare una pipeline tramite il sito web nf-core e recuperarne una copia locale del codice sorgente.

### Cosa c'è dopo?

Imparate come provare una pipeline nf-core con il minimo sforzo.

---

## 2. Provare la pipeline con il suo profilo di test

Comodamente, ogni pipeline nf-core viene fornita con un profilo di test.
Si tratta di un insieme minimo di impostazioni di configurazione per eseguire la pipeline utilizzando un piccolo dataset di test ospitato nel repository [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
È un ottimo modo per provare rapidamente una pipeline su piccola scala.

!!! tip "Suggerimento"

    Il sistema di profili di configurazione di Nextflow vi permette di passare facilmente tra diversi motori di container o ambienti di esecuzione.
    Per maggiori dettagli, consultate [Hello Nextflow Parte 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Esaminare il profilo di test

È buona pratica verificare cosa specifica il profilo di test di una pipeline prima di eseguirla.
Il profilo `test` per `nf-core/demo` si trova nel file di configurazione `conf/test.config`.
Potete trovarlo localmente all'interno del codice sorgente della pipeline scaricato da `nextflow pull`, tramite il symlink `pipelines` creato nella sezione 1.2.4:

```bash
code pipelines/nf-core/demo/conf/test.config
```

Ecco il contenuto di quel file:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Noterete subito che il blocco di commento in cima include un esempio d'uso che mostra come eseguire la pipeline con questo profilo di test.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Le uniche cose che dobbiamo fornire sono quelle mostrate tra parentesi angolari nel comando di esempio: `<docker/singularity>` e `<OUTDIR>`.

Come promemoria, `<docker/singularity>` si riferisce alla scelta del sistema di container. Tutte le pipeline nf-core sono progettate per essere utilizzabili con container (Docker, Singularity, ecc.) per garantire la riproducibilità ed eliminare i problemi di installazione del software.
Quindi dovremo specificare se vogliamo usare Docker o Singularity per testare la pipeline.

La parte `--outdir <OUTDIR>` si riferisce alla directory in cui Nextflow scriverà gli output della pipeline.
Dobbiamo fornire un nome per essa, che possiamo semplicemente inventare.
Se non esiste già, Nextflow la creerà per noi al momento dell'esecuzione.

Passando alla sezione dopo il blocco di commento, il profilo di test ci mostra cosa è stato pre-configurato per il testing: in particolare, il parametro `input` è già impostato per puntare a un dataset di test, quindi non dobbiamo fornire i nostri dati.
Se seguite il link all'input pre-configurato, vedrete che è un file csv contenente identificatori di campione e percorsi di file per diversi campioni sperimentali.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Questo si chiama samplesheet ed è la forma più comune di input per le pipeline nf-core.
Non preoccupatevi se non avete familiarità con i formati e i tipi di dati, non è importante per quello che segue.

Ora abbiamo tutto il necessario per provare la pipeline.

### 2.2. Eseguire la pipeline

Come indicato sopra, possiamo usare il comando di test di esempio quasi così com'è; dobbiamo solo specificare quale sistema di packaging software utilizzare e come nominare la directory di output.
Qui useremo Docker come sistema di container e `demo-results`, rispettivamente.

Con questo, possiamo eseguire il comando di test:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Se il vostro output corrisponde a quello mostrato, Congratulazioni! Avete appena eseguito la vostra prima pipeline nf-core.

Noterete che c'è molto più output nella console rispetto a quando eseguite una pipeline Nextflow di base.
C'è un'intestazione che include un riepilogo della versione della pipeline, degli input e output, e alcuni elementi di configurazione.

!!! info "Info"

    Il vostro output mostrerà timestamp, nomi di esecuzione e percorsi di file diversi, ma la struttura generale e l'esecuzione dei processi dovrebbero essere simili.

Notate la riga vicino alla cima dell'output:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Questo vi dice quale revisione della pipeline è stata utilizzata.
Poiché non abbiamo specificato una versione, Nextflow ha utilizzato l'ultimo commit su `master`.
Per esecuzioni riproducibili, dovreste fissare una release specifica usando il flag `-r`:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

Questo garantisce che lo stesso codice della pipeline venga utilizzato ogni volta, indipendentemente da nuovi commit o release.
Per questa formazione omettiamo `-r` per semplicità, ma in produzione dovreste sempre specificarlo.

Passando all'output di esecuzione, diamo un'occhiata alle righe che ci dicono quali processi sono stati eseguiti:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Questo ci dice che sono stati eseguiti quattro processi, corrispondenti ai quattro strumenti mostrati nella pagina di documentazione della pipeline sul sito web nf-core: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` e `COWPY`.

I nomi completi dei processi come mostrati qui, ad esempio `NFCORE_DEMO:DEMO:MULTIQC`, sono più lunghi di quelli che potreste aver visto nel materiale introduttivo Hello Nextflow.
Questi includono i nomi dei loro flussi di lavoro padre e riflettono la modularità del codice della pipeline.
Se volete imparare a sviluppare pipeline in stile nf-core voi stessi, consultate il corso [Build with nf-core](../nfcore_build/index.md).

### 2.3. Esaminare gli output della pipeline

Infine, diamo un'occhiata alla directory `demo-results` prodotta dalla pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contenuto della directory"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Potrebbe sembrare molto.
Per saperne di più sugli output della pipeline `nf-core/demo`, consultate la sua [pagina di documentazione](https://nf-co.re/demo/1.2.0/docs/output/).

In questa fase, ciò che è importante osservare è che i risultati sono organizzati per modulo, e c'è inoltre una directory chiamata `pipeline_info` contenente vari report con timestamp sull'esecuzione della pipeline.

Ad esempio, il file `execution_timeline_*` vi mostra quali processi sono stati eseguiti, in quale ordine e quanto tempo hanno impiegato:

![report della timeline di esecuzione](./img/execution_timeline.png)

!!! info "Info"

    Qui le attività non sono state eseguite in parallelo perché stiamo lavorando su una macchina minimalista in Github Codespaces.
    Per vedere queste esecuzioni in parallelo, provate ad aumentare l'allocazione di CPU del vostro codespace e i limiti di risorse nella configurazione di test.

Questi report vengono generati automaticamente per tutte le pipeline nf-core.

### Takeaway

Sapete come eseguire una pipeline nf-core usando il suo profilo di test integrato e dove trovare i suoi output.

### Cosa c'è dopo?

Passate alla [Parte 2](./02_configure_execution.md), dove imparerete come configurare l'esecuzione della pipeline.

---

## Riepilogo

In questa parte avete imparato a:

- Trovare e recuperare una pipeline nf-core ed esaminarne la struttura del codice
- Eseguire una pipeline usando il suo profilo di test integrato
