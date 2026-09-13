# Parte 1: Eseguire una pipeline demo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa prima parte del corso di formazione Build with nf-core, mostreremo come trovare e provare una pipeline nf-core, configurare e personalizzare la sua esecuzione in base alle proprie esigenze, e capire come la validazione dell'input protegge dagli errori più comuni.

Utilizzeremo una pipeline chiamata nf-core/demo che è mantenuta dal progetto nf-core come parte del suo inventario di pipeline per scopi dimostrativi e di formazione.

Assicuratevi che la vostra directory di lavoro sia impostata su `hello-nf-core/` come indicato nella pagina [Iniziare](./00_orientation.md).

---

## 1. Trovare e recuperare la pipeline nf-core/demo

Iniziamo localizzando la pipeline nf-core/demo sul sito web del progetto [nf-co.re](https://nf-co.re), che centralizza tutte le informazioni come: documentazione generale e articoli di aiuto, documentazione per ciascuna delle pipeline, post di blog, annunci di eventi e così via.

### 1.1. Trovare la pipeline sul sito web

Nel vostro browser web, andate su [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) e digitate `demo` nella barra di ricerca.

![risultati della ricerca](./img/search-results.png)

Cliccate sul nome della pipeline, `demo`, per accedere alla pagina di documentazione della pipeline.

Ogni pipeline rilasciata ha una pagina dedicata che include le seguenti sezioni di documentazione:

- **Introduction:** Un'introduzione e panoramica della pipeline
- **Usage:** Descrizioni di come eseguire la pipeline
- **Parameters:** Parametri della pipeline raggruppati con descrizioni
- **Output:** Descrizioni ed esempi dei file di output previsti
- **Results:** File di output di esempio generati dal dataset di test completo
- **Releases & Statistics:** Cronologia delle versioni della pipeline e statistiche

Quando state considerando di adottare una nuova pipeline, dovreste leggere attentamente la documentazione della pipeline prima per comprendere cosa fa e come dovrebbe essere configurata prima di tentare di eseguirla.

Date un'occhiata ora e vedete se riuscite a scoprire:

- Quali strumenti la pipeline eseguirà (Controllate la scheda: `Introduction`)
- Quali input e parametri la pipeline accetta o richiede (Controllate la scheda: `Parameters`)
- Quali sono gli output prodotti dalla pipeline (Controllate la scheda: `Output`)

#### 1.1.1. Panoramica della pipeline

La scheda `Introduction` fornisce una panoramica della pipeline, inclusa una rappresentazione visiva (chiamata mappa della metropolitana) e un elenco di strumenti che vengono eseguiti come parte della pipeline.

![mappa della metropolitana della pipeline](./img/nf-core-demo-subway-cropped.png)

1. Read QC ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Present QC for raw reads ([MULTIQC](http://multiqc.info/))
4. Generate a lighthearted text message from a cow ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Esempio di riga di comando

La documentazione fornisce anche un file di input di esempio (discusso ulteriormente più avanti) e un esempio di riga di comando.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Noterete che il comando di esempio NON specifica un file workflow, solo il riferimento al repository della pipeline, `nf-core/demo`.

Quando invocato in questo modo, Nextflow assumerà che il codice sia organizzato in un certo modo.
Recuperiamo il codice così possiamo esaminare questa struttura.

### 1.2. Recuperare il codice della pipeline

Una volta determinato che la pipeline sembra essere adatta ai nostri scopi, proviamola.
Fortunatamente Nextflow rende facile recuperare pipeline da repository formattati correttamente senza dover scaricare nulla manualmente.

#### 1.2.1. Usare `nextflow pull`

Torniamo al terminale ed eseguiamo quanto segue:

```bash
nextflow pull nf-core/demo
```

??? success "Output del comando"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow esegue un `pull` del codice della pipeline, cioè scarica il repository completo sul vostro disco locale.

Per essere chiari, potete farlo con qualsiasi pipeline Nextflow che sia configurata appropriatamente in GitHub, non solo le pipeline nf-core.
Tuttavia nf-core è la più grande collezione open-source di pipeline Nextflow.

#### 1.2.2. Usare `nextflow list`

Potete ottenere da Nextflow un elenco di quali pipeline avete recuperato in questo modo:

```bash
nextflow list
```

??? success "Output del comando"

    ```console
    nf-core/demo
    ```

Potete provare a fare il pull di alcune altre pipeline per vedere come vengono elencate quando ne avete più di una.

#### 1.2.3. Trovare dove è stata scaricata la pipeline

Noterete che i file non sono nella vostra directory di lavoro corrente.
Per impostazione predefinita, Nextflow salva le pipeline scaricate in `$NXF_HOME/assets`.

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
      1.0.0 [t]
      1.0.1 [t]
      1.0.2 [t]
      1.1.0 [t]
    > 1.2.0 [t]
    ```

!!! info "Info"

    Il percorso completo potrebbe differire sul vostro sistema se non state utilizzando il nostro ambiente di formazione.

Nextflow mantiene intenzionalmente il codice sorgente scaricato 'fuori mano' sul principio che queste pipeline dovrebbero essere utilizzate più come librerie che come codice con cui interagire direttamente.

Internamente, Nextflow memorizza ogni pipeline scaricata come repository git in `$NXF_HOME/assets/.repos/`, e fa il checkout del codice per ogni revisione in una sottodirectory `clones/<commit>/`.
Poiché `.repos` è una directory nascosta, un semplice `tree -L 2 $NXF_HOME/assets/` apparirà vuoto.

#### 1.2.4. Creare un symlink per accedere facilmente al codice sorgente

Non esamineremo il codice in dettaglio, ma diamo una rapida occhiata per farci un'idea di come appare l'organizzazione generale.

Per rendere più facile sfogliare il codice sorgente della pipeline, create un collegamento simbolico che punta alla copia estratta della pipeline:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Questo crea una scorciatoia che vi permette di esplorare il codice con `tree -L 2 pipelines/nf-core/demo` o di aprire i file direttamente.

#### 1.2.5. Panoramica dell'organizzazione del codice

Potete usare `tree` o utilizzare l'esploratore di file per trovare e aprire la directory `nf-core/demo`.

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

Come potete vedere, c'è molto in corso là dentro, ma la maggior parte non richiede la vostra attenzione.

In breve, notiamo che al livello superiore potete trovare un file README con informazioni di riepilogo, così come file accessori che riassumono informazioni sul progetto come licenza, linee guida per i contributi, citazioni e codice di condotta.
La documentazione dettagliata della pipeline si trova nella directory `docs`.
Tutto questo contenuto viene utilizzato per generare le pagine web sul sito web nf-core in modo programmatico, quindi sono sempre aggiornate con il codice.

Per il resto, possiamo distinguere tre gruppi funzionali di file di codice:

1. Componenti del codice della pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configurazione della pipeline
3. Parametri della pipeline / input e validazione

Non esamineremo i componenti del codice della pipeline in questa parte del corso, ma toccheremo gli elementi di configurazione e validazione che potrebbero essere rilevanti per voi come utenti finali delle pipeline nf-core.

!!! tip "Suggerimento"

    Potete anche sfogliare il codice sorgente di qualsiasi pipeline nf-core su GitHub, ad esempio [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Ogni pipeline nf-core segue la stessa struttura di directory, quindi una volta che conoscete la struttura, potete trovare file di configurazione, moduli e workflow per qualsiasi pipeline nello stesso modo.

Ma per ora, passiamo all'esecuzione della pipeline!

### Takeaway

Ora sapete come trovare una pipeline tramite il sito web nf-core e recuperare una copia locale del codice sorgente.

### Cosa c'è dopo?

Imparate come provare una pipeline nf-core con il minimo sforzo.

---

## 2. Provare la pipeline con il suo profilo di test

Convenientemente, ogni pipeline nf-core viene fornita con un profilo di test.
Questo è un set minimo di impostazioni di configurazione per l'esecuzione della pipeline utilizzando un piccolo dataset di test ospitato nel repository [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
È un ottimo modo per provare rapidamente una pipeline su piccola scala.

!!! tip "Suggerimento"

    Il sistema di profili di configurazione di Nextflow vi permette di passare facilmente tra diversi motori di container o ambienti di esecuzione.
    Per maggiori dettagli, vedete [Hello Nextflow Parte 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Esaminare il profilo di test

È buona pratica verificare cosa specifica il profilo di test di una pipeline prima di eseguirla.
Il profilo `test` per `nf-core/demo` risiede nel file di configurazione `conf/test.config`.
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

    // Dati di input
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Noterete subito che il blocco di commenti in alto include un esempio di utilizzo che mostra come eseguire la pipeline con questo profilo di test.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Le uniche cose che dobbiamo fornire sono quelle mostrate tra parentesi angolari nell'esempio di comando: `<docker/singularity>` e `<OUTDIR>`.

Come promemoria, `<docker/singularity>` si riferisce alla scelta del sistema di container. Tutte le pipeline nf-core sono progettate per essere utilizzabili con container (Docker, Singularity, ecc.) per garantire la riproducibilità ed eliminare problemi di installazione del software.
Quindi dovremo specificare se vogliamo usare Docker o Singularity per testare la pipeline.

La parte `--outdir <OUTDIR>` si riferisce alla directory in cui Nextflow scriverà gli output della pipeline.
Dobbiamo fornire un nome per essa, che possiamo semplicemente inventare.
Se non esiste già, Nextflow la creerà per noi durante l'esecuzione.

Passando alla sezione dopo il blocco di commenti, il profilo di test ci mostra cosa è stato preconfigurato per il test: in particolare, il parametro `input` è già impostato per puntare a un dataset di test, quindi non dobbiamo fornire i nostri dati.
Se seguite il link all'input preconfigurato, vedrete che è un file CSV contenente identificatori di campioni e percorsi di file per diversi campioni sperimentali.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Questo è chiamato samplesheet, ed è la forma più comune di input per le pipeline nf-core.
Non preoccupatevi se non avete familiarità con i formati e i tipi di dati, non è importante per quello che segue.

Ora abbiamo tutto ciò di cui abbiamo bisogno per provare la pipeline.

### 2.2. Eseguire la pipeline

Come indicato sopra, possiamo usare il comando di test di esempio quasi così com'è; dobbiamo solo specificare quale sistema di packaging del software usare e come chiamare la directory di output.
Qui useremo Docker per il sistema di container e `demo-results`, rispettivamente.

Con questo, possiamo eseguire il comando di test:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
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
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : docker,test
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
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     [100%] 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) [100%] 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   [100%] 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          [100%] 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Se il vostro output corrisponde a quello, congratulazioni! Avete appena eseguito la vostra prima pipeline nf-core.

Noterete che c'è molto più output sulla console rispetto a quando eseguite una pipeline Nextflow di base.
C'è un'intestazione che include un riepilogo della versione della pipeline, input e output, e alcuni elementi di configurazione.

!!! info "Info"

    Il vostro output mostrerà timestamp, nomi di esecuzione e percorsi di file diversi, ma la struttura complessiva e l'esecuzione dei processi dovrebbero essere simili.

Notate la riga vicino alla cima dell'output:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Questo vi dice quale revisione della pipeline è stata utilizzata.
Poiché non abbiamo specificato una versione, Nextflow ha usato l'ultimo commit su `master`.
Per esecuzioni riproducibili, dovreste fissare una release specifica usando il flag `-r`:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile docker,test --outdir demo-results
```

Questo garantisce che lo stesso codice della pipeline venga utilizzato ogni volta, indipendentemente da nuovi commit o release.
Per questa formazione omettiamo `-r` per semplicità, ma in produzione dovreste sempre specificarlo.

Passando all'output di esecuzione, diamo un'occhiata alle righe che ci dicono quali processi sono stati eseguiti:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     [100%] 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) [100%] 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   [100%] 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          [100%] 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Questo ci dice che sono stati eseguiti quattro processi, corrispondenti ai quattro strumenti mostrati nella pagina di documentazione della pipeline sul sito web nf-core: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` e `COWPY`.

I nomi completi dei processi come mostrati qui, come `NFCORE_DEMO:DEMO:MULTIQC`, sono più lunghi di quelli che potreste aver visto nel materiale introduttivo Hello Nextflow.
Questi includono i nomi dei loro workflow padre e riflettono la modularità del codice della pipeline.
Entreremo più nel dettaglio nella Parte 2 di questo corso.

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

Per esempio, il file `execution_timeline_*` vi mostra quali processi sono stati eseguiti, in quale ordine e quanto tempo hanno impiegato per essere eseguiti:

![report della timeline di esecuzione](./img/execution_timeline.png)

!!! info "Info"

    Qui le attività non sono state eseguite in parallelo perché stiamo eseguendo su una macchina minimalista in Github Codespaces.
    Per vedere queste eseguite in parallelo, provate ad aumentare l'allocazione CPU del vostro codespace e i limiti di risorse nella configurazione di test.

Questi report sono generati automaticamente per tutte le pipeline nf-core.

### Takeaway

Sapete come eseguire una pipeline nf-core utilizzando il suo profilo di test integrato e dove trovare i suoi output.

### Cosa c'è dopo?

Imparate come configurare la pipeline per personalizzare la sua esecuzione.

---

## 3. Configurare l'esecuzione della pipeline

Come spiegato in [Hello Config](../hello_nextflow/06_hello_config.md), vogliamo poter cambiare su quali dati la nostra pipeline verrà eseguita e come verrà eseguita senza modificare il codice della pipeline stesso.
A tal fine, Nextflow supporta diversi modi per controllare la configurazione della pipeline, il che può risultare un po' opprimente.

Il progetto nf-core specifica convenzioni per organizzare gli elementi di configurazione, distinguendo due tipi di configurazione al livello superiore: **parametri della pipeline** e **configurazione** in senso stretto.

- **Parametri della pipeline** (impostati tramite il sistema `params`) includono tipicamente cose come file di input, flag di comportamento degli strumenti e parametri di analisi.
- **Configurazione** in senso stretto si riferisce alla logistica di come la pipeline viene eseguita, cioè l'executor, le allocazioni di risorse di calcolo e così via.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/params_vs_config.excalidraw.svg"
</figure>

Iniziamo affrontando i parametri della pipeline, poi esamineremo la configurazione in senso stretto.

### 3.1. Parametri della pipeline

Per tutte le pipeline nf-core, potete ottenere un elenco completo dei parametri della pipeline direttamente dalla riga di comando usando il flag `--help`, che è esso stesso un parametro della pipeline.

#### 3.1.1. Ottenere l'elenco dei parametri con `--help`

Eseguite il comando di aiuto per la pipeline demo:

```bash
nextflow run nf-core/demo --help
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
    !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Come potete vedere, l'output raggruppa i parametri in categorie (opzioni di input/output, opzioni del genoma di riferimento, ecc.) con tipi e descrizioni per ciascuno.

Questa categorizzazione è determinata da un file di schema, di cui parleremo più avanti.
Nelle pipeline Nextflow semplici, `--help` funziona solo se lo sviluppatore lo ha implementato manualmente.

!!! tip "Suggerimento"

    Usate `--help --show_hidden` per vedere i parametri aggiuntivi che sono nascosti per impostazione predefinita, come `--publish_dir_mode` o `--monochrome_logs`.

#### 3.1.2. Impostare i valori dei parametri

Come trattato in [Hello Config](../hello_nextflow/06_hello_config.md), potete impostare i valori dei parametri dalla riga di comando con `--nome_parametro` o raccogliere un insieme di parametri in un file YAML e passarlo con `-params-file`.
Entrambi gli approcci funzionano allo stesso modo con le pipeline nf-core.

Per esempio, per saltare la fase di trimming, vogliamo impostare il parametro boolean `skip_trim` a `true`.
Nella vostra directory di lavoro è fornito un file di parametri chiamato `my_params.yml` con quel valore già impostato:

```yaml title="my_params.yml"
skip_trim: true
```

Passatelo con `-params-file`:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


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
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) [100%] 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               [100%] 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      [100%] 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Il processo `SEQTK_TRIM` non appare più nell'output.

!!! warning "Avviso: limitazioni importanti sull'input dei parametri"

    **Impostare parametri boolean dalla riga di comando**

    A partire da Nextflow versione 26.04, tutti i valori forniti dalla riga di comando sono tipizzati come stringhe.
    Per un parametro boolean come `skip_trim`, passarlo come flag semplice (`--skip_trim`) o come `--skip_trim true` viene valutato come la **stringa** `"true"`, che non supera la validazione dello schema:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Per impostare un parametro boolean a un valore genuino `true`/`false`, usate un `-params-file` come mostrato sopra, oppure impostatelo in un file di configurazione.
    I parametri di tipo string, integer e file-path non sono interessati e possono ancora essere impostati direttamente dalla riga di comando.
    Questo corso usa questo schema per tutti i parametri boolean.

    **Usare file di configurazione personalizzati**

    Sebbene sia tecnicamente possibile impostare i parametri della pipeline in un file di configurazione personalizzato passato con `-c`, questo potrebbe non sovrascrivere i valori predefiniti già impostati nel `nextflow.config` della pipeline, a seconda delle regole di precedenza della configurazione di Nextflow.
    Usare `--nome_parametro` dalla riga di comando o `-params-file` è più affidabile, poiché questi hanno sempre la precedenza.

    Come regola generale: se appare nell'output di `--help`, impostatelo tramite la riga di comando o un file di parametri piuttosto che un file di configurazione.

#### 3.1.3. Validazione dei parametri

Curiosità: il comando `--help` funziona per tutte le pipeline nf-core perché il progetto nf-core richiede agli sviluppatori di definire formalmente tutti i parametri della pipeline in un file di schema JSON (`nextflow_schema.json`).
Questo schema registra il tipo, la descrizione, il valore predefinito e il raggruppamento di ciascun parametro.

Oltre a supportare l'output di `--help`, il file di schema abilita anche la validazione automatizzata al momento del lancio.
Ciò significa che Nextflow può verificare che ogni parametro passato esista e abbia ricevuto un valore appropriato (del tipo appropriato, nell'intervallo di valori consentito, ecc.).

Tratteremo questo in maggior dettaglio nella [Parte 5: Input Validation](05_input_validation.md), ma potete già vederlo in azione fornendo alla pipeline demo un input di parametri non valido.

##### 3.1.3.1. Parametri non riconosciuti

Provate a passare un parametro che non esiste:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --foobar "invalid"
```

L'output della console include un avviso:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

La pipeline continua a essere eseguita, ma l'avviso vi avvisa immediatamente che `--foobar` non è un parametro riconosciuto.
Questo è pensato per attirare la vostra attenzione su errori di battitura non bloccanti, come `--outDir` usato invece di `--outdir`, che può aiutarvi a evitare di sprecare tempo e risorse di calcolo.

##### 3.1.3.2. Valori di parametri non validi

La validazione controlla anche i **valori** dei parametri.
Il parametro `--skip_trim` è un flag boolean, quindi passare un valore stringa causa il fallimento immediato della pipeline:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]
```

La pipeline si ferma prima che vengano eseguiti i processi, risparmiandovi un'esecuzione fallita o errata.
Come indicato nella sezione 3.1.2, i parametri boolean dovrebbero essere impostati a un valore genuino `true`/`false` in un file di parametri piuttosto che passati dalla riga di comando, poiché i valori da riga di comando sono tipizzati come stringhe.

#### 3.1.4. Validazione dell'input

La stessa logica di validazione può essere utilizzata anche per verificare la validità dei file di input.
Per esempio, se una pipeline si aspetta un samplesheet come input principale dei dati (il che è il caso di molte se non della maggior parte delle pipeline nf-core), lo sviluppatore può fornire uno schema di input (distinto dallo schema dei parametri) che descrive come il file di input dovrebbe essere strutturato.

Poi, durante l'esecuzione, Nextflow può verificare che il file di input fornito sia valido.

Tratteremo anche questo in maggior dettaglio nella [Parte 5: Input Validation](05_input_validation.md), ma potete già vederlo in azione fornendo alla pipeline demo un samplesheet di input non valido.

La pipeline `nf-core/demo` si aspetta un file CSV con le colonne `sample`, `fastq_1` e `fastq_2`.
Questo è definito in un file di schema (`assets/schema_input.json`) che specifica la struttura attesa, i tipi di colonna e i vincoli.

??? abstract "File di schema per gli input"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Lo schema specifica che `sample` e `fastq_1` sono obbligatori, mentre `fastq_2` è opzionale (supportando sia dati paired-end che single-end).
I percorsi dei file vengono validati per esistenza e pattern dell'estensione.

Per dimostrarlo, nella vostra directory di lavoro è fornito un samplesheet non valido chiamato `malformed_samplesheet.csv`:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Questo samplesheet manca della colonna obbligatoria `fastq_1` e ha un percorso di file inesistente in `fastq_2`.

Eseguite la pipeline demo usando `malformed_samplesheet.csv` come input:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory
       '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces
       and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1
```

Come potete vedere, la pipeline fallisce immediatamente e riporta **tutti** gli errori di validazione in una volta sola.
nf-schema non si ferma al primo errore — raccoglie ogni problema e li elenca insieme, così potete correggere tutto in una volta sola invece di scoprire i problemi uno alla volta.

Ogni errore identifica la voce e il campo esatti che hanno causato il problema, così potete correggere il vostro samplesheet e poi rilanciare la pipeline con la certezza che non fallirà in un momento successivo quando Nextflow andrà effettivamente ad accedere al percorso del file.

Per gli sviluppatori, tutto questo è trattato in maggior dettaglio nella [Parte 5](./05_input_validation.md) di questo corso.

### 3.2. Configurazione

La configurazione in senso stretto controlla **come** viene eseguita la pipeline: allocazione delle risorse, argomenti specifici degli strumenti, dove vengono eseguiti i job e quale sistema di packaging del software utilizzare.

Le pipeline nf-core includono la configurazione predefinita in `nextflow.config` e nella directory `conf/`.
Prima di sovrascrivere qualsiasi cosa, è utile sapere dove si trovano i valori predefiniti.

Avete già visto nella sezione 2.1 che il codice sorgente della pipeline si trova in `$NXF_HOME/assets`.
Usando il symlink `pipelines` dalla sezione 1.2.4, elencate i file di configurazione per vedere cosa è disponibile:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/nfcore_config_files.excalidraw.svg"
</figure>

I file di configurazione più importanti sono:

- **`conf/base.config`**: Definisce le etichette delle risorse (`process_low`, `process_medium`, `process_high`) che assegnano CPU, memoria e tempo ai processi. Quando vedete un processo che utilizza più risorse del previsto, è qui che si trovano quei valori predefiniti.
- **`conf/modules.config`**: Imposta gli argomenti degli strumenti per processo (`ext.args`) e le impostazioni di pubblicazione dell'output (`publishDir`). Aprite questo file per vedere quali argomenti riceve ogni strumento per impostazione predefinita.
- **`conf/test.config`**: Il profilo di test che avete usato nella sezione 2.1, che limita le risorse tramite `resourceLimits` e imposta un samplesheet di test. Attivato con `-profile test`.
  Esiste anche un `conf/test_full.config` per l'esecuzione con un dataset di test di dimensioni complete, utile per il benchmarking.

Il `nextflow.config` centrale carica tutti i file sopra indicati e imposta i valori predefiniti appropriati per tutto.

Se desiderate modificare qualsiasi impostazione specificata in questi file, non modificate nessuno di essi direttamente.
Create invece il vostro file di configurazione e passatelo con `-c`.
I valori che specificate sovrascriveranno i valori predefiniti impostati in quegli altri file.

Proviamo questo in pratica.

#### 3.2.1. Personalizzare le risorse dei processi e gli argomenti degli strumenti

I moduli nf-core supportano due tipi comuni di sovrascrittura della configurazione: **allocazione delle risorse** (CPU, memoria, tempo) e **argomenti degli strumenti** tramite `ext.args`.

Molti strumenti da riga di comando hanno argomenti che non sono abbastanza comunemente usati da essere esposti come parametri della pipeline.
La convenzione `ext.args` vi permette di passare questi argomenti allo strumento sottostante tramite un file di configurazione.

Il file `custom.config` fornito nella vostra directory di lavoro dimostra entrambe le sovrascritture:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

Il primo blocco sovrascrive l'allocazione delle risorse di `FASTQC`.
Per impostazione predefinita, `FASTQC` usa l'etichetta `process_medium` da `base.config`, che alloca 6 CPU e 36 GB di memoria; qui la limitiamo a 2 CPU e 4 GB.

Il secondo blocco passa un argomento aggiuntivo a `SEQTK_TRIM` tramite `ext.args`.
Il flag `-b 5` dice a `seqtk trimfq` di tagliare 5 basi dall'inizio di ogni read in aggiunta al trimming di qualità.

Eseguite la pipeline con questa configurazione:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-custom -c custom.config
```

??? success "Output del comando"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Il flag `-c` aggiunge la vostra configurazione sopra la configurazione integrata della pipeline.

Per verificare che la sovrascrittura di `ext.args` abbia avuto effetto, trovate l'hash della directory di lavoro di `SEQTK_TRIM` dall'output dell'esecuzione (ad es. `work/17/428668...`) e controllate il file `.command.sh` al suo interno:

```bash
cat work/17/428668/.command.sh
```

??? success "Output del comando"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Dovreste vedere `-b 5` nel comando `seqtk trimfq`.

Una cosa importante da sapere su `ext.args`: se un modulo ha già un valore predefinito impostato, il vostro valore lo **sostituirà completamente** invece di aggiungersi ad esso.
Per esempio, `FASTQC` ha `ext.args = '--quiet'` impostato per impostazione predefinita in `conf/modules.config`:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Se impostate `ext.args = '--kmers 8'` per `FASTQC`, il flag `--quiet` non verrà più applicato.
Per mantenere entrambi, impostate `ext.args = '--quiet --kmers 8'`.

Dovreste sempre verificare la configurazione predefinita di un modulo prima di sovrascrivere `ext.args`.

### Takeaway

Sapete come ottenere aiuto da una pipeline nf-core, impostare i parametri e capire come vengono validati, e personalizzare la configurazione tramite file di configurazione.

### Cosa c'è dopo?

Se volete solo eseguire pipeline nf-core, avete finito!

Se volete imparare a sviluppare le vostre pipeline secondo gli standard nf-core, prendetevi una pausa e passate alla Parte 2 quando siete pronti. Imparerete a creare la vostra pipeline compatibile con nf-core usando gli strumenti basati sul template nf-core.
