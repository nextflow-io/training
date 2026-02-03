# Parte 1: Eseguire una pipeline demo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa prima parte del corso di formazione Hello nf-core, mostreremo come trovare e provare una pipeline nf-core, comprendere come è organizzato il codice e riconoscere come differisce dal codice Nextflow semplice mostrato in [Hello Nextflow](../hello_nextflow/index.md).

Utilizzeremo una pipeline chiamata nf-core/demo che è mantenuta dal progetto nf-core come parte del suo inventario di pipeline per dimostrare la struttura del codice e le operazioni degli strumenti.

Assicuratevi che la vostra directory di lavoro sia impostata su `hello-nf-core/` come indicato nella pagina [Iniziare](./00_orientation.md).

---

## 1. Trovare e recuperare la pipeline nf-core/demo

Iniziamo localizzando la pipeline nf-core/demo sul sito web del progetto [nf-co.re](https://nf-co.re), che centralizza tutte le informazioni come: documentazione generale e articoli di aiuto, documentazione per ciascuna delle pipeline, post di blog, annunci di eventi e così via.

### 1.1. Trovare la pipeline sul sito web

Nel vostro browser web, vada su [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) e digiti `demo` nella barra di ricerca.

![risultati della ricerca](./img/search-results.png)

Clicchi sul nome della pipeline, `demo`, per accedere alla pagina di documentazione della pipeline.

Ogni pipeline rilasciata ha una pagina dedicata che include le seguenti sezioni di documentazione:

- **Introduction:** Un'introduzione e panoramica della pipeline
- **Usage:** Descrizioni di come eseguire la pipeline
- **Parameters:** Parametri della pipeline raggruppati con descrizioni
- **Output:** Descrizioni ed esempi dei file di output previsti
- **Results:** File di output di esempio generati dal dataset di test completo
- **Releases & Statistics:** Cronologia delle versioni della pipeline e statistiche

Quando sta considerando di adottare una nuova pipeline, dovrebbe leggere attentamente la documentazione della pipeline prima per comprendere cosa fa e come dovrebbe essere configurata prima di tentare di eseguirla.

Dia un'occhiata ora e veda se riesce a scoprire:

- Quali strumenti la pipeline eseguirà (Controlli la scheda: `Introduction`)
- Quali input e parametri la pipeline accetta o richiede (Controlli la scheda: `Parameters`)
- Quali sono gli output prodotti dalla pipeline (Controlli la scheda: `Output`)

#### 1.1.1. Panoramica della pipeline

La scheda `Introduction` fornisce una panoramica della pipeline, inclusa una rappresentazione visiva (chiamata mappa della metropolitana) e un elenco di strumenti che vengono eseguiti come parte della pipeline.

![mappa della metropolitana della pipeline](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Esempio di riga di comando

La documentazione fornisce anche un file di input di esempio (discusso ulteriormente più avanti) e un esempio di riga di comando.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Noterà che il comando di esempio NON specifica un file workflow, solo il riferimento al repository della pipeline, `nf-core/demo`.

Quando invocato in questo modo, Nextflow assumerà che il codice sia organizzato in un certo modo.
Recuperiamo il codice così possiamo esaminare questa struttura.

### 1.2. Recuperare il codice della pipeline

Una volta determinato che la pipeline sembra essere adatta ai nostri scopi, proviamola.
Fortunatamente Nextflow rende facile recuperare pipeline da repository formattati correttamente senza dover scaricare nulla manualmente.

Torniamo al terminale ed eseguiamo quanto segue:

```bash
nextflow pull nf-core/demo
```

??? success "Output del comando"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow esegue un `pull` del codice della pipeline, cioè scarica il repository completo sul vostro disco locale.

Per essere chiari, potete farlo con qualsiasi pipeline Nextflow che sia configurata appropriatamente in GitHub, non solo le pipeline nf-core.
Tuttavia nf-core è la più grande collezione open-source di pipeline Nextflow.

Potete ottenere da Nextflow un elenco di quali pipeline avete recuperato in questo modo:

```bash
nextflow list
```

??? success "Output del comando"

    ```console
    nf-core/demo
    ```

Noterete che i file non sono nella vostra directory di lavoro corrente.
Per impostazione predefinita, Nextflow li salva in `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Contenuto della directory"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    Il percorso completo potrebbe differire sul vostro sistema se non state utilizzando il nostro ambiente di formazione.

Nextflow mantiene intenzionalmente il codice sorgente scaricato 'fuori mano' sul principio che queste pipeline dovrebbero essere utilizzate più come librerie che come codice con cui interagire direttamente.

Tuttavia, per gli scopi di questa formazione, vogliamo essere in grado di esplorare e vedere cosa c'è dentro.
Quindi per rendere ciò più facile, creiamo un collegamento simbolico a quella posizione dalla nostra directory di lavoro corrente.

```bash
ln -s $NXF_HOME/assets pipelines
```

Questo crea una scorciatoia che rende più facile esplorare il codice appena scaricato.

```bash
tree -L 2 pipelines
```

```console title="Contenuto della directory"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Ora possiamo più facilmente sbirciare nel codice sorgente secondo necessità.

Ma prima, proviamo ad eseguire la nostra prima pipeline nf-core!

### Takeaway

Ora sapete come trovare una pipeline tramite il sito web nf-core e recuperare una copia locale del codice sorgente.

### Prossimi passi

Imparate come provare una pipeline nf-core con il minimo sforzo.

---

## 2. Provare la pipeline con il suo profilo di test

Convenientemente, ogni pipeline nf-core viene fornita con un profilo di test.
Questo è un set minimo di impostazioni di configurazione per l'esecuzione della pipeline utilizzando un piccolo dataset di test ospitato nel repository [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
È un ottimo modo per provare rapidamente una pipeline su piccola scala.

!!! note

    Il sistema di profili di configurazione di Nextflow Le permette di passare facilmente tra diversi motori di container o ambienti di esecuzione.
    Per maggiori dettagli, veda [Hello Nextflow Parte 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Esaminare il profilo di test

È buona pratica verificare cosa specifica il profilo di test di una pipeline prima di eseguirla.
Il profilo `test` per `nf-core/demo` risiede nel file di configurazione `conf/test.config` ed è mostrato di seguito.

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
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Dati di input
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Noterà subito che il blocco di commenti in alto include un esempio di utilizzo che mostra come eseguire la pipeline con questo profilo di test.

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
Se segue il link all'input preconfigurato, vedrà che è un file csv contenente identificatori di campioni e percorsi di file per diversi campioni sperimentali.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Questo è chiamato samplesheet, ed è la forma più comune di input per le pipeline nf-core.

!!! note

    Non si preoccupi se non ha familiarità con i formati e i tipi di dati, non è importante per quello che segue.

Quindi questo conferma che abbiamo tutto ciò di cui abbiamo bisogno per provare la pipeline.

### 2.2. Eseguire la pipeline

Decidiamo di usare Docker per il sistema di container e `demo-results` come directory di output, e siamo pronti per eseguire il comando di test:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Se il vostro output corrisponde a quello, congratulazioni! Avete appena eseguito la vostra prima pipeline nf-core.

Noterete che c'è molto più output sulla console rispetto a quando eseguite una pipeline Nextflow di base.
C'è un'intestazione che include un riepilogo della versione della pipeline, input e output, e alcuni elementi di configurazione.

!!! note

    Il vostro output mostrerà timestamp, nomi di esecuzione e percorsi di file diversi, ma la struttura complessiva e l'esecuzione dei processi dovrebbero essere simili.

Passando all'output di esecuzione, diamo un'occhiata alle righe che ci dicono quali processi sono stati eseguiti:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Questo ci dice che sono stati eseguiti tre processi, corrispondenti ai tre strumenti mostrati nella pagina di documentazione della pipeline sul sito web nf-core: FASTQC, SEQTK_TRIM e MULTIQC.

I nomi completi dei processi come mostrati qui, come `NFCORE_DEMO:DEMO:MULTIQC`, sono più lunghi di quelli che potreste aver visto nel materiale introduttivo Hello Nextflow.
Questi includono i nomi dei loro workflow padre e riflettono la modularità del codice della pipeline.
Entreremo più nel dettaglio tra poco.

### 2.3. Esaminare gli output della pipeline

Infine, diamo un'occhiata alla directory `demo-results` prodotta dalla pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contenuto della directory"

    ```console
    demo-results
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
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

Potrebbe sembrare molto.
Per saperne di più sugli output della pipeline `nf-core/demo`, consulti la sua [pagina di documentazione](https://nf-co.re/demo/1.0.2/docs/output/).

In questa fase, ciò che è importante osservare è che i risultati sono organizzati per modulo, e c'è inoltre una directory chiamata `pipeline_info` contenente vari report con timestamp sull'esecuzione della pipeline.

Per esempio, il file `execution_timeline_*` vi mostra quali processi sono stati eseguiti, in quale ordine e quanto tempo hanno impiegato per essere eseguiti:

![report della timeline di esecuzione](./img/execution_timeline.png)

!!! note

    Qui le attività non sono state eseguite in parallelo perché stiamo eseguendo su una macchina minimalista in Github Codespaces.
    Per vedere queste eseguite in parallelo, provate ad aumentare l'allocazione CPU del vostro codespace e i limiti di risorse nella configurazione di test.

Questi report sono generati automaticamente per tutte le pipeline nf-core.

### Takeaway

Sapete come eseguire una pipeline nf-core utilizzando il suo profilo di test integrato e dove trovare i suoi output.

### Prossimi passi

Imparate come è organizzato il codice della pipeline.

---

## 3. Esaminare la struttura del codice della pipeline

Ora che abbiamo eseguito con successo la pipeline come utenti, spostiamo la nostra prospettiva per guardare come le pipeline nf-core sono strutturate internamente.

Il progetto nf-core applica linee guida rigorose su come le pipeline sono strutturate e su come il codice è organizzato, configurato e documentato.
Comprendere come tutto questo è organizzato è il primo passo verso lo sviluppo delle proprie pipeline compatibili con nf-core, che affronteremo nella Parte 2 di questo corso.

Diamo un'occhiata a come il codice della pipeline è organizzato nel repository `nf-core/demo`, utilizzando il symlink `pipelines` che abbiamo creato in precedenza.

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
    ```

C'è molto in corso là dentro, quindi affronteremo questo passo per passo.

Prima, notiamo che al livello superiore, potete trovare un file README con informazioni di riepilogo, così come file accessori che riassumono informazioni sul progetto come licenza, linee guida per i contributi, citazioni e codice di condotta.
La documentazione dettagliata della pipeline si trova nella directory `docs`.
Tutto questo contenuto viene utilizzato per generare le pagine web sul sito web nf-core in modo programmatico, quindi sono sempre aggiornate con il codice.

Ora, per il resto, divideremo la nostra esplorazione in tre fasi:

1. Componenti del codice della pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configurazione della pipeline
3. Input e validazione

Iniziamo con i componenti del codice della pipeline.
Ci concentreremo sulla gerarchia dei file e sull'organizzazione strutturale, piuttosto che immergerci nel codice all'interno dei singoli file.

### 3.1. Componenti del codice della pipeline

L'organizzazione standard del codice della pipeline nf-core segue una struttura modulare progettata per massimizzare il riutilizzo del codice, come introdotto in [Hello Modules](../hello_nextflow/04_hello_modules.md), Parte 4 del corso [Hello Nextflow](../hello_nextflow/index.md), sebbene nel vero stile nf-core, questo sia implementato con un po' di complessità aggiuntiva.
Specificamente, le pipeline nf-core fanno uso abbondante di subworkflow, cioè script di workflow che sono importati da un workflow padre.

Questo potrebbe sembrare un po' astratto, quindi diamo un'occhiata a come viene utilizzato nella pratica nella pipeline `nf-core/demo`.

!!! note

    Non esamineremo il codice effettivo per _come_ questi componenti modulari sono connessi, perché c'è una certa complessità aggiuntiva associata all'uso dei subworkflow che può risultare confusa, e comprendere questo non è necessario in questa fase della formazione.
    Per ora, ci concentreremo sull'organizzazione generale e sulla logica.

#### 3.1.1. Panoramica generale

Ecco come appaiono le relazioni tra i componenti di codice rilevanti per la pipeline `nf-core/demo`:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

C'è un cosiddetto script _entrypoint_ chiamato `main.nf`, che funge da wrapper per due tipi di workflow nidificati: il workflow contenente la logica di analisi effettiva, situato sotto `workflows/` e chiamato `demo.nf`, e un set di workflow di gestione situati sotto `subworkflows/`.
Il workflow `demo.nf` richiama **moduli** situati sotto `modules/`; questi contengono i **processi** che eseguiranno le effettive fasi di analisi.

!!! note

    I subworkflow non sono limitati a funzioni di gestione e possono utilizzare moduli di processo.

    La pipeline `nf-core/demo` mostrata qui si trova sul lato più semplice dello spettro, ma altre pipeline nf-core (come `nf-core/rnaseq`) utilizzano subworkflow che sono coinvolti nell'analisi effettiva.

Ora, esaminiamo questi componenti a turno.

#### 3.1.2. Lo script entrypoint: `main.nf`

Lo script `main.nf` è l'entrypoint da cui parte Nextflow quando eseguiamo `nextflow run nf-core/demo`.
Ciò significa che quando esegue `nextflow run nf-core/demo` per eseguire la pipeline, Nextflow trova ed esegue automaticamente lo script `main.nf`.
Questo funziona per qualsiasi pipeline Nextflow che segua questa convenzione di denominazione e struttura, non solo le pipeline nf-core.

L'uso di uno script entrypoint rende facile eseguire subworkflow di 'gestione' standardizzati prima e dopo l'esecuzione dello script di analisi effettivo.
Esamineremo questi dopo aver esaminato il workflow di analisi effettivo e i suoi moduli.

#### 3.1.3. Lo script di analisi: `workflows/demo.nf`

Il workflow `workflows/demo.nf` è dove è memorizzata la logica centrale della pipeline.
È strutturato molto come un normale workflow Nextflow, tranne che è progettato per essere chiamato da un workflow padre, il che richiede alcune funzionalità extra.
Copriremo le differenze rilevanti nella prossima parte di questo corso, quando affronteremo la conversione della semplice pipeline Hello da Hello Nextflow in una forma compatibile con nf-core.

Il workflow `demo.nf` richiama **moduli** situati sotto `modules/`, che esamineremo successivamente.

!!! note

    Alcuni workflow di analisi nf-core mostrano livelli aggiuntivi di nidificazione richiamando subworkflow di livello inferiore.
    Questo è usato principalmente per incapsulare due o più moduli che sono comunemente usati insieme in segmenti di pipeline facilmente riutilizzabili.
    Può vedere alcuni esempi esplorando i [subworkflow nf-core](https://nf-co.re/subworkflows/) disponibili sul sito web nf-core.

    Quando lo script di analisi usa subworkflow, questi sono memorizzati sotto la directory `subworkflows/`.

#### 3.1.4. I moduli

I moduli sono dove risiede il codice del processo, come descritto nella [Parte 4 del corso di formazione Hello Nextflow](../hello_nextflow/04_hello_modules.md).

Nel progetto nf-core, i moduli sono organizzati utilizzando una struttura nidificata multi-livello che riflette sia la loro origine che il loro contenuto.
Al livello superiore, i moduli sono differenziati come `nf-core` o `local` (non parte del progetto nf-core), e poi ulteriormente posizionati in una directory denominata in base allo strumento/agli strumenti che incapsulano.
Se lo strumento appartiene a un toolkit (cioè un pacchetto contenente più strumenti) allora c'è un livello di directory intermedio denominato in base al toolkit.

Può vedere questo applicato nella pratica ai moduli della pipeline `nf-core/demo`:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Contenuto della directory"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Qui vede che i moduli `fastqc` e `multiqc` si trovano al livello superiore all'interno dei moduli `nf-core`, mentre il modulo `trim` si trova sotto il toolkit a cui appartiene, `seqtk`.
In questo caso non ci sono moduli `local`.

Il file del codice del modulo che descrive il processo è sempre chiamato `main.nf`, ed è accompagnato da test e file `.yml` che ignoreremo per ora.

Presi insieme, il workflow entrypoint, il workflow di analisi e i moduli sono sufficienti per eseguire le parti 'interessanti' della pipeline.
Tuttavia, sappiamo che ci sono anche subworkflow di gestione lì dentro, quindi guardiamoli ora.

#### 3.1.5. I subworkflow di gestione

Come i moduli, i subworkflow sono differenziati in directory `local` e `nf-core`, e ogni subworkflow ha la propria struttura di directory nidificata con il proprio script `main.nf`, test e file `.yml`.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Contenuto della directory"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Come notato sopra, la pipeline `nf-core/demo` non include subworkflow specifici per l'analisi, quindi tutti i subworkflow che vediamo qui sono cosiddetti workflow di 'gestione' o 'utility', come denotato dal prefisso `utils_` nei loro nomi.
Questi subworkflow sono ciò che produce l'intestazione nf-core elegante nell'output della console, tra altre funzioni accessorie.

!!! tip

    Oltre al loro pattern di denominazione, un'altra indicazione che questi subworkflow non eseguono alcuna funzione realmente correlata all'analisi è che non richiamano alcun processo.

Questo completa il riepilogo dei componenti di codice centrali che costituiscono la pipeline `nf-core/demo`.
Ora diamo un'occhiata agli elementi rimanenti che dovrebbe conoscere un po' prima di immergersi nello sviluppo: configurazione della pipeline e validazione dell'input.

### 3.2. Configurazione della pipeline

Ha appreso in precedenza che Nextflow offre molte opzioni per configurare l'esecuzione della pipeline, sia in termini di input e parametri, risorse di calcolo e altri aspetti dell'orchestrazione.
Il progetto nf-core applica linee guida altamente standardizzate per la configurazione della pipeline che mirano a costruire sulle opzioni di personalizzazione flessibili di Nextflow in un modo che fornisca maggiore coerenza e manutenibilità tra le pipeline.

Il file di configurazione centrale `nextflow.config` è utilizzato per impostare i valori predefiniti per i parametri e altre opzioni di configurazione.
La maggior parte di queste opzioni di configurazione sono applicate per impostazione predefinita mentre altre (ad esempio, profili di dipendenze software) sono incluse come profili opzionali.

Ci sono diversi file di configurazione aggiuntivi che sono memorizzati nella cartella `conf` e che possono essere aggiunti alla configurazione per impostazione predefinita o opzionalmente come profili:

- `base.config`: Un file di configurazione 'da zero', appropriato per l'uso generale nella maggior parte degli ambienti di calcolo ad alte prestazioni. Questo definisce ampie categorie di utilizzo delle risorse, per esempio, che sono convenienti da applicare ai moduli.
- `modules.config`: Direttive e argomenti aggiuntivi del modulo.
- `test.config`: Un profilo per eseguire la pipeline con dati di test minimi, che abbiamo usato quando abbiamo eseguito la pipeline demo.
- `test_full.config`: Un profilo per eseguire la pipeline con un dataset di test completo.

Toccheremo alcuni di questi file più avanti nel corso.

### 3.3. Input e validazione

Come abbiamo notato in precedenza, quando abbiamo esaminato il profilo di test della pipeline `nf-core/demo`, è progettata per prendere come input un samplesheet contenente percorsi di file e identificatori di campioni.
I percorsi di file collegati a dati reali situati nel repository `nf-core/test-datasets`.

Un samplesheet di esempio è fornito anche sotto la directory `assets`, sebbene i percorsi in questo non siano reali.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Questo particolare samplesheet è abbastanza semplice, ma alcune pipeline vengono eseguite su samplesheet più complessi, con molti più metadati associati agli input primari.

Sfortunatamente, poiché questi file possono essere difficili da verificare a occhio, la formattazione impropria dei dati di input è una fonte molto comune di fallimenti della pipeline.
Un problema correlato è quando i parametri sono forniti in modo errato.

La soluzione a questi problemi è eseguire controlli di validazione automatizzati su tutti i file di input per garantire che contengano i tipi di informazioni previsti, formattati correttamente, e sui parametri per garantire che siano del tipo previsto.
Questo è chiamato validazione dell'input, e dovrebbe idealmente essere fatto _prima_ di provare a eseguire una pipeline, piuttosto che aspettare che la pipeline fallisca per scoprire che c'era un problema con gli input.

Proprio come per la configurazione, il progetto nf-core ha opinioni molto forti sulla validazione dell'input, e raccomanda l'uso del [plugin nf-schema](https://nextflow-io.github.io/nf-schema/latest/), un plugin Nextflow che fornisce capacità di validazione complete per le pipeline Nextflow.

Copriremo questo argomento in maggior dettaglio nella Parte 5 di questo corso.
Per ora, sia consapevole che ci sono due file JSON forniti per quello scopo, `nextflow_schema.json` e `assets/schema_input.json`.

Il `nextflow_schema.json` è un file utilizzato per memorizzare informazioni sui parametri della pipeline inclusi tipo, descrizione e testo di aiuto in un formato leggibile dalle macchine.
Questo è utilizzato per vari scopi, inclusa la validazione automatizzata dei parametri, la generazione di testo di aiuto e il rendering di form di parametri interattivi nelle interfacce UI.

Il `schema_input.json` è un file utilizzato per definire la struttura del samplesheet di input.
Ogni colonna può avere un tipo, pattern, descrizione e testo di aiuto in un formato leggibile dalle macchine.

### Takeaway

Sapete quali sono i componenti principali di una pipeline nf-core e come il codice è organizzato; dove si trovano gli elementi principali di configurazione; e siete consapevoli di cosa serve la validazione dell'input.

### Prossimi passi

Si prenda una pausa! È stato molto. Quando si sente riposato e pronto, passi alla sezione successiva per applicare ciò che avete appreso per scrivere una pipeline compatibile con nf-core.

!!! tip

    Se desiderate imparare come comporre workflow con subworkflow prima di passare alla parte successiva, consulti la [Side Quest Workflows of Workflows](../side_quests/workflows_of_workflows.md).
