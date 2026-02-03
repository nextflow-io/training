# Parte 2: Riscrivere Hello per nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa seconda parte del corso di formazione Hello nf-core, Le mostriamo come creare una versione compatibile con nf-core della pipeline prodotta dal corso per principianti [Hello Nextflow](../hello_nextflow/index.md).

Avrà notato nella prima sezione della formazione che le pipeline nf-core seguono una struttura abbastanza elaborata con molti file accessori.
Creare tutto ciò da zero sarebbe molto tedioso, quindi la comunità nf-core ha sviluppato strumenti per farlo invece da un template, per avviare il processo.

Le mostreremo come utilizzare questi strumenti per creare uno scaffold della pipeline, quindi adattare il codice esistente della pipeline 'regolare' sullo scaffold nf-core.

Se non ha familiarità con la pipeline Hello o potrebbe aver bisogno di un ripasso, consulti [questa pagina informativa](../info/hello_pipeline.md).

---

## 1. Creare un nuovo progetto pipeline

Prima di tutto, creiamo lo scaffold per la nuova pipeline.

!!! note "Nota"

    Assicuratevi di trovarsi nella directory `hello-nf-core` nel suo terminale.

### 1.1. Eseguire lo strumento di creazione pipeline basato su template

Iniziamo creando una nuova pipeline con il comando `nf-core pipelines create`.
Questo creerà un nuovo scaffold di pipeline utilizzando il template base nf-core, personalizzato con un nome, una descrizione e un autore della pipeline.

```bash
nf-core pipelines create
```

L'esecuzione di questo comando aprirà un'interfaccia utente testuale (TUI) per la creazione della pipeline:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Questa TUI Le chiederà di fornire informazioni di base sulla vostra pipeline e Le offrirà una scelta di funzionalità da includere o escludere nello scaffold della pipeline.

- Nella schermata di benvenuto, cliccate su **Let's go!**.
- Nella schermata `Choose pipeline type`, cliccate su **Custom**.
- Inserite i dettagli della vostra pipeline come segue (sostituendo `< IL VOSTRO NOME >` con il vostro nome), quindi cliccate su **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < IL VOSTRO NOME >
```

- Nella schermata Template features, impostate `Toggle all features` su **off**, quindi **abilitate** selettivamente i seguenti. Controllate le vostre selezioni e cliccate su **Continue**.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- Nella schermata `Final details`, cliccate su **Finish**. Attendete che la pipeline venga creata, quindi cliccate su **Continue**.
- Nella schermata Create GitHub repository, cliccate su **Finish without creating a repo**. Questo mostrerà le istruzioni per creare successivamente un repository GitHub. Ignoratele e cliccate su **Close**.

Una volta chiusa la TUI, dovrebbe vedere il seguente output nella console.

??? success "Output del comando"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

Non c'è una conferma esplicita nell'output della console che la creazione della pipeline abbia funzionato, ma dovrebbe vedere una nuova directory chiamata `core-hello`.

Visualizzi i contenuti della nuova directory per vedere quanto lavoro si è risparmiato utilizzando il template.

```bash
tree core-hello
```

??? abstract "Contenuti della directory"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

Sono molti file!

Speriamo riconosca molti di essi come gli stessi che abbiamo incontrato quando abbiamo esplorato la struttura della pipeline `nf-core/demo`.
Ma non si preoccupi se si sente ancora un po' spaesato; percorreremo insieme le parti importanti nel corso di questa formazione.

!!! note "Nota"

    Una differenza importante rispetto alla pipeline `nf-core/demo` che abbiamo esaminato nella prima parte di questa formazione è che non c'è una directory `modules`.
    Questo perché non abbiamo scelto di includere nessuno dei moduli nf-core predefiniti.

### 1.2. Testare che lo scaffold sia funzionale

Che ci crediate o no, anche se non avete ancora aggiunto alcun modulo per farle svolgere un lavoro reale, lo scaffold della pipeline può effettivamente essere eseguito utilizzando il profilo test, nello stesso modo in cui abbiamo eseguito la pipeline `nf-core/demo`.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Questo Le mostra che tutto il cablaggio di base è a posto.
Quindi, dove sono gli output? Ce ne sono?

In effetti, è stata creata una nuova directory di risultati chiamata `core-hello-results` contenente i report di esecuzione standard:

```bash
tree core-hello-results
```

??? abstract "Contenuti della directory"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Può dare un'occhiata ai report per vedere cosa è stato eseguito, e la risposta è: niente del tutto!

![report timeline di esecuzione vuoto](./img/execution_timeline_empty.png)

Diamo un'occhiata a cosa c'è effettivamente nel codice.

### 1.3. Esaminare il workflow placeholder

Se guarda dentro il file `main.nf`, vedrà che importa un workflow chiamato `HELLO` da `workflows/hello`.

Questo è equivalente al workflow `workflows/demo.nf` che abbiamo incontrato nella Parte 1, e serve come workflow placeholder per il nostro workflow di interesse, con alcune funzionalità nf-core già in atto.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Rispetto a un workflow Nextflow di base come quello sviluppato in [Hello Nextflow](../hello_nextflow/index.md), noterà alcune cose nuove qui (righe evidenziate sopra):

- Il blocco workflow ha un nome
- Gli input del workflow sono dichiarati utilizzando la parola chiave `take:` e la costruzione del canale viene spostata al workflow genitore
- Il contenuto del workflow è posizionato all'interno di un blocco `main:`
- Gli output sono dichiarati utilizzando la parola chiave `emit:`

Queste sono funzionalità opzionali di Nextflow che rendono il workflow **componibile**, il che significa che può essere richiamato dall'interno di un altro workflow.

!!! note "Workflow componibili in profondità"

    La [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest esplora la composizione dei workflow in modo molto più approfondito, incluso come comporre più workflow insieme e gestire flussi di dati complessi tra di essi. Stiamo introducendo la componibilità qui perché è un requisito fondamentale dell'architettura del template nf-core, che utilizza workflow annidati per organizzare l'inizializzazione della pipeline, il workflow di analisi principale e le attività di completamento in componenti separati e riutilizzabili.

Dovremo collegare la logica pertinente dal nostro workflow di interesse in quella struttura.
Il primo passo per questo è rendere il nostro workflow originale componibile.

### Riepilogo

Ora sa come creare uno scaffold di pipeline utilizzando gli strumenti nf-core.

### Prossimi passi?

Imparare come rendere un workflow semplice componibile come preludio a renderlo compatibile con nf-core.

---

## 2. Rendere il workflow Hello Nextflow originale componibile

Ora è il momento di mettersi al lavoro per integrare il nostro workflow nello scaffold nf-core.
Come promemoria, stiamo lavorando con il workflow presentato nel nostro corso di formazione [Hello Nextflow](../hello_nextflow/index.md).

!!! tip "Suggerimento"

    Se non ha familiarità con quella pipeline o potrebbe aver bisogno di un ripasso, consulti [The Hello pipeline](../info/hello_pipeline.md).

Le forniamo una copia pulita e completamente funzionale del workflow Hello Nextflow completato nella directory `original-hello` insieme ai suoi moduli e al file CSV predefinito che si aspetta di utilizzare come input.

```bash
tree original-hello/
```

??? abstract "Contenuti della directory"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

Si senta libero di eseguirla per assicurarsi che funzioni:

```bash
nextflow run original-hello/hello.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Apriamo il file workflow `hello.nf` per ispezionare il codice, che è mostrato per intero di seguito (senza contare i processi, che sono nei moduli):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include i moduli
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // crea un canale per gli input da un file CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emette un saluto
  sayHello(greeting_ch)

  // converte il saluto in maiuscolo
  convertToUpper(sayHello.out)

  // raccoglie tutti i saluti in un file
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // genera arte ASCII dei saluti con cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Come potete vedere, questo workflow è stato scritto come un semplice workflow senza nome che può essere eseguito autonomamente.
Per renderlo eseguibile dall'interno di un workflow genitore come richiede il template nf-core, dobbiamo renderlo **componibile**.

Esaminiamo le modifiche necessarie una per una.

### 2.1. Nominare il workflow

Prima di tutto, diamo un nome al workflow così possiamo fare riferimento ad esso da un workflow genitore.

=== "Dopo"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Prima"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

Le stesse convenzioni si applicano ai nomi dei workflow come ai nomi dei moduli.

### 2.2. Sostituire la costruzione del canale con `take`

Ora, sostituisca la costruzione del canale con una semplice dichiarazione `take` che dichiara gli input attesi.

=== "Dopo"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // channel of greetings
        greeting_ch
    ```

=== "Prima"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Questo lascia i dettagli di come vengono forniti gli input al workflow genitore.

Mentre ci siamo, possiamo anche commentare la riga `params.greeting = 'greetings.csv'`

=== "Dopo"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Prima"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Nota"

    Se ha installato l'estensione del language server di Nextflow, il controllo della sintassi evidenzierà il suo codice con sottolineature rosse ondulate.
    Questo perché se inserisce una dichiarazione `take:`, deve anche avere un `main:`.

    Lo aggiungeremo nel prossimo passaggio.

### 2.3. Prefare le operazioni del workflow con la dichiarazione `main`

Successivamente, aggiunga una dichiarazione `main` prima del resto delle operazioni chiamate nel corpo del workflow.

=== "Dopo"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // emette un saluto
        sayHello(greeting_ch)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emette un saluto
        sayHello(greeting_ch)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Questo sostanzialmente dice 'questo è ciò che questo workflow _fa_'.

### 2.4. Aggiungere la dichiarazione `emit`

Infine, aggiunga una dichiarazione `emit` che dichiara quali sono gli output finali del workflow.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Questa è un'aggiunta completamente nuova al codice rispetto al workflow originale.

### 2.5. Riepilogo delle modifiche completate

Se ha effettuato tutte le modifiche come descritto, il suo workflow dovrebbe ora apparire così:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include i moduli
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // channel of greetings
    greeting_ch

    main:

    // emette un saluto
    sayHello(greeting_ch)

    // converte il saluto in maiuscolo
    convertToUpper(sayHello.out)

    // raccoglie tutti i saluti in un file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // genera arte ASCII dei saluti con cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Questo descrive tutto ciò di cui Nextflow ha bisogno TRANNE cosa alimentare nel canale di input.
Ciò sarà definito nel workflow genitore, chiamato anche workflow **entrypoint**.

### 2.6. Creare un workflow entrypoint fittizio

Prima di integrare il nostro workflow componibile nello scaffold complesso nf-core, verifichiamo che funzioni correttamente.
Possiamo creare un semplice workflow entrypoint fittizio per testare il workflow componibile in isolamento.

Crei un file vuoto chiamato `main.nf` nella stessa directory `original-hello`.

```bash
touch original-hello/main.nf
```

Copi il seguente codice nel file `main.nf`.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// import the workflow code from the hello.nf file
include { HELLO } from './hello.nf'

// declare input parameter
params.greeting = 'greetings.csv'

workflow {
  // crea un canale per gli input da un file CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // call the imported workflow on the channel of greetings
  HELLO(greeting_ch)

  // view the outputs emitted by the workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

Ci sono due osservazioni importanti da fare qui:

- La sintassi per chiamare il workflow importato è essenzialmente la stessa della sintassi per chiamare i moduli.
- Tutto ciò che è correlato al trasferimento degli input nel workflow (parametro di input e costruzione del canale) è ora dichiarato in questo workflow genitore.

!!! note "Nota"

    Nominare il file del workflow entrypoint `main.nf` è una convenzione, non un requisito.

    Se seguite questa convenzione, potete omettere di specificare il nome del file del workflow nel vostro comando `nextflow run`.
    Nextflow cercherà automaticamente un file chiamato `main.nf` nella directory di esecuzione.

    Tuttavia, potete nominare il file del workflow entrypoint in altro modo se preferite.
    In tal caso, assicuratevi di specificare il nome del file del workflow nel suo comando `nextflow run`.

### 2.7. Testare che il workflow venga eseguito

Abbiamo finalmente tutti i pezzi di cui abbiamo bisogno per verificare che il workflow componibile funzioni.
Eseguiamolo!

```bash
nextflow run ./original-hello
```

Qui vede il vantaggio di utilizzare la convenzione di denominazione `main.nf`.
Se avessimo nominato il workflow entrypoint `something_else.nf`, avremmo dovuto fare `nextflow run original-hello/something_else.nf`.

Se ha effettuato tutte le modifiche correttamente, questo dovrebbe essere eseguito fino al completamento.

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Questo significa che abbiamo aggiornato con successo il nostro workflow HELLO per essere componibile.

### Riepilogo

Sa come rendere un workflow componibile dandogli un nome e aggiungendo dichiarazioni `take`, `main` ed `emit`, e come chiamarlo da un workflow entrypoint.

### Prossimi passi?

Imparare come innestare un workflow componibile di base sullo scaffold nf-core.

---

## 3. Adattare la logica del workflow aggiornato nel workflow placeholder

Ora che abbiamo verificato che il nostro workflow componibile funziona correttamente, torniamo allo scaffold della pipeline nf-core che abbiamo creato nella sezione 1.
Vogliamo integrare il workflow componibile che abbiamo appena sviluppato nella struttura del template nf-core, quindi il risultato finale dovrebbe apparire così.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Quindi come facciamo a farlo accadere? Diamo un'occhiata al contenuto attuale del workflow `HELLO` in `core-hello/workflows/hello.nf` (lo scaffold nf-core).

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Nel complesso, questo codice fa molto poco oltre a qualche housekeeping che ha a che fare con la cattura della versione di qualsiasi strumento software che viene eseguito nella pipeline.

Dobbiamo aggiungere il codice pertinente dalla versione componibile del workflow originale che abbiamo sviluppato nella sezione 2.

Affronteremo questo nelle seguenti fasi:

1. Copiare i moduli e configurare le importazioni dei moduli
2. Lasciare la dichiarazione `take` così com'è
3. Aggiungere la logica del workflow al blocco `main`
4. Aggiornare il blocco `emit`

!!! note "Nota"

    Ignoreremo la cattura della versione per questo primo passaggio e vedremo come collegarla in una parte successiva di questa formazione.

### 3.1. Copiare i moduli e configurare le importazioni dei moduli

I quattro processi del nostro workflow Hello Nextflow sono memorizzati come moduli in `original-hello/modules/`.
Dobbiamo copiare quei moduli nella struttura del progetto nf-core (sotto `core-hello/modules/local/`) e aggiungere dichiarazioni di importazione al file del workflow nf-core.

Prima copiamo i file dei moduli da `original-hello/` a `core-hello/`:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Ora dovrebbe vedere la directory dei moduli elencata sotto `core-hello/`.

```bash
tree core-hello/modules
```

??? abstract "Contenuti della directory"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Ora configuriamo le dichiarazioni di importazione dei moduli.

Queste erano le dichiarazioni di importazione nel workflow `original-hello/hello.nf`:

```groovy title="original-hello/hello.nf" linenums="9"
// Include i moduli
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Aprite il file `core-hello/workflows/hello.nf` e trasponga quelle dichiarazioni di importazione in esso come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Altre due osservazioni interessanti qui:

- Abbiamo adattato la formattazione delle dichiarazioni di importazione per seguire la convenzione di stile nf-core.
- Abbiamo aggiornato i percorsi relativi ai moduli per riflettere che ora sono memorizzati a un livello diverso di annidamento.

### 3.2. Lasciare la dichiarazione `take` così com'è

Il progetto nf-core ha molte funzionalità pre-costruite intorno al concetto di samplesheet, che è tipicamente un file CSV contenente dati in colonne.
Poiché è essenzialmente ciò che è il nostro file `greetings.csv`, manterremo l'attuale dichiarazione `take` così com'è, e aggiorneremo semplicemente il nome del canale di input nel prossimo passaggio.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

La gestione dell'input sarà fatta a monte di questo workflow (non in questo file di codice).

### 3.3. Aggiungere la logica del workflow al blocco `main`

Ora che i nostri moduli sono disponibili per il workflow, possiamo collegare la logica del workflow nel blocco `main`.

Come promemoria, questo è il codice pertinente nel workflow originale, che non è cambiato molto quando l'abbiamo reso componibile (abbiamo solo aggiunto la riga `main:`):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // emette un saluto
    sayHello(greeting_ch)

    // converte il saluto in maiuscolo
    convertToUpper(sayHello.out)

    // raccoglie tutti i saluti in un file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // genera arte ASCII dei saluti con cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

Dobbiamo copiare il codice che viene dopo `main:` nella nuova versione del workflow.

C'è già del codice lì che ha a che fare con la cattura delle versioni degli strumenti che vengono eseguiti dal workflow. Lo lasceremo in pace per ora (ci occuperemo delle versioni degli strumenti più tardi).
Manterremo l'inizializzazione `ch_versions = channel.empty()` in alto, quindi inseriremo la nostra logica del workflow, mantenendo il codice di raccolta delle versioni alla fine.
Questo ordinamento ha senso perché in una pipeline reale, i processi emetterebbero informazioni sulla versione che verrebbero aggiunte al canale `ch_versions` mentre il workflow viene eseguito.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // emette un saluto
        sayHello(greeting_ch)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

Noterà che abbiamo anche aggiunto una riga vuota prima di `main:` per rendere il codice più leggibile.

Sembra ottimo, ma dobbiamo ancora aggiornare il nome del canale che stiamo passando al processo `sayHello()` da `greeting_ch` a `ch_samplesheet` come mostrato di seguito, per corrispondere a ciò che è scritto sotto la parola chiave `take:`.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emette un saluto (updated to use the nf-core convention for samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emette un saluto
        sayHello(greeting_ch)
    ```

Ora la logica del workflow è correttamente collegata.

### 3.4. Aggiornare il blocco `emit`

Infine, dobbiamo aggiornare il blocco `emit` per includere la dichiarazione degli output finali del workflow.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Questo conclude le modifiche che dobbiamo apportare al workflow HELLO stesso.
A questo punto, abbiamo raggiunto la struttura complessiva del codice che ci eravamo proposti di implementare.

### Riepilogo

Sa come adattare i pezzi principali di un workflow componibile in un workflow placeholder nf-core.

### Prossimi passi?

Imparare come adattare la gestione degli input nello scaffold della pipeline nf-core.

---

## 4. Adattare la gestione degli input

Ora che abbiamo integrato con successo la nostra logica del workflow nello scaffold nf-core, dobbiamo affrontare un altro pezzo critico: assicurarci che i nostri dati di input siano elaborati correttamente.
Il template nf-core viene fornito con una gestione degli input sofisticata progettata per dataset genomici complessi, quindi dobbiamo adattarla per funzionare con il nostro file `greetings.csv` più semplice.

### 4.1. Identificare dove vengono gestiti gli input

Il primo passo è capire dove viene eseguita la gestione degli input.

Potrà ricordare che quando abbiamo riscritto il workflow Hello Nextflow per renderlo componibile, abbiamo spostato la dichiarazione del parametro di input di un livello verso l'alto, nel workflow entrypoint `main.nf`.
Quindi diamo un'occhiata al workflow entrypoint `main.nf` di livello superiore che è stato creato come parte dello scaffold della pipeline:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Il progetto nf-core fa un uso intensivo di subworkflow annidati, quindi questa parte può risultare un po' confusa al primo approccio.

Ciò che conta qui è che ci sono due workflow definiti:

- `CORE_HELLO` è un wrapper sottile per l'esecuzione del workflow HELLO che abbiamo appena finito di adattare in `core-hello/workflows/hello.nf`.
- Un workflow senza nome che chiama `CORE_HELLO` così come altri due subworkflow, `PIPELINE_INITIALISATION` e `PIPELINE_COMPLETION`.

Ecco un diagramma di come si relazionano tra loro:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Importante, non possiamo trovare alcun codice che costruisce un canale di input a questo livello, solo riferimenti a un samplesheet fornito tramite il parametro `--input`.

Un po' di ricerca rivela che la gestione degli input è eseguita dal subworkflow `PIPELINE_INITIALISATION`, appropriatamente, che è importato da `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`.

Se apriamo quel file e scorriamo verso il basso, arriviamo a questo blocco di codice:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

Questa è la factory del canale che analizza il samplesheet e lo passa in una forma pronta per essere consumata dal workflow HELLO.

!!! note "Nota"

    La sintassi sopra è un po' diversa da quella che abbiamo usato in precedenza, ma fondamentalmente questo:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    è equivalente a questo:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Questo codice coinvolge alcuni passaggi di analisi e validazione che sono altamente specifici per il samplesheet di esempio incluso con il template della pipeline nf-core, che al momento della scrittura è molto specifico del dominio e non adatto per il nostro progetto di pipeline semplice.

### 4.2. Sostituire il codice del canale di input del template

La buona notizia è che le esigenze della nostra pipeline sono molto più semplici, quindi possiamo sostituire tutto ciò con il codice di costruzione del canale che abbiamo sviluppato nel workflow Hello Nextflow originale.

Come promemoria, ecco come appariva la costruzione del canale (come visto nella directory delle soluzioni):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // crea un canale per gli input da un file CSV
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Quindi dobbiamo solo collegarlo nel workflow di inizializzazione, con modifiche minori: aggiorniamo il nome del canale da `greeting_ch` a `ch_samplesheet`, e il nome del parametro da `params.greeting` a `params.input` (vedi riga evidenziata).

=== "Dopo"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Create channel from input file provided through params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Prima"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Create channel from input file provided through params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Questo completa le modifiche di cui abbiamo bisogno per far funzionare l'elaborazione degli input.

Nella sua forma attuale, questo non ci permetterà di sfruttare le capacità integrate di nf-core per la validazione dello schema, ma possiamo aggiungere ciò in seguito.
Per ora, ci stiamo concentrando nel mantenerlo il più semplice possibile per arrivare a qualcosa che possiamo eseguire con successo sui dati di test.

### 4.3. Aggiornare il profilo test

Parlando di dati e parametri di test, aggiorniamo il profilo test per questa pipeline per utilizzare il mini-samplesheet `greetings.csv` invece del samplesheet di esempio fornito nel template.

Sotto `core-hello/conf`, troviamo due profili test del template: `test.config` e `test_full.config`, che sono pensati per testare un piccolo campione di dati e uno di dimensioni complete.
Dato lo scopo della nostra pipeline, non c'è davvero un punto nell'impostare un profilo test di dimensioni complete, quindi si senta libero di ignorare o eliminare `test_full.config`.
Ci concentreremo sulla configurazione di `test.config` per essere eseguito sul nostro file `greetings.csv` con alcuni parametri predefiniti.

#### 4.3.1. Copiare il file `greetings.csv`

Prima dobbiamo copiare il file `greetings.csv` in un posto appropriato nel nostro progetto pipeline.
Tipicamente i piccoli file di test sono memorizzati nella directory `assets`, quindi copiamo il file dalla nostra directory di lavoro.

```bash
cp greetings.csv core-hello/assets/.
```

Ora il file `greetings.csv` è pronto per essere utilizzato come input di test.

#### 4.3.2. Aggiornare il file `test.config`

Ora possiamo aggiornare il file `test.config` come segue:

=== "Dopo"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Dati di input
        input  = "${projectDir}/assets/greetings.csv"

        // Other parameters
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Prima"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Dati di input
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Punti chiave:

- **Utilizzo di `${projectDir}`**: Questa è una variabile implicita di Nextflow che punta alla directory dove si trova lo script del workflow principale (la radice della pipeline). Utilizzarla garantisce che il percorso funzioni indipendentemente da dove viene eseguita la pipeline.
- **Percorsi assoluti**: Utilizzando `${projectDir}`, creiamo un percorso assoluto, che è importante per i dati di test che vengono forniti con la pipeline.
- **Posizione dei dati di test**: Le pipeline nf-core tipicamente memorizzano i dati di test nella directory `assets/` all'interno del repository della pipeline per i piccoli file di test, o fanno riferimento a dataset di test esterni per i file più grandi.

E mentre ci siamo, stringiamo i limiti di risorse predefiniti per assicurarci che questo venga eseguito su macchine molto basilari (come le VM minimali in Github Codespaces):

=== "Dopo"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Prima"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Questo completa le modifiche del codice che dobbiamo fare.

### 4.4. Eseguire la pipeline con il profilo test

È stato molto, ma possiamo finalmente provare a eseguire la pipeline!
Noti che dobbiamo aggiungere `--validate_params false` alla riga di comando perché non abbiamo ancora configurato la validazione (che arriverà più tardi).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Se ha effettuato tutte le modifiche correttamente, dovrebbe essere eseguita fino al completamento.

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Come potete vedere, questo ha prodotto il tipico riepilogo nf-core all'inizio grazie al subworkflow di inizializzazione, e le righe per ogni modulo ora mostrano i nomi completi PIPELINE:WORKFLOW:module.

### 4.5. Trovare gli output della pipeline

La domanda ora è: dove sono gli output della pipeline?
E la risposta è abbastanza interessante: ci sono ora due posti diversi dove cercare i risultati.

Come potrà ricordare da prima, la nostra prima esecuzione del workflow appena creato ha prodotto una directory chiamata `core-hello-results/` che conteneva vari report di esecuzione e metadati.

```bash
tree core-hello-results
```

??? abstract "Contenuti della directory"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

Vede che abbiamo ottenuto un altro set di report di esecuzione oltre a quelli che abbiamo ottenuto dalla prima esecuzione, quando il workflow era ancora solo un placeholder.
Questa volta vede tutte le attività che sono state eseguite come previsto.

![report timeline di esecuzione per la pipeline Hello](./img/execution_timeline_hello.png)

!!! note "Nota"

    Ancora una volta le attività non sono state eseguite in parallelo perché stiamo eseguendo su una macchina minimalista in Github Codespaces.
    Per vederle eseguire in parallelo, provi ad aumentare l'allocazione della CPU del suo codespace e i limiti di risorse nella configurazione di test.

È fantastico, ma i nostri risultati effettivi della pipeline non sono lì!

Ecco cosa è successo: non abbiamo cambiato nulla ai moduli stessi, quindi gli output gestiti dalle direttive `publishDir` a livello di modulo vanno ancora in una directory `results` come specificato nella pipeline originale.

```bash
tree results
```

??? abstract "Contenuti della directory"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, eccoli, mescolati con gli output delle esecuzioni precedenti della pipeline Hello originale.

Se vogliamo che siano organizzati ordinatamente come gli output della pipeline demo, dovremo cambiare il modo in cui impostiamo la pubblicazione degli output.
Le mostreremo come farlo più tardi in questo corso di formazione.

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

Ed eccolo! Può sembrare molto lavoro per ottenere lo stesso risultato della pipeline originale, ma ottiene tutti quei bei report generati automaticamente, e ora ha una solida base per sfruttare le funzionalità aggiuntive di nf-core, inclusa la validazione degli input e alcune interessanti capacità di gestione dei metadati che tratteremo in una sezione successiva.

---

### Riepilogo

Sa come convertire una pipeline Nextflow normale in una pipeline in stile nf-core utilizzando il template nf-core.
Come parte di ciò, ha imparato come rendere un workflow componibile e come identificare gli elementi del template nf-core che più comunemente necessitano di essere adattati quando si sviluppa una pipeline personalizzata in stile nf-core.

### Prossimi passi?

Si prenda una pausa, è stato un lavoro duro! Quando sarà pronto, passi a [Part 3: Use an nf-core module](./03_use_module.md) per imparare come sfruttare i moduli mantenuti dalla comunità dal repository nf-core/modules.
