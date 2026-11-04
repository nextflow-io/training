# Parte 2: Configurare l'esecuzione della pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella [Parte 1](./01_run_demo.md), hai trovato ed eseguito la pipeline nf-core/demo usando il suo profilo di test.
Ora vediamo come configurare l'esecuzione della pipeline: impostare i parametri, capire la validazione e personalizzare l'allocazione delle risorse e gli argomenti degli strumenti.

Come spiegato in [Hello Config](../hello_nextflow/06_hello_config.md), vogliamo poter cambiare su quali dati la nostra pipeline verrà eseguita e come verrà eseguita senza modificare il codice della pipeline stessa.
A tal fine, Nextflow supporta diversi modi per controllare la configurazione della pipeline, il che può risultare un po' disorientante.

Il progetto nf-core specifica delle convenzioni per organizzare gli elementi di configurazione, distinguendo due tipi di configurazione al livello superiore: **parametri della pipeline** e **configurazione** in senso stretto.

- I **parametri della pipeline** (impostati tramite il sistema `params`) includono tipicamente cose come i file di input, i flag di comportamento degli strumenti e i parametri di analisi.
- La **configurazione** in senso stretto si riferisce alla logistica di come la pipeline viene eseguita, ovvero l'executor, le allocazioni di risorse di calcolo e così via.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Iniziamo affrontando i parametri della pipeline, poi esamineremo la configurazione in senso stretto.

---

## 1. Parametri della pipeline

Per tutte le pipeline nf-core, è possibile ottenere un elenco completo dei parametri della pipeline direttamente dalla riga di comando usando il flag `--help`, che è esso stesso un parametro della pipeline.

### 1.1. Ottenere l'elenco dei parametri con `--help`

Eseguiamo il comando help per la pipeline demo:

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

Come si può vedere, l'output raggruppa i parametri in categorie (opzioni di input/output, opzioni del genoma di riferimento, ecc.) con tipi e descrizioni per ciascuno.

Questa categorizzazione è determinata da un file di schema, di cui parleremo più avanti.
Nelle pipeline Nextflow semplici, `--help` funziona solo se lo sviluppatore lo ha implementato manualmente.

!!! tip "Suggerimento"

    Usa `--help --show_hidden` per vedere i parametri aggiuntivi che sono nascosti per impostazione predefinita, come `--publish_dir_mode` o `--monochrome_logs`.

### 1.2. Impostare i valori dei parametri

Come illustrato in [Hello Config](../hello_nextflow/06_hello_config.md), è possibile impostare i valori dei parametri dalla riga di comando con `--nome_parametro` oppure raccogliere un insieme di parametri in un file YAML e passarlo con `-params-file`.
Entrambi gli approcci funzionano allo stesso modo con le pipeline nf-core.

Ad esempio, per saltare la fase di trimming, vogliamo impostare il parametro booleano `skip_trim` a `true`.
Nella tua directory di lavoro è disponibile un file params chiamato `my_params.yml` con quel valore già impostato:

```yaml title="my_params.yml"
skip_trim: true
```

Passiamolo con `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
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

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Il processo `SEQTK_TRIM` non compare più nell'output.

!!! warning "Limitazioni importanti riguardo agli input dei parametri"

    **Impostare parametri booleani dalla riga di comando**

    A partire dalla versione 26.04 di Nextflow, tutti i valori forniti dalla riga di comando vengono tipizzati come stringhe.
    Per un parametro booleano come `skip_trim`, passarlo come flag semplice (`--skip_trim`) o come `--skip_trim true` viene valutato come la **stringa** `"true"`, che non supera la validazione dello schema:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Per impostare un parametro booleano a un valore genuino `true`/`false`, usa un `-params-file` come mostrato sopra, oppure impostalo in un file di configurazione.
    I parametri di tipo stringa, intero e percorso file non sono interessati e possono ancora essere impostati direttamente dalla riga di comando.
    Questo corso usa questo schema per tutti i parametri booleani.

    **Usare file di configurazione personalizzati**

    Sebbene sia tecnicamente possibile impostare i parametri della pipeline in un file di configurazione personalizzato passato con `-c`, questo potrebbe non sovrascrivere i valori predefiniti già impostati nel `nextflow.config` della pipeline, a seconda delle regole di precedenza della configurazione di Nextflow.
    Usare `--nome_parametro` dalla riga di comando o `-params-file` è più affidabile, poiché questi hanno sempre la precedenza.

    Come regola generale: se appare nell'output di `--help`, impostalo tramite la riga di comando o un file params piuttosto che un file di configurazione.

### 1.3. Validazione dei parametri

Curiosità: il comando `--help` funziona per tutte le pipeline nf-core perché il progetto nf-core richiede agli sviluppatori di definire formalmente tutti i parametri della pipeline in un file di schema JSON (`nextflow_schema.json`).
Questo schema registra il tipo, la descrizione, il valore predefinito e il raggruppamento di ciascun parametro.

Oltre a generare l'output di `--help`, il file di schema abilita anche la validazione automatica al momento dell'avvio.
Ciò significa che Nextflow può verificare che ogni parametro passato esista e abbia ricevuto un valore appropriato (del tipo corretto, nell'intervallo di valori consentito, ecc.).

Approfondiremo questo argomento nella [sezione sulla validazione dell'input](../nfcore_build/04_input_validation.md), ma puoi già vederlo in azione fornendo alla pipeline demo un input di parametri non valido.

#### 1.3.1. Parametri non riconosciuti

Proviamo a passare un parametro che non esiste:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

L'output della console include un avviso:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

La pipeline viene comunque eseguita, ma l'avviso ti segnala immediatamente che `--foobar` non è un parametro riconosciuto.
Lo scopo è attirare la tua attenzione su errori di battitura non bloccanti, come l'uso di `--outDir` invece di `--outdir`, che può aiutarti a evitare sprechi di tempo e risorse di calcolo.

#### 1.3.2. Valori di parametri non validi

La validazione controlla anche i **valori** dei parametri.
Il parametro `--skip_trim` è un flag booleano, quindi passare un valore stringa causa il fallimento immediato della pipeline:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

La pipeline si interrompe prima che vengano eseguiti i processi, evitando un'esecuzione fallita o errata.
Come indicato nella [sezione 1.2](#12-set-parameter-values), i parametri booleani dovrebbero essere impostati a un valore genuino `true`/`false` in un file params piuttosto che passati dalla riga di comando, poiché i valori da riga di comando vengono tipizzati come stringhe.

### 1.4. Validazione dell'input

La stessa logica di validazione può essere usata anche per verificare la validità dei file di input.
Ad esempio, se una pipeline si aspetta un samplesheet come input principale dei dati (il che è il caso di molte, se non della maggior parte, delle pipeline nf-core), lo sviluppatore può fornire uno schema di input (distinto dallo schema dei parametri) che descrive come deve essere strutturato il file di input.

In fase di esecuzione, Nextflow può quindi verificare che il file di input fornito sia valido.

Approfondiremo anche questo nella [sezione sulla validazione dell'input](../nfcore_build/04_input_validation.md), ma puoi già vederlo in azione fornendo alla pipeline demo un samplesheet di input non valido.

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

Per dimostrarlo, nella tua directory di lavoro è disponibile un samplesheet malformato chiamato `malformed_samplesheet.csv`:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Questo samplesheet manca della colonna obbligatoria `fastq_1` e contiene un percorso file inesistente in `fastq_2`.

Eseguiamo la pipeline demo usando `malformed_samplesheet.csv` come input:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Come si può vedere, la pipeline fallisce immediatamente e riporta **tutti** gli errori di validazione in una volta sola.
nf-schema non si ferma al primo errore — raccoglie ogni problema e li elenca insieme, così puoi correggere tutto in un'unica volta invece di scoprire i problemi uno alla volta.

Ogni errore identifica la voce e il campo esatti che hanno causato il problema, così puoi correggere il tuo samplesheet e rilanciare la pipeline con la certezza che non fallirà in un momento successivo quando Nextflow andrà effettivamente ad accedere al percorso del file.

Per gli sviluppatori, tutto questo è trattato in modo più dettagliato nella [Parte 4 di Build with nf-core](../nfcore_build/04_input_validation.md).

### Takeaway

Sai come ottenere un elenco completo dei parametri di una pipeline con `--help`, come impostarli tramite la riga di comando o un file params, e come Nextflow valida sia i valori dei parametri che i file di input rispetto agli schemi della pipeline.

### Cosa c'è dopo?

Scopriamo l'altro tipo di configurazione: come viene eseguita la pipeline, che riguarda l'allocazione delle risorse e gli argomenti degli strumenti.

---

## 2. Configurazione

La configurazione in senso stretto controlla **come** viene eseguita la pipeline: allocazione delle risorse, argomenti specifici degli strumenti, dove vengono eseguiti i job e quale sistema di packaging software utilizzare.

Le pipeline nf-core includono la configurazione predefinita in `nextflow.config` e nella directory `conf/`.
Prima di sovrascrivere qualsiasi cosa, è utile sapere dove si trovano i valori predefiniti.

### 2.1. Esplorare i file di configurazione

Hai già visto nella [Parte 1](./01_run_demo.md) che il codice sorgente della pipeline si trova sotto `$NXF_HOME/assets`.
Usando il symlink `pipelines` che hai creato nella [Parte 1](./01_run_demo.md), elenchiamo i file di configurazione per vedere cosa è disponibile:

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
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

I file di configurazione più importanti sono:

- **`conf/base.config`**: Definisce le etichette delle risorse (`process_low`, `process_medium`, `process_high`) che assegnano CPU, memoria e tempo ai processi. Quando vedi un processo che usa più risorse del previsto, è qui che si trovano quei valori predefiniti.
- **`conf/modules.config`**: Imposta gli argomenti degli strumenti per processo (`ext.args`) e le impostazioni di pubblicazione dell'output (`publishDir`). Apri questo file per vedere quali argomenti riceve ciascuno strumento per impostazione predefinita.
- **`conf/test.config`**: Il profilo di test che hai usato nella [Parte 1](./01_run_demo.md), che limita le risorse tramite `resourceLimits` e imposta un samplesheet di test. Viene attivato con `-profile test`.
  Esiste anche un `conf/test_full.config` per l'esecuzione con un dataset di test di dimensioni complete, utile per il benchmarking.

Il `nextflow.config` centrale carica tutti i file sopra indicati e imposta i valori predefiniti appropriati per tutto.

Se desideri modificare qualsiasi impostazione specificata in questi file, non modificare direttamente nessuno di essi.
Crea invece il tuo file di configurazione e passalo con `-c`.
I valori che specifichi sovrascriveranno i valori predefiniti impostati in quegli altri file.

Proviamo questo in pratica.

### 2.2. Personalizzare le risorse dei processi e gli argomenti degli strumenti

I moduli nf-core supportano due tipi comuni di override della configurazione: **allocazione delle risorse** (CPU, memoria, tempo) e **argomenti degli strumenti** tramite `ext.args`.

Molti strumenti da riga di comando hanno argomenti che non sono abbastanza comunemente usati da essere esposti come parametri della pipeline.
La convenzione `ext.args` ti permette di passare questi argomenti allo strumento sottostante tramite un file di configurazione.

Il file `custom.config` fornito nella tua directory di lavoro dimostra entrambi gli override:

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
Il flag `-b 5` dice a `seqtk trimfq` di rimuovere 5 basi dall'inizio di ogni read in aggiunta al trimming per qualità.

Eseguiamo la pipeline con questa configurazione:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
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

Il flag `-c` aggiunge la tua configurazione sopra la configurazione integrata della pipeline.

Per verificare che l'override di `ext.args` abbia avuto effetto, trova l'hash della directory di lavoro di `SEQTK_TRIM` dall'output dell'esecuzione (ad es. `work/17/428668...`) e controlla il file `.command.sh` al suo interno:

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

Dovresti vedere `-b 5` nel comando `seqtk trimfq`.

Una cosa importante da sapere su `ext.args`: se un modulo ha già un valore predefinito impostato, il tuo valore lo **sostituirà completamente** invece di aggiungersi ad esso.
Ad esempio, `FASTQC` ha `ext.args = '--quiet'` impostato per impostazione predefinita in `conf/modules.config`:

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

Se imposti `ext.args = '--kmers 8'` per `FASTQC`, il flag `--quiet` non verrà più applicato.
Per mantenere entrambi, imposta `ext.args = '--quiet --kmers 8'`.

Dovresti sempre controllare la configurazione predefinita di un modulo prima di sovrascrivere `ext.args`.

### Takeaway

Sai dove si trovano i valori predefiniti della configurazione delle pipeline nf-core e come sovrascrivere le allocazioni delle risorse e gli argomenti degli strumenti con un file di configurazione personalizzato.

### Cosa c'è dopo?

Passa alla [Parte 3](./03_run_production_pipeline.md), dove applicherai ciò che hai imparato a una vera pipeline di produzione.

---

## Riepilogo

In questa parte hai imparato a:

- Ottenere aiuto, impostare parametri e comprendere la validazione dei parametri e dell'input
- Personalizzare l'allocazione delle risorse e gli argomenti degli strumenti tramite file di configurazione
