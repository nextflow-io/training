# Parte 4: Creare un modulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa quarta parte del corso di formazione Hello nf-core, vi mostreremo come creare un modulo nf-core applicando le convenzioni chiave che rendono i moduli portabili e manutenibili.

Il progetto nf-core fornisce un comando (`nf-core modules create`) che genera automaticamente template di moduli correttamente strutturati, simile a quello che abbiamo utilizzato per il workflow nella Parte 2.
Tuttavia, per scopi didattici, inizieremo facendolo manualmente: trasformando il modulo locale `cowpy` nel vostro pipeline `core-hello` in un modulo in stile nf-core passo dopo passo.
Successivamente, vi mostreremo come utilizzare la creazione di moduli basata su template per lavorare in modo più efficiente in futuro.

??? info "Come iniziare da questa sezione"

    Questa sezione presuppone che abbiate completato la [Parte 3: Utilizzare un modulo nf-core](./03_use_module.md) e abbiate integrato il modulo `CAT_CAT` nel vostro pipeline.

    Se non avete completato la Parte 3 o desiderate iniziare da zero per questa parte, potete utilizzare la soluzione `core-hello-part3` come punto di partenza.
    Eseguite questi comandi dall'interno della directory `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Questo vi fornisce un pipeline con il modulo `CAT_CAT` già integrato.
    Potete verificare che funzioni correttamente eseguendo il seguente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Trasformare `cowpy` in un modulo nf-core

In questa sezione, applicheremo le convenzioni nf-core al modulo locale `cowpy` nel vostro pipeline `core-hello`, trasformandolo in un modulo che segue gli standard della community nf-core.

Questo è il codice attuale per il modulo di processo `cowpy`:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Applicheremo le seguenti convenzioni nf-core in modo incrementale:

1. **Mettere in maiuscolo il nome del processo in `COWPY`** per seguire la convenzione.
2. **Aggiornare `COWPY` per utilizzare tuple di metadati** per propagare i metadati dei campioni attraverso il workflow.
3. **Centralizzare la configurazione degli argomenti dello strumento con `ext.args`** per aumentare la versatilità del modulo mantenendo l'interfaccia minimale.
4. **Standardizzare la denominazione dell'output con `ext.prefix`** per promuovere la coerenza.
5. **Centralizzare la configurazione di pubblicazione** per promuovere la coerenza.

Dopo ogni passaggio, eseguiremo il pipeline per verificare che tutto funzioni come previsto.

!!! warning "Directory di lavoro"

    Assicuratevi di trovarvi nella directory `core-hello` (la radice del vostro pipeline) per tutte le modifiche ai file e le esecuzioni di comandi in questa sezione.

    ```bash
    cd core-hello
    ```

### 1.1. Mettere in maiuscolo il nome del processo

Questa è puramente una convenzione stilistica (non c'è alcuna giustificazione tecnica) ma poiché è la norma per i moduli nf-core, conformiamoci.

Dobbiamo effettuare tre serie di modifiche:

1. Aggiornare il nome del processo nel modulo
2. Aggiornare l'istruzione di import del modulo nell'intestazione del workflow
3. Aggiornare la chiamata del processo e la dichiarazione emit nel corpo del workflow

Iniziamo!

#### 1.1.1. Aggiornare il nome del processo nel modulo

Aprite il file del modulo `cowpy.nf` (sotto `core-hello/modules/local/`) e modificate il nome del processo in maiuscolo:

=== "Dopo"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Prima"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

In questo caso la conversione in maiuscolo è completamente diretta.

Se il nome del processo fosse composto da più parole, ad esempio se avessimo un processo chiamato MyCowpyTool originariamente in camel case, la convenzione nf-core sarebbe utilizzare underscore per separarle, ottenendo MY_COWPY_TOOL.

#### 1.1.2. Aggiornare l'istruzione di import del modulo

I nomi dei processi sono case-sensitive, quindi ora che abbiamo cambiato il nome del processo, dobbiamo aggiornare di conseguenza l'istruzione di import del modulo nell'intestazione del workflow di `hello.nf`:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Potremmo utilizzare un alias nell'istruzione di import per evitare di dover aggiornare le chiamate al processo, ma ciò vanificherebbe in qualche modo lo scopo di adottare la convenzione delle maiuscole.

#### 1.1.3. Aggiornare la chiamata del processo e la dichiarazione emit

Quindi ora aggiorniamo i due riferimenti al processo nel blocco workflow di `hello.nf`:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // genera arte ASCII dei saluti con cowpy
    COWPY(CAT_CAT.out.file_out)

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
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // genera arte ASCII dei saluti con cowpy
    cowpy(CAT_CAT.out.file_out)

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
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Assicuratevi di effettuare **entrambe** le modifiche, altrimenti riceverà un errore quando eseguirà questo.

#### 1.1.4. Eseguire il pipeline per testarlo

Eseguiamo il workflow per verificare che tutto funzioni correttamente dopo queste modifiche.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Bene, funziona! Ora passiamo a effettuare modifiche più sostanziali.

### 1.2. Aggiornare `COWPY` per utilizzare tuple di metadati

Nella versione attuale del pipeline `core-hello`, stiamo estraendo il file dalla tupla di output di `CAT_CAT` per passarlo a `COWPY`, come mostrato nella metà superiore del diagramma sottostante.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Sarebbe meglio avere `COWPY` che accetta direttamente tuple di metadati, permettendo ai metadati di fluire attraverso il workflow, come mostrato nella metà inferiore del diagramma.

A tal fine, dovremo effettuare le seguenti modifiche:

1. Aggiornare le definizioni di input e output
2. Aggiornare la chiamata del processo nel workflow
3. Aggiornare il blocco emit nel workflow

Una volta fatto tutto ciò, eseguiremo il pipeline per verificare che tutto funzioni ancora come prima.

#### 1.2.1. Aggiornare le definizioni di input e output

Torni al file del modulo `cowpy.nf` e lo modifichi per accettare tuple di metadati come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Prima"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Come potete vedere, abbiamo modificato sia l'**input principale** che l'**output** in una tupla che segue il pattern `tuple val(meta), path(input_file)` introdotto nella Parte 3 di questa formazione.
Per l'output, abbiamo anche colto questa opportunità per aggiungere `emit: cowpy_output` al fine di dare al canale di output un nome descrittivo.

Ora che abbiamo cambiato ciò che il processo si aspetta, dobbiamo aggiornare ciò che gli forniamo nella chiamata del processo.

#### 1.2.2. Aggiornare la chiamata del processo nel workflow

La buona notizia è che questa modifica semplificherà la chiamata del processo.
Ora che l'output di `CAT_CAT` e l'input di `COWPY` hanno la stessa 'forma', cioè entrambi consistono in una struttura `tuple val(meta), path(input_file)`, possiamo semplicemente connetterli direttamente invece di dover estrarre esplicitamente il file dall'output del processo `CAT_CAT`.

Aprite il file del workflow `hello.nf` (sotto `core-hello/workflows/`) e aggiorni la chiamata a `COWPY` come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // genera arte ASCII dei saluti con cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // genera arte ASCII dei saluti con cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Ora chiamiamo `COWPY` su `CAT_CAT.out.file_out` direttamente.

Di conseguenza, non abbiamo più bisogno di costruire il canale `ch_for_cowpy`, quindi quella riga (e la sua riga di commento) può essere rimossa completamente.

#### 1.2.3. Aggiornare il blocco emit nel workflow

Poiché `COWPY` ora emette un output nominato, `cowpy_output`, possiamo aggiornare il blocco `emit:` del workflow `hello.nf` per utilizzarlo.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Questo tecnicamente non è richiesto, ma è buona pratica fare riferimento a output nominati quando possibile.

#### 1.2.4. Eseguire il pipeline per testarlo

Eseguiamo il workflow per verificare che tutto funzioni correttamente dopo queste modifiche.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Il pipeline dovrebbe essere eseguito con successo, con i metadati che ora fluiscono da `CAT_CAT` attraverso `COWPY`.

Questo completa ciò che dovevamo fare per far gestire a `COWPY` le tuple di metadati.
Ora, vediamo cos'altro possiamo fare per sfruttare i pattern dei moduli nf-core.

### 1.3. Centralizzare la configurazione degli argomenti dello strumento con `ext.args`

Nel suo stato attuale, il processo `COWPY` si aspetta di ricevere un valore per il parametro `character`.
Di conseguenza, dobbiamo fornire un valore ogni volta che chiamiamo il processo, anche se saremmo contenti con i valori predefiniti impostati dallo strumento.
Per `COWPY` questo non è certamente un grande problema, ma per strumenti con molti parametri opzionali, può diventare piuttosto oneroso.

Il progetto nf-core raccomanda di utilizzare una funzionalità Nextflow chiamata [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) per gestire gli argomenti degli strumenti in modo più conveniente tramite file di configurazione.

Invece di dichiarare input di processo per ogni opzione dello strumento, si scrive il modulo per riferirsi a `ext.args` nella costruzione della sua riga di comando.
Quindi è solo questione di configurare la variabile `ext.args` per contenere gli argomenti e i valori che si desidera utilizzare nel file `modules.config`, che consolida i dettagli di configurazione per tutti i moduli.
Nextflow aggiungerà quegli argomenti con i loro valori nella riga di comando dello strumento a runtime.

Applichiamo questo approccio al modulo `COWPY`.
Dovremo effettuare le seguenti modifiche:

1. Aggiornare il modulo `COWPY`
2. Configurare `ext.args` nel file `modules.config`
3. Aggiornare il workflow `hello.nf`

Una volta fatto tutto ciò, eseguiremo il pipeline per verificare che tutto funzioni ancora come prima.

#### 1.3.1. Aggiornare il modulo `COWPY`

Facciamolo.
Aprite il file del modulo `cowpy.nf` (sotto `core-hello/modules/local/`) e lo modifichi per riferirsi a `ext.args` come mostrato di seguito.

=== "Dopo"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Prima"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Può vedere che abbiamo effettuato tre modifiche.

1. **Nel blocco `input:`, abbiamo rimosso l'input `val character`.**
   D'ora in poi, forniremo quell'argomento tramite la configurazione `ext.args` come descritto più avanti.

2. **Nel blocco `script:`, abbiamo aggiunto la riga `def args = task.ext.args ?: ''`.**
   Quella riga utilizza l'operatore `?:` per determinare il valore della variabile `args`: il contenuto di `task.ext.args` se non è vuoto, o una stringa vuota se lo è.
   Notate che mentre generalmente ci riferiamo a `ext.args`, questo codice deve riferirsi a `task.ext.args` per estrarre la configurazione `ext.args` a livello di modulo.

3. **Nella riga di comando, abbiamo sostituito `-c "$character"` con `$args`.**
   Qui è dove Nextflow inietterà eventuali argomenti dello strumento impostati in `ext.args` nel file `modules.config`.

Di conseguenza, l'interfaccia del modulo è ora più semplice: si aspetta solo gli input essenziali di metadati e file.

!!! note

    L'operatore `?:` è spesso chiamato 'operatore Elvis' perché assomiglia a un viso di Elvis Presley laterale, con il carattere `?` che simboleggia l'onda nei suoi capelli.

#### 1.3.2. Configurare `ext.args` nel file `modules.config`

Ora che abbiamo tolto la dichiarazione `character` dal modulo, dobbiamo aggiungerla a `ext.args` nel file di configurazione `modules.config`.

Nello specifico, aggiungeremo questo piccolo pezzo di codice al blocco `process {}`:

```groovy title="Codice da aggiungere"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

La sintassi `withName:` assegna questa configurazione solo al processo `COWPY`, e `ext.args = { "-c ${params.character}" }` semplicemente compone una stringa che includerà il valore del parametro `character`.
Notate l'uso delle parentesi graffe, che dicono a Nextflow di valutare il valore del parametro a runtime.

Ha senso? Aggiungiamolo.

Aprite `conf/modules.config` e aggiungete il codice di configurazione all'interno del blocco `process {}` come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Prima"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Speriamo possiate immaginare di avere tutti i moduli in un pipeline con i loro `ext.args` specificati in questo file, con i seguenti vantaggi:

- L'**interfaccia del modulo rimane semplice** - Accetta solo gli input essenziali di metadati e file
- Il **pipeline espone ancora `params.character`** - Gli utenti finali possono ancora configurarlo come prima
- Il **modulo è ora portabile** - Può essere riutilizzato in altri pipeline senza dover aspettarsi un nome di parametro specifico
- La configurazione è **centralizzata** in `modules.config`, mantenendo pulita la logica del workflow

Utilizzando il file `modules.config` come luogo dove tutti i pipeline centralizzano la configurazione per-modulo, rendiamo i nostri moduli più riutilizzabili attraverso diversi pipeline.

#### 1.3.3. Aggiornare il workflow `hello.nf`

Poiché il modulo `COWPY` non richiede più il parametro `character` come input, dobbiamo aggiornare di conseguenza la chiamata del workflow.

Aprite il file del workflow `hello.nf` (sotto `core-hello/workflows/`) e aggiorni la chiamata a `COWPY` come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // genera arte ASCII dei saluti con cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // genera arte ASCII dei saluti con cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

Il codice del workflow è ora più pulito: non abbiamo bisogno di passare `params.character` direttamente al processo.
L'interfaccia del modulo è mantenuta minimale, rendendola più portabile, mentre il pipeline fornisce ancora l'opzione esplicita tramite configurazione.

#### 1.3.4. Eseguire il pipeline per testarlo

Verifichiamo che il workflow funzioni ancora come previsto, specificando un personaggio diverso per verificare che la configurazione `ext.args` funzioni.

Eseguite questo comando usando `kosh`, una delle opzioni più... enigmatiche:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Questo dovrebbe essere eseguito con successo come in precedenza.

Verifichiamo che la configurazione `ext.args` abbia funzionato controllando l'output.
Trovi l'output nel file browser o utilizzi l'hash dell'attività (la parte `38/eb29ea` nell'esempio sopra) per guardare il file di output:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Output del comando"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

Dovrebbe vedere l'arte ASCII visualizzata con il personaggio `kosh`, confermando che la configurazione `ext.args` ha funzionato!

??? info "(Opzionale) Ispezionare il file di comando"

    Se desideratete vedere esattamente come è stata applicata la configurazione, potete ispezionare il file `.command.sh`:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Vedrà il comando `cowpy` con l'argomento `-c kosh`:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Questo mostra che il file `.command.sh` è stato generato correttamente in base alla configurazione `ext.args`.

Si prenda un momento per riflettere su ciò che abbiamo ottenuto qui.
Questo approccio mantiene l'interfaccia del modulo focalizzata sui dati essenziali (file, metadati e eventuali parametri obbligatori per campione), mentre le opzioni che controllano il comportamento dello strumento sono gestite separatamente tramite configurazione.

Questo può sembrare superfluo per uno strumento semplice come `cowpy`, ma potete fare una grande differenza per gli strumenti di analisi dati che hanno molti argomenti opzionali.

Per riassumere i vantaggi di questo approccio:

- **Interfaccia pulita**: Il modulo si concentra sugli input di dati essenziali (metadati e file)
- **Flessibilità**: Gli utenti possono specificare argomenti dello strumento tramite configurazione, inclusi valori specifici per campione
- **Coerenza**: Tutti i moduli nf-core seguono questo pattern
- **Portabilità**: I moduli possono essere riutilizzati senza opzioni dello strumento hardcoded
- **Nessuna modifica al workflow**: Aggiungere o modificare opzioni dello strumento non richiede l'aggiornamento del codice del workflow

!!! note

    Il sistema `ext.args` ha potenti capacità aggiuntive non trattate qui, inclusa la commutazione dinamica dei valori degli argomenti in base ai metadati. Veda le [specifiche dei moduli nf-core](https://nf-co.re/docs/guidelines/components/modules) per maggiori dettagli.

### 1.4. Standardizzare la denominazione dell'output con `ext.prefix`

Ora che abbiamo dato al processo `COWPY` accesso alla metamap, possiamo iniziare a sfruttare un altro utile pattern nf-core: denominare i file di output in base ai metadati.

Qui utilizzeremo una funzionalità Nextflow chiamata `ext.prefix` che ci permetterà di standardizzare la denominazione dei file di output attraverso i moduli usando `meta.id` (l'identificatore incluso nella metamap), pur essendo ancora in grado di configurare i moduli individualmente se desiderateto.

Questo sarà simile a ciò che abbiamo fatto con `ext.args`, con alcune differenze che dettaglieremo man mano che procediamo.

Applichiamo questo approccio al modulo `COWPY`.
Dovremo effettuare le seguenti modifiche:

1. Aggiornare il modulo `COWPY`
2. Configurare `ext.prefix` nel file `modules.config`

(Nessuna modifica necessaria al workflow.)

Una volta fatto ciò, eseguiremo il pipeline per verificare che tutto funzioni ancora come prima.

#### 1.4.1. Aggiornare il modulo `COWPY`

Aprite il file del modulo `cowpy.nf` (sotto `core-hello/modules/local/`) e lo modifichi per riferirsi a `ext.prefix` come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Prima"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Potete vedere che abbiamo effettuato tre modifiche.

1. **Nel blocco `script:`, abbiamo aggiunto la riga `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Quella riga utilizza l'operatore `?:` per determinare il valore della variabile `prefix`: il contenuto di `task.ext.prefix` se non è vuoto, o l'identificatore dalla metamap (`meta.id`) se lo è.
   Notate che mentre generalmente ci riferiamo a `ext.prefix`, questo codice deve riferirsi a `task.ext.prefix` per estrarre la configurazione `ext.prefix` a livello di modulo.

2. **Nella riga di comando, abbiamo sostituito `cowpy-${input_file}` con `${prefix}.txt`.**
   Qui è dove Nextflow inietterà il valore di `prefix` determinato dalla riga sopra.

3. **Nel blocco `output:`, abbiamo sostituito `path("cowpy-${input_file}")` con `path("${prefix}.txt")`.**
   Questo semplicemente ribadisce quale sarà il percorso del file secondo ciò che è scritto nella riga di comando.

Di conseguenza, il nome del file di output è ora costruito utilizzando un valore predefinito sensato (l'identificatore dalla metamap) combinato con l'estensione del formato file appropriata.

#### 1.4.2. Configurare `ext.prefix` nel file `modules.config`

In questo caso il valore predefinito sensato non è sufficientemente espressivo per il nostro gusto; vogliamo utilizzare un pattern di denominazione personalizzato che includa il nome dello strumento, `cowpy-<id>.txt`, come avevamo prima.

Lo faremo configurando `ext.prefix` in `modules.config`, proprio come abbiamo fatto per il parametro `character` con `ext.args`, tranne che questa volta il blocco `withName: 'COWPY' {}` esiste già, e dobbiamo solo aggiungere la seguente riga:

```groovy title="Codice da aggiungere"
ext.prefix = { "cowpy-${meta.id}" }
```

Questo comporrà la stringa che vogliamo.
Notate che ancora una volta utilizziamo le parentesi graffe, questa volta per dire a Nextflow di valutare il valore di `meta.id` a runtime.

Aggiungiamolo.

Aprite `conf/modules.config` e aggiungete il codice di configurazione all'interno del blocco `process {}` come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Prima"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Nel caso si stesse chiedendo, la closure `ext.prefix` ha accesso al pezzo corretto di metadati perché la configurazione viene valutata nel contesto dell'esecuzione del processo, dove i metadati sono disponibili.

#### 1.4.3. Eseguire il pipeline per testarlo

Verifichiamo che il workflow funzioni ancora come previsto.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dia un'occhiata all'output nella directory dei risultati.
Dovrebbe vedere il file di output cowpy con la stessa denominazione di prima: `cowpy-test.txt`, basato sul nome batch predefinito.

??? abstract "Contenuto della directory"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Sentitevi liberi di cambiare la configurazione `ext.prefix` in `conf/modules.config` per convincervi che potete cambiare il pattern di denominazione senza dover apportare modifiche al codice del modulo o del workflow.

In alternativa, potete anche provare a eseguirlo di nuovo con un parametro `--batch` diverso specificato sulla riga di comando per convincervi che quella parte sia ancora personalizzabile al volo.

Questo dimostra come `ext.prefix` vi permetta di mantenere la vostra convenzione di denominazione preferita mantenendo flessibile l'interfaccia del modulo.

Per riassumere i vantaggi di questo approccio:

- **Denominazione standardizzata**: I file di output sono tipicamente denominati utilizzando gli ID dei campioni dai metadati
- **Configurabile**: Gli utenti possono sovrascrivere la denominazione predefinita se necessario
- **Coerente**: Tutti i moduli nf-core seguono questo pattern
- **Prevedibile**: Facile sapere come saranno chiamati i file di output

Abbastanza bene, vero?
Beh, c'è un'altra modifica importante che dobbiamo fare per migliorare il nostro modulo per adattarlo alle linee guida nf-core.

### 1.5. Centralizzare la configurazione di pubblicazione

Potrebbe aver notato che abbiamo pubblicato output in due directory diverse:

- **`results`** — La directory di output originale che abbiamo utilizzato dall'inizio per i nostri moduli locali, impostata individualmente utilizzando direttive `publishDir` per-modulo;
- **`core-hello-results`** — La directory di output impostata con `--outdir` sulla riga di comando, che ha ricevuto i log nf-core e i risultati pubblicati da `CAT_CAT`.

Questo è disordinato e subottimale; sarebbe meglio avere una posizione per tutto.
Ovviamente, potremmo andare in ciascuno dei nostri moduli locali e aggiornare manualmente la direttiva `publishDir` per utilizzare la directory `core-hello-results`, ma che dire della prossima volta che decidiamo di cambiare la directory di output?

Avere moduli individuali che prendono decisioni di pubblicazione non è chiaramente la strada da percorrere, specialmente in un mondo dove lo stesso modulo potrebbe essere utilizzato in molti pipeline diversi, da persone che hanno esigenze o preferenze diverse.
Vogliamo poter controllare dove vengono pubblicati gli output a livello di configurazione del workflow.

"Ehi," potrebbe dire, "`CAT_CAT` sta inviando i suoi output a `--outdir`. Forse dovremmo copiare la sua direttiva `publishDir`?"

Sì, è un'ottima idea.

Tranne che non ha una direttiva `publishDir`. (Vada avanti, guardi il codice del modulo.)

Questo perché i pipeline nf-core centralizzano il controllo a livello di workflow configurando `publishDir` in `conf/modules.config` piuttosto che nei singoli moduli.
Nello specifico, il template nf-core dichiara una direttiva `publishDir` predefinita (con una struttura di directory predefinita) che si applica a tutti i moduli a meno che non venga fornita una direttiva sovrascrivente.

Non suona fantastico? Potrebbe essere che per sfruttare questa direttiva predefinita, tutto ciò che dobbiamo fare è rimuovere la direttiva `publishDir` attuale dai nostri moduli locali?

Proviamolo su `COWPY` per vedere cosa succede, poi guarderemo il codice per la configurazione predefinita per capire come funziona.

Infine, dimostreremo come sovrascrivere il comportamento predefinito se desiderateto.

#### 1.5.1. Rimuovere la direttiva `publishDir` da `COWPY`

Facciamolo.
Aprite il file del modulo `cowpy.nf` (sotto `core-hello/modules/local/`) e rimuova la direttiva `publishDir` come mostrato di seguito.

=== "Dopo"

    ```groovy title="core-hello/modules/local/cowpy.nf (estratto)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Prima"

    ```groovy title="core-hello/modules/local/cowpy.nf (estratto)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

È tutto!

#### 1.5.2. Eseguire il pipeline per testarlo

Diamo un'occhiata a cosa succede se eseguiamo il pipeline ora.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dia un'occhiata alla vostra directory di lavoro corrente.
Ora `core-hello-results` contiene anche gli output del modulo `COWPY`.

??? abstract "Contenuto della directory"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Potete vedere che Nextflow ha creato questa gerarchia di directory basata sui nomi del workflow e del modulo.

Il codice responsabile si trova nel file `conf/modules.config`.
Questa è la configurazione `publishDir` predefinita che fa parte del template nf-core e si applica a tutti i processi:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Questo può sembrare complicato, quindi esaminiamo ciascuno dei tre componenti:

- **`path:`** Determina la directory di output in base al nome del processo.
  Il nome completo di un processo contenuto in `task.process` include la gerarchia di import di workflow e moduli (come `CORE_HELLO:HELLO:CAT_CAT`).
  Le operazioni `tokenize` eliminano quella gerarchia per ottenere solo il nome del processo, poi prendono la prima parte prima di qualsiasi underscore (se applicabile), e la convertono in minuscolo.
  Questo è ciò che determina che i risultati di `CAT_CAT` vengano pubblicati in `${params.outdir}/cat/`.
- **`mode:`** Controlla come vengono pubblicati i file (copy, symlink, ecc.).
  Questo è configurabile tramite il parametro `params.publish_dir_mode`.
- **`saveAs:`** Filtra quali file pubblicare.
  Questo esempio esclude i file `versions.yml` restituendo `null` per loro, impedendo che vengano pubblicati.

Questo fornisce una logica coerente per organizzare gli output.

L'output appare ancora meglio quando tutti i moduli in un pipeline adottano questa convenzione, quindi sentitevi liberi di andare a eliminare le direttive `publishDir` dagli altri moduli nel vostro pipeline.
Questo valore predefinito verrà applicato anche ai moduli che non abbiamo esplicitamente modificato per seguire le linee guida nf-core.

Detto questo, potreste decidere di voler organizzare i vostri input in modo diverso, e la buona notizia è che è facile farlo.

#### 1.5.3. Sovrascrivere il valore predefinito

Per sovrascrivere la direttiva `publishDir` predefinita, potete semplicemente aggiungere le vostre direttive al file `conf/modules.config`.

Ad esempio, potreste sovrascrivere il valore predefinito per un singolo processo utilizzando il selettore `withName:`, come in questo esempio dove aggiungiamo una direttiva `publishDir` personalizzata per il processo 'COWPY'.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

Non effettueremo effettivamente quella modifica, ma sentitevi liberi di sperimentare con questo e vedere quale logica potete implementare.

Il punto è che questo sistema vi dà il meglio di entrambi i mondi: coerenza per default e la flessibilità di personalizzare la configurazione su richiesta.

Per riassumere, ottenete:

- **Unica fonte di verità**: Tutta la configurazione di pubblicazione risiede in `modules.config`
- **Default utile**: I processi funzionano out-of-the-box senza configurazione per-modulo
- **Personalizzazione facile**: Sovrascrivete il comportamento di pubblicazione nella configurazione, non nel codice del modulo
- **Moduli portabili**: I moduli non hardcodeano le posizioni di output

Questo completa l'insieme di funzionalità dei moduli nf-core che dovreste assolutamente imparare a utilizzare, ma ce ne sono altre che potete leggere nelle [specifiche dei moduli nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Takeaway

Ora sa come adattare i moduli locali per seguire le convenzioni nf-core:

- Progettare i vostri moduli per accettare e propagare tuple di metadati;
- Utilizzare `ext.args` per mantenere le interfacce dei moduli minimali e portabili;
- Utilizzare `ext.prefix` per una denominazione dei file di output configurabile e standardizzata;
- Adottare la direttiva `publishDir` centralizzata predefinita per una struttura della directory dei risultati coerente.

### Prossimi passi

Imparate come utilizzare gli strumenti integrati nf-core basati su template per creare moduli nel modo più semplice.

---

## 2. Creare un modulo con gli strumenti nf-core

Ora che ha imparato i pattern dei moduli nf-core applicandoli manualmente, vediamo come creerebbe i moduli nella pratica.

### 2.1. Generare una struttura di modulo da un template

Simile a ciò che esiste per la creazione di pipeline, il progetto nf-core fornisce strumenti per generare moduli correttamente strutturati basati su un template, con tutti questi pattern integrati dall'inizio.

#### 2.1.1. Eseguire il comando di creazione del modulo

Il comando `nf-core modules create` genera un template di modulo che segue già tutte le convenzioni che ha imparato.

Creiamo una nuova versione del modulo `COWPY` con un template minimale eseguendo questo comando:

```bash
nf-core modules create --empty-template COWPY
```

Il flag `--empty-template` crea un template iniziale pulito senza codice extra, rendendo più facile vedere la struttura essenziale.

Il comando viene eseguito in modo interattivo, guidandoLa attraverso la configurazione.
Cerca automaticamente le informazioni sullo strumento da repository di pacchetti come Bioconda e bio.tools per pre-popolare i metadati.

Le verrà richiesto di configurare diverse opzioni:

- **Informazioni sull'autore**: Il vostro nome utente GitHub per l'attribuzione
- **Etichetta di risorsa**: Un insieme predefinito di requisiti computazionali.
  Il progetto nf-core fornisce etichette standard come `process_single` per strumenti leggeri e `process_high` per quelli impegnativi.
  Queste etichette aiutano a gestire l'allocazione delle risorse in diversi ambienti di esecuzione.
- **Requisito di metadati**: Se il modulo necessita di informazioni specifiche per campione tramite una mappa `meta` (di solito sì per i moduli di elaborazione dati).

Lo strumento gestisce la complessità di trovare informazioni sui pacchetti e configurare la struttura, permettendoLe di concentrarSi sull'implementazione della logica specifica dello strumento.

#### 2.1.2. Esaminare la struttura del modulo

Lo strumento crea una struttura di modulo completa in `modules/local/` (o `modules/nf-core/` se si trova nel repository nf-core/modules):

??? abstract "Contenuto della directory"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Ogni file serve uno scopo specifico:

- **`main.nf`**: Definizione del processo con tutti i pattern nf-core integrati
- **`meta.yml`**: Documentazione del modulo che descrive input, output e lo strumento
- **`environment.yml`**: Specifica dell'ambiente Conda per le dipendenze
- **`tests/main.nf.test`**: Casi di test nf-test per validare il funzionamento del modulo

!!! tip "Impari di più sui test"

    Il file di test generato utilizza nf-test, un framework di testing per pipeline e moduli Nextflow. Per imparare come scrivere ed eseguire questi test, veda la [side quest nf-test](../side_quests/nf-test.md).

Il `main.nf` generato include tutti i pattern che ha appena imparato, più alcune funzionalità aggiuntive:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.
```
