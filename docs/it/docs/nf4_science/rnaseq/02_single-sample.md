# Parte 2: Implementazione per singolo campione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa parte del corso, scriveremo il workflow più semplice possibile che racchiude tutti i comandi eseguiti nella Parte 1 per automatizzarne l'esecuzione, e ci limiteremo a processare un campione alla volta.

!!! warning "Prerequisito"

    È necessario completare la [Parte 1: Panoramica del metodo](./01_method.md) prima di iniziare questa lezione.
    Nello specifico, il completamento della sezione 1.2.3 crea il file di indice del genoma (`data/genome_index.tar.gz`) richiesto per il passaggio di allineamento in questa lezione.

## Obiettivo

In questa parte del corso, svilupperemo un workflow che esegue le seguenti operazioni:

1. Eseguire il controllo di qualità (FastQC) sulle letture di input
2. Trimmare gli adapter ed eseguire il QC post-trimming (Trim Galore)
3. Allineare le letture trimmate a un genoma di riferimento (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

Questo automatizza i passaggi della prima sezione della [Parte 1: Panoramica del metodo](./01_method.md#1-single-sample-processing), dove avete eseguito questi comandi manualmente nei loro container.

Come punto di partenza, vi forniamo un file di workflow, `rnaseq.nf`, che delinea le parti principali del workflow, insieme a quattro file modulo nella directory `modules/` (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf` e `multiqc.nf`) che delineano la struttura di ciascun processo.

??? full-code "File di base"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Istruzioni INCLUDE dei moduli

    /*
     * Parametri della pipeline
     */

    // Input primario

    workflow {

        main:
        // Crea canale di input

        // Chiama i processi

        publish:
        // Dichiara gli output da pubblicare
    }

    output {
        // Configura i target di pubblicazione
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Esegue FastQC sulle letture di input
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trimma gli adapter ed esegue il QC post-trimming
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Allinea le letture a un genoma di riferimento
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggrega i report QC con MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

Questi file non sono funzionali; il loro scopo è solo quello di servire come scheletri da riempire con le parti interessanti del codice.

## Piano della lezione

Per rendere il processo di sviluppo più educativo, abbiamo suddiviso questo in tre fasi:

1. **Scrivere un workflow a singolo stadio che esegue il passaggio di QC iniziale.**
   Questo copre la configurazione di un parametro CLI, la creazione di un canale di input, la scrittura di un modulo di processo e la configurazione della pubblicazione degli output.
2. **Aggiungere il trimming degli adapter e il QC post-trimming.**
   Questo introduce il concatenamento dei processi collegando l'output di un processo all'input di un altro.
3. **Aggiungere l'allineamento al genoma di riferimento.**
   Questo copre la gestione di input di riferimento aggiuntivi e il lavoro con archivi compressi.

Ogni passaggio si concentra su un aspetto specifico dello sviluppo del workflow.

!!! tip "Suggerimento"

     Assicuratevi di essere nella directory di lavoro corretta:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Scrivere un workflow a singolo stadio che esegue il QC iniziale

Questo primo passaggio si concentra sulle basi: caricare un file FASTQ ed eseguire il controllo di qualità su di esso.

Ricordate il comando `fastqc` dalla [Parte 1](01_method.md):

```bash
fastqc <reads>
```

Il comando prende un file FASTQ come input e produce un report di controllo qualità come archivio `.zip` e un riepilogo `.html`.
L'URI del container era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Prenderemo queste informazioni e le racchiuderemo in Nextflow in tre fasi:

1. Configurare l'input
2. Scrivere il processo QC e chiamarlo nel workflow
3. Configurare la gestione dell'output

### 1.1. Configurare l'input

Dobbiamo dichiarare un parametro di input, creare un profilo di test per fornire un valore predefinito conveniente e creare un canale di input.

#### 1.1.1. Aggiungere una dichiarazione di parametro di input

In `rnaseq.nf`, sotto la sezione `Pipeline parameters`, dichiarate un parametro chiamato `reads` con il tipo `Path`.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Parametri della pipeline
     */
    params {
        // Input primario
        input: Path
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Parametri della pipeline
     */

    // Input primario
    ```

Questo configura il parametro CLI, ma non vogliamo digitare il percorso del file ogni volta che eseguiamo il workflow durante lo sviluppo.
Ci sono diverse opzioni per fornire un valore predefinito; qui usiamo un profilo di test.

#### 1.1.2. Creare un profilo di test con un valore predefinito in `nextflow.config`

Un profilo di test fornisce valori predefiniti convenienti per provare un workflow senza specificare input dalla riga di comando.
Questa è una convenzione comune nell'ecosistema Nextflow (vedi [Hello Config](../../hello_nextflow/06_hello_config.md) per maggiori dettagli).

Aggiungete un blocco `profiles` a `nextflow.config` con un profilo `test` che imposta il parametro `reads` su uno dei file FASTQ di test.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Qui stiamo usando `#!groovy ${projectDir}`, una variabile integrata di Nextflow che punta alla directory dove si trova lo script del workflow.
Questo rende facile fare riferimento a file di dati e altre risorse senza codificare percorsi assoluti.

Il parametro ora ha un valore predefinito conveniente. Successivamente, dobbiamo creare un canale da esso.

#### 1.1.3. Configurare il canale di input

Nel blocco workflow, create un canale di input dal valore del parametro usando la fabbrica di canali `.fromPath` (come usato in [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Chiama i processi

        publish:
        // Dichiara gli output da pubblicare
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Crea canale di input

        // Chiama i processi

        publish:
        // Dichiara gli output da pubblicare
    }
    ```

Successivamente, dovremo creare il processo per eseguire il QC su questo input.

### 1.2. Scrivere il processo QC e chiamarlo nel workflow

Dobbiamo riempire la definizione del processo nel file modulo, importarlo nel workflow usando un'istruzione include e chiamarlo sull'input.

#### 1.2.1. Riempire il modulo per il processo QC

Aprite `modules/fastqc.nf` ed esaminate la struttura della definizione del processo.
Dovreste riconoscere gli elementi strutturali principali; in caso contrario, considerate di leggere [Hello Nextflow](../../hello_nextflow/01_hello_world.md) per un ripasso.

Procedete e riempite la definizione del processo da soli usando le informazioni fornite sopra, quindi verificate il vostro lavoro confrontandolo con la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Esegue FastQC sulle letture di input
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Dopo"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Esegue FastQC sulle letture di input
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

L'accessor `simpleName` rimuove tutte le estensioni dal nome del file, quindi `ENCSR000COQ1_1.fastq.gz` diventa `ENCSR000COQ1_1`.
Usiamo la sintassi `emit:` per assegnare nomi a ciascun canale di output, il che sarà utile per collegare gli output nel blocco publish.

Una volta completato questo, il processo è completo.
Per usarlo nel workflow, dovrete importare il modulo e aggiungere una chiamata al processo.

#### 1.2.2. Includere il modulo

In `rnaseq.nf`, aggiungete un'istruzione `include` per rendere il processo disponibile al workflow:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="3"
    // Istruzioni INCLUDE dei moduli
    ```

Il processo è ora disponibile nello scope del workflow.

#### 1.2.3. Chiamare il processo QC sull'input

Aggiungete una chiamata a `FASTQC` nel blocco workflow, passando il canale di input come argomento.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Controllo di qualità iniziale
        FASTQC(read_ch)

        publish:
        // Dichiara gli output da pubblicare
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Chiama i processi

        publish:
        // Dichiara gli output da pubblicare
    }
    ```

Il workflow ora carica l'input ed esegue il processo QC su di esso.
Successivamente, dobbiamo configurare come viene pubblicato l'output.

### 1.3. Configurare la gestione dell'output

Dobbiamo dichiarare quali output del processo pubblicare e specificare dove devono andare.

#### 1.3.1. Dichiarare gli output nella sezione `publish:`

La sezione `publish:` all'interno del blocco workflow dichiara quali output del processo devono essere pubblicati.
Assegnate gli output di `FASTQC` a target nominati.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Dichiara gli output da pubblicare
    }
    ```

Successivamente, dovremo dire a Nextflow dove mettere gli output pubblicati.

#### 1.3.2. Configurare i target di output nel blocco `output {}`

Il blocco `output {}` si trova all'esterno del workflow e specifica dove viene pubblicato ciascun target nominato.
Configurate entrambi i target per pubblicare in una sottodirectory `fastqc/`.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configura i target di pubblicazione
    }
    ```

!!! note "Nota"

    Per impostazione predefinita, Nextflow pubblica i file di output come link simbolici, il che evita duplicazioni non necessarie.
    Anche se i file di dati che stiamo usando qui sono molto piccoli, in genomica possono diventare molto grandi.
    I link simbolici si romperanno quando pulite la vostra directory `work`, quindi per i workflow di produzione potreste voler sovrascrivere la modalità di pubblicazione predefinita con `'copy'`.

### 1.4. Eseguire il workflow

A questo punto, abbiamo un workflow QC a un passaggio che dovrebbe essere completamente funzionale.

Eseguiamo con `-profile test` per usare il valore predefinito configurato nel profilo di test, evitando la necessità di scrivere il percorso sulla riga di comando.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

Questo dovrebbe essere eseguito molto rapidamente se avete completato la Parte 1 e avete già scaricato il container.
Se l'avete saltata, Nextflow scaricherà il container per voi; non dovete fare nulla perché accada, ma potreste dover attendere fino a un minuto.

Potete verificare gli output nella directory results.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

I report QC per il campione sono ora pubblicati nella sottodirectory `fastqc/`.

### Takeaway

Sapete come creare un modulo contenente un processo, importarlo in un workflow, chiamarlo con un canale di input e pubblicare i risultati usando il blocco output a livello di workflow.

### Cosa c'è dopo?

Aggiungete il trimming degli adapter con QC post-trimming come secondo passaggio nel workflow.

---

## 2. Aggiungere il trimming degli adapter e il QC post-trimming

Ora che abbiamo il QC iniziale in atto, possiamo aggiungere il passaggio di trimming degli adapter con il suo QC post-trimming integrato.

Ricordate il comando `trim_galore` dalla [Parte 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

Il comando trimma gli adapter da un file FASTQ ed esegue FastQC sull'output trimmato.
Produce letture trimmate, un report di trimming e report FastQC per le letture trimmate.
L'URI del container era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Dobbiamo solo scrivere la definizione del processo, importarlo, chiamarlo nel workflow e aggiornare la gestione dell'output.

### 2.1. Scrivere il processo di trimming e chiamarlo nel workflow

Come prima, dobbiamo riempire la definizione del processo, importare il modulo e aggiungere la chiamata al processo.

#### 2.1.1. Riempire il modulo per il processo di trimming

Aprite `modules/trim_galore.nf` ed esaminate la struttura della definizione del processo.

Procedete e riempite la definizione del processo da soli usando le informazioni fornite sopra, quindi verificate il vostro lavoro confrontandolo con la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trimma gli adapter ed esegue il QC post-trimming
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Dopo"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trimma gli adapter ed esegue il QC post-trimming
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    }
    ```

Questo processo ha tre output nominati: le letture trimmate che alimentano il passaggio di allineamento, il report di trimming e i report FastQC post-trimming.
Il flag `--fastqc` dice a Trim Galore di eseguire automaticamente FastQC sull'output trimmato.

#### 2.1.2. Includere il modulo

Aggiornate `rnaseq.nf` per importare il nuovo modulo:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="3"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    ```

Successivamente, aggiungeremo la chiamata al processo nel workflow.

#### 2.1.3. Chiamare il processo di trimming sull'input

Aggiungete la chiamata al processo nel blocco workflow:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Controllo di qualità iniziale
        FASTQC(read_ch)

        // Trimming degli adapter e QC post-trimming
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Controllo di qualità iniziale
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Il processo di trimming è ora collegato al workflow.

### 2.2. Aggiornare la gestione dell'output

Dobbiamo aggiungere gli output di trimming alla dichiarazione publish e configurare dove vanno.

#### 2.2.1. Aggiungere target di pubblicazione per gli output di trimming

Aggiungete gli output di trimming alla sezione `publish:`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Successivamente, dovremo dire a Nextflow dove mettere questi output.

#### 2.2.2. Configurare i nuovi target di output

Aggiungete voci per i target di trimming nel blocco `output {}`, pubblicandoli in una sottodirectory `trimming/`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

La configurazione dell'output è completa.

### 2.3. Eseguire il workflow

Il workflow ora include sia il QC iniziale che il trimming degli adapter.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

Anche questo dovrebbe essere eseguito molto rapidamente, dato che stiamo lavorando su un file di input così piccolo.

Potete trovare gli output di trimming nella directory results.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

Gli output di trimming e i report QC post-trimming sono ora nella sottodirectory `trimming/`.

### Takeaway

Sapete come aggiungere un secondo passaggio di elaborazione che viene eseguito indipendentemente sullo stesso input, producendo più output nominati.

### Cosa c'è dopo?

Aggiungete il passaggio di allineamento che si concatena dall'output delle letture trimmate.

---

## 3. Aggiungere l'allineamento al genoma di riferimento

Infine possiamo aggiungere il passaggio di allineamento del genoma usando HISAT2.

Ricordate il comando di allineamento dalla [Parte 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

Il comando allinea le letture a un genoma di riferimento e converte l'output in formato BAM.
Richiede un archivio di indice del genoma pre-costruito e produce un file BAM e un log di riepilogo dell'allineamento.
L'URI del container era `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

Questo processo richiede un input aggiuntivo (l'archivio dell'indice del genoma), quindi dobbiamo configurarlo prima, poi scrivere e collegare il processo.

### 3.1. Configurare gli input

Dobbiamo dichiarare un parametro per l'archivio dell'indice del genoma.

#### 3.1.1. Aggiungere un parametro per l'indice del genoma

Aggiungete una dichiarazione di parametro per l'archivio dell'indice del genoma in `rnaseq.nf`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Input primario
        input: Path

        // Archivio del genoma di riferimento
        hisat2_index_zip: Path
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Input primario
        input: Path
    }
    ```

#### 3.1.2. Aggiungere il valore predefinito dell'indice del genoma al profilo di test

Proprio come abbiamo fatto per `reads` nella sezione 1.1.2, aggiungete un valore predefinito per l'indice del genoma al profilo di test in `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

Il parametro è pronto; ora possiamo creare il processo di allineamento.

### 3.2. Scrivere il processo di allineamento e chiamarlo nel workflow

Come prima, dobbiamo riempire la definizione del processo, importare il modulo e aggiungere la chiamata al processo.

#### 3.2.1. Riempire il modulo per il processo di allineamento

Aprite `modules/hisat2_align.nf` ed esaminate la struttura della definizione del processo.

Procedete e riempite la definizione del processo da soli usando le informazioni fornite sopra, quindi verificate il vostro lavoro confrontandolo con la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Allinea le letture a un genoma di riferimento
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Dopo"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Allinea le letture a un genoma di riferimento
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    }
    ```

Questo processo prende due input: le letture e l'archivio dell'indice del genoma.
Il blocco script prima estrae l'indice dall'archivio, poi esegue l'allineamento HISAT2 inviato tramite pipe a `samtools view` per convertire l'output in formato BAM.
L'accessor `simpleName` su `index_zip` estrae il nome base dell'archivio (`genome_index`) da usare come prefisso dell'indice.

#### 3.2.2. Includere il modulo

Aggiornate `rnaseq.nf` per importare il nuovo modulo:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="3"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Successivamente, aggiungeremo la chiamata al processo nel workflow.

#### 3.2.3. Chiamare il processo di allineamento

Le letture trimmate si trovano nel canale `TRIM_GALORE.out.trimmed_reads` prodotto dal passaggio precedente.
Usiamo `#!groovy file(params.hisat2_index_zip)` per fornire l'archivio dell'indice del genoma.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Controllo di qualità iniziale
        FASTQC(read_ch)

        // Trimming degli adapter e QC post-trimming
        TRIM_GALORE(read_ch)

        // Allineamento a un genoma di riferimento
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Crea canale di input da un percorso file
        read_ch = channel.fromPath(params.input)

        // Controllo di qualità iniziale
        FASTQC(read_ch)

        // Trimming degli adapter e QC post-trimming
        TRIM_GALORE(read_ch)
    ```

Il processo di allineamento è ora collegato al workflow.

### 3.3. Aggiornare la gestione dell'output

Dobbiamo aggiungere gli output di allineamento alla dichiarazione publish e configurare dove vanno.

#### 3.3.1. Aggiungere target di pubblicazione per gli output di allineamento

Aggiungete gli output di allineamento alla sezione `publish:`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

Successivamente, dovremo dire a Nextflow dove mettere questi output.

#### 3.3.2. Configurare i nuovi target di output

Aggiungete voci per i target di allineamento nel blocco `output {}`, pubblicandoli in una sottodirectory `align/`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

La configurazione dell'output è completa.

### 3.4. Eseguire il workflow

Il workflow ora include tutti e tre i passaggi di elaborazione: QC, trimming e allineamento.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

Potete trovare gli output di allineamento nella directory results.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Questo completa l'elaborazione di base che dobbiamo applicare a ciascun campione.

_Aggiungeremo l'aggregazione dei report MultiQC nella Parte 3, dopo aver modificato il workflow per accettare più campioni contemporaneamente._

---

### Takeaway

Sapete come racchiudere tutti i passaggi principali per processare campioni RNAseq single-end individualmente.

### Cosa c'è dopo?

Prendetevi una pausa! È stato molto.

Quando vi sentite riposati, passate alla [Parte 3](./03_multi-sample.md), dove imparerete come modificare il workflow per processare più campioni in parallelo, aggregare i report QC attraverso tutti i passaggi per tutti i campioni e abilitare l'esecuzione del workflow su dati RNAseq paired-end.
