# Parte 2: Variant calling per campione

Nella Parte 1, avete testato manualmente i comandi di Samtools e GATK nei rispettivi container.
Ora andremo ad integrare quegli stessi comandi in un flusso di lavoro Nextflow.

## Obiettivo

In questa parte del corso, svilupperemo un flusso di lavoro che esegue le seguenti operazioni:

1. Generare un file indice per ciascun file BAM di input utilizzando [Samtools](https://www.htslib.org/)
2. Eseguire GATK HaplotypeCaller su ciascun file BAM di input per generare chiamate di varianti per campione in formato VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Questo replica i passaggi della Parte 1, dove avete eseguito questi comandi manualmente nei container.

Come punto di partenza, vi forniamo un file di flusso di lavoro, `genomics.nf`, che delinea le parti principali del flusso di lavoro, insieme a due file modulo, samtools_index.nf e gatk_haplotypecaller.nf, che descrivono la struttura dei moduli.
Questi file non sono funzionali; il loro scopo è solo quello di servire da struttura per voi da riempire con le parti interessanti del codice.

## Piano della lezione

Per rendere il processo di sviluppo più educativo, abbiamo suddiviso il tutto in quattro passaggi:

1. **Scrivere un flusso di lavoro a singolo stadio che esegue Samtools index su un file BAM.**
   Questo copre la creazione di un modulo, l'importazione e la chiamata in un flusso di lavoro.
2. **Aggiungere un secondo processo per eseguire GATK HaplotypeCaller sul file BAM indicizzato.**
   Questo introduce il concatenamento degli output dei processi agli input e la gestione dei file accessori.
3. **Adattare il flusso di lavoro per essere eseguito su un lotto di campioni.**
   Questo copre l'esecuzione parallela e introduce le tuple per mantenere insieme i file associati.
4. **Far accettare al flusso di lavoro un file di testo contenente un lotto di file di input.**
   Questo dimostra un pattern comune per fornire input in blocco.

Ogni passaggio si concentra su un aspetto specifico dello sviluppo del flusso di lavoro.

---

## 1. Scrivere un flusso di lavoro a singolo stadio che esegue Samtools index su un file BAM

Questo primo passaggio si concentra sulle basi: caricare un file BAM e generare un indice per esso.

Ricordiamo il comando `samtools index` dalla [Parte 1](01_method.md):

```bash
samtools index '<input_bam>'
```

Il comando prende un file BAM come input e produce un file indice `.bai` accanto ad esso.
L'URI del container era `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Andremo ad integrare queste informazioni in Nextflow in tre fasi:

1. Configurare l'input
2. Scrivere il processo di indicizzazione e chiamarlo nel flusso di lavoro
3. Configurare la gestione dell'output

### 1.1. Configurare l'input

Dobbiamo dichiarare un parametro di input, creare un profilo di test per fornire un valore predefinito conveniente e creare un canale di input.

#### 1.1.1. Aggiungere una dichiarazione di parametro di input

Nel file principale del flusso di lavoro `genomics.nf`, sotto la sezione `Pipeline parameters`, dichiarate un parametro CLI chiamato `reads_bam`.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Questo configura il parametro CLI, ma non vogliamo digitare il percorso del file ogni volta che eseguiamo il flusso di lavoro durante lo sviluppo.
Ci sono più opzioni per fornire un valore predefinito; qui usiamo un profilo di test.

#### 1.1.2. Creare un profilo di test con un valore predefinito in `nextflow.config`

Un profilo di test fornisce valori predefiniti convenienti per provare un flusso di lavoro senza specificare input sulla riga di comando.
Questa è una convenzione comune nell'ecosistema Nextflow (vedi [Hello Config](../../hello_nextflow/06_hello_config.md) per maggiori dettagli).

Aggiungete un blocco `profiles` a `nextflow.config` con un profilo `test` che imposta il parametro `reads_bam` a uno dei file BAM di test.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Qui stiamo usando `${projectDir}`, una variabile Nextflow integrata che punta alla directory dove si trova lo script del flusso di lavoro.
Questo facilita il riferimento ai file di dati e ad altre risorse senza codificare percorsi assoluti.

#### 1.1.3. Configurare il canale di input

Nel blocco workflow, create un canale di input dal valore del parametro usando la fabbrica di canali `.fromPath` (come usato in [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Dopo"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Ora dobbiamo creare il processo per eseguire l'indicizzazione su questo input.

### 1.2. Scrivere il processo di indicizzazione e chiamarlo nel flusso di lavoro

Dobbiamo scrivere la definizione del processo nel file modulo, importarlo nel flusso di lavoro usando un'istruzione include e chiamarlo sull'input.

#### 1.2.1. Completare il modulo per il processo di indicizzazione

Aprite `modules/samtools_index.nf` ed esaminate lo schema della definizione del processo.
Dovreste riconoscere gli elementi strutturali principali; in caso contrario, considerate di leggere [Hello Nextflow](../../hello_nextflow/01_hello_world.md) per un ripasso.

Andate avanti e completate la definizione del processo da soli usando le informazioni fornite sopra, poi verificate il vostro lavoro contro la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Dopo"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Genera file indice BAM
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

Una volta completato questo, il processo è completo.
Per usarlo nel flusso di lavoro, dovrete importare il modulo e aggiungere una chiamata al processo.

#### 1.2.2. Includere il modulo

In `genomics.nf`, aggiungete un'istruzione `include` per rendere il processo disponibile al flusso di lavoro:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

Il processo è ora disponibile nell'ambito del flusso di lavoro.

#### 1.2.3. Chiamare il processo di indicizzazione sull'input

Ora aggiungiamo una chiamata a `SAMTOOLS_INDEX` nel blocco workflow, passando il canale di input come argomento.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

Il flusso di lavoro ora carica l'input ed esegue il processo di indicizzazione su di esso.
Successivamente, dobbiamo configurare come viene pubblicato l'output.

### 1.3. Configurare la gestione dell'output

Dobbiamo dichiarare quali output del processo pubblicare e specificare dove devono andare.

#### 1.3.1. Dichiarare un output nella sezione `publish:`

La sezione `publish:` all'interno del blocco workflow dichiara quali output del processo devono essere pubblicati.
Assegnate l'output di `SAMTOOLS_INDEX` a un target denominato `bam_index`.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Ora dobbiamo dire a Nextflow dove mettere l'output pubblicato.

#### 1.3.2. Configurare il target di output nel blocco `output {}`

Il blocco `output {}` si trova fuori dal flusso di lavoro e specifica dove viene pubblicato ciascun target denominato.
Aggiungiamo un target per `bam_index` che pubblica in una sottodirectory `bam/`.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note

    Di default, Nextflow pubblica i file di output come link simbolici, il che evita duplicazioni non necessarie.
    Anche se i file di dati che stiamo usando qui sono molto piccoli, nella genomica possono diventare molto grandi.
    I link simbolici si romperanno quando pulite la vostra directory `work`, quindi per flussi di lavoro di produzione potreste voler sovrascrivere la modalità di pubblicazione predefinita a `'copy'`.

### 1.4. Eseguire il flusso di lavoro

A questo punto, abbiamo un flusso di lavoro di indicizzazione a un passaggio che dovrebbe essere completamente funzionale. Proviamo che funzioni!

Possiamo eseguirlo con `-profile test` per usare il valore predefinito configurato nel profilo di test ed evitare di dover scrivere il percorso sulla riga di comando.

```bash
nextflow run genomics.nf -profile test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Potete verificare che il file indice sia stato generato correttamente guardando nella directory di lavoro o nella directory results.

??? abstract "Directory di lavoro"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Directory results"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

Eccolo!

### Takeaway

Sapete come creare un modulo contenente un processo, importarlo in un flusso di lavoro, chiamarlo con un canale di input e pubblicare i risultati.

### Cosa c'è dopo?

Aggiungere un secondo passaggio che prende l'output del processo di indicizzazione e lo usa per eseguire il variant calling.

---

## 2. Aggiungere un secondo processo per eseguire GATK HaplotypeCaller sul file BAM indicizzato

Ora che abbiamo un indice per il nostro file di input, possiamo passare alla configurazione del passaggio di variant calling.

Ricordiamo il comando `gatk HaplotypeCaller` dalla [Parte 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

Il comando prende un file BAM (`-I`), un genoma di riferimento (`-R`) e un file di intervalli (`-L`), e produce un file VCF (`-O`) insieme al suo indice.
Lo strumento si aspetta anche che l'indice BAM, l'indice di riferimento e il dizionario di riferimento siano co-localizzati con i rispettivi file.
L'URI del container era `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Seguiamo le stesse tre fasi di prima:

1. Configurare gli input
2. Scrivere il processo di variant calling e chiamarlo nel flusso di lavoro
3. Configurare la gestione dell'output

### 2.1. Configurare gli input

Il passaggio di variant calling richiede diversi file di input aggiuntivi.
Dobbiamo dichiarare parametri per essi, aggiungere valori predefiniti al profilo di test e creare variabili per caricarli.

#### 2.1.1. Aggiungere dichiarazioni di parametri per gli input accessori

Poiché il nostro nuovo processo si aspetta una manciata di file aggiuntivi da fornire, aggiungete dichiarazioni di parametri per essi in `genomics.nf` sotto la sezione `Pipeline parameters`:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

Come prima, forniamo valori predefiniti attraverso il profilo di test piuttosto che inline.

#### 2.1.2. Aggiungere valori predefiniti per i file accessori al profilo di test

Proprio come abbiamo fatto per `reads_bam` nella sezione 1.1.2, aggiungete valori predefiniti per i file accessori al profilo di test in `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Ora dobbiamo creare variabili che carichino questi percorsi di file per l'uso nel flusso di lavoro.

#### 2.1.3. Creare variabili per i file accessori

Aggiungete variabili per i percorsi dei file accessori all'interno del blocco workflow:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

La sintassi `file()` dice esplicitamente a Nextflow di gestire questi input come percorsi di file.
Potete saperne di più su questo nella Side Quest [Working with files](../../side_quests/working_with_files.md).

### 2.2. Scrivere il processo di variant calling e chiamarlo nel flusso di lavoro

Dobbiamo scrivere la definizione del processo nel file modulo, importarlo nel flusso di lavoro usando un'istruzione include e chiamarlo sull'input reads più l'output del passaggio di indicizzazione e i file accessori.

#### 2.2.1. Completare il modulo per il processo di variant calling

Aprite `modules/gatk_haplotypecaller.nf` ed esaminate lo schema della definizione del processo.

Andate avanti e completate la definizione del processo da soli usando le informazioni fornite sopra, poi verificate il vostro lavoro contro la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Dopo"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Chiama varianti con GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

Noterete che questo processo ha più input di quanti ne richieda effettivamente il comando GATK.
GATK sa cercare il file indice BAM e i file accessori del genoma di riferimento basandosi su convenzioni di denominazione, ma Nextflow è indipendente dal dominio e non conosce queste convenzioni.
Dobbiamo elencarli esplicitamente in modo che Nextflow li metta in scena nella directory di lavoro a runtime; altrimenti GATK genererà un errore sui file mancanti.

Similmente, elenchiamo esplicitamente il file indice dell'output VCF (`"${input_bam}.vcf.idx"`) in modo che Nextflow ne tenga traccia per i passaggi successivi.
Usiamo la sintassi `emit:` per assegnare un nome a ciascun canale di output, il che diventerà utile quando collegheremo gli output nel blocco publish.

Una volta completato questo, il processo è completo.
Per usarlo nel flusso di lavoro, dovrete importare il modulo e aggiungere una chiamata al processo.

#### 2.2.2. Importare il nuovo modulo

Aggiornate `genomics.nf` per importare il nuovo modulo:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

Il processo è ora disponibile nell'ambito del flusso di lavoro.

#### 2.2.3. Aggiungere la chiamata al processo

Aggiungete la chiamata al processo nel corpo del flusso di lavoro, sotto `main:`:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

Dovreste riconoscere la sintassi `*.out` dalla serie di formazione Hello Nextflow; stiamo dicendo a Nextflow di prendere il canale in output da `SAMTOOLS_INDEX` e collegarlo alla chiamata del processo `GATK_HAPLOTYPECALLER`.

!!! note

    Notate che gli input sono forniti esattamente nello stesso ordine nella chiamata al processo come sono elencati nel blocco input del processo.
    In Nextflow, gli input sono posizionali, il che significa che _dovete_ seguire lo stesso ordine; e naturalmente ci deve essere lo stesso numero di elementi.

### 2.3. Configurare la gestione dell'output

Dobbiamo aggiungere i nuovi output alla dichiarazione publish e configurare dove vanno.

#### 2.3.1. Aggiungere target di pubblicazione per gli output del variant calling

Aggiungete gli output VCF e indice alla sezione `publish:`:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Ora dobbiamo dire a Nextflow dove mettere i nuovi output.

#### 2.3.2. Configurare i nuovi target di output

Aggiungete voci per i target `vcf` e `vcf_idx` nel blocco `output {}`, pubblicando entrambi in una sottodirectory `vcf/`:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

Il VCF e il suo indice sono pubblicati come target separati che entrambi vanno nella sottodirectory `vcf/`.

### 2.4. Eseguire il flusso di lavoro

Eseguite il flusso di lavoro espanso, aggiungendo `-resume` questa volta in modo da non dover eseguire nuovamente il passaggio di indicizzazione.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Ora se guardiamo l'output della console, vediamo i due processi elencati.

Il primo processo è stato saltato grazie alla cache, come previsto, mentre il secondo processo è stato eseguito poiché è completamente nuovo.

Troverete i file di output nella directory results (come link simbolici alla directory di lavoro).

??? abstract "Directory contents"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Se aprite il file VCF, dovreste vedere gli stessi contenuti del file che avete generato eseguendo il comando GATK direttamente nel container.

??? abstract "File contents"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Questo è l'output che ci interessa generare per ciascun campione nel nostro studio.

### Takeaway

Sapete come creare un flusso di lavoro modulare a due passaggi che fa un lavoro di analisi reale ed è capace di gestire idiosincrasie dei formati di file genomici come i file accessori.

### Cosa c'è dopo?

Far gestire al flusso di lavoro più campioni in blocco.

---

## 3. Adattare il flusso di lavoro per essere eseguito su un lotto di campioni

Va tutto bene avere un flusso di lavoro che può automatizzare l'elaborazione su un singolo campione, ma cosa succede se avete 1000 campioni?
Dovete scrivere uno script bash che fa un ciclo attraverso tutti i vostri campioni?

No, per fortuna! Basta fare una piccola modifica al codice e Nextflow gestirà anche quello per voi.

### 3.1. Aggiornare l'input per elencare tre campioni

Per eseguire su più campioni, aggiornate il profilo di test per fornire un array di percorsi di file invece di uno singolo.
Questo è un modo rapido per testare l'esecuzione multi-campione; nel prossimo passaggio passeremo a un approccio più scalabile usando un file di input.

Per prima cosa, commentate l'annotazione di tipo nella dichiarazione del parametro, poiché gli array non possono usare dichiarazioni tipizzate:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Poi aggiornate il profilo di test per elencare tutti e tre i campioni:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

La fabbrica di canali nel corpo del flusso di lavoro (`.fromPath`) accetta più percorsi di file così come uno singolo, quindi non sono necessarie altre modifiche.

### 3.2. Eseguire il flusso di lavoro

Provate a eseguire il flusso di lavoro ora che l'infrastruttura è configurata per essere eseguita su tutti e tre i campioni di test.

```bash
nextflow run genomics.nf -profile test -resume
```

Cosa curiosa: questo _potrebbe funzionare_, OPPURE _potrebbe fallire_. Ad esempio, ecco un'esecuzione che è riuscita:

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Se l'esecuzione del vostro flusso di lavoro è riuscita, eseguitelo nuovamente finché non ottenete un errore come questo:

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Se guardate l'output di errore del comando GATK, ci sarà una riga come questa:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Beh, è strano, considerando che abbiamo esplicitamente indicizzato i file BAM nel primo passaggio del flusso di lavoro. Potrebbe esserci qualcosa di sbagliato nell'infrastruttura?

### 3.3. Risolvere il problema

Ispezioniamo le directory di lavoro e usiamo l'operatore `view()` per capire cosa è andato storto.

#### 3.3.1. Verificare le directory di lavoro per le chiamate pertinenti

Date un'occhiata all'interno della directory di lavoro per la chiamata al processo `GATK_HAPLOTYPECALLER` fallita elencata nell'output della console.

??? abstract "Directory contents"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Prestate particolare attenzione ai nomi del file BAM e dell'indice BAM che sono elencati in questa directory: `reads_son.bam` e `reads_father.bam.bai`.

Ma che diavolo? Nextflow ha messo in scena un file indice nella directory di lavoro di questa chiamata al processo, ma è quello sbagliato. Come potrebbe essere successo?

#### 3.3.2. Usare l'[operatore view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) per ispezionare i contenuti del canale

Aggiungete queste due righe nel corpo del flusso di lavoro prima della chiamata al processo `GATK_HAPLOTYPECALLER` per visualizzare i contenuti del canale:

=== "Dopo"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "Prima"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

Poi eseguite nuovamente il comando del flusso di lavoro.

```bash
nextflow run genomics.nf -profile test
```

Ancora una volta, questo può riuscire o fallire. Ecco come appare l'output delle due chiamate `.view()` per un'esecuzione fallita:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Le prime tre righe corrispondono al canale di input e le seconde, al canale di output.
Potete vedere che i file BAM e i file indice per i tre campioni non sono elencati nello stesso ordine!

!!! note

    Quando chiamate un processo Nextflow su un canale contenente più elementi, Nextflow cercherà di parallelizzare l'esecuzione il più possibile e raccoglierà gli output in qualsiasi ordine diventino disponibili.
    La conseguenza è che gli output corrispondenti possono essere raccolti in un ordine diverso da quello in cui sono stati forniti gli input originali.

Come scritto attualmente, il nostro script del flusso di lavoro presuppone che i file indice escano dal passaggio di indicizzazione elencati nello stesso ordine madre/padre/figlio in cui sono stati dati gli input.
Ma questo non è garantito essere il caso, motivo per cui a volte (anche se non sempre) i file sbagliati vengono abbinati nel secondo passaggio.

Per risolvere questo, dobbiamo assicurarci che i file BAM e i loro file indice viaggino insieme attraverso i canali.

!!! tip

    Le istruzioni `view()` nel codice del flusso di lavoro non fanno nulla, quindi non è un problema lasciarle.
    Tuttavia riempiranno il vostro output della console, quindi raccomandiamo di rimuoverle quando avete finito di risolvere il problema.

### 3.4. Aggiornare il flusso di lavoro per gestire correttamente i file indice

La soluzione è raggruppare ciascun file BAM con il suo indice in una tupla, poi aggiornare il processo downstream e l'infrastruttura del flusso di lavoro per corrispondere.

#### 3.4.1. Cambiare l'output del modulo SAMTOOLS_INDEX in una tupla

Il modo più semplice per garantire che un file BAM e il suo indice rimangano strettamente associati è impacchettarli insieme in una tupla in uscita dall'attività di indicizzazione.

!!! note

    Una **tupla** è una lista finita e ordinata di elementi che è comunemente usata per restituire più valori da una funzione. Le tuple sono particolarmente utili per passare più input o output tra processi preservando la loro associazione e ordine.

Aggiornate l'output in `modules/samtools_index.nf` per includere il file BAM:

=== "Dopo"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Prima"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

In questo modo, ciascun file indice sarà strettamente accoppiato con il suo file BAM originale, e l'output complessivo del passaggio di indicizzazione sarà un singolo canale contenente coppie di file.

#### 3.4.2. Cambiare l'input del modulo GATK_HAPLOTYPECALLER per accettare una tupla

Poiché abbiamo cambiato la 'forma' dell'output del primo processo, dobbiamo aggiornare la definizione dell'input del secondo processo per corrispondere.

Aggiornate `modules/gatk_haplotypecaller.nf`:

=== "Dopo"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Prima"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Ora dobbiamo aggiornare il flusso di lavoro per riflettere la nuova struttura della tupla nella chiamata al processo e nei target di pubblicazione.

#### 3.4.3. Aggiornare la chiamata a GATK_HAPLOTYPECALLER nel flusso di lavoro

Non abbiamo più bisogno di fornire l'originale `reads_ch` al processo `GATK_HAPLOTYPECALLER`, poiché il file BAM è ora raggruppato nel canale di output da `SAMTOOLS_INDEX`.

Aggiornate la chiamata in `genomics.nf`:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Infine, dobbiamo aggiornare i target di pubblicazione per riflettere la nuova struttura di output.

#### 3.4.4. Aggiornare il target di pubblicazione per l'output del BAM indicizzato

Poiché l'output di SAMTOOLS_INDEX è ora una tupla contenente sia il file BAM che il suo indice, rinominate il target di pubblicazione da `bam_index` a `indexed_bam` per riflettere meglio i suoi contenuti:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Con queste modifiche, il BAM e il suo indice sono garantiti viaggiare insieme, quindi l'abbinamento sarà sempre corretto.

### 3.5. Eseguire il flusso di lavoro corretto

Eseguite nuovamente il flusso di lavoro per assicurarvi che questo funzionerà in modo affidabile andando avanti.

```bash
nextflow run genomics.nf -profile test
```

Questa volta (e ogni volta) tutto dovrebbe funzionare correttamente:

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

La directory results ora contiene sia file BAM che BAI per ciascun campione (dalla tupla), insieme agli output VCF:

??? abstract "Results directory contents"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

Raggruppando i file associati in tuple, abbiamo garantito che i file corretti viaggino sempre insieme attraverso il flusso di lavoro.
Il flusso di lavoro ora elabora qualsiasi numero di campioni in modo affidabile, ma elencarli individualmente nella configurazione non è molto scalabile.
Nel prossimo passaggio, passeremo alla lettura degli input da un file.

### Takeaway

Sapete come far eseguire il vostro flusso di lavoro su più campioni (indipendentemente).

### Cosa c'è dopo?

Rendere più facile gestire i campioni in blocco.

---

## 4. Far accettare al flusso di lavoro un file di testo contenente un lotto di file di input

Un modo molto comune per fornire più file di dati di input a un flusso di lavoro è farlo con un file di testo contenente i percorsi dei file.
Può essere semplice come un file di testo che elenca un percorso di file per riga e nient'altro, oppure il file può contenere metadati aggiuntivi, nel qual caso è spesso chiamato samplesheet.

Qui vi mostreremo come fare il caso semplice.

### 4.1. Esaminare il file di testo fornito che elenca i percorsi dei file di input

Abbiamo già creato un file di testo che elenca i percorsi dei file di input, chiamato `sample_bams.txt`, che potete trovare nella directory `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Come potete vedere, abbiamo elencato un percorso di file per riga, e sono percorsi assoluti.

!!! note

    I file che stiamo usando qui sono solo sul filesystem locale dei vostri GitHub Codespaces, ma potremmo anche puntare a file nello storage cloud.
    Se non state usando l'ambiente Codespaces fornito, potrebbe essere necessario adattare i percorsi dei file per corrispondere alla vostra configurazione locale.

### 4.2. Aggiornare il parametro e il profilo di test

Cambiate il parametro `reads_bam` per puntare al file `sample_bams.txt` invece di elencare campioni individuali.

Ripristinate l'annotazione di tipo nel blocco params (poiché è di nuovo un percorso singolo):

=== "Dopo"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

Poi aggiornate il profilo di test per puntare al file di testo:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

L'elenco dei file non vive più nel codice, il che è un grande passo nella direzione giusta.

### 4.3. Aggiornare la fabbrica di canali per leggere righe da un file

Attualmente, la nostra fabbrica di canali di input tratta tutti i file che le diamo come gli input di dati che vogliamo fornire al processo di indicizzazione.
Poiché ora gli stiamo dando un file che elenca i percorsi dei file di input, dobbiamo cambiare il suo comportamento per analizzare il file e trattare i percorsi dei file che contiene come gli input di dati.

Possiamo farlo usando lo stesso pattern che abbiamo usato nella [Parte 2 di Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): applicando l'operatore [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) per analizzare il file, poi un'operazione `map` per selezionare il primo campo di ciascuna riga.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Tecnicamente potremmo farlo più semplicemente usando l'operatore [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext), poiché il nostro file di input attualmente contiene solo percorsi di file.
Tuttavia, usando il più versatile operatore `splitCsv` (integrato da `map`), possiamo rendere il nostro flusso di lavoro a prova di futuro nel caso decidessimo di aggiungere metadati al file contenente i percorsi dei file.

!!! tip

    Se non siete sicuri di aver capito cosa stanno facendo gli operatori qui, questa è un'altra grande opportunità per usare l'operatore `.view()` per guardare come appaiono i contenuti del canale prima e dopo averli applicati.

### 4.4. Eseguire il flusso di lavoro

Eseguite il flusso di lavoro ancora una volta. Questo dovrebbe produrre lo stesso risultato di prima, giusto?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Sì! Infatti, Nextflow rileva correttamente che le chiamate ai processi sono esattamente le stesse, e non si preoccupa nemmeno di rieseguire tutto, poiché stavamo eseguendo con `-resume`.

E questo è tutto! Il nostro semplice flusso di lavoro di variant calling ha tutte le caratteristiche di base che volevamo.

### Takeaway

Sapete come creare un flusso di lavoro modulare multi-passaggio per indicizzare un file BAM e applicare variant calling per campione usando GATK.

Più in generale, avete imparato come usare componenti e logica essenziali di Nextflow per costruire una semplice pipeline genomica che fa un lavoro reale, tenendo conto delle idiosincrasie dei formati di file genomici e dei requisiti degli strumenti.

### Cosa c'è dopo?

Celebrate il vostro successo e prendetevi una pausa extra lunga!

Nella prossima parte di questo corso, imparerete come trasformare questo semplice flusso di lavoro di variant calling per campione per applicare il joint variant calling ai dati.
