# Parte 1: Chiamata di varianti per campione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella prima parte di questo corso, Le mostreremo come costruire una semplice pipeline di chiamata di varianti che applica la chiamata di varianti GATK a singoli campioni di sequenziamento.

### Panoramica del metodo

La chiamata di varianti è un metodo di analisi genomica che mira a identificare variazioni in una sequenza genomica rispetto a un genoma di riferimento.
Qui utilizzeremo strumenti e metodi progettati per chiamare varianti brevi, _cioè_ SNP e indel.

![GATK pipeline](img/gatk-pipeline.png)

Una pipeline completa di chiamata di varianti coinvolge tipicamente molti passaggi, incluso il mapping al riferimento (a volte chiamato allineamento del genoma) e il filtraggio e la prioritizzazione delle varianti.
Per semplicità, in questa parte del corso ci concentreremo solo sulla parte di chiamata delle varianti.

### Dataset

Forniamo i seguenti dati e risorse correlate:

- **Un genoma di riferimento** costituito da una piccola regione del cromosoma umano 20 (da hg19/b37) e i suoi file accessori (indice e dizionario di sequenza).
- **Tre campioni di sequenziamento dell'intero genoma** corrispondenti a un trio familiare (madre, padre e figlio), che sono stati ridotti a una piccola porzione di dati sul cromosoma 20 per mantenere le dimensioni dei file ridotte.
  Si tratta di dati di sequenziamento Illumina a lettura breve che sono già stati mappati al genoma di riferimento, forniti in formato [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, una versione compressa di SAM, Sequence Alignment Map).
- **Un elenco di intervalli genomici**, cioè coordinate sul genoma dove i nostri campioni hanno dati adatti per chiamare varianti, fornito in formato BED.

### Workflow

In questa parte del corso, svilupperemo un workflow che fa quanto segue:

1. Generare un file indice per ciascun file BAM di input utilizzando [Samtools](https://www.htslib.org/)
2. Eseguire GATK HaplotypeCaller su ciascun file BAM di input per generare chiamate di varianti per campione in VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note "Nota"

    I file indice sono una caratteristica comune dei formati di file bioinformatici; contengono informazioni sulla struttura del file principale che consente a strumenti come GATK di accedere a un sottoinsieme dei dati senza dover leggere l'intero file.
    Questo è importante a causa delle dimensioni che questi file possono raggiungere.

---

## 0. Riscaldamento: Testare i comandi Samtools e GATK in modo interattivo

Prima vogliamo provare i comandi manualmente prima di tentare di includerli in un workflow.
Gli strumenti di cui abbiamo bisogno (Samtools e GATK) non sono installati nell'ambiente GitHub Codespaces, quindi li useremo tramite container (vedere [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Assicuratevi di trovarsi nella directory `nf4-science/genomics` in modo che l'ultima parte del percorso mostrato quando digita `pwd` sia `genomics`.

### 0.1. Indicizzare un file BAM di input con Samtools

Scaricheremo un container Samtools, lo avvieremo in modo interattivo ed eseguiremo il comando `samtools index` su uno dei file BAM.

#### 0.1.1. Scaricare il container Samtools

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

#### 0.1.2. Avviare il container Samtools in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

#### 0.1.3. Eseguire il comando di indicizzazione

La [documentazione di Samtools](https://www.htslib.org/doc/samtools-index.html) ci fornisce la riga di comando da eseguire per indicizzare un file BAM.

Dobbiamo solo fornire il file di input; lo strumento genererà automaticamente un nome per l'output aggiungendo `.bai` al nome del file di input.

```bash
samtools index /data/bam/reads_mother.bam
```

Questo dovrebbe completarsi immediatamente, e ora dovrebbe vedere un file chiamato `reads_mother.bam.bai` nella stessa directory del file BAM di input originale.

??? abstract "Contenuto della directory"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Uscire dal container Samtools

```bash
exit
```

### 0.2. Chiamare varianti con GATK HaplotypeCaller

Scaricheremo un container GATK, lo avvieremo in modo interattivo ed eseguiremo il comando `gatk HaplotypeCaller` sul file BAM che abbiamo appena indicizzato.

#### 0.2.1. Scaricare il container GATK

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

#### 0.2.2. Avviare il container GATK in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

#### 0.2.3. Eseguire il comando di chiamata di varianti

La [documentazione di GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) ci fornisce la riga di comando da eseguire per effettuare la chiamata di varianti su un file BAM.

Dobbiamo fornire il file BAM di input (`-I`) così come il genoma di riferimento (`-R`), un nome per il file di output (`-O`) e un elenco di intervalli genomici da analizzare (`-L`).

Tuttavia, non è necessario specificare il percorso del file indice; lo strumento lo cercherà automaticamente nella stessa directory, basandosi sulla convenzione stabilita di denominazione e co-localizzazione.
Lo stesso vale per i file accessori del genoma di riferimento (file indice e dizionario di sequenza, `*.fai` e `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

Il file di output `reads_mother.vcf` viene creato all'interno della vostra directory di lavoro nel container, quindi non lo vedrà nell'explorer di file di VS Code a meno che non modifichi il percorso del file di output.
Tuttavia, è un piccolo file di test, quindi potete usare `cat` per aprirlo e visualizzarne il contenuto.
Se scorrete fino all'inizio del file, troverete un'intestazione composta da molte righe di metadati, seguita da un elenco di chiamate di varianti, una per riga.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Ogni riga descrive una possibile variante identificata nei dati di sequenziamento del campione. Per una guida sull'interpretazione del formato VCF, consulti [questo articolo utile](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Il file VCF di output è accompagnato da un file indice chiamato `reads_mother.vcf.idx` che è stato creato automaticamente da GATK.
Ha la stessa funzione del file indice BAM, per consentire agli strumenti di cercare e recuperare sottoinsiemi di dati senza caricare l'intero file.

#### 0.2.4. Uscire dal container GATK

```bash
exit
```

### Takeaway

Sa come testare i comandi di indicizzazione Samtools e di chiamata di varianti GATK nei rispettivi container.

### Cosa c'è dopo?

Imparare come includere quegli stessi comandi in un workflow a due passaggi che utilizza container per eseguire il lavoro.

---

## 1. Scrivere un workflow a singolo stadio che esegue Samtools index su un file BAM

Le forniamo un file workflow, `genomics-1.nf`, che delinea le parti principali del workflow.
Non è funzionale; il suo scopo è solo quello di servire come scheletro che utilizzerà per scrivere il workflow effettivo.

### 1.1. Definire il processo di indicizzazione

Iniziamo scrivendo un processo, che chiameremo `SAMTOOLS_INDEX`, che descrive l'operazione di indicizzazione.

```groovy title="genomics-1.nf" linenums="9"
/*
 * Generate BAM index file
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

Dovrebbe riconoscere tutti i pezzi da quello che ha imparato nella Parte 1 e Parte 2 di questa serie di formazione.

Questo processo richiederà che passiamo un percorso di file tramite l'input `input_bam`, quindi configuriamo questo aspetto successivamente.

### 1.2. Aggiungere una dichiarazione di parametro di input

All'inizio del file, sotto la sezione `Pipeline parameters`, dichiariamo un parametro CLI chiamato `reads_bam` e gli diamo un valore predefinito.
In questo modo, possiamo essere pigri e non specificare l'input quando digitiamo il comando per avviare la pipeline (a scopo di sviluppo).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Pipeline parameters
 */
params {
    // Input primario
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Ora abbiamo un processo pronto, così come un parametro per dargli un input su cui eseguire, quindi colleghiamo queste cose insieme.

!!! note "Nota"

    `${projectDir}` è una variabile Nextflow integrata che punta alla directory dove si trova lo script workflow Nextflow corrente (`genomics-1.nf`).

    Questo facilita il riferimento a file, directory di dati e altre risorse incluse nel repository del workflow senza codificare percorsi assoluti.

### 1.3. Aggiungere il blocco workflow per eseguire SAMTOOLS_INDEX

Nel blocco `workflow`, dobbiamo configurare un **channel** per alimentare l'input al processo `SAMTOOLS_INDEX`; quindi possiamo chiamare il processo stesso per eseguirlo sul contenuto di quel channel.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Crea canale di input (singolo file tramite parametro CLI)
    reads_ch = channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

Il blocco workflow ha due sezioni:

- `main:` contiene le operazioni sui channel e le chiamate ai processi
- `publish:` dichiara quali output devono essere pubblicati, assegnandoli a target nominati

Noterà che stiamo usando lo stesso channel factory `.fromPath` che abbiamo usato in [Hello Channels](../../hello_nextflow/02_hello_channels.md).
Infatti, stiamo facendo qualcosa di molto simile.
La differenza è che stiamo dicendo a Nextflow di caricare semplicemente il percorso del file stesso nel channel come elemento di input, piuttosto che leggerne il contenuto.

### 1.4. Aggiungere un blocco output per definire dove vengono pubblicati i risultati

Dopo il blocco workflow, aggiungiamo un blocco `output` che specifica dove pubblicare gli output del workflow.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

Ogni target nominato dalla sezione `publish:` (come `bam_index`) ottiene il proprio blocco dove è possibile configurare il percorso di output relativo alla directory di output base.

!!! note "Nota"

    Anche se i file di dati che stiamo usando qui sono molto piccoli, nella genomica possono diventare molto grandi.
    Per impostazione predefinita, Nextflow crea collegamenti simbolici ai file di output nella directory di pubblicazione, il che evita copie di file non necessarie.
    È possibile modificare questo comportamento utilizzando l'opzione `mode` (ad esempio, `mode 'copy'`) per creare copie effettive invece.
    Sia consapevole che i collegamenti simbolici si romperanno quando pulirà la vostra directory `work`, quindi per i workflow di produzione potrebbe voler usare `mode 'copy'`.

### 1.5. Configurare la directory di output

La directory di output base viene impostata tramite l'opzione di configurazione `outputDir`. La aggiunga a `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Prima"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Eseguire il workflow per verificare che il passaggio di indicizzazione funzioni

Eseguiamo il workflow! Come promemoria, non è necessario specificare un input nella riga di comando perché abbiamo impostato un valore predefinito per l'input quando abbiamo dichiarato il parametro di input.

```bash
nextflow run genomics-1.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Può verificare che il file indice sia stato generato correttamente guardando nella directory di lavoro o nella directory dei risultati.

??? abstract "Contenuto della directory di lavoro"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Contenuto della directory dei risultati"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

Eccolo!

### Takeaway

Sa come includere uno strumento genomico in un workflow Nextflow a passaggio singolo e farlo eseguire utilizzando un container.

### Cosa c'è dopo?

Aggiungere un secondo passaggio che consuma l'output del primo.

---

## 2. Aggiungere un secondo processo per eseguire GATK HaplotypeCaller sul file BAM indicizzato

Ora che abbiamo un indice per il nostro file di input, possiamo passare alla configurazione del passaggio di chiamata di varianti, che è la parte interessante del workflow.

### 2.1. Definire il processo di chiamata di varianti

Scriviamo un processo, che chiameremo `GATK_HAPLOTYPECALLER`, che descrive l'operazione di chiamata di varianti.

```groovy title="genomics-1.nf" linenums="44"
/*
 * Call variants with GATK HaplotypeCaller
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

Noterà che abbiamo introdotto una nuova sintassi qui (`emit:`) per nominare in modo univoco ciascuno dei nostri channel di output, e le ragioni di questo diventeranno presto chiare.

Questo comando richiede molti più input, perché GATK ha bisogno di più informazioni per eseguire l'analisi rispetto a un semplice lavoro di indicizzazione.
Ma noterà che ci sono ancora più input definiti nel blocco degli input rispetto a quelli elencati nel comando GATK. Perché?

!!! note "Nota"

    GATK sa cercare il file indice BAM e i file accessori del genoma di riferimento perché è consapevole delle convenzioni relative a quei file.
    Tuttavia, Nextflow è progettato per essere indipendente dal dominio e non sa nulla dei requisiti dei formati di file bioinformatici.

Dobbiamo dire esplicitamente a Nextflow che deve preparare quei file nella directory di lavoro in fase di esecuzione; altrimenti non lo farà, e GATK (correttamente) genererà un errore riguardo alla mancanza dei file indice.

Allo stesso modo, dobbiamo elencare esplicitamente il file indice del VCF di output (il file `"${input_bam}.vcf.idx"`) in modo che Nextflow sappia di tenere traccia di quel file nel caso sia necessario nei passaggi successivi.

### 2.2. Aggiungere definizioni per gli input accessori

Poiché il nostro nuovo processo si aspetta che vengano forniti alcuni file aggiuntivi, configuriamo alcuni parametri CLI per essi nella sezione `Pipeline parameters`, insieme ad alcuni valori predefiniti (per le stesse ragioni di prima).

```groovy title="genomics-1.nf" linenums="8"
    // Accessory files
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Creare variabili per contenere i percorsi dei file accessori

Mentre gli input di dati principali vengono trasmessi dinamicamente attraverso i channel, ci sono due approcci per gestire i file accessori. L'approccio raccomandato è creare channel espliciti, il che rende il flusso di dati più chiaro e coerente. In alternativa, la funzione file() per creare variabili può essere utilizzata per casi più semplici, in particolare quando è necessario fare riferimento allo stesso file in più processi - anche se sia consapevole che questo crea comunque channel implicitamente.

Aggiunga questo al blocco workflow (dopo la creazione di `reads_ch`, all'interno della sezione `main:`):

```groovy title="genomics-1.nf" linenums="79"
    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

Questo renderà disponibili i percorsi dei file accessori per essere forniti come input a qualsiasi processo che ne abbia bisogno.

### 2.4. Aggiungere una chiamata al blocco workflow per eseguire GATK_HAPLOTYPECALLER

Ora che abbiamo configurato il nostro secondo processo e tutti gli input e i file accessori sono pronti e disponibili, possiamo aggiungere una chiamata al processo `GATK_HAPLOTYPECALLER` nel corpo del workflow.

```groovy title="genomics-1.nf" linenums="88"
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

Dovrebbe riconoscere la sintassi `*.out` dalla Parte 1 di questa serie di formazione; stiamo dicendo a Nextflow di prendere l'output del channel da `SAMTOOLS_INDEX` e collegarlo alla chiamata del processo `GATK_HAPLOTYPECALLER`.

!!! note "Nota"

    Noterà che gli input vengono forniti esattamente nello stesso ordine nella chiamata al processo come sono elencati nel blocco input del processo.
    In Nextflow, gli input sono posizionali, il che significa che _deve_ seguire lo stesso ordine; e ovviamente ci devono essere lo stesso numero di elementi.

### 2.5. Aggiornare la sezione publish e il blocco output

Dobbiamo aggiornare la sezione `publish:` per includere gli output VCF, e aggiungere target corrispondenti nel blocco `output`.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Eseguire il workflow per verificare che il passaggio di chiamata di varianti funzioni

Eseguiamo il workflow espanso con `-resume` in modo da non dover eseguire nuovamente il passaggio di indicizzazione.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Ora se guardiamo l'output della console, vediamo i due processi elencati.

Il primo processo è stato saltato grazie alla memorizzazione nella cache, come previsto, mentre il secondo processo è stato eseguito poiché è completamente nuovo.

Troverà i file di output nella directory dei risultati (come collegamenti simbolici alla directory di lavoro).

??? abstract "Contenuto della directory"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

Se apre il file VCF, dovrebbe vedere lo stesso contenuto del file che ha generato eseguendo il comando GATK direttamente nel container.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Questo è l'output che ci interessa generare per ciascun campione nel nostro studio.

### Takeaway

Sa come creare un workflow a due passaggi molto semplice che fa un lavoro di analisi reale ed è in grado di gestire le idiosincrasie dei formati di file genomici come i file accessori.

### Cosa c'è dopo?

Fare in modo che il workflow gestisca più campioni in blocco.

---

## 3. Adattare il workflow per eseguirlo su un batch di campioni

È tutto bello avere un workflow che può automatizzare l'elaborazione su un singolo campione, ma cosa succede se avete 1000 campioni?
Ha bisogno di scrivere uno script bash che scorra attraverso tutti i vostri campioni?

No, per fortuna! Basta fare una piccola modifica al codice e Nextflow gestirà anche questo per voi.

### 3.1. Trasformare la dichiarazione del parametro di input in un array che elenca i tre campioni

Trasformiamo quel percorso di file predefinito nella dichiarazione del file BAM di input in un array che elenca i percorsi dei file per i nostri tre campioni di test, nella sezione `Pipeline parameters`.

=== "Dopo"

    ```groovy title="genomics-1.nf" linenums="7"
    // Input primario (array di tre campioni)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Prima"

    ```groovy title="genomics-1.nf" linenums="7"
        // Input primario
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note "Nota"

    Quando si utilizzano dichiarazioni di parametri tipizzate (come `reads_bam: Path`), non è possibile assegnare un valore array.
    Per gli array, ometta l'annotazione del tipo.

E questo è in realtà tutto ciò che dobbiamo fare, perché il channel factory che usiamo nel corpo del workflow (`.fromPath`) è altrettanto felice di accettare più percorsi di file da caricare nel channel di input quanto lo era di caricarne uno solo.

!!! note "Nota"

    Normalmente, non vorrebbe codificare l'elenco dei campioni nel vostro file workflow, ma lo stiamo facendo qui per mantenere le cose semplici.
    Presenteremo modi più eleganti per gestire gli input più avanti in questa serie di formazione.

### 3.2. Eseguire il workflow per verificare che venga eseguito su tutti e tre i campioni

Proviamo a eseguire il workflow ora che l'impianto è configurato per eseguire tutti e tre i campioni di test.

```bash
nextflow run genomics-1.nf -resume
```

Cosa divertente: questo _potrebbe funzionare_, OPPURE _potrebbe fallire_. Ad esempio, ecco un'esecuzione riuscita:

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Se l'esecuzione del vostro workflow è riuscita, eseguitela di nuovo finché non ottenete un errore come questo:

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

Se guarda l'output di errore del comando GATK, ci sarà una riga come questa:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Beh, è strano, considerando che abbiamo esplicitamente indicizzato i file BAM nel primo passaggio del workflow. Potrebbe esserci qualcosa di sbagliato nell'impianto?

#### 3.2.1. Controllare le directory di lavoro per le chiamate rilevanti

Diamo un'occhiata all'interno della directory di lavoro per la chiamata al processo `GATK_HAPLOTYPECALLER` fallita elencata nell'output della console.

??? abstract "Contenuto della directory"

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

Presti particolare attenzione ai nomi del file BAM e dell'indice BAM elencati in questa directory: `reads_son.bam` e `reads_father.bam.bai`.

Che diavolo? Nextflow ha preparato un file indice nella directory di lavoro di questa chiamata al processo, ma è quello sbagliato. Come è potuto succedere?

#### 3.2.2. Usare l'[operatore view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) per ispezionare il contenuto del channel

Aggiunga queste due righe nel corpo del workflow prima della chiamata al processo `GATK_HAPLOTYPER`:

```groovy title="genomics-1.nf" linenums="84"
    // temporary diagnostics
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Quindi eseguite di nuovo il comando del workflow.

```bash
nextflow run genomics-1.nf
```

Ancora una volta, questo può avere successo o fallire. Ecco un'esecuzione riuscita:

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

Ed ecco una fallita:

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Potrebbe dover eseguire più volte perché fallisca di nuovo.
Questo errore non si riprodurrà in modo coerente perché dipende da una certa variabilità nei tempi di esecuzione delle singole chiamate ai processi.

Questo è l'aspetto dell'output delle due chiamate `.view()` che abbiamo aggiunto per un'esecuzione fallita:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Le prime tre righe corrispondono al channel di input e la seconda, al channel di output.
Può vedere che i file BAM e i file indice per i tre campioni non sono elencati nello stesso ordine!

!!! note "Nota"

    Quando chiama un processo Nextflow su un channel contenente più elementi, Nextflow tenterà di parallelizzare l'esecuzione il più possibile e raccoglierà gli output in qualsiasi ordine diventino disponibili.
    La conseguenza è che gli output corrispondenti possono essere raccolti in un ordine diverso rispetto a quello in cui sono stati forniti gli input originali.

Come attualmente scritto, il nostro script workflow presume che i file indice escano dal passaggio di indicizzazione elencati nello stesso ordine madre/padre/figlio degli input forniti.
Ma questo non è garantito, motivo per cui a volte (anche se non sempre) i file sbagliati vengono accoppiati nel secondo passaggio.

Per risolvere questo problema, dobbiamo assicurarci che i file BAM e i loro file indice viaggino insieme attraverso i channel.

!!! tip "Suggerimento"

    Le istruzioni `view()` nel codice del workflow non fanno nulla, quindi non è un problema lasciarle.
    Tuttavia, ingombreranno l'output della vostra console, quindi vi consigliamo di rimuoverle quando ha finito di risolvere il problema.

### 3.3. Modificare l'output del processo SAMTOOLS_INDEX in una tupla che mantiene insieme il file di input e il suo indice

Il modo più semplice per garantire che un file BAM e il suo indice rimangano strettamente associati è confezionarli insieme in una tupla in uscita dall'attività di indicizzazione.

!!! note "Nota"

    Una **tupla** è una lista ordinata finita di elementi che viene comunemente utilizzata per restituire più valori da una funzione. Le tuple sono particolarmente utili per passare più input o output tra processi preservandone l'associazione e l'ordine.

Prima, modifichiamo l'output del processo `SAMTOOLS_INDEX` per includere il file BAM nella sua dichiarazione di output.

=== "Dopo"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Prima"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

In questo modo, ciascun file indice sarà strettamente accoppiato con il suo file BAM originale, e l'output complessivo del passaggio di indicizzazione sarà un singolo channel contenente coppie di file.

### 3.4. Modificare l'input del processo GATK_HAPLOTYPECALLER per essere una tupla

Poiché abbiamo modificato la 'forma' dell'output del primo processo nel workflow, dobbiamo aggiornare la definizione dell'input del secondo processo per corrispondere.

Specificamente, dove in precedenza dichiaravamo due percorsi di input separati nel blocco input del processo `GATK_HAPLOTYPECALLER`, ora dichiariamo un singolo input che corrisponde alla struttura della tupla emessa da `SAMTOOLS_INDEX`.

=== "Dopo"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Prima"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

Naturalmente, poiché ora abbiamo modificato la forma degli input che `GATK_HAPLOTYPECALLER` si aspetta, dobbiamo aggiornare di conseguenza la chiamata al processo nel corpo del workflow.

### 3.5. Aggiornare la chiamata a GATK_HAPLOTYPECALLER nel blocco workflow

Non è più necessario fornire il `reads_ch` originale al processo `GATK_HAPLOTYPECALLER`, poiché il file BAM è ora incluso nel channel di output da `SAMTOOLS_INDEX`.

Di conseguenza, possiamo semplicemente eliminare quella riga.

=== "Dopo"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Prima"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

Questa è tutta la riconfigurazione necessaria per risolvere il problema di mancata corrispondenza dell'indice.

### 3.6. Aggiornare la sezione publish e il blocco output per la tupla

Poiché `SAMTOOLS_INDEX.out` è ora una tupla contenente sia il BAM che il suo indice, entrambi i file verranno pubblicati insieme.
Rinominiamo il target da `bam_index` a `indexed_bam` per riflettere che ora contiene entrambi i file.

=== "Dopo"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Prima"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

Dobbiamo anche aggiornare il blocco output per usare il nuovo nome del target:

=== "Dopo"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Prima"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Eseguire il workflow per verificare che funzioni correttamente su tutti e tre i campioni ogni volta

Naturalmente, la prova è nel budino, quindi eseguiamo il workflow di nuovo alcune volte per assicurarci che funzionerà in modo affidabile andando avanti.

```bash
nextflow run genomics-1.nf
```

Questa volta (e ogni volta) tutto dovrebbe funzionare correttamente:

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

La directory dei risultati ora contiene sia file BAM che BAI per ciascun campione (dalla tupla), insieme agli output VCF:

??? abstract "Contenuto della directory dei risultati"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

Se desiderate, potete usare `.view()` di nuovo per dare un'occhiata a come appare il contenuto del channel di output di `SAMTOOLS_INDEX`:

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Vedrà che il channel contiene le tre tuple previste (percorsi dei file troncati per leggibilità).

```console title="Output"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

Questo sarà molto più sicuro, andando avanti.

### Takeaway

Sa come fare in modo che il vostro workflow venga eseguito su più campioni (in modo indipendente).

### Cosa c'è dopo?

Rendere più facile gestire campioni in blocco.

---

## 4. Fare in modo che il workflow accetti un file di testo contenente un batch di file di input

Un modo molto comune per fornire più file di dati di input a un workflow è farlo con un file di testo contenente i percorsi dei file.
Può essere semplice come un file di testo che elenca un percorso di file per riga e nient'altro, oppure il file può contenere metadati aggiuntivi, nel qual caso è spesso chiamato samplesheet.

Qui Le mostreremo come fare nel caso semplice.

### 4.1. Esaminare il file di testo fornito che elenca i percorsi dei file di input

Abbiamo già creato un file di testo che elenca i percorsi dei file di input, chiamato `sample_bams.txt`, che potete trovare nella directory `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Come potete vedere, abbiamo elencato un percorso di file per riga, e sono percorsi assoluti.

!!! note "Nota"

    I file che stiamo usando qui sono semplicemente sul filesystem locale del vostro GitHub Codespaces, ma potremmo anche puntare a file nello storage cloud.

### 4.2. Aggiornare il valore predefinito del parametro

Cambiamo il valore predefinito per il nostro parametro di input `reads_bam` per puntare al file `sample_bams.txt`.

=== "Dopo"

    ```groovy title="genomics-1.nf" linenums="7"
        // Input primario (file di file di input, uno per riga)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Prima"

    ```groovy title="genomics-1.nf" linenums="7"
    // Input primario (array di tre campioni)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

In questo modo possiamo continuare a essere pigri, ma l'elenco dei file non vive più nel codice del workflow stesso, il che è un grande passo nella giusta direzione.

### 4.3. Aggiornare il channel factory per leggere righe da un file

Attualmente, il nostro channel factory di input tratta qualsiasi file che gli forniamo come input di dati che vogliamo alimentare al processo di indicizzazione.
Poiché ora gli stiamo dando un file che elenca i percorsi dei file di input, dobbiamo cambiare il suo comportamento per analizzare il file e trattare i percorsi dei file che contiene come input di dati.

Fortunatamente possiamo farlo molto semplicemente, semplicemente aggiungendo l'[operatore `.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) al passaggio di costruzione del channel.

=== "Dopo"

    ```groovy title="genomics-1.nf" linenums="68"
        // Crea canale di input da un file di testo che elenca i percorsi dei file di input
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Prima"

    ```groovy title="genomics-1.nf" linenums="68"
        // Crea canale di input (singolo file tramite parametro CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip "Suggerimento"

    Questa è un'altra grande opportunità per usare l'operatore `.view()` per vedere come appaiono i contenuti del channel prima e dopo l'applicazione di un operatore.

### 4.4. Eseguire il workflow per verificare che funzioni correttamente

Eseguiamo il workflow un'ultima volta. Questo dovrebbe produrre lo stesso risultato di prima, giusto?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Sì! Infatti, Nextflow rileva correttamente che le chiamate ai processi sono esattamente le stesse, e non si preoccupa nemmeno di rieseguire tutto, poiché stavamo eseguendo con `-resume`.

E questo è tutto! Il nostro semplice workflow di chiamata di varianti ha tutte le funzionalità di base che volevamo.

### Takeaway

Sa come creare un workflow lineare a più passaggi per indicizzare un file BAM e applicare la chiamata di varianti per campione utilizzando GATK.

Più in generale, ha imparato come utilizzare componenti e logica essenziali di Nextflow per costruire una semplice pipeline genomica che fa un lavoro reale, tenendo conto delle idiosincrasie dei formati di file genomici e dei requisiti degli strumenti.

### Cosa c'è dopo?

Celebrate il vostro successo e prendetevi una pausa extra lunga!

Nella prossima parte di questo corso, imparerà come utilizzare alcune funzionalità aggiuntive di Nextflow (inclusi più operatori di channel) per applicare la chiamata congiunta di varianti ai dati.
