# Parte 2: Implementazione per singolo campione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa parte del corso, scriveremo il workflow più semplice possibile che racchiude tutti i comandi eseguiti nella Parte 1 per automatizzarne l'esecuzione, e ci limiteremo a processare un campione alla volta.

Lo faremo in tre fasi:

1. Scrivere un workflow a singolo stadio che esegue il passaggio iniziale di QC
2. Aggiungere il trimming degli adapter e il QC post-trimming
3. Aggiungere l'allineamento al genoma di riferimento

!!! warning "Prerequisito"

    È necessario completare la Parte 1 del corso prima di iniziare questa lezione.
    Nello specifico, il completamento delle sezioni 2.1-3 crea il file di indice del genoma (`data/genome_index.tar.gz`) richiesto per il passaggio di allineamento in questa lezione.

---

## 1. Scrivere un workflow a singolo stadio che esegue il QC iniziale

Iniziamo scrivendo un semplice workflow che esegue lo strumento FastQC su un file FASTQ contenente letture RNAseq single-end.

Forniamo un file di workflow, `rnaseq.nf`, che delinea le parti principali del workflow.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Istruzioni INCLUDE dei moduli

/*
 * Parametri della pipeline
 */

// Input primario

workflow {

    // Crea canale di input

    // Chiama i processi

}
```

Tenere presente che questo codice del workflow è corretto ma non è funzionale; il suo scopo è solo quello di servire come scheletro da utilizzare per scrivere il workflow effettivo.

### 1.1. Creare una directory per memorizzare i moduli

Creeremo moduli standalone per ciascun processo per facilitarne la gestione e il riutilizzo, quindi creiamo una directory per memorizzarli.

```bash
mkdir modules
```

### 1.2. Creare un modulo per il processo di raccolta delle metriche QC

Creiamo un file modulo chiamato `modules/fastqc.nf` per contenere il processo `FASTQC`:

```bash
touch modules/fastqc.nf
```

Aprire il file nell'editor di codice e copiarvi il seguente codice:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

Dovrebbe riconoscere tutti gli elementi da quanto appreso nella Parte 1 e nella Parte 2 di questa serie di formazione; l'unica modifica degna di nota è che questa volta stiamo usando `mode: symlink` per la direttiva `publishDir`, e stiamo usando un parametro per definire la `publishDir`.

!!! note "Nota"

    Anche se i file di dati che stiamo utilizzando qui sono molto piccoli, in genomica possono diventare molto grandi. Ai fini della dimostrazione nell'ambiente di formazione, stiamo usando la modalità di pubblicazione 'symlink' per evitare copie di file non necessarie. Non dovrebbe farlo nei workflow finali, poiché perderà i risultati quando pulisce la directory `work`.

### 1.3. Importare il modulo nel file del workflow

Aggiungere l'istruzione `include { FASTQC } from './modules/fastqc.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Istruzioni INCLUDE dei moduli
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Aggiungere una dichiarazione di input

Dichiarare un parametro di input con un valore predefinito:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Input primario
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Creare un canale di input nel blocco workflow

Utilizzare un factory di canale `.fromPath()` di base per creare il canale di input:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Crea canale di input da un percorso file
    read_ch = channel.fromPath(params.reads)

    // Chiama i processi

}
```

### 1.6. Chiamare il processo `FASTQC` sul canale di input

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Crea canale di input da un percorso file
    read_ch = channel.fromPath(params.reads)

    // Controllo di qualità iniziale
    FASTQC(read_ch)

}
```

### 1.7. Eseguire il workflow per verificare che funzioni

Potremmo usare il parametro `--reads` per specificare un input dalla riga di comando, ma durante lo sviluppo possiamo essere pigri e semplicemente usare il test predefinito che abbiamo configurato.

```bash
nextflow run rnaseq.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Questo dovrebbe essere eseguito molto rapidamente se avete completato la Parte 1 e ha già scaricato il container.
Se lo avete saltato, Nextflow scaricherà il container per voi; non deve fare nulla perché accada, ma potrebbe dover attendere fino a un minuto.

Può trovare gli output in `results/fastqc` come specificato nel processo `FASTQC` dalla direttiva `publishDir`.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Aggiungere il trimming degli adapter e il controllo di qualità post-trimming

Utilizzeremo il wrapper Trim_Galore, che raggruppa Cutadapt per il trimming stesso e FastQC per il controllo di qualità post-trimming.

### 2.1. Creare un modulo per il processo di trimming e QC

Creiamo un file modulo chiamato `modules/trim_galore.nf` per contenere il processo `TRIM_GALORE`:

```bash
touch modules/trim_galore.nf
```

Aprire il file nell'editor di codice e copiarvi il seguente codice:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Importare il modulo nel file del workflow

Aggiungere l'istruzione `include { TRIM_GALORE } from './modules/trim_galore.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Istruzioni INCLUDE dei moduli
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Chiamare il processo sul canale di input

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Crea canale di input da un percorso file
    read_ch = channel.fromPath(params.reads)

    // Controllo di qualità iniziale
    FASTQC(read_ch)

    // Trimming degli adapter e QC post-trimming
    TRIM_GALORE(read_ch)
}
```

### 2.4. Eseguire il workflow per verificare che funzioni

```bash
nextflow run rnaseq.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

Anche questo dovrebbe essere eseguito molto rapidamente, dato che stiamo lavorando su un file di input così piccolo.

Può trovare gli output in `results/trimming` come specificato nel processo `TRIM_GALORE` dalla direttiva `publishDir`.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Allineare le letture al genoma di riferimento

Infine possiamo eseguire il passaggio di allineamento del genoma utilizzando Hisat2, che emetterà anche metriche di controllo della qualità in stile FastQC.

### 3.1. Creare un modulo per il processo HiSat2

Creiamo un file modulo chiamato `modules/hisat2_align.nf` per contenere il processo `HISAT2_ALIGN`:

```bash
touch modules/hisat2_align.nf
```

Aprire il file nell'editor di codice e copiarvi il seguente codice:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Importare il modulo nel file del workflow

Aggiungere l'istruzione `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Istruzioni INCLUDE dei moduli
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Aggiungere una dichiarazione di parametro per fornire l'indice del genoma

Dichiarare un parametro di input con un valore predefinito:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Input primario
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Archivio del genoma di riferimento
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Chiamare il processo `HISAT2_ALIGN` sulle letture trimmate prodotte da `TRIM_GALORE`

Le letture trimmate si trovano nel canale `TRIM_GALORE.out.trimmed_reads` prodotto dal passaggio precedente.

Inoltre, utilizziamo `file (params.hisat2_index_zip)` per fornire allo strumento Hisat2 il tarball compresso dell'indice del genoma.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Crea canale di input da un percorso file
    read_ch = channel.fromPath(params.reads)

    // Controllo di qualità iniziale
    FASTQC(read_ch)

    // Trimming degli adapter e QC post-trimming
    TRIM_GALORE(read_ch)

    // Allineamento a un genoma di riferimento
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Eseguire il workflow per verificare che funzioni

```bash
nextflow run rnaseq.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Può trovare gli output in `results/align` come specificato nel processo `HISAT2_ALIGN` dalla direttiva `publishDir`.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Questo completa l'elaborazione di base che dobbiamo applicare a ciascun campione.

_Aggiungeremo l'aggregazione dei report MultiQC nella Parte 2, dopo aver modificato il workflow per accettare più campioni contemporaneamente._

---

### Takeaway

Sa come racchiudere tutti i passaggi principali per processare campioni RNAseq single-end individualmente.

### Qual è il prossimo passo?

Imparate come modificare il workflow per processare più campioni in parallelo, aggregare i report QC attraverso tutti i passaggi per tutti i campioni e abilitare l'esecuzione del workflow su dati RNAseq paired-end.
