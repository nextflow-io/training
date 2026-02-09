# Parte 2: Implementazione per singolo campione

In questa parte del corso, scriveremo il flusso di lavoro più semplice possibile che racchiude tutti i comandi eseguiti nella Parte 1 per automatizzarne l'esecuzione, e ci concentreremo sull'elaborazione di un campione alla volta.

Lo faremo in tre fasi:

1. Scrivere un flusso di lavoro a singolo stadio che esegue il controllo qualità iniziale
2. Aggiungere il trimming degli adattatori e il controllo qualità post-trimming
3. Aggiungere l'allineamento al genoma di riferimento

!!! warning "Prerequisito"

    Dovete completare la Parte 1 del corso prima di iniziare questa lezione.
    In particolare, completando le sezioni 2.1-3 si crea il file di indice del genoma (`data/genome_index.tar.gz`) necessario per la fase di allineamento in questa lezione.

---

## 1. Scrivere un flusso di lavoro a singolo stadio che esegue il controllo qualità iniziale

Iniziamo scrivendo un semplice flusso di lavoro che esegue lo strumento FastQC su un file FASTQ contenente letture RNAseq single-end.

Vi forniamo un file di flusso di lavoro, `rnaseq.nf`, che delinea le parti principali del flusso di lavoro.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Module INCLUDE statements

/*
 * Pipeline parameters
 */

// Primary input

workflow {

    // Create input channel

    // Call processes

}
```

Tenete presente che questo codice del flusso di lavoro è corretto ma non è funzionale; il suo scopo è solo quello di servire come scheletro che userete per scrivere il flusso di lavoro effettivo.

### 1.1. Creare una directory per memorizzare i moduli

Creeremo moduli autonomi per ogni processo per facilitarne la gestione e il riutilizzo, quindi creiamo una directory per memorizzarli.

```bash
mkdir modules
```

### 1.2. Creare un modulo per il processo di raccolta delle metriche QC

Creiamo un file modulo chiamato `modules/fastqc.nf` per contenere il processo `FASTQC`:

```bash
touch modules/fastqc.nf
```

Aprite il file nell'editor di codice e copiatevi il seguente codice:

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

Dovreste riconoscere tutti i pezzi da ciò che avete imparato nella Parte 1 e Parte 2 di questa serie di formazione; l'unica modifica degna di nota è che questa volta stiamo usando `mode: symlink` per la direttiva `publishDir`, e stiamo usando un parametro per definire il `publishDir`.

!!! note

    Anche se i file di dati che stiamo usando qui sono molto piccoli, in genomica possono diventare molto grandi. Per scopi dimostrativi nell'ambiente di formazione, stiamo usando la modalità di pubblicazione 'symlink' per evitare copie di file non necessarie. Non dovreste farlo nei vostri flussi di lavoro finali, poiché perderete i risultati quando pulirete la vostra directory `work`.

### 1.3. Importare il modulo nel file del flusso di lavoro

Aggiungete l'istruzione `include { FASTQC } from './modules/fastqc.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Aggiungere una dichiarazione di input

Dichiarate un parametro di input con un valore predefinito:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Creare un canale di input nel blocco workflow

Usate una semplice fabbrica di canali `.fromPath()` per creare il canale di input:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Call processes

}
```

### 1.6. Chiamare il processo `FASTQC` sul canale di input

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
```

### 1.7. Eseguire il flusso di lavoro per verificare che funzioni

Potremmo usare il parametro `--reads` per specificare un input da riga di comando, ma durante lo sviluppo possiamo essere pigri e usare semplicemente il valore predefinito di test che abbiamo impostato.

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

Questo dovrebbe essere eseguito molto rapidamente se avete completato la Parte 1 e avete già scaricato il container.
Se l'avete saltata, Nextflow scaricherà il container per voi; non dovete fare nulla perché ciò avvenga, ma potreste dover attendere fino a un minuto.

Potete trovare gli output sotto `results/fastqc` come specificato nel processo `FASTQC` dalla direttiva `publishDir`.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Aggiungere il trimming degli adattatori e il controllo qualità post-trimming

Useremo il wrapper Trim_Galore, che raggruppa Cutadapt per il trimming stesso e FastQC per il controllo qualità post-trimming.

### 2.1. Creare un modulo per il processo di trimming e QC

Creiamo un file modulo chiamato `modules/trim_galore.nf` per contenere il processo `TRIM_GALORE`:

```bash
touch modules/trim_galore.nf
```

Aprite il file nell'editor di codice e copiatevi il seguente codice:

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

### 2.2. Importare il modulo nel file del flusso di lavoro

Aggiungete l'istruzione `include { TRIM_GALORE } from './modules/trim_galore.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Chiamare il processo sul canale di input

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. Eseguire il flusso di lavoro per verificare che funzioni

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

Potete trovare gli output sotto `results/trimming` come specificato nel processo `TRIM_GALORE` dalla direttiva `publishDir`.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Allineare le letture al genoma di riferimento

Infine possiamo eseguire la fase di allineamento del genoma usando Hisat2, che emetterà anche metriche di controllo qualità in stile FastQC.

### 3.1. Creare un modulo per il processo HiSat2

Creiamo un file modulo chiamato `modules/hisat2_align.nf` per contenere il processo `HISAT2_ALIGN`:

```bash
touch modules/hisat2_align.nf
```

Aprite il file nell'editor di codice e copiatevi il seguente codice:

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

### 3.2. Importare il modulo nel file del flusso di lavoro

Aggiungete l'istruzione `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Aggiungere una dichiarazione di parametro per fornire l'indice del genoma

Dichiarate un parametro di input con un valore predefinito:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Chiamare il processo `HISAT2_ALIGN` sulle letture trimmate in output da `TRIM_GALORE`

Le letture trimmate si trovano nel canale `TRIM_GALORE.out.trimmed_reads` in output dalla fase precedente.

Inoltre, usiamo `file (params.hisat2_index_zip)` per fornire allo strumento Hisat2 il tarball compresso dell'indice del genoma.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Eseguire il flusso di lavoro per verificare che funzioni

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

Potete trovare gli output sotto `results/align` come specificato nel processo `HISAT2_ALIGN` dalla direttiva `publishDir`.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Questo completa l'elaborazione di base che dobbiamo applicare a ciascun campione.

_Aggiungeremo l'aggregazione dei report MultiQC nella Parte 2, dopo aver modificato il flusso di lavoro per accettare più campioni alla volta._

---

### Takeaway

Sapete come racchiudere tutte le fasi principali per elaborare campioni RNAseq single-end individualmente.

### Cosa c'è dopo?

Imparate come modificare il flusso di lavoro per elaborare più campioni in parallelo, aggregare i report QC attraverso tutte le fasi per tutti i campioni e abilitare l'esecuzione del flusso di lavoro su dati RNAseq paired-end.
