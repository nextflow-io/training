# Parte 1: Panoramica del metodo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Esistono molteplici metodi validi per elaborare e analizzare dati RNAseq in bulk.
Per questo corso, seguiamo il metodo descritto [qui](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) dai Dott. Simon Andrews e Laura Biggins presso il [Babraham Institute](https://www.babraham.ac.uk/).

Il nostro obiettivo è sviluppare un flusso di lavoro che implementi le seguenti fasi di elaborazione: eseguire il controllo qualità iniziale sulle read in un campione RNAseq in bulk, rimuovere le sequenze degli adapter dalle read, allineare le read a un genoma di riferimento e produrre un report completo di controllo qualità (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Eseguire il QC sui dati delle read prima del trimming utilizzando FastQC
- **TRIM_GALORE:** Rimuovere le sequenze degli adapter ed eseguire il QC dopo il trimming utilizzando Trim Galore (che include Cutadapt e FastQC)
- **HISAT2_ALIGN:** Allineare le read al genoma di riferimento utilizzando Hisat2
- **MULTIQC:** Generare un report QC completo utilizzando MultiQC

### Metodi

Vi mostreremo come applicare queste fasi di elaborazione in due approcci.
Prima inizieremo con l'**elaborazione di un singolo campione** che esegue gli strumenti di QC, trimming e allineamento su un campione.
Poi estenderemo all'**elaborazione multi-campione** che esegue gli stessi strumenti su più campioni e genera un report di controllo qualità aggregato.

Prima di iniziare a scrivere qualsiasi codice del flusso di lavoro per entrambi gli approcci, proveremo i comandi manualmente su alcuni dati di test.

### Dataset

Forniamo i seguenti dati e risorse correlate:

- **Dati RNAseq** (`reads/`): file FASTQ da sei campioni, ridotti a una piccola regione per mantenere le dimensioni dei file contenute. Ogni campione ha read paired-end (due file per campione), anche se iniziamo lavorando solo con read single-end.
- **Un genoma di riferimento** (`genome.fa`): una piccola regione del cromosoma umano 20 (da hg19/b37).
- **Samplesheet CSV** (`single-end.csv` e `paired-end.csv`): file che elencano gli ID e i percorsi dei file di dati di esempio.

### Software

I quattro strumenti principali coinvolti sono [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) per la raccolta delle metriche di controllo qualità, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) per il trimming degli adapter (include Cutadapt e FastQC per il QC post-trimming), [HISAT2](http://daehwankimlab.github.io/hisat2/) per l'allineamento spliced a un genoma di riferimento, e [MultiQC](https://multiqc.info/) per la generazione di report QC aggregati.

Questi strumenti non sono installati nell'ambiente GitHub Codespaces, quindi li utilizzeremo tramite container recuperati tramite il servizio Seqera Containers (vedere [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Suggerimento"

     Verificate di trovarvi nella directory `nf4-science/rnaseq`. L'ultima parte del percorso mostrato quando digitate `pwd` dovrebbe essere `rnaseq`.

---

## 1. Elaborazione di un singolo campione

In questa sezione testiamo i comandi che elaborano un singolo campione RNAseq: controllo qualità, trimming degli adapter e allineamento a un genoma di riferimento.
Questi sono i comandi che inseriremo in un flusso di lavoro Nextflow nella Parte 2 di questo corso.

1. Eseguire il QC iniziale su un file FASTQ utilizzando FastQC
2. Rimuovere le sequenze degli adapter ed eseguire il QC post-trimming utilizzando Trim Galore
3. Allineare le read trimmati al genoma di riferimento utilizzando HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Iniziamo testando questi comandi su un solo campione.

### 1.1. QC e trimming degli adapter

Prima, vogliamo eseguire i comandi di QC e trimming su uno dei file di dati di esempio.

#### 1.1.1. Ottenere il container

Otteniamo un'immagine container che ha sia `fastqc` che `trim_galore` installati:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Output del comando"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

Se non avete scaricato questa immagine in precedenza, potrebbe richiedere un minuto per completare.
Una volta terminato, avrete una copia locale dell'immagine container.

#### 1.1.2. Avviare il container in modo interattivo

Per eseguire il container in modo interattivo, utilizzate `docker run` con i flag `-it`.
L'opzione `-v ./data:/data` monta la nostra directory locale `data/` in modo da poter accedere ai file di input dall'interno del container.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Output del comando"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Il prompt cambierà in qualcosa come `(base) root@b645838b3314:/tmp#`, il che indica che vi trovate all'interno del container.

Verificate di poter vedere i file di dati delle sequenze sotto `/data/reads`:

```bash
ls /data/reads
```

??? abstract "Contenuto della directory"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Con questo, siete pronti per provare il vostro primo comando.

#### 1.1.3. Eseguire il comando FastQC

Il metodo citato sopra ci fornisce la riga di comando per eseguire il QC su un singolo file.
Dobbiamo solo fornire il file di input; lo strumento genererà automaticamente i file di output nella stessa directory dei dati originali.

Eseguite il comando `fastqc` su un file di dati:

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Output del comando"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Questo dovrebbe essere eseguito molto rapidamente.
Potete trovare i file di output nella stessa directory dei dati originali:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "Contenuto della directory"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Dovreste vedere un report HTML e un archivio ZIP contenente le metriche QC.
Questo completa il test della prima fase.

#### 1.1.4. Rimuovere le sequenze degli adapter con Trim Galore

Ora eseguiamo `trim_galore`, che include Cutadapt e FastQC, per rimuovere le sequenze degli adapter e raccogliere le metriche QC post-trimming.
Come notato sopra, il software è incluso nello stesso container, quindi non è necessario alcun cambiamento.

Il comando è semplice; dobbiamo solo aggiungere il flag `--fastqc` per eseguire automaticamente una fase di raccolta QC dopo il completamento del trimming.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Output del comando"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

L'output è molto dettagliato, quindi abbiamo evidenziato le righe più rilevanti nell'esempio sopra.
Potete trovare i file di output nella directory di lavoro:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contenuto della directory"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Questo include le read trimmati, il report di trimming e i file QC post-trimming.

#### 1.1.5. Spostare i file di output

Tutto ciò che rimane all'interno del container sarà inaccessibile per il lavoro futuro, quindi dobbiamo spostare questi file in una directory sul filesystem montato.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Contenuto della directory"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

I file sono ora accessibili nel vostro filesystem normale.

#### 1.1.6. Uscire dal container

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il prompt dovrebbe tornare normale; questo completa il test delle prime due fasi.

### 1.2. Allineare le read al genoma di riferimento

Successivamente, vogliamo eseguire il comando di allineamento per allineare le read RNAseq trimmati a un genoma di riferimento.

#### 1.2.1. Ottenere il container

Otteniamo un'immagine container che ha `hisat2` e `samtools` installati:

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Output del comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

Noterete che alcuni layer mostrano `Already exists` perché sono condivisi con l'immagine container Trim Galore che abbiamo ottenuto in precedenza.
Di conseguenza, questo download dovrebbe essere più veloce del primo.

#### 1.2.2. Avviare il container in modo interattivo

Avviate il container in modo interattivo, utilizzando lo stesso approccio di prima con l'URI del container pertinente sostituito.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Il prompt cambierà di nuovo per indicare che vi trovate all'interno del container.

#### 1.2.3. Creare i file di indice del genoma

HISAT2 richiede che il riferimento del genoma sia fornito in un formato molto specifico e non può semplicemente utilizzare il file FASTA `genome.fa` che forniamo, quindi coglieremo questa opportunità per creare le risorse pertinenti.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Output del comando"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

L'output è molto dettagliato, quindi abbiamo evidenziato alcune righe rilevanti nell'esempio sopra.

Questo crea molteplici file di indice del genoma, che potete trovare nella directory di lavoro.

```bash
ls genome_index.*
```

??? abstract "Contenuto della directory"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Avremo bisogno di questi file più avanti, e generarli non è tipicamente qualcosa che vogliamo fare come parte di un flusso di lavoro, quindi genereremo un tarball compresso contenente i file di indice del genoma che possiamo facilmente passare secondo necessità.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Output del comando"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

Sposteremo il tarball risultante `genome_index.tar.gz` contenente i file di indice del genoma nella directory `data/` sul nostro filesystem tra pochi minuti.
Questo tornerà utile nella Parte 2 di questo corso.

#### 1.2.4. Eseguire il comando di allineamento

Ora possiamo eseguire il comando di allineamento, che esegue la fase di allineamento con `hisat2` e poi invia l'output a `samtools` per scrivere l'output come file BAM.

L'input dei dati delle read è il file `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` che abbiamo generato con `trim_galore` nella fase precedente.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Output del comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

Questo viene eseguito quasi istantaneamente perché è un file di test molto piccolo.
Su scala reale potrebbe richiedere molto più tempo.

Ancora una volta potete trovare i file di output nella directory di lavoro:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contenuto della directory"

    ```console title="Output"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

L'allineamento ha prodotto un file BAM e un file di log con le statistiche di allineamento.

#### 1.2.5. Spostare i file di output

Come prima, spostate i file di output in una directory sul filesystem montato in modo che rimangano accessibili dopo l'uscita dal container.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Con questo fatto, abbiamo tutto ciò di cui abbiamo bisogno.

#### 1.2.6. Uscire dal container

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il prompt dovrebbe tornare normale.
Questo conclude l'esecuzione di test dell'elaborazione di un singolo campione.

!!! example "Scrivetelo come un flusso di lavoro!"

    Sentitevi liberi di passare direttamente alla [Parte 2](./02_single-sample.md) se volete iniziare a implementare questa analisi come un flusso di lavoro Nextflow.
    Dovrete solo tornare per completare il secondo round di test prima di passare alla Parte 3.

---

## 2. Aggregazione QC multi-campione

I comandi che abbiamo appena testato elaborano un campione alla volta.
In pratica, tipicamente dobbiamo elaborare molti campioni e poi aggregare i risultati QC su tutti loro per valutare la qualità del dataset complessivo.

[MultiQC](https://multiqc.info/) è uno strumento che cerca nelle directory i report QC di molti strumenti bioinformatici comuni e li aggrega in un singolo report HTML completo.
Può riconoscere l'output di FastQC, Cutadapt (tramite Trim Galore) e HISAT2, tra molti altri.

Qui elaboriamo due campioni aggiuntivi attraverso gli stessi strumenti per-campione, poi utilizziamo MultiQC per aggregare i report QC su tutti e tre i campioni.
Questi sono i comandi che inseriremo in un flusso di lavoro Nextflow nella Parte 3 di questo corso.

1. Eseguire QC e trimming su campioni aggiuntivi utilizzando Trim Galore
2. Eseguire l'allineamento su campioni aggiuntivi utilizzando HISAT2
3. Aggregare tutti i report QC in un report completo utilizzando MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC e trimming di campioni aggiuntivi

I comandi di QC e trimming per-campione sono identici a quelli che abbiamo eseguito nella sezione 1.1.
Abbiamo già ottenuto l'immagine container, quindi possiamo avviarla direttamente.

#### 2.1.1. Avviare il container

Abbiamo già ottenuto questa immagine container nella sezione 1.1, quindi possiamo avviarla direttamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Il prompt cambia per indicare che vi trovate all'interno del container.

#### 2.1.2. Eseguire QC e trimming su campioni aggiuntivi

Eseguite FastQC e Trim Galore su altri due campioni, uno dopo l'altro.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Una volta completato, dovreste avere i file di output di Trim Galore per entrambi i campioni nella directory di lavoro.

#### 2.1.3. Spostare i file di output

Spostate i file di output di Trim Galore nella stessa directory che abbiamo utilizzato nella sezione 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Contenuto della directory"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

I file sono ora accessibili nel vostro filesystem normale.

#### 2.1.4. Uscire dal container

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il prompt dovrebbe tornare normale.

### 2.2. Allineare campioni aggiuntivi

I comandi di allineamento sono identici a quelli che abbiamo eseguito nella sezione 1.2.
Dobbiamo estrarre l'indice del genoma dal tarball che abbiamo salvato in precedenza, poiché i file di indice originali sono stati creati all'interno di un container che non esiste più.

#### 2.2.1. Avviare il container

Abbiamo già ottenuto questa immagine container nella sezione 1.2, quindi possiamo avviarla direttamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Il prompt cambia per indicare che vi trovate all'interno del container.

#### 2.2.2. Estrarre l'indice del genoma

Estraete i file di indice del genoma dal tarball che abbiamo salvato sul filesystem montato:

```bash
tar -xzf /data/genome_index.tar.gz
```

Questo ripristina i file `genome_index.*` nella directory di lavoro.

#### 2.2.3. Eseguire l'allineamento su campioni aggiuntivi

Eseguite l'allineamento HISAT2 sui due campioni appena trimmati, uno dopo l'altro.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Output del comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "Output del comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Una volta completato, dovreste avere i file BAM e log per entrambi i campioni nella directory di lavoro.

#### 2.2.4. Spostare i file di output

Spostate i file di output dell'allineamento nella stessa directory che abbiamo utilizzato nella sezione 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Contenuto della directory"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

I file sono ora accessibili nel vostro filesystem normale.

#### 2.2.5. Uscire dal container

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il prompt dovrebbe tornare normale.

### 2.3. Generare un report QC completo

Ora che abbiamo l'output di QC, trimming e allineamento per tre campioni, possiamo utilizzare MultiQC per aggregarli in un singolo report.
MultiQC cerca nelle directory i report QC compatibili e aggrega tutto ciò che trova.

#### 2.3.1. Ottenere il container

Otteniamo un'immagine container che ha `multiqc` installato:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Output del comando"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

Noterete che alcuni layer mostrano `Already exists` perché sono condivisi con le immagini container che abbiamo ottenuto in precedenza.
Di conseguenza, questo download dovrebbe essere più veloce dei precedenti.

#### 2.3.2. Avviare il container in modo interattivo

Avviate il container in modo interattivo con la directory data montata, come prima.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Il prompt cambierà per indicare che vi trovate all'interno del container.

#### 2.3.3. Eseguire il comando MultiQC

Eseguite `multiqc`, puntandolo alle directory dove abbiamo memorizzato i file di output relativi al QC per tutti e tre i campioni.
Il flag `-n` imposta il nome del report di output.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Output del comando"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

Qui vediamo che lo strumento ha trovato i report QC per tutti e tre i campioni: il QC iniziale da `fastqc`, i report post-trimming da `cutadapt` (tramite `trim_galore`) e i riepiloghi di allineamento da `hisat2`.

I file di output sono nella directory di lavoro:

```bash
ls all_samples_QC*
```

??? abstract "Contenuto della directory"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

L'output principale è il report `all_samples_QC.html`, accompagnato da una directory di dati contenente le metriche sottostanti.

#### 2.3.4. Spostare i file di output

Spostate il report e la sua directory di dati sul filesystem montato.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

I file sono ora accessibili nel vostro filesystem normale.

#### 2.3.5. Uscire dal container

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il prompt dovrebbe tornare normale.
Questo conclude il test di tutti i comandi di elaborazione RNAseq.

---

### Takeaway

Sapete come eseguire i comandi FastQC, Trim Galore, HISAT2 e MultiQC nei rispettivi container, incluso come elaborare più campioni e aggregare i report QC.

### Cosa c'è dopo?

Prendetevi una pausa, poi passate alla [Parte 2](./02_single-sample.md) per imparare come inserire gli stessi comandi in flussi di lavoro che utilizzano container per eseguire il lavoro.
