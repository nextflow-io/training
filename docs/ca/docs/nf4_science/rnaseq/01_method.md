# Part 1: Visió general del mètode

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hi ha múltiples mètodes vàlids per processar i analitzar dades de RNAseq en bulk.
Per a aquest curs, seguim el mètode descrit [aquí](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) pels Drs. Simon Andrews i Laura Biggins al [Babraham Institute](https://www.babraham.ac.uk/).

El nostre objectiu és desenvolupar un workflow que implementi els següents passos de processament: executar un control de qualitat inicial sobre les lectures d'una mostra de RNAseq en bulk, retallar les seqüències adaptadores de les lectures, alinear les lectures a un genoma de referència i produir un informe de control de qualitat (QC) complet.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Realitzar QC sobre les dades de lectura abans del retallat utilitzant FastQC
- **TRIM_GALORE:** Retallar seqüències adaptadores i realitzar QC després del retallat utilitzant Trim Galore (agrupa Cutadapt i FastQC)
- **HISAT2_ALIGN:** Alinear lectures al genoma de referència utilitzant Hisat2
- **MULTIQC:** Generar un informe de QC complet utilitzant MultiQC

### Mètodes

Us mostrarem com aplicar aquests passos de processament en dues fases.
Primer començarem amb el **processament d'una sola mostra** que executa les eines de QC, retallat i alineament sobre una mostra.
Després ampliarem al **processament de múltiples mostres** que executa les mateixes eines sobre múltiples mostres i genera un informe de control de qualitat agregat.

Abans d'endinsar-nos en escriure cap codi de workflow per a qualsevol dels dos enfocaments, provarem les comandes manualment sobre algunes dades de prova.

### Conjunt de dades

Proporcionem les següents dades i recursos relacionats:

- **Dades RNAseq** (`reads/`): fitxers FASTQ de sis mostres, reduïts a una petita regió per mantenir les mides de fitxer baixes. Cada mostra té lectures paired-end (dos fitxers per mostra), tot i que comencem treballant només amb lectures single-end.
- **Un genoma de referència** (`genome.fa`): una petita regió del cromosoma humà 20 (de hg19/b37).
- **Fulls de càlcul CSV** (`single-end.csv` i `paired-end.csv`): fitxers que llisten els IDs i camins dels fitxers de dades d'exemple.

### Programari

Les quatre eines principals implicades són [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) per a la recollida de mètriques de control de qualitat, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) per al retallat d'adaptadors (agrupa Cutadapt i FastQC per a QC post-retallat), [HISAT2](http://daehwankimlab.github.io/hisat2/) per a l'alineament amb splicing a un genoma de referència, i [MultiQC](https://multiqc.info/) per a la generació d'informes de QC agregats.

Aquestes eines no estan instal·lades a l'entorn de GitHub Codespaces, així que les utilitzarem via contenidors recuperats mitjançant el servei Seqera Containers (vegeu [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Consell"

     Assegureu-vos que esteu al directori `nf4-science/rnaseq`. L'última part del camí mostrat quan escriviu `pwd` hauria de ser `rnaseq`.

---

## 1. Processament d'una sola mostra

En aquesta secció provem les comandes que processen una sola mostra de RNAseq: control de qualitat, retallat d'adaptadors i alineament a un genoma de referència.
Aquestes són les comandes que encapsularem en un workflow de Nextflow a la Part 2 d'aquest curs.

1. Executar QC inicial sobre un fitxer FASTQ utilitzant FastQC
2. Retallar seqüències adaptadores i executar QC post-retallat utilitzant Trim Galore
3. Alinear les lectures retallades al genoma de referència utilitzant HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Comencem provant aquestes comandes sobre només una mostra.

### 1.1. QC i retallat d'adaptadors

Primer, volem executar les comandes de QC i retallat sobre un dels fitxers de dades d'exemple.

#### 1.1.1. Descarregar el contenidor

Descarreguem una imatge de contenidor que té tant `fastqc` com `trim_galore` instal·lats:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Sortida de la comanda"

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

Si no heu descarregat aquesta imatge abans, pot trigar un minut a completar-se.
Un cop finalitzi, tindreu una còpia local de la imatge del contenidor.

#### 1.1.2. Iniciar el contenidor interactivament

Per executar el contenidor interactivament, utilitzeu `docker run` amb les opcions `-it`.
L'opció `-v ./data:/data` munta el nostre directori local `data/` perquè puguem accedir als fitxers d'entrada des de dins del contenidor.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Sortida de la comanda"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

El vostre prompt canviarà a quelcom com `(base) root@b645838b3314:/tmp#`, que indica que ara esteu dins del contenidor.

Verifiqueu que podeu veure els fitxers de dades de seqüenciació sota `/data/reads`:

```bash
ls /data/reads
```

??? abstract "Contingut del directori"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Amb això, esteu preparats per provar la vostra primera comanda.

#### 1.1.3. Executar la comanda FastQC

El mètode referenciat anteriorment ens dóna la línia de comandes per executar QC sobre un sol fitxer.
Només necessitem proporcionar el fitxer d'entrada; l'eina generarà automàticament fitxers de sortida al mateix directori que les dades originals.

Executeu la comanda `fastqc` sobre un fitxer de dades:

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Sortida de la comanda"

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

Això hauria d'executar-se molt ràpidament.
Podeu trobar els fitxers de sortida al mateix directori que les dades originals:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "Contingut del directori"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Hauríeu de veure un informe HTML i un arxiu ZIP que conté les mètriques de QC.
Això completa la prova del primer pas.

#### 1.1.4. Retallar seqüències adaptadores amb Trim Galore

Ara executem `trim_galore`, que agrupa Cutadapt i FastQC, per retallar les seqüències adaptadores i recollir mètriques de QC post-retallat.
Com s'ha indicat anteriorment, el programari està inclòs al mateix contenidor, així que no cal cap canvi.

La comanda és senzilla; simplement necessitem afegir l'opció `--fastqc` per executar automàticament un pas de recollida de QC després que el retallat estigui complet.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Sortida de la comanda"

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

La sortida és molt detallada, així que hem destacat les línies més rellevants a l'exemple anterior.
Podeu trobar els fitxers de sortida al directori de treball:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contingut del directori"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Això inclou les lectures retallades, l'informe de retallat i els fitxers de QC post-retallat.

#### 1.1.5. Moure els fitxers de sortida

Qualsevol cosa que romangui dins del contenidor serà inaccessible per a treballs futurs, així que necessitem moure aquests fitxers a un directori al sistema de fitxers muntat.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Contingut del directori"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Els fitxers ara són accessibles al vostre sistema de fitxers normal.

#### 1.1.6. Sortir del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre prompt hauria de tornar a la normalitat; això completa la prova dels dos primers passos.

### 1.2. Alinear lectures al genoma de referència

A continuació, volem executar la comanda d'alineament per alinear les lectures de RNAseq retallades a un genoma de referència.

#### 1.2.1. Descarregar el contenidor

Descarreguem una imatge de contenidor que té `hisat2` i `samtools` instal·lats:

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Sortida de la comanda"

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

Notareu que algunes capes mostren `Already exists` perquè es comparteixen amb la imatge del contenidor Trim Galore que vam descarregar anteriorment.
Com a resultat, aquesta descàrrega hauria d'anar més ràpida que la primera.

#### 1.2.2. Iniciar el contenidor interactivament

Inicieu el contenidor interactivament, utilitzant el mateix enfocament que abans amb l'URI del contenidor rellevant intercanviat.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

El vostre prompt canviarà de nou per indicar que esteu dins del contenidor.

#### 1.2.3. Crear els fitxers d'índex del genoma

HISAT2 requereix que la referència del genoma es proporcioni en un format molt específic, i no pot simplement consumir el fitxer FASTA `genome.fa` que proporcionem, així que aprofitarem aquesta oportunitat per crear els recursos rellevants.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Sortida de la comanda"

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

La sortida és molt detallada, així que hem destacat algunes línies rellevants a l'exemple anterior.

Això crea múltiples fitxers d'índex del genoma, que podeu trobar al directori de treball.

```bash
ls genome_index.*
```

??? abstract "Contingut del directori"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Necessitarem aquests fitxers més endavant, i generar-los no és típicament quelcom que vulguem fer com a part d'un workflow, així que generarem un tarball comprimit amb gzip que contingui els fitxers d'índex del genoma que podem passar fàcilment segons sigui necessari.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Sortida de la comanda"

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

Mourem el tarball resultant `genome_index.tar.gz` que conté els fitxers d'índex del genoma al directori `data/` del nostre sistema de fitxers en uns minuts.
Això serà útil a la Part 2 d'aquest curs.

#### 1.2.4. Executar la comanda d'alineament

Ara podem executar la comanda d'alineament, que realitza el pas d'alineament amb `hisat2` i després canalitza la sortida a `samtools` per escriure la sortida com un fitxer BAM.

L'entrada de dades de lectura és el fitxer `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que vam generar amb `trim_galore` al pas anterior.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

Això s'executa gairebé instantàniament perquè és un fitxer de prova molt petit.
A escala completa, això podria trigar molt més.

Un cop més podeu trobar els fitxers de sortida al directori de treball:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contingut del directori"

    ```console title="Output"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

L'alineament va produir un fitxer BAM i un fitxer de registre amb estadístiques d'alineament.

#### 1.2.5. Moure els fitxers de sortida

Com abans, moveu els fitxers de sortida a un directori al sistema de fitxers muntat perquè romanguin accessibles després de sortir del contenidor.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Amb això fet, tenim tot el que necessitem.

#### 1.2.6. Sortir del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre prompt hauria de tornar a la normalitat.
Això conclou l'execució de prova del processament d'una sola mostra.

!!! example "Escriviu-ho com un workflow!"

    Sentiu-vos lliures de passar directament a la [Part 2](./02_single-sample.md) si voleu començar a implementar aquesta anàlisi com un workflow de Nextflow.
    Només haureu de tornar per completar la segona ronda de proves abans de passar a la Part 3.

---

## 2. Agregació de QC de múltiples mostres

Les comandes que acabem de provar processen una mostra cada vegada.
A la pràctica, típicament necessitem processar moltes mostres i després agregar els resultats de QC a través de totes elles per avaluar la qualitat del conjunt de dades global.

[MultiQC](https://multiqc.info/) és una eina que cerca a través de directoris informes de QC de moltes eines bioinformàtiques comunes i els agrega en un sol informe HTML complet.
Pot reconèixer sortides de FastQC, Cutadapt (via Trim Galore) i HISAT2, entre moltes altres.

Aquí processem dues mostres addicionals a través de les mateixes eines per mostra, després utilitzem MultiQC per agregar informes de QC a través de les tres mostres.
Aquestes són les comandes que encapsularem en un workflow de Nextflow a la Part 3 d'aquest curs.

1. Executar QC i retallat sobre mostres addicionals utilitzant Trim Galore
2. Executar alineament sobre mostres addicionals utilitzant HISAT2
3. Agregar tots els informes de QC en un informe complet utilitzant MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC i retallat de mostres addicionals

Les comandes de QC i retallat per mostra són idèntiques al que vam executar a la secció 1.1.
Ja vam descarregar la imatge del contenidor, així que podem iniciar-la directament.

#### 2.1.1. Iniciar el contenidor

Ja vam descarregar aquesta imatge de contenidor a la secció 1.1, així que podem iniciar-la directament:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

El vostre prompt canvia per indicar que esteu dins del contenidor.

#### 2.1.2. Executar QC i retallat sobre mostres addicionals

Executeu FastQC i Trim Galore sobre dues mostres més, una després de l'altra.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Un cop això es completi, hauríeu de tenir fitxers de sortida de Trim Galore per a ambdues mostres al directori de treball.

#### 2.1.3. Moure els fitxers de sortida

Moveu els fitxers de sortida de Trim Galore al mateix directori que vam utilitzar a la secció 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Contingut del directori"

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

Els fitxers ara són accessibles al vostre sistema de fitxers normal.

#### 2.1.4. Sortir del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre prompt hauria de tornar a la normalitat.

### 2.2. Alinear mostres addicionals

Les comandes d'alineament són idèntiques al que vam executar a la secció 1.2.
Necessitem extreure l'índex del genoma del tarball que vam desar anteriorment, ja que els fitxers d'índex originals es van crear dins d'un contenidor que ja no existeix.

#### 2.2.1. Iniciar el contenidor

Ja vam descarregar aquesta imatge de contenidor a la secció 1.2, així que podem iniciar-la directament:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

El vostre prompt canvia per indicar que esteu dins del contenidor.

#### 2.2.2. Extreure l'índex del genoma

Extraieu els fitxers d'índex del genoma del tarball que vam desar al sistema de fitxers muntat:

```bash
tar -xzf /data/genome_index.tar.gz
```

Això restaura els fitxers `genome_index.*` al directori de treball.

#### 2.2.3. Executar alineament sobre mostres addicionals

Executeu l'alineament HISAT2 sobre les dues mostres retallades recentment, una després de l'altra.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Sortida de la comanda"

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

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Un cop això es completi, hauríeu de tenir fitxers BAM i de registre per a ambdues mostres al directori de treball.

#### 2.2.4. Moure els fitxers de sortida

Moveu els fitxers de sortida d'alineament al mateix directori que vam utilitzar a la secció 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Contingut del directori"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Els fitxers ara són accessibles al vostre sistema de fitxers normal.

#### 2.2.5. Sortir del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre prompt hauria de tornar a la normalitat.

### 2.3. Generar un informe de QC complet

Ara que tenim sortides de QC, retallat i alineament per a tres mostres, podem utilitzar MultiQC per agregar-les en un sol informe.
MultiQC cerca a través de directoris informes de QC compatibles i agrega tot el que troba.

#### 2.3.1. Descarregar el contenidor

Descarreguem una imatge de contenidor que té `multiqc` instal·lat:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Sortida de la comanda"

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

Notareu que algunes capes mostren `Already exists` perquè es comparteixen amb les imatges de contenidor que vam descarregar anteriorment.
Com a resultat, aquesta descàrrega hauria d'anar més ràpida que les anteriors.

#### 2.3.2. Iniciar el contenidor interactivament

Inicieu el contenidor interactivament amb el directori de dades muntat, com abans.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

El vostre prompt canviarà per indicar que esteu dins del contenidor.

#### 2.3.3. Executar la comanda MultiQC

Executeu `multiqc`, apuntant-lo als directoris on vam emmagatzemar fitxers de sortida relacionats amb QC per a les tres mostres.
L'opció `-n` estableix el nom de l'informe de sortida.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Sortida de la comanda"

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

Aquí veiem que l'eina va trobar informes de QC per a les tres mostres: el QC inicial de `fastqc`, els informes post-retallat de `cutadapt` (via `trim_galore`) i els resums d'alineament de `hisat2`.

Els fitxers de sortida estan al directori de treball:

```bash
ls all_samples_QC*
```

??? abstract "Contingut del directori"

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

La sortida principal és l'informe `all_samples_QC.html`, acompanyat d'un directori de dades que conté les mètriques subjacents.

#### 2.3.4. Moure els fitxers de sortida

Moveu l'informe i el seu directori de dades al sistema de fitxers muntat.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Els fitxers ara són accessibles al vostre sistema de fitxers normal.

#### 2.3.5. Sortir del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre prompt hauria de tornar a la normalitat.
Això conclou la prova de totes les comandes de processament de RNAseq.

---

### Conclusió

Sabeu com executar les comandes FastQC, Trim Galore, HISAT2 i MultiQC als seus respectius contenidors, incloent com processar múltiples mostres i agregar informes de QC.

### Què segueix?

Feu una pausa i després aneu a la [Part 2](./02_single-sample.md) per aprendre com encapsular aquestes mateixes comandes en workflows que utilitzen contenidors per executar el treball.
