# Teil 1: Methodenübersicht

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Es gibt mehrere gültige Methoden zur Verarbeitung und Analyse von Bulk-RNAseq-Daten.
Für diesen Kurs folgen wir der Methode, die [hier](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) von Dr. Simon Andrews und Dr. Laura Biggins am [Babraham Institute](https://www.babraham.ac.uk/) beschrieben wird.

Unser Ziel ist es, einen Workflow zu entwickeln, der die folgenden Verarbeitungsschritte implementiert: initiale Qualitätskontrolle der Reads in einer Bulk-RNAseq-Probe durchführen, Adaptersequenzen aus den Reads trimmen, die Reads auf ein Referenzgenom alignieren und einen umfassenden Qualitätskontroll-Bericht (QC) erstellen.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** QC der Read-Daten vor dem Trimmen mit FastQC durchführen
- **TRIM_GALORE:** Adaptersequenzen trimmen und QC nach dem Trimmen mit Trim Galore durchführen (bündelt Cutadapt und FastQC)
- **HISAT2_ALIGN:** Reads mit Hisat2 auf das Referenzgenom alignieren
- **MULTIQC:** Einen umfassenden QC-Bericht mit MultiQC generieren

### Methoden

Wir zeigen dir, wie du diese Verarbeitungsschritte in zwei Phasen anwendest.
Zuerst beginnen wir mit der **Verarbeitung einzelner Proben**, die die QC-, Trimming- und Alignment-Tools auf eine Probe anwendet.
Dann erweitern wir auf die **Verarbeitung mehrerer Proben**, die dieselben Tools auf mehrere Proben anwendet und einen aggregierten Qualitätskontroll-Bericht generiert.

Bevor wir mit dem Schreiben von Workflow-Code für einen der beiden Ansätze beginnen, werden wir die Befehle manuell an einigen Testdaten ausprobieren.

### Datensatz

Wir stellen die folgenden Daten und zugehörigen Ressourcen bereit:

- **RNAseq-Daten** (`reads/`): FASTQ-Dateien von sechs Proben, auf eine kleine Region reduziert, um die Dateigrößen klein zu halten. Jede Probe hat Paired-End-Reads (zwei Dateien pro Probe), obwohl wir zunächst nur mit Single-End-Reads arbeiten.
- **Ein Referenzgenom** (`genome.fa`): eine kleine Region des menschlichen Chromosoms 20 (aus hg19/b37).
- **CSV-Samplesheets** (`single-end.csv` und `paired-end.csv`): Dateien, die die IDs und Pfade der Beispieldateien auflisten.

### Software

Die vier Haupttools sind [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) für die Sammlung von Qualitätskontroll-Metriken, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) für das Adapter-Trimming (bündelt Cutadapt und FastQC für Post-Trimming-QC), [HISAT2](http://daehwankimlab.github.io/hisat2/) für das gespleißte Alignment auf ein Referenzgenom und [MultiQC](https://multiqc.info/) für die Generierung aggregierter QC-Berichte.

Diese Tools sind in der GitHub Codespaces-Umgebung nicht installiert, daher werden wir sie über Container verwenden, die über den Seqera Containers-Service abgerufen werden (siehe [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Tipp"

     Stelle sicher, dass du dich im Verzeichnis `nf4-science/rnaseq` befindest. Der letzte Teil des Pfads, der angezeigt wird, wenn du `pwd` eingibst, sollte `rnaseq` sein.

---

## 1. Verarbeitung einzelner Proben

In diesem Abschnitt testen wir die Befehle, die eine einzelne RNAseq-Probe verarbeiten: Qualitätskontrolle, Adapter-Trimming und Alignment auf ein Referenzgenom.
Dies sind die Befehle, die wir in Teil 2 dieses Kurses in einen Nextflow-Workflow verpacken werden.

1. Initiale QC auf einer FASTQ-Datei mit FastQC durchführen
2. Adaptersequenzen trimmen und Post-Trimming-QC mit Trim Galore durchführen
3. Die getrimmten Reads mit HISAT2 auf das Referenzgenom alignieren

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Wir beginnen damit, diese Befehle an nur einer Probe zu testen.

### 1.1. QC und Adapter-Trimming

Zuerst möchten wir die QC- und Trimming-Befehle auf einer der Beispieldateien ausführen.

#### 1.1.1. Container herunterladen

Lass uns ein Container-Image herunterladen, das sowohl `fastqc` als auch `trim_galore` installiert hat:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Befehlsausgabe"

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

Falls du dieses Image noch nicht heruntergeladen hast, kann es eine Minute dauern, bis es fertig ist.
Sobald es abgeschlossen ist, hast du eine lokale Kopie des Container-Images.

#### 1.1.2. Container interaktiv starten

Um den Container interaktiv auszuführen, verwende `docker run` mit den Flags `-it`.
Die Option `-v ./data:/data` mountet unser lokales Verzeichnis `data/`, sodass wir von innerhalb des Containers auf die Eingabedateien zugreifen können.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Befehlsausgabe"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Dein Prompt ändert sich zu etwas wie `(base) root@b645838b3314:/tmp#`, was anzeigt, dass du dich jetzt innerhalb des Containers befindest.

Überprüfe, dass du die Sequenzdateien unter `/data/reads` sehen kannst:

```bash
ls /data/reads
```

??? abstract "Verzeichnisinhalt"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Damit bist du bereit, deinen ersten Befehl auszuprobieren.

#### 1.1.3. FastQC-Befehl ausführen

Die oben referenzierte Methode gibt uns die Befehlszeile, um QC auf einer einzelnen Datei auszuführen.
Wir müssen nur die Eingabedatei angeben; das Tool generiert automatisch Ausgabedateien im selben Verzeichnis wie die ursprünglichen Daten.

Führe den `fastqc`-Befehl auf einer Datendatei aus:

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Befehlsausgabe"

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

Dies sollte sehr schnell ablaufen.
Du findest die Ausgabedateien im selben Verzeichnis wie die ursprünglichen Daten:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "Verzeichnisinhalt"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Du solltest einen HTML-Bericht und ein ZIP-Archiv mit den QC-Metriken sehen.
Damit ist das Testen des ersten Schritts abgeschlossen.

#### 1.1.4. Adaptersequenzen mit Trim Galore trimmen

Jetzt führen wir `trim_galore` aus, das Cutadapt und FastQC bündelt, um die Adaptersequenzen zu trimmen und Post-Trimming-QC-Metriken zu sammeln.
Wie oben erwähnt, ist die Software im selben Container enthalten, daher ist keine Änderung erforderlich.

Der Befehl ist unkompliziert; wir müssen nur das Flag `--fastqc` hinzufügen, um automatisch einen QC-Sammelschritt auszuführen, nachdem das Trimmen abgeschlossen ist.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Befehlsausgabe"

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

Die Ausgabe ist sehr ausführlich, daher haben wir die relevantesten Zeilen im obigen Beispiel hervorgehoben.
Du findest die Ausgabedateien im Arbeitsverzeichnis:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Verzeichnisinhalt"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Dies umfasst die getrimmten Reads, den Trimming-Bericht und die Post-Trimming-QC-Dateien.

#### 1.1.5. Ausgabedateien verschieben

Alles, was im Container verbleibt, wird für zukünftige Arbeiten nicht zugänglich sein, daher müssen wir diese Dateien in ein Verzeichnis auf dem gemounteten Dateisystem verschieben.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Verzeichnisinhalt"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Die Dateien sind jetzt in deinem normalen Dateisystem zugänglich.

#### 1.1.6. Container beenden

Um den Container zu beenden, gib `exit` ein.

```bash
exit
```

Dein Prompt sollte wieder normal werden; damit ist das Testen der ersten beiden Schritte abgeschlossen.

### 1.2. Reads auf das Referenzgenom alignieren

Als Nächstes möchten wir den Alignment-Befehl ausführen, um die getrimmten RNAseq-Reads auf ein Referenzgenom zu alignieren.

#### 1.2.1. Container herunterladen

Lass uns ein Container-Image herunterladen, das `hisat2` und `samtools` installiert hat:

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Befehlsausgabe"

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

Du wirst bemerken, dass einige Layer `Already exists` anzeigen, weil sie mit dem Trim Galore-Container-Image geteilt werden, das wir zuvor heruntergeladen haben.
Daher sollte dieser Download schneller gehen als der erste.

#### 1.2.2. Container interaktiv starten

Starte den Container interaktiv, mit demselben Ansatz wie zuvor, wobei die entsprechende Container-URI ausgetauscht wird.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Dein Prompt ändert sich erneut, um anzuzeigen, dass du dich im Container befindest.

#### 1.2.3. Genom-Indexdateien erstellen

HISAT2 benötigt die Genomreferenz in einem sehr spezifischen Format und kann nicht einfach die Datei `genome.fa` im FASTA-Format verarbeiten, die wir bereitstellen. Daher nutzen wir diese Gelegenheit, um die entsprechenden Ressourcen zu erstellen.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Befehlsausgabe"

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

Die Ausgabe ist sehr ausführlich, daher haben wir einige relevante Zeilen im obigen Beispiel hervorgehoben.

Dies erstellt mehrere Genom-Indexdateien, die du im Arbeitsverzeichnis finden kannst.

```bash
ls genome_index.*
```

??? abstract "Verzeichnisinhalt"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Wir werden diese Dateien später benötigen, und das Generieren dieser Dateien ist normalerweise nichts, was wir als Teil eines Workflows tun möchten. Daher werden wir einen gzippten Tarball erstellen, der die Genom-Indexdateien enthält, die wir bei Bedarf einfach weitergeben können.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Befehlsausgabe"

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

Wir werden den resultierenden Tarball `genome_index.tar.gz`, der die Genom-Indexdateien enthält, in wenigen Minuten in das Verzeichnis `data/` auf unserem Dateisystem verschieben.
Das wird in Teil 2 dieses Kurses nützlich sein.

#### 1.2.4. Alignment-Befehl ausführen

Jetzt können wir den Alignment-Befehl ausführen, der den Alignment-Schritt mit `hisat2` durchführt und dann die Ausgabe an `samtools` weiterleitet, um die Ausgabe als BAM-Datei zu schreiben.

Die Read-Dateneingabe ist die Datei `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz`, die wir im vorherigen Schritt mit `trim_galore` generiert haben.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

Dies läuft fast sofort ab, weil es sich um eine sehr kleine Testdatei handelt.
In voller Größe könnte dies viel länger dauern.

Wieder kannst du die Ausgabedateien im Arbeitsverzeichnis finden:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Verzeichnisinhalt"

    ```console title="Ausgabe"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

Das Alignment hat eine BAM-Datei und eine Log-Datei mit Alignment-Statistiken erstellt.

#### 1.2.5. Ausgabedateien verschieben

Wie zuvor verschieben wir die Ausgabedateien in ein Verzeichnis auf dem gemounteten Dateisystem, damit sie nach dem Beenden des Containers zugänglich bleiben.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Damit haben wir alles, was wir brauchen.

#### 1.2.6. Container beenden

Um den Container zu beenden, gib `exit` ein.

```bash
exit
```

Dein Prompt sollte wieder normal werden.
Damit ist der Testlauf für die Verarbeitung einzelner Proben abgeschlossen.

!!! example "Schreibe es als Workflow!"

    Du kannst direkt zu [Teil 2](./02_single-sample.md) übergehen, wenn du sofort mit der Implementierung dieser Analyse als Nextflow-Workflow beginnen möchtest.
    Du musst nur zurückkommen, um die zweite Testrunde abzuschließen, bevor du zu Teil 3 übergehst.

---

## 2. QC-Aggregation über mehrere Proben

Die Befehle, die wir gerade getestet haben, verarbeiten jeweils eine Probe.
In der Praxis müssen wir normalerweise viele Proben verarbeiten und dann die QC-Ergebnisse über alle Proben hinweg aggregieren, um die Qualität des gesamten Datensatzes zu bewerten.

[MultiQC](https://multiqc.info/) ist ein Tool, das Verzeichnisse nach QC-Berichten vieler gängiger Bioinformatik-Tools durchsucht und sie in einem einzigen umfassenden HTML-Bericht aggregiert.
Es kann Ausgaben von FastQC, Cutadapt (über Trim Galore) und HISAT2 sowie vielen anderen erkennen.

Hier verarbeiten wir zwei zusätzliche Proben durch dieselben Tools pro Probe und verwenden dann MultiQC, um QC-Berichte über alle drei Proben hinweg zu aggregieren.
Dies sind die Befehle, die wir in Teil 3 dieses Kurses in einen Nextflow-Workflow verpacken werden.

1. QC und Trimming auf zusätzlichen Proben mit Trim Galore durchführen
2. Alignment auf zusätzlichen Proben mit HISAT2 durchführen
3. Alle QC-Berichte mit MultiQC in einen umfassenden Bericht aggregieren

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC und Trimming zusätzlicher Proben

Die QC- und Trimming-Befehle pro Probe sind identisch mit denen, die wir in Abschnitt 1.1 ausgeführt haben.
Wir haben das Container-Image bereits heruntergeladen, daher können wir es direkt starten.

#### 2.1.1. Container starten

Wir haben dieses Container-Image bereits in Abschnitt 1.1 heruntergeladen, daher können wir es direkt starten:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Dein Prompt ändert sich, um anzuzeigen, dass du dich im Container befindest.

#### 2.1.2. QC und Trimming auf zusätzlichen Proben durchführen

Führe FastQC und Trim Galore auf zwei weiteren Proben aus, eine nach der anderen.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Sobald dies abgeschlossen ist, solltest du Trim Galore-Ausgabedateien für beide Proben im Arbeitsverzeichnis haben.

#### 2.1.3. Ausgabedateien verschieben

Verschiebe die Trim Galore-Ausgabedateien in dasselbe Verzeichnis, das wir in Abschnitt 1 verwendet haben.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Verzeichnisinhalt"

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

Die Dateien sind jetzt in deinem normalen Dateisystem zugänglich.

#### 2.1.4. Container beenden

Um den Container zu beenden, gib `exit` ein.

```bash
exit
```

Dein Prompt sollte wieder normal werden.

### 2.2. Zusätzliche Proben alignieren

Die Alignment-Befehle sind identisch mit denen, die wir in Abschnitt 1.2 ausgeführt haben.
Wir müssen den Genom-Index aus dem Tarball extrahieren, den wir zuvor gespeichert haben, da die ursprünglichen Indexdateien in einem Container erstellt wurden, der nicht mehr existiert.

#### 2.2.1. Container starten

Wir haben dieses Container-Image bereits in Abschnitt 1.2 heruntergeladen, daher können wir es direkt starten:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Dein Prompt ändert sich, um anzuzeigen, dass du dich im Container befindest.

#### 2.2.2. Genom-Index extrahieren

Extrahiere die Genom-Indexdateien aus dem Tarball, den wir auf dem gemounteten Dateisystem gespeichert haben:

```bash
tar -xzf /data/genome_index.tar.gz
```

Dies stellt die Dateien `genome_index.*` im Arbeitsverzeichnis wieder her.

#### 2.2.3. Alignment auf zusätzlichen Proben durchführen

Führe das HISAT2-Alignment auf den beiden neu getrimmten Proben aus, eine nach der anderen.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Befehlsausgabe"

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

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Sobald dies abgeschlossen ist, solltest du BAM- und Log-Dateien für beide Proben im Arbeitsverzeichnis haben.

#### 2.2.4. Ausgabedateien verschieben

Verschiebe die Alignment-Ausgabedateien in dasselbe Verzeichnis, das wir in Abschnitt 1 verwendet haben.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Verzeichnisinhalt"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Die Dateien sind jetzt in deinem normalen Dateisystem zugänglich.

#### 2.2.5. Container beenden

Um den Container zu beenden, gib `exit` ein.

```bash
exit
```

Dein Prompt sollte wieder normal werden.

### 2.3. Umfassenden QC-Bericht generieren

Jetzt, da wir QC-, Trimming- und Alignment-Ausgaben für drei Proben haben, können wir MultiQC verwenden, um sie in einem einzigen Bericht zu aggregieren.
MultiQC durchsucht Verzeichnisse nach kompatiblen QC-Berichten und aggregiert alles, was es findet.

#### 2.3.1. Container herunterladen

Lass uns ein Container-Image herunterladen, das `multiqc` installiert hat:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Befehlsausgabe"

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

Du wirst bemerken, dass einige Layer `Already exists` anzeigen, weil sie mit den Container-Images geteilt werden, die wir zuvor heruntergeladen haben.
Daher sollte dieser Download schneller gehen als die vorherigen.

#### 2.3.2. Container interaktiv starten

Starte den Container interaktiv mit dem gemounteten Datenverzeichnis, wie zuvor.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Dein Prompt ändert sich, um anzuzeigen, dass du dich im Container befindest.

#### 2.3.3. MultiQC-Befehl ausführen

Führe `multiqc` aus und zeige auf die Verzeichnisse, in denen wir QC-bezogene Ausgabedateien für alle drei Proben gespeichert haben.
Das Flag `-n` legt den Namen des Ausgabeberichts fest.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Befehlsausgabe"

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

Hier sehen wir, dass das Tool QC-Berichte für alle drei Proben gefunden hat: die initiale QC von `fastqc`, die Post-Trimming-Berichte von `cutadapt` (über `trim_galore`) und die Alignment-Zusammenfassungen von `hisat2`.

Die Ausgabedateien befinden sich im Arbeitsverzeichnis:

```bash
ls all_samples_QC*
```

??? abstract "Verzeichnisinhalt"

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

Die Hauptausgabe ist der Bericht `all_samples_QC.html`, begleitet von einem Datenverzeichnis, das die zugrunde liegenden Metriken enthält.

#### 2.3.4. Ausgabedateien verschieben

Verschiebe den Bericht und sein Datenverzeichnis auf das gemountete Dateisystem.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Die Dateien sind jetzt in deinem normalen Dateisystem zugänglich.

#### 2.3.5. Container beenden

Um den Container zu beenden, gib `exit` ein.

```bash
exit
```

Dein Prompt sollte wieder normal werden.
Damit ist das Testen aller RNAseq-Verarbeitungsbefehle abgeschlossen.

---

### Fazit

Du weißt, wie du die Befehle FastQC, Trim Galore, HISAT2 und MultiQC in ihren jeweiligen Containern ausführst, einschließlich der Verarbeitung mehrerer Proben und der Aggregation von QC-Berichten.

### Wie geht es weiter?

Mach eine Pause und gehe dann zu [Teil 2](./02_single-sample.md), um zu lernen, wie du dieselben Befehle in Workflows verpackst, die Container zur Ausführung der Arbeit verwenden.
