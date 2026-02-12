# Część 1: Przegląd metody

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Istnieje wiele prawidłowych metod przetwarzania i analizy danych bulk RNAseq.
W tym szkoleniu stosujemy metodę opisaną [tutaj](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) przez dr. Simona Andrewsa i dr. Laurę Biggins z [Babraham Institute](https://www.babraham.ac.uk/).

Naszym celem jest opracowanie workflow'u implementującego następujące etapy przetwarzania: przeprowadzenie wstępnej kontroli jakości odczytów w próbce bulk RNAseq, przycięcie sekwencji adapterów z odczytów, dopasowanie odczytów do genomu referencyjnego i wygenerowanie kompleksowego raportu kontroli jakości (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Przeprowadzenie QC na danych odczytów przed przycięciem przy użyciu FastQC
- **TRIM_GALORE:** Przycięcie sekwencji adapterów i przeprowadzenie QC po przycięciu przy użyciu Trim Galore (łączy Cutadapt i FastQC)
- **HISAT2_ALIGN:** Dopasowanie odczytów do genomu referencyjnego przy użyciu Hisat2
- **MULTIQC:** Wygenerowanie kompleksowego raportu QC przy użyciu MultiQC

### Metody

Pokażemy Ci, jak zastosować te etapy przetwarzania w dwóch fazach.
Najpierw zaczniemy od **przetwarzania pojedynczej próbki**, które uruchamia narzędzia QC, przycinania i dopasowania na jednej próbce.
Następnie rozszerzymy to do **przetwarzania wielu próbek**, które uruchamia te same narzędzia na wielu próbkach i generuje zagregowany raport kontroli jakości.

Zanim jednak przejdziemy do pisania jakiegokolwiek kodu workflow'u dla któregokolwiek z tych podejść, wypróbujemy polecenia ręcznie na danych testowych.

### Zbiór danych

Dostarczamy następujące dane i powiązane zasoby:

- **Dane RNAseq** (`reads/`): pliki FASTQ z sześciu próbek, ograniczone do małego regionu w celu zmniejszenia rozmiaru plików. Każda próbka ma odczyty sparowane (dwa pliki na próbkę), choć zaczynamy od pracy tylko z odczytami jednostronnymi.
- **Genom referencyjny** (`genome.fa`): mały region ludzkiego chromosomu 20 (z hg19/b37).
- **Arkusze CSV z próbkami** (`single-end.csv` i `paired-end.csv`): pliki zawierające identyfikatory i ścieżki przykładowych plików danych.

### Oprogramowanie

Cztery główne narzędzia to [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) do zbierania metryk kontroli jakości, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) do przycinania adapterów (łączy Cutadapt i FastQC dla QC po przycięciu), [HISAT2](http://daehwankimlab.github.io/hisat2/) do dopasowania ze splajsowaniem do genomu referencyjnego oraz [MultiQC](https://multiqc.info/) do generowania zagregowanych raportów QC.

Narzędzia te nie są zainstalowane w środowisku GitHub Codespaces, więc użyjemy ich za pośrednictwem kontenerów pobranych z usługi Seqera Containers (zobacz [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Wskazówka"

     Upewnij się, że jesteś w katalogu `nf4-science/rnaseq`. Ostatnia część ścieżki wyświetlana po wpisaniu `pwd` powinna brzmieć `rnaseq`.

---

## 1. Przetwarzanie pojedynczej próbki

W tej sekcji testujemy polecenia przetwarzające pojedynczą próbkę RNAseq: kontrolę jakości, przycinanie adapterów i dopasowanie do genomu referencyjnego.
To są polecenia, które opakujemy w workflow Nextflow'a w Części 2 tego szkolenia.

1. Uruchomienie wstępnego QC na pliku FASTQ przy użyciu FastQC
2. Przycięcie sekwencji adapterów i uruchomienie QC po przycięciu przy użyciu Trim Galore
3. Dopasowanie przyciętych odczytów do genomu referencyjnego przy użyciu HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Zaczynamy od przetestowania tych poleceń na tylko jednej próbce.

### 1.1. QC i przycinanie adapterów

Najpierw chcemy uruchomić polecenia QC i przycinania na jednym z przykładowych plików danych.

#### 1.1.1. Pobranie kontenera

Pobierzmy obraz kontenera z zainstalowanymi narzędziami `fastqc` i `trim_galore`:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Wyjście polecenia"

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

Jeśli nie pobierałeś wcześniej tego obrazu, może to potrwać minutę.
Po zakończeniu będziesz mieć lokalną kopię obrazu kontenera.

#### 1.1.2. Uruchomienie kontenera interaktywnie

Aby uruchomić kontener interaktywnie, użyj `docker run` z flagami `-it`.
Opcja `-v ./data:/data` montuje nasz lokalny katalog `data/`, dzięki czemu możemy uzyskać dostęp do plików wejściowych z wnętrza kontenera.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Wyjście polecenia"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Twój prompt zmieni się na coś w rodzaju `(base) root@b645838b3314:/tmp#`, co oznacza, że jesteś teraz wewnątrz kontenera.

Sprawdź, czy widzisz pliki danych sekwencji w `/data/reads`:

```bash
ls /data/reads
```

??? abstract "Zawartość katalogu"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Teraz jesteś gotowy do wypróbowania pierwszego polecenia.

#### 1.1.3. Uruchomienie polecenia FastQC

Metoda, do której odwołujemy się powyżej, podaje nam wiersz poleceń do uruchomienia QC na pojedynczym pliku.
Musimy tylko podać plik wejściowy; narzędzie automatycznie wygeneruje pliki wyjściowe w tym samym katalogu co oryginalne dane.

Uruchom polecenie `fastqc` na jednym pliku danych:

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Wyjście polecenia"

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

Powinno to działać bardzo szybko.
Pliki wyjściowe znajdziesz w tym samym katalogu co oryginalne dane:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "Zawartość katalogu"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Powinieneś zobaczyć raport HTML i archiwum ZIP zawierające metryki QC.
To kończy testowanie pierwszego etapu.

#### 1.1.4. Przycięcie sekwencji adapterów za pomocą Trim Galore

Teraz uruchommy `trim_galore`, który łączy Cutadapt i FastQC, aby przyciąć sekwencje adapterów i zebrać metryki QC po przycięciu.
Jak wspomniano powyżej, oprogramowanie jest zawarte w tym samym kontenerze, więc nie trzeba nic zmieniać.

Polecenie jest proste; wystarczy dodać flagę `--fastqc`, aby automatycznie uruchomić etap zbierania QC po zakończeniu przycinania.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Wyjście polecenia"

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

Wynik jest bardzo obszerny, więc podświetliliśmy najbardziej istotne linie w powyższym przykładzie.
Pliki wyjściowe znajdziesz w katalogu roboczym:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Zawartość katalogu"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Obejmuje to przycięte odczyty, raport przycinania i pliki QC po przycięciu.

#### 1.1.5. Przeniesienie plików wyjściowych

Wszystko, co pozostanie wewnątrz kontenera, będzie niedostępne dla przyszłej pracy, więc musimy przenieść te pliki do katalogu w zamontowanym systemie plików.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Zawartość katalogu"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Pliki są teraz dostępne w Twoim normalnym systemie plików.

#### 1.1.6. Wyjście z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normy; to kończy testowanie pierwszych dwóch etapów.

### 1.2. Dopasowanie odczytów do genomu referencyjnego

Następnie chcemy uruchomić polecenie dopasowania, aby dopasować przycięte odczyty RNAseq do genomu referencyjnego.

#### 1.2.1. Pobranie kontenera

Pobierzmy obraz kontenera z zainstalowanymi narzędziami `hisat2` i `samtools`:

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Wyjście polecenia"

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

Zauważysz, że niektóre warstwy pokazują `Already exists`, ponieważ są współdzielone z obrazem kontenera Trim Galore, który pobraliśmy wcześniej.
W rezultacie to pobieranie powinno przebiec szybciej niż pierwsze.

#### 1.2.2. Uruchomienie kontenera interaktywnie

Uruchom kontener interaktywnie, używając tego samego podejścia co wcześniej, z odpowiednim URI kontenera.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Twój prompt ponownie się zmieni, wskazując, że jesteś wewnątrz kontenera.

#### 1.2.3. Utworzenie plików indeksu genomu

HISAT2 wymaga dostarczenia genomu referencyjnego w bardzo konkretnym formacie i nie może po prostu korzystać z pliku FASTA `genome.fa`, który dostarczamy, więc wykorzystamy tę okazję do utworzenia odpowiednich zasobów.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Wyjście polecenia"

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

Wynik jest bardzo obszerny, więc podświetliliśmy niektóre istotne linie w powyższym przykładzie.

To tworzy wiele plików indeksu genomu, które znajdziesz w katalogu roboczym.

```bash
ls genome_index.*
```

??? abstract "Zawartość katalogu"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Będziemy potrzebować tych plików później, a ich generowanie nie jest zazwyczaj czymś, co chcemy robić w ramach workflow'u, więc wygenerujemy skompresowany tarball zawierający pliki indeksu genomu, który możemy łatwo przekazywać w razie potrzeby.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Wyjście polecenia"

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

Za kilka minut przeniesiemy powstały tarball `genome_index.tar.gz` zawierający pliki indeksu genomu do katalogu `data/` w naszym systemie plików.
Przyda się to w Części 2 tego szkolenia.

#### 1.2.4. Uruchomienie polecenia dopasowania

Teraz możemy uruchomić polecenie dopasowania, które wykonuje etap dopasowania za pomocą `hisat2`, a następnie przekazuje wynik do `samtools`, aby zapisać wyjście jako plik BAM.

Wejściem danych odczytu jest plik `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz`, który wygenerowaliśmy za pomocą `trim_galore` w poprzednim kroku.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

Działa to niemal natychmiast, ponieważ jest to bardzo mały plik testowy.
W pełnej skali może to potrwać znacznie dłużej.

Ponownie pliki wyjściowe znajdziesz w katalogu roboczym:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Zawartość katalogu"

    ```console title="Wyjście"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

Dopasowanie wygenerowało plik BAM i plik dziennika ze statystykami dopasowania.

#### 1.2.5. Przeniesienie plików wyjściowych

Jak wcześniej, przenieś pliki wyjściowe do katalogu w zamontowanym systemie plików, aby pozostały dostępne po wyjściu z kontenera.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Mając to za sobą, mamy wszystko, czego potrzebujemy.

#### 1.2.6. Wyjście z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normy.
To kończy testowanie przetwarzania pojedynczej próbki.

!!! example "Napisz to jako workflow!"

    Możesz od razu przejść do [Części 2](./02_single-sample.md), jeśli chcesz zacząć implementować tę analizę jako workflow Nextflow'a.
    Będziesz musiał tylko wrócić, aby ukończyć drugą rundę testowania przed przejściem do Części 3.

---

## 2. Agregacja QC dla wielu próbek

Polecenia, które właśnie przetestowaliśmy, przetwarzają jedną próbkę na raz.
W praktyce zazwyczaj musimy przetworzyć wiele próbek, a następnie zagregować wyniki QC dla wszystkich z nich, aby ocenić jakość całego zbioru danych.

[MultiQC](https://multiqc.info/) to narzędzie przeszukujące katalogi w poszukiwaniu raportów QC z wielu popularnych narzędzi bioinformatycznych i agregujące je w jeden kompleksowy raport HTML.
Potrafi rozpoznać wyjście z FastQC, Cutadapt (za pośrednictwem Trim Galore) i HISAT2, między innymi.

Tutaj przetwarzamy dwie dodatkowe próbki przez te same narzędzia dla pojedynczych próbek, a następnie używamy MultiQC do agregacji raportów QC dla wszystkich trzech próbek.
To są polecenia, które opakujemy w workflow Nextflow'a w Części 3 tego szkolenia.

1. Uruchomienie QC i przycinania na dodatkowych próbkach przy użyciu Trim Galore
2. Uruchomienie dopasowania na dodatkowych próbkach przy użyciu HISAT2
3. Agregacja wszystkich raportów QC w kompleksowy raport przy użyciu MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC i przycinanie dodatkowych próbek

Polecenia QC i przycinania dla pojedynczych próbek są identyczne z tymi, które uruchomiliśmy w sekcji 1.1.
Już pobraliśmy obraz kontenera, więc możemy go uruchomić bezpośrednio.

#### 2.1.1. Uruchomienie kontenera

Już pobraliśmy ten obraz kontenera w sekcji 1.1, więc możemy go uruchomić bezpośrednio:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Twój prompt zmienia się, wskazując, że jesteś wewnątrz kontenera.

#### 2.1.2. Uruchomienie QC i przycinania na dodatkowych próbkach

Uruchom FastQC i Trim Galore na dwóch kolejnych próbkach, jednej po drugiej.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Po zakończeniu powinieneś mieć pliki wyjściowe Trim Galore dla obu próbek w katalogu roboczym.

#### 2.1.3. Przeniesienie plików wyjściowych

Przenieś pliki wyjściowe Trim Galore do tego samego katalogu, którego użyliśmy w sekcji 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Zawartość katalogu"

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

Pliki są teraz dostępne w Twoim normalnym systemie plików.

#### 2.1.4. Wyjście z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normy.

### 2.2. Dopasowanie dodatkowych próbek

Polecenia dopasowania są identyczne z tymi, które uruchomiliśmy w sekcji 1.2.
Musimy wyodrębnić indeks genomu z tarball'a, który zapisaliśmy wcześniej, ponieważ oryginalne pliki indeksu zostały utworzone wewnątrz kontenera, który już nie istnieje.

#### 2.2.1. Uruchomienie kontenera

Już pobraliśmy ten obraz kontenera w sekcji 1.2, więc możemy go uruchomić bezpośrednio:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Twój prompt zmienia się, wskazując, że jesteś wewnątrz kontenera.

#### 2.2.2. Wyodrębnienie indeksu genomu

Wyodrębnij pliki indeksu genomu z tarball'a, który zapisaliśmy w zamontowanym systemie plików:

```bash
tar -xzf /data/genome_index.tar.gz
```

To przywraca pliki `genome_index.*` w katalogu roboczym.

#### 2.2.3. Uruchomienie dopasowania na dodatkowych próbkach

Uruchom dopasowanie HISAT2 na dwóch nowo przyciętych próbkach, jednej po drugiej.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Wyjście polecenia"

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

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Po zakończeniu powinieneś mieć pliki BAM i dziennika dla obu próbek w katalogu roboczym.

#### 2.2.4. Przeniesienie plików wyjściowych

Przenieś pliki wyjściowe dopasowania do tego samego katalogu, którego użyliśmy w sekcji 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Zawartość katalogu"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Pliki są teraz dostępne w Twoim normalnym systemie plików.

#### 2.2.5. Wyjście z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normy.

### 2.3. Wygenerowanie kompleksowego raportu QC

Teraz, gdy mamy wyjście QC, przycinania i dopasowania dla trzech próbek, możemy użyć MultiQC do zagregowania ich w jeden raport.
MultiQC przeszukuje katalogi w poszukiwaniu kompatybilnych raportów QC i agreguje wszystko, co znajdzie.

#### 2.3.1. Pobranie kontenera

Pobierzmy obraz kontenera z zainstalowanym `multiqc`:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Wyjście polecenia"

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

Zauważysz, że niektóre warstwy pokazują `Already exists`, ponieważ są współdzielone z obrazami kontenerów, które pobraliśmy wcześniej.
W rezultacie to pobieranie powinno przebiec szybciej niż poprzednie.

#### 2.3.2. Uruchomienie kontenera interaktywnie

Uruchom kontener interaktywnie z zamontowanym katalogiem danych, jak wcześniej.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Twój prompt zmieni się, wskazując, że jesteś wewnątrz kontenera.

#### 2.3.3. Uruchomienie polecenia MultiQC

Uruchom `multiqc`, wskazując mu katalogi, w których przechowywaliśmy pliki wyjściowe związane z QC dla wszystkich trzech próbek.
Flaga `-n` ustawia nazwę raportu wyjściowego.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Wyjście polecenia"

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

Widzimy tu, że narzędzie znalazło raporty QC dla wszystkich trzech próbek: wstępny QC z `fastqc`, raporty po przycięciu z `cutadapt` (za pośrednictwem `trim_galore`) oraz podsumowania dopasowania z `hisat2`.

Pliki wyjściowe są w katalogu roboczym:

```bash
ls all_samples_QC*
```

??? abstract "Zawartość katalogu"

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

Głównym wyjściem jest raport `all_samples_QC.html`, któremu towarzyszy katalog danych zawierający podstawowe metryki.

#### 2.3.4. Przeniesienie plików wyjściowych

Przenieś raport i jego katalog danych do zamontowanego systemu plików.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Pliki są teraz dostępne w Twoim normalnym systemie plików.

#### 2.3.5. Wyjście z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normy.
To kończy testowanie wszystkich poleceń przetwarzania RNAseq.

---

### Podsumowanie

Wiesz, jak uruchamiać polecenia FastQC, Trim Galore, HISAT2 i MultiQC w ich odpowiednich kontenerach, w tym jak przetwarzać wiele próbek i agregować raporty QC.

### Co dalej?

Zrób sobie przerwę, a następnie przejdź do [Części 2](./02_single-sample.md), aby dowiedzieć się, jak opakować te same polecenia w workflow'y wykorzystujące kontenery do wykonywania pracy.
