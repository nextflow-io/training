# Część 1: Przegląd metody i manualne testowanie

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Istnieje wiele prawidłowych metod przetwarzania i analizy danych bulk RNAseq.
W tym szkoleniu stosujemy metodę opisaną [tutaj](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) przez dr. Simona Andrewsa i dr. Laurę Biggins z [Babraham Institute](https://www.babraham.ac.uk/).

Naszym celem jest opracowanie workflow'u implementującego następujące etapy przetwarzania: przeprowadzenie wstępnej kontroli jakości odczytów w próbce bulk RNAseq, przycięcie sekwencji adapterów z odczytów, dopasowanie odczytów do genomu referencyjnego i wygenerowanie kompleksowego raportu kontroli jakości (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** Przeprowadzenie QC na danych odczytów przed przycięciem przy użyciu FastQC
- **TRIM_GALORE:** Przycięcie sekwencji adapterów i przeprowadzenie QC po przycięciu przy użyciu Trim Galore (łączy Cutadapt i FastQC)
- **HISAT2_ALIGN:** Dopasowanie odczytów do genomu referencyjnego przy użyciu Hisat2
- **MULTIQC:** Wygenerowanie kompleksowego raportu QC przy użyciu MultiQC

Zanim jednak przejdziemy do pisania jakiegokolwiek kodu workflow'u, wypróbujemy polecenia ręcznie na danych testowych.
Narzędzia, których potrzebujemy, nie są zainstalowane w środowisku GitHub Codespaces, więc użyjemy ich za pośrednictwem kontenerów (zobacz [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Uwaga"

     Upewnij się, że jesteś w katalogu `nf4-science/rnaseq`. Ostatnia część ścieżki wyświetlana po wpisaniu `pwd` powinna brzmieć `rnaseq`.

---

## 1. Wstępne QC i przycinanie adapterów

Pobierzemy obraz kontenera z zainstalowanymi narzędziami `fastqc` i `trim_galore`, uruchomimy go interaktywnie i wykonamy polecenia przycinania oraz QC na jednym z przykładowych plików danych.

### 1.1. Pobranie kontenera

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Otrzymasz następujący wynik w konsoli podczas pobierania obrazu:

??? success "Wynik polecenia"

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

### 1.2. Uruchomienie kontenera interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Wynik polecenia"

    ```console

    ```
-->

Twój prompt zmieni się na coś w rodzaju `(base) root@b645838b3314:/tmp#`, co oznacza, że jesteś teraz wewnątrz kontenera.

Część `-v ./data:/data` polecenia umożliwi nam dostęp do zawartości katalogu `data/` z wnętrza kontenera.

```bash
ls /data/reads
```

??? success "Wynik polecenia"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. Uruchomienie pierwszego polecenia `fastqc`

Uruchommy `fastqc`, aby zebrać metryki kontroli jakości danych odczytów.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Wynik polecenia"

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

<!-- switch to tree -->

```console title="Wyjście"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Przycięcie sekwencji adapterów za pomocą `trim_galore`

Teraz uruchommy `trim_galore`, który łączy Cutadapt i FastQC, aby przyciąć sekwencje adapterów i zebrać metryki QC po przycięciu.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

Flaga `--fastqc` powoduje automatyczne uruchomienie etapu zbierania QC po zakończeniu przycinania.

_Wynik jest bardzo obszerny, więc poniżej przedstawiono wersję skróconą._

??? success "Wynik polecenia"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Pliki wyjściowe znajdziesz w katalogu roboczym:

```bash
ls ENCSR000COQ1_1*
```

```console title="Wyjście"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Przeniesienie plików wyjściowych do systemu plików poza kontenerem

Wszystko, co pozostanie wewnątrz kontenera, będzie niedostępne w przyszłości, więc przenieśmy te pliki do nowego katalogu.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Wyjście z kontenera

```bash
exit
```

---

## 2. Dopasowanie odczytów do genomu referencyjnego

Pobierzemy obraz kontenera z zainstalowanym `hisat2`, uruchomimy go interaktywnie i wykonamy polecenie dopasowania danych RNAseq do genomu referencyjnego.

### 2.1. Pobranie kontenera `hisat2`

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Wynik polecenia"

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

### 2.2. Uruchomienie kontenera `hisat2` interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Polecenie jest takie samo jak wcześniej, z odpowiednim URI kontenera.

### 2.3. Utworzenie plików indeksu genomu Hisat2

Hisat2 wymaga dostarczenia genomu referencyjnego w bardzo konkretnym formacie i nie może po prostu korzystać z pliku FASTA `genome.fa`, który dostarczamy, więc wykorzystamy tę okazję do utworzenia odpowiednich zasobów.

```bash
hisat2-build /data/genome.fa genome_index
```

Wynik jest bardzo obszerny, więc poniżej przedstawiono wersję skróconą:

<!-- TODO: switch to full output -->

??? success "Wynik polecenia"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

To tworzy wiele plików indeksu genomu, które znajdziesz w katalogu roboczym.

```bash
ls genome_index.*
```

```console title="Wyjście"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Użyjemy ich za chwilę, ale najpierw wygenerujmy skompresowany tarball zawierający te pliki indeksu genomu; będziemy ich potrzebować później, a ich generowanie nie jest zazwyczaj czymś, co chcemy robić w ramach workflow'u.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

To zapisuje tarball `genome_index.tar.gz` zawierający pliki indeksu genomu w katalogu `data/` w naszym systemie plików, co przyda się w Części 2 tego szkolenia.

### 2.4. Uruchomienie polecenia `hisat2`

Teraz możemy uruchomić polecenie dopasowania, które wykonuje etap dopasowania za pomocą `hisat2`, a następnie przekazuje wynik do `samtools`, aby zapisać wyjście jako plik BAM.

Wejściem danych odczytu jest plik `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz`, który wygenerowaliśmy za pomocą `trim_galore` w poprzednim kroku.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Wynik polecenia"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Działa to niemal natychmiast, ponieważ jest to bardzo mały plik testowy.
W rzeczywistej skali może to potrwać znacznie dłużej.

Ponownie pliki wyjściowe znajdziesz w katalogu roboczym:

```bash
ls ENCSR000COQ1_1*
```

```console title="Wyjście"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Przeniesienie plików wyjściowych do systemu plików poza kontenerem

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Wyjście z kontenera

```bash
exit
```

---

## 3. Wygenerowanie kompleksowego raportu QC

Pobierzemy obraz kontenera z zainstalowanym `multiqc`, uruchomimy go interaktywnie i wykonamy polecenie generowania raportu na podstawie plików raportów FastQC przed i po przycięciu.

### 3.1. Pobranie kontenera `multiqc`

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Wynik polecenia"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
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
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. Uruchomienie kontenera `multiqc` interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. Uruchomienie polecenia `multiqc`

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Wynik polecenia"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC potrafi przeszukiwać katalogi w poszukiwaniu kompatybilnych raportów QC i zagreguje wszystko, co znajdzie.

Widzimy tu, że narzędzie znalazło wszystkie trzy raporty QC, które wygenerowaliśmy: wstępny QC wykonany za pomocą `fastqc`, raport po przycięciu z `cutadapt` (wykonany za pośrednictwem `trim_galore`) oraz QC po dopasowaniu wygenerowany przez `hisat2`.

Pliki wyjściowe są ponownie w katalogu roboczym:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Wyjście"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Przeniesienie plików wyjściowych do systemu plików poza kontenerem

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Wyjście z kontenera

```bash
exit
```

---

### Podsumowanie

Przetestowałeś wszystkie poszczególne polecenia interaktywnie w odpowiednich kontenerach.

### Co dalej?

Dowiedz się, jak opakować te same polecenia w wieloetapowy workflow, który wykorzystuje kontenery do wykonywania zadań.
