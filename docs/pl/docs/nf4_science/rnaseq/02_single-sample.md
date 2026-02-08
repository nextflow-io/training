# Część 2: Implementacja dla pojedynczej próbki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej części kursu napiszemy najprostszy możliwy workflow, który obejmie wszystkie polecenia uruchomione w Części 1, aby zautomatyzować ich wykonywanie. Będziemy dążyć do przetwarzania jednej próbki na raz.

Zrobimy to w trzech etapach:

1. Napisanie jednoetapowego workflow'u, który uruchamia początkowy krok kontroli jakości
2. Dodanie przycinania adapterów i kontroli jakości po przycięciu
3. Dodanie dopasowania do genomu referencyjnego

!!! warning "Wymaganie wstępne"

    Musisz przejść przez Część 1 kursu przed rozpoczęciem tej lekcji.
    W szczególności, praca przez sekcje 2.1-3 tworzy plik indeksu genomu (`data/genome_index.tar.gz`) wymagany do kroku dopasowania w tej lekcji.

---

## 1. Napisanie jednoetapowego workflow'u, który uruchamia początkową kontrolę jakości

Zacznijmy od napisania prostego workflow'u, który uruchamia narzędzie FastQC na pliku FASTQ zawierającym odczyty RNAseq pojedynczego końca.

Udostępniamy plik workflow'u, `rnaseq.nf`, który określa główne części workflow'u.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Instrukcje INCLUDE modułów

/*
 * Parametry pipeline'u
 */

// Główne wejście

workflow {

    // Utwórz kanał wejściowy

    // Wywołaj procesy

}
```

Pamiętaj, że ten kod workflow'u jest poprawny, ale nie jest funkcjonalny; jego celem jest jedynie służenie jako szkielet, którego użyjesz do napisania rzeczywistego workflow'u.

### 1.1. Utworzenie katalogu do przechowywania modułów

Utworzymy samodzielne moduły dla każdego procesu, aby ułatwić zarządzanie nimi i ponowne wykorzystanie, więc stwórzmy katalog do ich przechowywania.

```bash
mkdir modules
```

### 1.2. Utworzenie modułu dla procesu zbierania metryk kontroli jakości

Utwórzmy plik modułu o nazwie `modules/fastqc.nf` do przechowywania procesu `FASTQC`:

```bash
touch modules/fastqc.nf
```

Otwórz plik w edytorze kodu i skopiuj do niego następujący kod:

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

Powinieneś rozpoznać wszystkie elementy z tego, czego nauczyłeś się w Części 1 i Części 2 tej serii szkoleń; jedyną godną uwagi zmianą jest to, że tym razem używamy `mode: symlink` dla dyrektywy `publishDir` i używamy parametru do zdefiniowania `publishDir`.

!!! note "Uwaga"

    Mimo że pliki danych, których tutaj używamy, są bardzo małe, w genomice mogą być bardzo duże. W celach demonstracyjnych w środowisku szkoleniowym używamy trybu publikowania 'symlink', aby uniknąć niepotrzebnych kopii plików. Nie powinieneś tego robić w Swoich finalnych workflow'ach, ponieważ stracisz wyniki podczas czyszczenia katalogu `work`.

### 1.3. Importowanie modułu do pliku workflow'u

Dodaj instrukcję `include { FASTQC } from './modules/fastqc.nf'` do pliku `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Instrukcje INCLUDE modułów
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Dodanie deklaracji danych wejściowych

Zadeklaruj parametr wejściowy z wartością domyślną:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Główne wejście
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Utworzenie kanału wejściowego w bloku workflow

Użyj podstawowej fabryki kanałów `.fromPath()`, aby utworzyć kanał wejściowy:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Utwórz kanał wejściowy ze ścieżki pliku
    read_ch = channel.fromPath(params.reads)

    // Wywołaj procesy

}
```

### 1.6. Wywołanie procesu `FASTQC` na kanale wejściowym

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Utwórz kanał wejściowy ze ścieżki pliku
    read_ch = channel.fromPath(params.reads)

    // Początkowa kontrola jakości
    FASTQC(read_ch)

}
```

### 1.7. Uruchomienie workflow'u w celu przetestowania, czy działa

Moglibyśmy użyć parametru `--reads`, aby określić wejście z linii poleceń, ale podczas programowania możemy być leniwi i po prostu użyć domyślnej wartości testowej, którą ustawiliśmy.

```bash
nextflow run rnaseq.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Powinno to działać bardzo szybko, jeśli pracowałeś przez Część 1 i już pobrałeś kontener.
Jeśli ją pominąłeś, Nextflow pobierze kontener za Ciebie; nie musisz nic robić, aby to się stało, ale możesz potrzebować poczekać do minuty.

Możesz znaleźć wyniki w `results/fastqc`, jak określono w procesie `FASTQC` przez dyrektywę `publishDir`.

```bash
ls results/fastqc
```

```console title="Wyjście"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Dodanie przycinania adapterów i kontroli jakości po przycięciu

Zamierzamy użyć wrappera Trim_Galore, który łączy Cutadapt do samego przycinania i FastQC do kontroli jakości po przycięciu.

### 2.1. Utworzenie modułu dla procesu przycinania i kontroli jakości

Utwórzmy plik modułu o nazwie `modules/trim_galore.nf` do przechowywania procesu `TRIM_GALORE`:

```bash
touch modules/trim_galore.nf
```

Otwórz plik w edytorze kodu i skopiuj do niego następujący kod:

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

### 2.2. Importowanie modułu do pliku workflow'u

Dodaj instrukcję `include { TRIM_GALORE } from './modules/trim_galore.nf'` do pliku `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Instrukcje INCLUDE modułów
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Wywołanie procesu na kanale wejściowym

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Utwórz kanał wejściowy ze ścieżki pliku
    read_ch = channel.fromPath(params.reads)

    // Początkowa kontrola jakości
    FASTQC(read_ch)

    // Przycinanie adapterów i kontrola jakości po przycięciu
    TRIM_GALORE(read_ch)
}
```

### 2.4. Uruchomienie workflow'u w celu przetestowania, czy działa

```bash
nextflow run rnaseq.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

To również powinno działać bardzo szybko, ponieważ uruchamiamy na tak małym pliku wejściowym.

Możesz znaleźć wyniki w `results/trimming`, jak określono w procesie `TRIM_GALORE` przez dyrektywę `publishDir`.

```bash
ls results/trimming
```

```console title="Wyjście"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Dopasowanie odczytów do genomu referencyjnego

Na koniec możemy uruchomić krok dopasowania genomu przy użyciu Hisat2, który również wyemituje metryki kontroli jakości w stylu FastQC.

### 3.1. Utworzenie modułu dla procesu HiSat2

Utwórzmy plik modułu o nazwie `modules/hisat2_align.nf` do przechowywania procesu `HISAT2_ALIGN`:

```bash
touch modules/hisat2_align.nf
```

Otwórz plik w edytorze kodu i skopiuj do niego następujący kod:

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

### 3.2. Importowanie modułu do pliku workflow'u

Dodaj instrukcję `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` do pliku `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Instrukcje INCLUDE modułów
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Dodanie deklaracji parametru do dostarczenia indeksu genomu

Zadeklaruj parametr wejściowy z wartością domyślną:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Główne wejście
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Archiwum genomu referencyjnego
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Wywołanie procesu `HISAT2_ALIGN` na przyciętych odczytach wyjściowych z `TRIM_GALORE`

Przycięte odczyty znajdują się w kanale `TRIM_GALORE.out.trimmed_reads` wyjściowym z poprzedniego kroku.

Dodatkowo używamy `file (params.hisat2_index_zip)`, aby dostarczyć narzędziu Hisat2 spakowane archiwum tar indeksu genomu.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Utwórz kanał wejściowy ze ścieżki pliku
    read_ch = channel.fromPath(params.reads)

    // Początkowa kontrola jakości
    FASTQC(read_ch)

    // Przycinanie adapterów i kontrola jakości po przycięciu
    TRIM_GALORE(read_ch)

    // Dopasowanie do genomu referencyjnego
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Uruchomienie workflow'u w celu przetestowania, czy działa

```bash
nextflow run rnaseq.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Możesz znaleźć wyniki w `results/align`, jak określono w procesie `HISAT2_ALIGN` przez dyrektywę `publishDir`.

```bash
ls results/align
```

```console title="Wyjście"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

To kończy podstawowe przetwarzanie, które musimy zastosować do każdej próbki.

_Dodamy agregację raportów MultiQC w Części 3, po tym, jak sprawimy, że workflow będzie akceptować wiele próbek naraz._

---

### Podsumowanie

Wiesz, jak opakować wszystkie podstawowe kroki do przetwarzania próbek RNAseq single-end indywidualnie.

### Co dalej?

Dowiedz się, jak zmodyfikować workflow, aby przetwarzać wiele próbek równolegle, agregować raporty kontroli jakości we wszystkich krokach dla wszystkich próbek i umożliwić uruchamianie workflow'u na danych RNAseq paired-end.
