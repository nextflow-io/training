# Część 3: Implementacja dla wielu próbek z danymi paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej ostatniej części kursu podniesiemy nasz prosty workflow na wyższy poziom, zamieniając go w potężne narzędzie automatyzacji wsadowej obsługujące dowolną liczbę próbek.
Przy okazji zmienimy go również tak, aby oczekiwał danych paired-end, które są bardziej powszechne we współczesnych badaniach.

Zrobimy to w trzech etapach:

1. Dostosujemy workflow do akceptowania wielu próbek wejściowych i paralelizacji wykonania
2. Dodamy kompleksowe generowanie raportów QC
3. Przełączymy się na dane RNAseq paired-end

---

## 1. Dostosowanie workflow'u do akceptowania wielu próbek wejściowych i paralelizacji wykonania

Musimy zmienić sposób zarządzania danymi wejściowymi.

### 1.1. Zmiana głównego wejścia na plik CSV ze ścieżkami plików zamiast pojedynczego pliku

W katalogu `data/` udostępniamy plik CSV zawierający identyfikatory próbek oraz ścieżki do plików FASTQ.
Ten plik CSV zawiera linię nagłówkową.
Zwróć uwagę, że ścieżki do plików FASTQ są ścieżkami bezwzględnymi.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Zmieńmy nazwę głównego parametru wejściowego na `input_csv` i zmieńmy wartość domyślną na ścieżkę do pliku `single-end.csv`.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Główne dane wejściowe
    input_csv: Path = "data/single-end.csv"

    // Archiwum genomu referencyjnego
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Aktualizacja fabryki kanału wejściowego w celu obsługi pliku CSV jako danych wejściowych

Chcemy wczytać zawartość pliku do kanału zamiast samej ścieżki pliku, więc używamy operatora `.splitCsv()` do parsowania formatu CSV, a następnie operatora `.map()` do pobrania konkretnej informacji, której potrzebujemy (ścieżki do pliku FASTQ).

```groovy title="rnaseq.nf" linenums="16"
    // Utwórz kanał wejściowy z zawartości pliku CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Uruchomienie workflow'u w celu sprawdzenia, czy działa

```bash
nextflow run rnaseq.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Tym razem widzimy, że każdy krok jest uruchamiany 6 razy – na każdym z 6 dostarczonych plików danych.

To wszystko, czego potrzeba, aby workflow uruchamiał się na wielu plikach!
Nextflow obsługuje całą paralelizację za nas.

---

## 2. Agregacja metryk QC wstępnego przetwarzania w pojedynczy raport MultiQC

Generuje to wiele raportów QC i nie chcemy musieć przeglądać poszczególnych raportów.
To idealny moment, aby dodać krok agregacji raportów MultiQC!

### 2.1. Utworzenie modułu dla procesu agregacji QC

Stwórzmy plik modułu o nazwie `modules/multiqc.nf`, który będzie zawierał **process** `MULTIQC`:

```bash
touch modules/multiqc.nf
```

Otwórz plik w edytorze kodu i skopiuj do niego następujący kod:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

    input:
    path '*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}
```

### 2.2. Zaimportowanie modułu do pliku workflow'u

Dodaj instrukcję `include { MULTIQC } from './modules/multiqc.nf'` do pliku `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Instrukcje INCLUDE modułów
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Dodanie parametru `report_id` z rozsądną wartością domyślną

```groovy title="rnaseq.nf" linenums="9"
params {
    // Główne dane wejściowe
    input_csv: Path = "data/single-end.csv"

    // Archiwum genomu referencyjnego
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID raportu
    report_id: String = "all_single-end"
}
```

### 2.4. Wywołanie procesu na wyjściach poprzednich kroków

Musimy przekazać procesowi `MULTIQC` wszystkie wyjścia związane z QC z poprzednich kroków.

W tym celu użyjemy operatora `.mix()`, który agreguje wiele kanałów w jeden.

Gdybyśmy mieli cztery procesy nazwane A, B, C i D, każdy z prostym kanałem `.out`, składnia wyglądałaby następująco: `A.out.mix( B.out, C.out, D.out )`. Jak widać, stosujesz go do pierwszego z kanałów, które chcesz połączyć (nie ma znaczenia którego) i po prostu dodajesz wszystkie pozostałe, oddzielone przecinkami, w nawiasie, który następuje.

W przypadku naszego workflow'u mamy następujące wyjścia do agregacji:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Więc przykładowa składnia staje się:

```groovy title="Zastosowanie .mix() w wywołaniu MULTIQC"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

To zbierze raporty QC dla każdej próbki.
Ponieważ jednak chcemy je zagregować dla wszystkich próbek, musimy dodać operator `collect()`, aby zebrać raporty dla wszystkich próbek w jedno wywołanie `MULTIQC`.
Musimy również przekazać parametr `report_id`.

Daje nam to następujący kod:

```groovy title="Ukończone wywołanie MULTIQC" linenums="33"
    // Generowanie kompleksowego raportu QC
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

W kontekście pełnego bloku workflow'u wygląda to tak:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Utwórz kanał wejściowy z zawartości pliku CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Początkowa kontrola jakości
    FASTQC(read_ch)

    // Przycinanie adapterów i kontrola jakości po przycięciu
    TRIM_GALORE(read_ch)

    // Dopasowanie do genomu referencyjnego
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Generowanie kompleksowego raportu QC
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Uruchomienie workflow'u w celu sprawdzenia, czy działa

```bash
nextflow run rnaseq.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Tym razem widzimy pojedyncze wywołanie MULTIQC dodane po zbuforowanych wywołaniach procesów:

Wyniki możesz znaleźć w katalogu `results/multiqc`, jak określono w procesie `MULTIQC` przez dyrektywę `publishDir`.

```bash
tree -L 2 results/multiqc
```

```console title="Wyjście"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Ten ostatni plik `all_single-end.html` to pełny zagregowany raport, wygodnie zapakowany w jeden łatwy do przeglądania plik HTML.

---

## 3. Umożliwienie przetwarzania danych RNAseq paired-end

Obecnie nasz workflow obsługuje tylko dane RNAseq single-end.
Coraz częściej spotyka się dane RNAseq paired-end, więc chcemy móc je obsługiwać.

Uczynienie workflow'u całkowicie niezależnym od typu danych wymagałoby użycia nieco bardziej zaawansowanych funkcji języka Nextflow, więc nie zrobimy tego tutaj, ale możemy stworzyć wersję do przetwarzania paired-end, aby zademonstrować, co należy dostosować.

### 3.1. Utworzenie kopii workflow'u o nazwie `rnaseq_pe.nf`

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Modyfikacja domyślnej wartości `input_csv`, aby wskazywała na dane paired-end

Udostępniamy drugi plik CSV zawierający identyfikatory próbek oraz ścieżki do sparowanych plików FASTQ w katalogu `data/`

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Zmieńmy wartość domyślną `input_csv` na ścieżkę do pliku `paired-end.csv`.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Główne dane wejściowe
    input_csv: Path = "data/paired-end.csv"

    // Archiwum genomu referencyjnego
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID raportu
    report_id: String = "all_single-end"
}
```

### 3.3. Aktualizacja fabryki kanału

Musimy teraz powiedzieć operatorowi `.map()`, aby pobierał obie ścieżki do plików FASTQ.

Więc `row -> file(row.fastq_path)` staje się `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Utwórz kanał wejściowy z zawartości pliku CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Utworzenie wersji paired-end procesu FASTQC

Utwórzmy kopię modułu, aby zachować obie wersje.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Otwórz nowy plik modułu `fastqc_pe.nf` w edytorze kodu i wprowadź następujące zmiany:

- Zmień `fastqc $reads` na `fastqc ${reads}` w bloku `script` (linia 17), aby wejście `reads` zostało rozpakowane, ponieważ jest teraz krotką dwóch ścieżek zamiast pojedynczej ścieżki.
- Zastąp `${reads.simpleName}` symbolem wieloznacznym (`*`), aby uniknąć konieczności indywidualnego obsługiwania plików wyjściowych.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Technicznie uogólnia to **process** `FASTQC` w sposób, który umożliwia mu obsługę zarówno danych RNAseq single-end, jak i paired-end.

Na koniec zaktualizuj instrukcję importu modułu, aby używała wersji paired-end modułu.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Utworzenie wersji paired-end procesu TRIM_GALORE

Utwórz kopię modułu, aby zachować obie wersje.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Otwórz nowy plik modułu `trim_galore_pe.nf` w edytorze kodu i wprowadź następujące zmiany:

- Zmień deklarację wejścia z `path reads` na `tuple path(read1), path(read2)`
- Zaktualizuj polecenie w bloku `script`, zastępując `$reads` przez `--paired ${read1} ${read2}`
- Zaktualizuj deklaracje wyjścia, aby odzwierciedlały dodane pliki i różne konwencje nazewnictwa, używając symboli wieloznacznych.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)

    output:
    tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
    path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

    script:
    """
    trim_galore --fastqc --paired ${read1} ${read2}
    """
```

Na koniec zaktualizuj instrukcję importu modułu, aby używała wersji paired-end modułu.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Aktualizacja wywołania procesu MULTIQC, aby oczekiwał dwóch raportów z TRIM_GALORE

**Process** `TRIM_GALORE` generuje teraz dodatkowy kanał wyjściowy, więc musimy przekazać go do MultiQC.

Zastąp `TRIM_GALORE.out.fastqc_reports,` przez `TRIM_GALORE.out.fastqc_reports_1,` plus `TRIM_GALORE.out.fastqc_reports_2,`:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Generowanie kompleksowego raportu QC
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Skoro już przy MultiQC, zaktualizujmy również wartość domyślną parametru `report_id` z `"all_single-end"` na `"all_paired-end"`.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Główne dane wejściowe
    input_csv: Path = "data/paired-end.csv"

    // Archiwum genomu referencyjnego
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID raportu
    report_id: String = "all_paired-end"
}
```

### 3.7. Utworzenie wersji paired-end procesu HISAT2_ALIGN

Stwórz kopię modułu, aby zachować obie wersje.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Otwórz nowy plik modułu `hisat2_align_pe.nf` w edytorze kodu i wprowadź następujące zmiany:

- Zmień deklarację wejścia z `path reads` na `tuple path(read1), path(read2)`
- Zaktualizuj polecenie w bloku `script`, zastępując `-U $reads` przez `-1 ${read1} -2 ${read2}`
- Zastąp wszystkie wystąpienia `${reads.simpleName}` przez `${read1.simpleName}` w poleceniu w bloku `script`, a także w deklaracjach wyjścia.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Na koniec zaktualizuj instrukcję importu modułu, aby używała wersji paired-end modułu.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Uruchomienie workflow'u w celu sprawdzenia, czy działa

Nie używamy `-resume`, ponieważ to nie zostałoby zbuforowane, a jest dwa razy więcej danych do przetworzenia niż wcześniej, ale i tak powinno zakończyć się w mniej niż minutę.

```bash
nextflow run rnaseq_pe.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

I to wszystko! Teraz mamy dwie nieco rozbieżne wersje naszego workflow'u, jedną dla danych single-end i jedną dla danych paired-end.
Następnym logicznym krokiem byłoby sprawienie, aby workflow akceptował którykolwiek typ danych w locie, co wykracza poza zakres tego kursu, ale możemy się tym zająć w kontynuacji.

---

### Podsumowanie

Wiesz, jak dostosować workflow dla pojedynczej próbki, aby sparalelizować przetwarzanie wielu próbek, wygenerować kompleksowy raport QC i dostosować workflow do używania danych paired-end, jeśli jest to potrzebne.

### Co dalej?

Gratulacje, ukończyłeś mini-kurs Nextflow For RNAseq! Świętuj Swój sukces i weź zasłużoną przerwę!

Następnie prosimy o wypełnienie bardzo krótkiej ankiety dotyczącej Twoich doświadczeń z tym kursem szkoleniowym, a następnie przekierujemy Cię na stronę z linkami do dalszych materiałów szkoleniowych i pomocnych odnośników.
