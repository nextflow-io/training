# Część 3: Implementacja dla wielu próbek z danymi paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wcześniej zbudowałeś pipeline do wywoływania wariantów dla pojedynczych próbek, który przetwarzał dane każdej próbki niezależnie.
W tej części kursu podniesiemy nasz prosty workflow na wyższy poziom, zamieniając go w potężne narzędzie automatyzacji wsadowej obsługujące dowolną liczbę próbek.
Przy okazji zaktualizujemy go również tak, aby oczekiwał danych paired-end, które są bardziej powszechne w nowszych badaniach.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś [Część 1: Przegląd metody](./01_method.md), [Część 2: Implementacja dla pojedynczej próbki](./02_single-sample.md) i masz działający pipeline `rnaseq.nf` z wypełnionymi plikami modułów.

    Jeśli nie ukończyłeś Części 2 lub chcesz zacząć od nowa w tej części, możesz użyć rozwiązania z Części 2 jako punktu wyjścia.
    Uruchom te polecenia z katalogu `nf4-science/rnaseq/`:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    To da Ci kompletny workflow przetwarzania pojedynczej próbki.
    Możesz sprawdzić, czy działa poprawnie:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Zadanie

W tej części kursu rozszerzymy workflow, aby wykonywał następujące czynności:

1. Odczytywał informacje o próbkach z arkusza CSV
2. Uruchamiał QC, przycinanie i dopasowanie dla każdej próbki równolegle
3. Agregował wszystkie raporty QC w kompleksowy raport MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

To automatyzuje kroki z drugiej sekcji [Części 1: Przegląd metody](./01_method.md#2-multi-sample-qc-aggregation), gdzie uruchamiałeś te polecenia ręcznie w ich kontenerach.

## Plan lekcji

Podzieliliśmy to na trzy etapy:

1. **Dostosowanie workflow'u do akceptowania wielu próbek wejściowych.**
   Obejmuje to przejście z pojedynczej ścieżki pliku na arkusz CSV, parsowanie go za pomocą `splitCsv()` i uruchamianie wszystkich istniejących procesów na wielu próbkach.
2. **Dodanie kompleksowego generowania raportów QC.**
   Wprowadza to operator `collect()` do agregacji wyjść z różnych próbek i dodaje proces MultiQC do wygenerowania połączonego raportu.
3. **Przełączenie na dane RNAseq paired-end.**
   Obejmuje to dostosowanie procesów do wejść paired-end (używając krotek), utworzenie modułów paired-end i skonfigurowanie osobnego profilu testowego.

To implementuje metodę opisaną w [Części 1: Przegląd metody](./01_method.md) (druga sekcja dotycząca przypadku wielu próbek) i bazuje bezpośrednio na workflow'ie stworzonym w Części 2.

!!! tip "Wskazówka"

     Upewnij się, że jesteś we właściwym katalogu roboczym:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Dostosowanie workflow'u do akceptowania wielu próbek wejściowych

Aby uruchomić workflow na wielu próbkach, musimy zmienić sposób zarządzania danymi wejściowymi: zamiast podawać pojedynczą ścieżkę pliku, odczytamy informacje o próbkach z pliku CSV.

Udostępniamy plik CSV zawierający identyfikatory próbek i ścieżki do plików FASTQ w katalogu `data/`.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Ten plik CSV zawiera linię nagłówkową, która nazywa kolumny.

Zwróć uwagę, że to wciąż dane single-end.

!!! warning "Ostrzeżenie"

    Ścieżki plików w CSV są ścieżkami bezwzględnymi, które muszą pasować do Twojego środowiska.
    Jeśli nie uruchamiasz tego w środowisku szkoleniowym, które udostępniamy, będziesz musiał zaktualizować ścieżki, aby pasowały do Twojego systemu.

### 1.1. Zmiana głównego wejścia na CSV ze ścieżkami plików w profilu testowym

Najpierw musimy zaktualizować profil testowy w `nextflow.config`, aby podawał ścieżkę do pliku CSV zamiast pojedynczej ścieżki FASTQ.

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Następnie musimy zaktualizować tworzenie kanału, aby odczytywał z tego CSV.

### 1.2. Aktualizacja fabryki kanału do parsowania wejścia CSV

Musimy wczytać zawartość pliku do kanału zamiast samej ścieżki pliku.

Możemy to zrobić używając tego samego wzorca, którego użyliśmy w [Części 2 Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): stosując operator [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) do parsowania pliku, a następnie operację `map` do wyodrębnienia ścieżki pliku FASTQ z każdego wiersza.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Utwórz kanał wejściowy z zawartości pliku CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)
    ```

Jedną rzeczą, która jest nowa w porównaniu z tym, co napotkałeś w kursie Hello Nextflow, jest to, że ten CSV ma linię nagłówkową, więc dodajemy `#!groovy header: true` do wywołania `splitCsv()`.
To pozwala nam odwoływać się do kolumn po nazwie w operacji `map`: `#!groovy row.fastq_path` wyodrębnia ścieżkę pliku z kolumny `fastq_path` każdego wiersza.

Obsługa wejścia jest zaktualizowana i workflow jest gotowy do testowania.

### 1.3. Uruchomienie workflow'u

Workflow teraz odczytuje informacje o próbkach z pliku CSV i przetwarza wszystkie próbki równolegle.

```bash
nextflow run rnaseq.nf -profile test
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

Tym razem każdy krok jest uruchamiany 6 razy, raz dla każdej próbki w pliku CSV.

To wszystko, czego potrzeba, aby workflow uruchamiał się na wielu plikach.
Nextflow obsługuje całą paralelizację za nas.

### Podsumowanie

Wiesz, jak przełączyć się z wejścia pojedynczego pliku na wejście wielu próbek oparte na CSV, które Nextflow przetwarza równolegle.

### Co dalej?

Dodaj krok agregacji raportów QC, który łączy metryki ze wszystkich próbek.

---

## 2. Agregacja metryk QC wstępnego przetwarzania w pojedynczy raport MultiQC

To wszystko generuje wiele raportów QC i nie chcemy musieć przeglądać poszczególnych raportów.
To idealny moment, aby dodać krok agregacji raportów MultiQC.

Przypomnij sobie polecenie `multiqc` z [Części 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

Polecenie skanuje bieżący katalog w poszukiwaniu rozpoznanych plików wyjściowych QC i agreguje je w pojedynczy raport HTML.
URI kontenera to `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Musimy skonfigurować dodatkowy parametr, przygotować wejścia, napisać proces, podłączyć go i zaktualizować obsługę wyjść.

### 2.1. Konfiguracja wejść

Proces MultiQC potrzebuje parametru nazwy raportu i zebranych wyjść QC ze wszystkich poprzednich kroków połączonych razem.

#### 2.1.1. Dodanie parametru `report_id`

Dodaj parametr do nazwania raportu wyjściowego.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Główne wejście
        input: Path

        // Archiwum genomu referencyjnego
        hisat2_index_zip: Path

        // ID raportu
        report_id: String
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Główne wejście
        input: Path

        // Archiwum genomu referencyjnego
        hisat2_index_zip: Path
    }
    ```

Dodaj domyślną wartość ID raportu do profilu testowego:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Następnie musimy przygotować wejścia dla procesu MultiQC.

#### 2.1.2. Zebranie i połączenie wyjść QC z poprzednich kroków

Musimy przekazać procesowi `MULTIQC` wszystkie wyjścia związane z QC z poprzednich kroków połączone razem.

W tym celu używamy operatora `.mix()`, który agreguje wiele kanałów w jeden.
Zaczynamy od `channel.empty()` i mieszamy wszystkie kanały wyjściowe, które chcemy połączyć.
Jest to czystsze niż łączenie `.mix()` bezpośrednio do jednego z kanałów wyjściowych, ponieważ traktuje wszystkie wejścia symetrycznie.

W naszym workflow'ie wyjścia związane z QC do agregacji to:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Mieszamy je w jeden kanał, a następnie używamy `.collect()` do agregacji raportów ze wszystkich próbek w jedną listę.

Dodaj te linie do ciała workflow'u po wywołaniu `HISAT2_ALIGN`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Dopasowanie do genomu referencyjnego
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Generowanie kompleksowego raportu QC
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="38"
        // Dopasowanie do genomu referencyjnego
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

Użycie zmiennych pośrednich sprawia, że każdy krok jest jasny: `multiqc_files_ch` zawiera wszystkie pojedyncze pliki QC zmieszane w jeden kanał, a `multiqc_files_list` to zebrana paczka gotowa do przekazania do MultiQC.

### 2.2. Napisanie procesu agregacji QC i wywołanie go w workflow'ie

Jak poprzednio, musimy wypełnić definicję procesu, zaimportować moduł i dodać wywołanie procesu.

#### 2.2.1. Wypełnienie modułu dla procesu agregacji QC

Otwórz `modules/multiqc.nf` i przejrzyj zarys definicji procesu.

Wypełnij definicję procesu samodzielnie, używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Agreguj raporty QC za pomocą MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Po"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Agreguj raporty QC za pomocą MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

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

Ten proces używa `#!groovy path '*'` jako kwalifikatora wejścia dla plików QC.
Symbol wieloznaczny `'*'` mówi Nextflow'owi, aby umieścił wszystkie zebrane pliki w katalogu roboczym bez wymagania konkretnych nazw.
Wejście `val output_name` to ciąg znaków, który kontroluje nazwę pliku raportu.

Polecenie `multiqc .` skanuje bieżący katalog (gdzie znajdują się wszystkie umieszczone pliki QC) i generuje raport.

Po ukończeniu tego proces jest gotowy do użycia.

#### 2.2.2. Zaimportowanie modułu

Dodaj instrukcję importu do `rnaseq.nf`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="3"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Teraz dodaj wywołanie procesu do workflow'u.

#### 2.2.3. Dodanie wywołania procesu

Przekaż zebrane pliki QC i ID raportu do procesu `MULTIQC`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

Proces MultiQC jest teraz podłączony do workflow'u.

### 2.3. Aktualizacja obsługi wyjść

Musimy dodać wyjścia MultiQC do deklaracji publikowania i skonfigurować, gdzie mają trafić.

#### 2.3.1. Dodanie celów publikowania dla wyjść MultiQC

Dodaj wyjścia MultiQC do sekcji `publish:`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

Następnie musimy powiedzieć Nextflow'owi, gdzie umieścić te wyjścia.

#### 2.3.2. Konfiguracja nowych celów wyjściowych

Dodaj wpisy dla celów MultiQC w bloku `output {}`, publikując je w podkatalogu `multiqc/`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="4-9"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

Konfiguracja wyjść jest ukończona.

### 2.4. Uruchomienie workflow'u

Używamy `-resume`, aby poprzednie kroki przetwarzania były buforowane i uruchamiał się tylko nowy krok MultiQC.

```bash
nextflow run rnaseq.nf -profile test -resume
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

Pojedyncze wywołanie MULTIQC zostało dodane po zbuforowanych wywołaniach procesów.

Wyjścia MultiQC możesz znaleźć w katalogu wyników.

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

### Podsumowanie

Wiesz, jak zbierać wyjścia z wielu kanałów, łączyć je za pomocą `.mix()` i `.collect()` oraz przekazywać je do procesu agregacji.

### Co dalej?

Dostosuj workflow do obsługi danych RNAseq paired-end.

---

## 3. Umożliwienie przetwarzania danych RNAseq paired-end

Obecnie nasz workflow obsługuje tylko dane RNAseq single-end.
Coraz częściej spotyka się dane RNAseq paired-end, więc chcemy móc je obsługiwać.

Uczynienie workflow'u całkowicie niezależnym od typu danych wymagałoby użycia nieco bardziej zaawansowanych funkcji języka Nextflow, więc nie zrobimy tego tutaj, ale możemy stworzyć wersję do przetwarzania paired-end, aby zademonstrować, co należy dostosować.

### 3.1. Skopiowanie workflow'u i aktualizacja wejść

Zaczynamy od skopiowania pliku workflow'u single-end i zaktualizowania go dla danych paired-end.

#### 3.1.1. Skopiowanie pliku workflow'u

Utwórz kopię pliku workflow'u, aby użyć jej jako punktu wyjścia dla wersji paired-end.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Teraz zaktualizuj parametry i obsługę wejść w nowym pliku.

#### 3.1.2. Dodanie profilu testowego paired-end

Udostępniamy drugi plik CSV zawierający identyfikatory próbek i ścieżki do sparowanych plików FASTQ w katalogu `data/`.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Dodaj profil `test_pe` do `nextflow.config`, który wskazuje na ten plik i używa ID raportu paired-end.

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

Profil testowy dla danych paired-end jest gotowy.

#### 3.1.3. Aktualizacja fabryki kanału

Operator `.map()` musi teraz pobrać obie ścieżki do plików FASTQ i zwrócić je jako listę.

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Utwórz kanał wejściowy z zawartości pliku CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Utwórz kanał wejściowy z zawartości pliku CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

Obsługa wejścia jest skonfigurowana dla danych paired-end.

### 3.2. Dostosowanie modułu FASTQC dla danych paired-end

Skopiuj moduł, aby utworzyć wersję paired-end:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Wejście procesu FASTQC nie wymaga zmian — gdy Nextflow otrzymuje listę dwóch plików, umieszcza oba, a `reads` rozwija się do obu nazw plików.
Jedyna potrzebna zmiana jest w bloku wyjścia: ponieważ teraz otrzymujemy dwa raporty FastQC na próbkę, przechodzimy ze wzorców opartych na `simpleName` na symbole wieloznaczne.

=== "Po"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Przed"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

To uogólnia proces w sposób, który umożliwia mu obsługę zarówno danych single-end, jak i paired-end.

Zaktualizuj import w `rnaseq_pe.nf`, aby używał wersji paired-end:

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

Moduł FASTQC i jego import są zaktualizowane dla danych paired-end.

### 3.3. Dostosowanie modułu TRIM_GALORE dla danych paired-end

Skopiuj moduł, aby utworzyć wersję paired-end:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Ten moduł wymaga bardziej znaczących zmian:

- Wejście zmienia się z pojedynczej ścieżki na krotkę dwóch ścieżek
- Polecenie dodaje flagę `--paired` i przyjmuje oba pliki odczytów
- Wyjście zmienia się, aby odzwierciedlić konwencje nazewnictwa Trim Galore dla paired-end, generując osobne raporty FastQC dla każdego pliku odczytu

=== "Po"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
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

=== "Przed"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

Zaktualizuj import w `rnaseq_pe.nf`:

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Moduł TRIM_GALORE i jego import są zaktualizowane dla danych paired-end.

### 3.4. Dostosowanie modułu HISAT2_ALIGN dla danych paired-end

Skopiuj moduł, aby utworzyć wersję paired-end:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Ten moduł wymaga podobnych zmian:

- Wejście zmienia się z pojedynczej ścieżki na krotkę dwóch ścieżek
- Polecenie HISAT2 zmienia się z `-U` (nieparzysty) na argumenty odczytów `-1` i `-2` (parzyste)
- Wszystkie użycia `reads.simpleName` zmieniają się na `read1.simpleName`, ponieważ teraz odwołujemy się do konkretnego członka pary

=== "Po"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Przed"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

Zaktualizuj import w `rnaseq_pe.nf`:

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Moduł HISAT2_ALIGN i jego import są zaktualizowane dla danych paired-end.

### 3.5. Aktualizacja agregacji MultiQC dla wyjść paired-end

Proces `TRIM_GALORE` paired-end generuje teraz dwa osobne kanały raportów FastQC (`fastqc_reports_1` i `fastqc_reports_2`) zamiast jednego.
Zaktualizuj blok `.mix()` w `rnaseq_pe.nf`, aby uwzględnić oba:

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

Agregacja MultiQC teraz uwzględnia oba zestawy raportów FastQC paired-end.

### 3.6. Aktualizacja obsługi wyjść dla wyjść paired-end

Sekcja `publish:` i blok `output {}` również muszą odzwierciedlać dwa osobne kanały raportów FastQC z procesu `TRIM_GALORE` paired-end.

Zaktualizuj sekcję `publish:` w `rnaseq_pe.nf`:

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Zaktualizuj odpowiednie wpisy w bloku `output {}`:

=== "Po"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Przed"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

Workflow paired-end jest teraz w pełni zaktualizowany i gotowy do uruchomienia.

### 3.7. Uruchomienie workflow'u

Nie używamy `-resume`, ponieważ to nie zostałoby zbuforowane, a jest dwa razy więcej danych do przetworzenia niż wcześniej, ale i tak powinno zakończyć się w mniej niż minutę.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
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

Teraz mamy dwie nieco rozbieżne wersje naszego workflow'u, jedną dla danych single-end i jedną dla danych paired-end.
Następnym logicznym krokiem byłoby sprawienie, aby workflow akceptował którykolwiek typ danych w locie, co wykracza poza zakres tego kursu, ale możemy się tym zająć w kontynuacji.

---

### Podsumowanie

Wiesz, jak dostosować workflow dla pojedynczej próbki, aby sparalelizować przetwarzanie wielu próbek, wygenerować kompleksowy raport QC i dostosować workflow do używania danych paired-end.

### Co dalej?

Pogratuluj sobie! Ukończyłeś kurs Nextflow dla RNAseq.

Przejdź do końcowego [podsumowania kursu](./next_steps.md), aby przejrzeć to, czego się nauczyłeś i dowiedzieć się, co dalej.
