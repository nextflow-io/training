# Część 2: Implementacja dla pojedynczej próbki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej części kursu napiszemy workflow, który obejmie wszystkie polecenia uruchomione w Części 1, aby zautomatyzować ich wykonywanie. Będziemy dążyć do przetwarzania jednej próbki na raz.

!!! warning "Wymaganie wstępne"

    Musisz przejść przez [Część 1: Przegląd metody](./01_method.md) przed rozpoczęciem tej lekcji.
    W szczególności, praca przez sekcję 1.2.3 tworzy plik indeksu genomu (`data/genome_index.tar.gz`) wymagany do kroku dopasowania w tej lekcji.

## Zadanie

W tej części kursu opracujemy workflow, który wykonuje następujące czynności:

1. Uruchamia kontrolę jakości (FastQC) na odczytach wejściowych
2. Przycina adaptery i uruchamia kontrolę jakości po przycięciu (Trim Galore)
3. Dopasowuje przycięte odczyty do genomu referencyjnego (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

To automatyzuje kroki z pierwszej sekcji [Części 1: Przegląd metody](./01_method.md#1-single-sample-processing), gdzie uruchamiałeś te polecenia ręcznie w ich kontenerach.

Jako punkt wyjścia udostępniamy plik workflow'u, `rnaseq.nf`, który określa główne części workflow'u, oraz cztery pliki modułów w katalogu `modules/` (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf` i `multiqc.nf`), które określają strukturę każdego procesu.

??? full-code "Pliki szkieletowe"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Instrukcje INCLUDE modułów

    /*
     * Parametry pipeline'u
     */

    // Główne wejście

    workflow {

        main:
        // Utwórz kanał wejściowy

        // Wywołaj procesy

        publish:
        // Zadeklaruj wyjścia do publikacji
    }

    output {
        // Skonfiguruj cele publikacji
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Uruchom FastQC na odczytach wejściowych
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Przytnij adaptery i uruchom kontrolę jakości po przycięciu
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Dopasuj odczyty do genomu referencyjnego
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Agreguj raporty kontroli jakości za pomocą MultiQC
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

Te pliki nie są funkcjonalne; ich celem jest jedynie służenie jako szkielety, które wypełnisz interesującymi częściami kodu.

## Plan lekcji

Aby proces rozwoju był bardziej edukacyjny, podzieliliśmy to na trzy etapy:

1. **Napisanie jednoetapowego workflow'u, który uruchamia początkowy krok kontroli jakości.**
   Obejmuje to konfigurację parametru CLI, utworzenie kanału wejściowego, napisanie modułu procesu i skonfigurowanie publikacji wyjścia.
2. **Dodanie przycinania adapterów i kontroli jakości po przycięciu.**
   Wprowadza to łączenie procesów poprzez połączenie wyjścia jednego procesu z wejściem drugiego.
3. **Dodanie dopasowania do genomu referencyjnego.**
   Obejmuje to obsługę dodatkowych wejść referencyjnych i pracę ze skompresowanymi archiwami.

Każdy krok koncentruje się na konkretnym aspekcie rozwoju workflow'u.

!!! tip "Wskazówka"

     Upewnij się, że jesteś we właściwym katalogu roboczym:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Napisanie jednoetapowego workflow'u, który uruchamia początkową kontrolę jakości

Ten pierwszy krok koncentruje się na podstawach: wczytaniu pliku FASTQ i uruchomieniu na nim kontroli jakości.

Przypomnij sobie polecenie `fastqc` z [Części 1](01_method.md):

```bash
fastqc <reads>
```

Polecenie przyjmuje plik FASTQ jako wejście i tworzy raport kontroli jakości jako archiwum `.zip` oraz podsumowanie `.html`.
URI kontenera to `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Opakujemy te informacje w Nextflow w trzech etapach:

1. Konfiguracja wejścia
2. Napisanie procesu kontroli jakości i wywołanie go w workflow
3. Skonfigurowanie obsługi wyjścia

### 1.1. Konfiguracja wejścia

Musimy zadeklarować parametr wejściowy, utworzyć profil testowy, aby zapewnić wygodną wartość domyślną, i utworzyć kanał wejściowy.

#### 1.1.1. Dodanie deklaracji parametru wejściowego

W `rnaseq.nf`, w sekcji `Pipeline parameters`, zadeklaruj parametr o nazwie `reads` z typem `Path`.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Parametry pipeline'u
     */
    params {
        // Główne wejście
        input: Path
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Parametry pipeline'u
     */

    // Główne wejście
    ```

To konfiguruje parametr CLI, ale nie chcemy wpisywać ścieżki pliku za każdym razem, gdy uruchamiamy workflow podczas rozwoju.
Istnieje wiele opcji dostarczenia wartości domyślnej; tutaj używamy profilu testowego.

#### 1.1.2. Utworzenie profilu testowego z wartością domyślną w `nextflow.config`

Profil testowy zapewnia wygodne wartości domyślne do wypróbowania workflow'u bez określania wejść w linii poleceń.
Jest to powszechna konwencja w ekosystemie Nextflow (zobacz [Hello Config](../../hello_nextflow/06_hello_config.md) po więcej szczegółów).

Dodaj blok `profiles` do `nextflow.config` z profilem `test`, który ustawia parametr `reads` na jeden z testowych plików FASTQ.

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Tutaj używamy `#!groovy ${projectDir}`, wbudowanej zmiennej Nextflow, która wskazuje na katalog, w którym znajduje się skrypt workflow'u.
Ułatwia to odwoływanie się do plików danych i innych zasobów bez kodowania bezwzględnych ścieżek.

Parametr ma teraz wygodną wartość domyślną. Następnie musimy utworzyć z niego kanał.

#### 1.1.3. Konfiguracja kanału wejściowego

W bloku workflow utwórz kanał wejściowy z wartości parametru używając fabryki kanałów `.fromPath` (jak w [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Po"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Wywołaj procesy

        publish:
        // Zadeklaruj wyjścia do publikacji
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Utwórz kanał wejściowy

        // Wywołaj procesy

        publish:
        // Zadeklaruj wyjścia do publikacji
    }
    ```

Następnie musimy utworzyć proces do uruchomienia kontroli jakości na tym wejściu.

### 1.2. Napisanie procesu kontroli jakości i wywołanie go w workflow

Musimy wypełnić definicję procesu w pliku modułu, zaimportować go do workflow'u za pomocą instrukcji include i wywołać go na wejściu.

#### 1.2.1. Wypełnienie modułu dla procesu kontroli jakości

Otwórz `modules/fastqc.nf` i przeanalizuj zarys definicji procesu.
Powinieneś rozpoznać główne elementy strukturalne; jeśli nie, rozważ przeczytanie [Hello Nextflow](../../hello_nextflow/01_hello_world.md) dla przypomnienia.

Wypełnij definicję procesu samodzielnie, używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Uruchom FastQC na odczytach wejściowych
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Po"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Uruchom FastQC na odczytach wejściowych
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

Akcesor `simpleName` usuwa wszystkie rozszerzenia z nazwy pliku, więc `ENCSR000COQ1_1.fastq.gz` staje się `ENCSR000COQ1_1`.
Używamy składni `emit:`, aby przypisać nazwy do każdego kanału wyjściowego, co będzie przydatne do podłączenia wyjść do bloku publish.

Po zakończeniu tego proces jest kompletny.
Aby użyć go w workflow, musisz zaimportować moduł i dodać wywołanie procesu.

#### 1.2.2. Dołączenie modułu

W `rnaseq.nf` dodaj instrukcję `include`, aby udostępnić proces workflow'owi:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="3"
    // Instrukcje INCLUDE modułów
    ```

Proces jest teraz dostępny w zakresie workflow'u.

#### 1.2.3. Wywołanie procesu kontroli jakości na wejściu

Dodaj wywołanie `FASTQC` w bloku workflow, przekazując kanał wejściowy jako argument.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Początkowa kontrola jakości
        FASTQC(read_ch)

        publish:
        // Zadeklaruj wyjścia do publikacji
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Wywołaj procesy

        publish:
        // Zadeklaruj wyjścia do publikacji
    }
    ```

Workflow teraz wczytuje wejście i uruchamia na nim proces kontroli jakości.
Następnie musimy skonfigurować sposób publikacji wyjścia.

### 1.3. Skonfigurowanie obsługi wyjścia

Musimy zadeklarować, które wyjścia procesu publikować i określić, gdzie powinny trafić.

#### 1.3.1. Zadeklarowanie wyjść w sekcji `publish:`

Sekcja `publish:` wewnątrz bloku workflow deklaruje, które wyjścia procesu powinny być publikowane.
Przypisz wyjścia `FASTQC` do nazwanych celów.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Zadeklaruj wyjścia do publikacji
    }
    ```

Następnie musimy powiedzieć Nextflow, gdzie umieścić opublikowane wyjścia.

#### 1.3.2. Skonfigurowanie celów wyjściowych w bloku `output {}`

Blok `output {}` znajduje się poza workflow i określa, gdzie każdy nazwany cel jest publikowany.
Skonfiguruj oba cele do publikacji w podkatalogu `fastqc/`.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Skonfiguruj cele publikacji
    }
    ```

!!! note "Uwaga"

    Domyślnie Nextflow publikuje pliki wyjściowe jako dowiązania symboliczne, co unika niepotrzebnego duplikowania.
    Mimo że pliki danych, których tutaj używamy, są bardzo małe, w genomice mogą być bardzo duże.
    Dowiązania symboliczne przestaną działać po wyczyszczeniu katalogu `work`, więc w finalnych workflow'ach możesz chcieć nadpisać domyślny tryb publikacji na `'copy'`.

### 1.4. Uruchomienie workflow'u

W tym momencie mamy jednoetapowy workflow kontroli jakości, który powinien być w pełni funkcjonalny.

Uruchamiamy z `-profile test`, aby użyć wartości domyślnej ustawionej w profilu testowym, unikając konieczności wpisywania ścieżki w linii poleceń.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

Powinno to działać bardzo szybko, jeśli pracowałeś przez Część 1 i już pobrałeś kontener.
Jeśli ją pominąłeś, Nextflow pobierze kontener za Ciebie; nie musisz nic robić, aby to się stało, ale możesz potrzebować poczekać do minuty.

Możesz sprawdzić wyjścia w katalogu results.

```bash
ls results/fastqc
```

```console title="Wyjście"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

Raporty kontroli jakości dla próbki są teraz opublikowane w podkatalogu `fastqc/`.

### Podsumowanie

Wiesz, jak utworzyć moduł zawierający proces, zaimportować go do workflow'u, wywołać go z kanałem wejściowym i opublikować wyniki używając bloku output na poziomie workflow'u.

### Co dalej?

Dodaj przycinanie adapterów z kontrolą jakości po przycięciu jako drugi krok w workflow.

---

## 2. Dodanie przycinania adapterów i kontroli jakości po przycięciu

Teraz, gdy mamy początkową kontrolę jakości, możemy dodać krok przycinania adapterów z wbudowaną kontrolą jakości po przycięciu.

Przypomnij sobie polecenie `trim_galore` z [Części 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

Polecenie przycina adaptery z pliku FASTQ i uruchamia FastQC na przyciętym wyjściu.
Tworzy przycięte odczyty, raport przycinania i raporty FastQC dla przyciętych odczytów.
URI kontenera to `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Musimy tylko napisać definicję procesu, zaimportować go, wywołać w workflow i zaktualizować obsługę wyjścia.

### 2.1. Napisanie procesu przycinania i wywołanie go w workflow

Jak poprzednio, musimy wypełnić definicję procesu, zaimportować moduł i dodać wywołanie procesu.

#### 2.1.1. Wypełnienie modułu dla procesu przycinania

Otwórz `modules/trim_galore.nf` i przeanalizuj zarys definicji procesu.

Wypełnij definicję procesu samodzielnie, używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Przytnij adaptery i uruchom kontrolę jakości po przycięciu
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Po"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Przytnij adaptery i uruchom kontrolę jakości po przycięciu
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

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
    }
    ```

Ten proces ma trzy nazwane wyjścia: przycięte odczyty, które trafiają do kroku dopasowania, raport przycinania i raporty FastQC po przycięciu.
Flaga `--fastqc` mówi Trim Galore, aby automatycznie uruchomił FastQC na przyciętym wyjściu.

#### 2.1.2. Dołączenie modułu

Zaktualizuj `rnaseq.nf`, aby zaimportować nowy moduł:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="3"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    ```

Następnie dodamy wywołanie procesu do workflow'u.

#### 2.1.3. Wywołanie procesu przycinania na wejściu

Dodaj wywołanie procesu w bloku workflow:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Początkowa kontrola jakości
        FASTQC(read_ch)

        // Przycinanie adapterów i kontrola jakości po przycięciu
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Początkowa kontrola jakości
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Proces przycinania jest teraz podłączony do workflow'u.

### 2.2. Aktualizacja obsługi wyjścia

Musimy dodać wyjścia przycinania do deklaracji publish i skonfigurować, gdzie mają trafić.

#### 2.2.1. Dodanie celów publikacji dla wyjść przycinania

Dodaj wyjścia przycinania do sekcji `publish:`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Następnie musimy powiedzieć Nextflow, gdzie umieścić te wyjścia.

#### 2.2.2. Skonfigurowanie nowych celów wyjściowych

Dodaj wpisy dla celów przycinania w bloku `output {}`, publikując je w podkatalogu `trimming/`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

Konfiguracja wyjścia jest kompletna.

### 2.3. Uruchomienie workflow'u

Workflow teraz obejmuje zarówno początkową kontrolę jakości, jak i przycinanie adapterów.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

To również powinno działać bardzo szybko, ponieważ uruchamiamy na tak małym pliku wejściowym.

Możesz znaleźć wyjścia przycinania w katalogu results.

```bash
ls results/trimming
```

```console title="Wyjście"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

Wyjścia przycinania i raporty kontroli jakości po przycięciu są teraz w podkatalogu `trimming/`.

### Podsumowanie

Wiesz, jak dodać drugi krok przetwarzania, który działa niezależnie na tym samym wejściu, tworząc wiele nazwanych wyjść.

### Co dalej?

Dodaj krok dopasowania, który łączy się z wyjściem przyciętych odczytów.

---

## 3. Dodanie dopasowania do genomu referencyjnego

Na koniec możemy dodać krok dopasowania genomu przy użyciu HISAT2.

Przypomnij sobie polecenie dopasowania z [Części 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

Polecenie dopasowuje odczyty do genomu referencyjnego i konwertuje wyjście do formatu BAM.
Wymaga wcześniej zbudowanego archiwum indeksu genomu i tworzy plik BAM oraz log podsumowania dopasowania.
URI kontenera to `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

Ten proces wymaga dodatkowego wejścia (archiwum indeksu genomu), więc musimy najpierw to skonfigurować, a następnie napisać i podłączyć proces.

### 3.1. Konfiguracja wejść

Musimy zadeklarować parametr dla archiwum indeksu genomu.

#### 3.1.1. Dodanie parametru dla indeksu genomu

Dodaj deklarację parametru dla archiwum indeksu genomu w `rnaseq.nf`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Główne wejście
        input: Path

        // Archiwum genomu referencyjnego
        hisat2_index_zip: Path
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Główne wejście
        input: Path
    }
    ```

#### 3.1.2. Dodanie wartości domyślnej indeksu genomu do profilu testowego

Tak jak zrobiliśmy dla `reads` w sekcji 1.1.2, dodaj wartość domyślną dla indeksu genomu do profilu testowego w `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
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
        }
    }
    ```

Parametr jest gotowy; teraz możemy utworzyć proces dopasowania.

### 3.2. Napisanie procesu dopasowania i wywołanie go w workflow

Jak poprzednio, musimy wypełnić definicję procesu, zaimportować moduł i dodać wywołanie procesu.

#### 3.2.1. Wypełnienie modułu dla procesu dopasowania

Otwórz `modules/hisat2_align.nf` i przeanalizuj zarys definicji procesu.

Wypełnij definicję procesu samodzielnie, używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Dopasuj odczyty do genomu referencyjnego
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Po"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Dopasuj odczyty do genomu referencyjnego
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

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
    }
    ```

Ten proces przyjmuje dwa wejścia: odczyty i archiwum indeksu genomu.
Blok script najpierw wypakowuje indeks z archiwum, a następnie uruchamia dopasowanie HISAT2 przekierowane do `samtools view`, aby przekonwertować wyjście do formatu BAM.
Akcesor `simpleName` na `index_zip` wyodrębnia nazwę bazową archiwum (`genome_index`) do użycia jako prefiks indeksu.

#### 3.2.2. Dołączenie modułu

Zaktualizuj `rnaseq.nf`, aby zaimportować nowy moduł:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="3"
    // Instrukcje INCLUDE modułów
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Następnie dodamy wywołanie procesu do workflow'u.

#### 3.2.3. Wywołanie procesu dopasowania

Przycięte odczyty znajdują się w kanale wyjściowym `TRIM_GALORE.out.trimmed_reads` z poprzedniego kroku.
Używamy `#!groovy file(params.hisat2_index_zip)`, aby dostarczyć archiwum indeksu genomu.

=== "Po"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Początkowa kontrola jakości
        FASTQC(read_ch)

        // Przycinanie adapterów i kontrola jakości po przycięciu
        TRIM_GALORE(read_ch)

        // Dopasowanie do genomu referencyjnego
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Utwórz kanał wejściowy ze ścieżki pliku
        read_ch = channel.fromPath(params.input)

        // Początkowa kontrola jakości
        FASTQC(read_ch)

        // Przycinanie adapterów i kontrola jakości po przycięciu
        TRIM_GALORE(read_ch)
    ```

Proces dopasowania jest teraz podłączony do workflow'u.

### 3.3. Aktualizacja obsługi wyjścia

Musimy dodać wyjścia dopasowania do deklaracji publish i skonfigurować, gdzie mają trafić.

#### 3.3.1. Dodanie celów publikacji dla wyjść dopasowania

Dodaj wyjścia dopasowania do sekcji `publish:`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
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

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

Następnie musimy powiedzieć Nextflow, gdzie umieścić te wyjścia.

#### 3.3.2. Skonfigurowanie nowych celów wyjściowych

Dodaj wpisy dla celów dopasowania w bloku `output {}`, publikując je w podkatalogu `align/`:

=== "Po"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Przed"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

Konfiguracja wyjścia jest kompletna.

### 3.4. Uruchomienie workflow'u

Workflow teraz obejmuje wszystkie trzy kroki przetwarzania: kontrolę jakości, przycinanie i dopasowanie.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

Możesz znaleźć wyjścia dopasowania w katalogu results.

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

Zrób sobie przerwę! To było dużo.

Kiedy poczujesz się odświeżony, przejdź do [Części 3](./03_multi-sample.md), gdzie dowiesz się, jak zmodyfikować workflow, aby przetwarzać wiele próbek równolegle, agregować raporty kontroli jakości we wszystkich krokach dla wszystkich próbek i umożliwić uruchamianie workflow'u na danych RNAseq paired-end.
