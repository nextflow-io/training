# Część 2: Wywoływanie wariantów dla poszczególnych próbek

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W Części 1 przetestowałeś polecenia Samtools i GATK ręcznie w ich odpowiednich kontenerach.
Teraz opakujemy te same polecenia w workflow Nextflow'a.

## Zadanie

W tej części kursu stworzymy workflow'a, który wykonuje następujące czynności:

1. Generuje plik indeksu dla każdego wejściowego pliku BAM przy użyciu [Samtools](https://www.htslib.org/)
2. Uruchamia GATK HaplotypeCaller na każdym wejściowym pliku BAM, aby wygenerować wywołania wariantów dla poszczególnych próbek w formacie VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

To odtwarza kroki z Części 1, gdzie uruchamiałeś te polecenia ręcznie w kontenerach.

Jako punkt wyjścia dostarczamy Ci plik workflow'a, `genomics.nf`, który przedstawia główne części workflow'a, oraz dwa pliki modułów, samtools_index.nf i gatk_haplotypecaller.nf, które przedstawiają strukturę modułów.
Te pliki nie są funkcjonalne; ich celem jest jedynie służenie jako szkielety, które wypełnisz interesującymi fragmentami kodu.

## Plan lekcji

Aby proces tworzenia był bardziej edukacyjny, podzieliliśmy to na cztery kroki:

1. **Napisz jednokrokowy workflow'a, który uruchamia Samtools index na pliku BAM.**
   Obejmuje to tworzenie modułu, importowanie go i wywoływanie w workflow.
2. **Dodaj drugi proces, aby uruchomić GATK HaplotypeCaller na zaindeksowanym pliku BAM.**
   Wprowadza to łączenie wyjść procesów z wejściami oraz obsługę plików pomocniczych.
3. **Dostosuj workflow'a do uruchamiania na zestawie próbek.**
   Obejmuje to równoległe wykonywanie i wprowadza krotki, aby utrzymać powiązane pliki razem.
4. **Spraw, aby workflow akceptował plik tekstowy zawierający zestaw plików wejściowych.**
   Demonstruje to powszechny wzorzec dostarczania wejść zbiorczo.

Każdy krok koncentruje się na konkretnym aspekcie tworzenia workflow'a.

---

## 1. Napisz jednokrokowy workflow'a, który uruchamia Samtools index na pliku BAM

Ten pierwszy krok koncentruje się na podstawach: wczytywaniu pliku BAM i generowaniu dla niego indeksu.

Przypomnij sobie polecenie `samtools index` z [Części 1](01_method.md):

```bash
samtools index '<input_bam>'
```

Polecenie przyjmuje plik BAM jako wejście i tworzy plik indeksu `.bai` obok niego.
URI kontenera to `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Opakujemy te informacje w Nextflow'a w trzech etapach:

1. Skonfiguruj wejście
2. Napisz proces indeksowania i wywołaj go w workflow
3. Skonfiguruj obsługę wyjścia

### 1.1. Skonfiguruj wejście

Musimy zadeklarować parametr wejściowy, utworzyć profil testowy, aby zapewnić wygodną wartość domyślną, oraz utworzyć kanał wejściowy.

#### 1.1.1. Dodaj deklarację parametru wejściowego

W głównym pliku workflow'a `genomics.nf`, w sekcji `Pipeline parameters`, zadeklaruj parametr CLI o nazwie `reads_bam`.

=== "Po"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

To konfiguruje parametr CLI, ale nie chcemy wpisywać ścieżki pliku za każdym razem, gdy uruchamiamy workflow'a podczas tworzenia.
Istnieje wiele opcji dostarczania wartości domyślnej; tutaj używamy profilu testowego.

#### 1.1.2. Utwórz profil testowy z wartością domyślną w `nextflow.config`

Profil testowy zapewnia wygodne wartości domyślne do wypróbowania workflow'a bez określania wejść w wierszu poleceń.
Jest to powszechna konwencja w ekosystemie Nextflow'a (zobacz [Hello Config](../../hello_nextflow/06_hello_config.md) po więcej szczegółów).

Dodaj blok `profiles` do `nextflow.config` z profilem `test`, który ustawia parametr `reads_bam` na jeden z testowych plików BAM.

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Tutaj używamy `${projectDir}`, wbudowanej zmiennej Nextflow'a, która wskazuje na katalog, w którym znajduje się skrypt workflow'a.
Ułatwia to odwoływanie się do plików danych i innych zasobów bez kodowania bezwzględnych ścieżek na stałe.

#### 1.1.3. Skonfiguruj kanał wejściowy

W bloku workflow utwórz kanał wejściowy z wartości parametru, używając fabryki kanałów `.fromPath` (jak w [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Po"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Teraz musimy utworzyć proces do uruchomienia indeksowania na tym wejściu.

### 1.2. Napisz proces indeksowania i wywołaj go w workflow

Musimy napisać definicję procesu w pliku modułu, zaimportować go do workflow'a za pomocą instrukcji include oraz wywołać go na wejściu.

#### 1.2.1. Wypełnij moduł dla procesu indeksowania

Otwórz `modules/samtools_index.nf` i przeanalizuj zarys definicji procesu.
Powinieneś rozpoznać główne elementy strukturalne; jeśli nie, rozważ przeczytanie [Hello Nextflow](../../hello_nextflow/01_hello_world.md) dla przypomnienia.

Wypełnij definicję procesu samodzielnie, używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Po"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Wygeneruj plik indeksu BAM
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

Po ukończeniu tego proces jest gotowy.
Aby użyć go w workflow, musisz zaimportować moduł i dodać wywołanie procesu.

#### 1.2.2. Dołącz moduł

W `genomics.nf` dodaj instrukcję `include`, aby udostępnić proces workflow'owi:

=== "Po"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

Proces jest teraz dostępny w zakresie workflow'a.

#### 1.2.3. Wywołaj proces indeksowania na wejściu

Teraz dodajmy wywołanie `SAMTOOLS_INDEX` w bloku workflow, przekazując kanał wejściowy jako argument.

=== "Po"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Utwórz plik indeksu dla wejściowego pliku BAM
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

Workflow teraz wczytuje wejście i uruchamia na nim proces indeksowania.
Następnie musimy skonfigurować sposób publikowania wyjścia.

### 1.3. Skonfiguruj obsługę wyjścia

Musimy zadeklarować, które wyjścia procesów publikować i określić, gdzie powinny trafić.

#### 1.3.1. Zadeklaruj wyjście w sekcji `publish:`

Sekcja `publish:` wewnątrz bloku workflow deklaruje, które wyjścia procesów powinny być publikowane.
Przypisz wyjście `SAMTOOLS_INDEX` do nazwanego celu o nazwie `bam_index`.

=== "Po"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Teraz musimy powiedzieć Nextflow'owi, gdzie umieścić opublikowane wyjście.

#### 1.3.2. Skonfiguruj cel wyjściowy w bloku `output {}`

Blok `output {}` znajduje się poza workflow'em i określa, gdzie publikowany jest każdy nazwany cel.
Dodajmy cel dla `bam_index`, który publikuje do podkatalogu `bam/`.

=== "Po"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note "Uwaga"

    Domyślnie Nextflow publikuje pliki wyjściowe jako dowiązania symboliczne, co pozwala uniknąć niepotrzebnego duplikowania.
    Mimo że pliki danych, których tutaj używamy, są bardzo małe, w genomice mogą być bardzo duże.
    Dowiązania symboliczne przestaną działać po wyczyszczeniu katalogu `work`, więc w przypadku workflow'ów produkcyjnych możesz chcieć nadpisać domyślny tryb publikowania na `'copy'`.

### 1.4. Uruchom workflow'a

W tym momencie mamy jednokrokowy workflow indeksowania, który powinien być w pełni funkcjonalny. Sprawdźmy, czy działa!

Możemy go uruchomić z `-profile test`, aby użyć wartości domyślnej ustawionej w profilu testowym i uniknąć konieczności wpisywania ścieżki w wierszu poleceń.

```bash
nextflow run genomics.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Możesz sprawdzić, czy plik indeksu został wygenerowany poprawnie, zaglądając do katalogu roboczego lub do katalogu wyników.

??? abstract "Zawartość katalogu roboczego"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Zawartość katalogu wyników"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

Oto on!

### Podsumowanie

Wiesz, jak utworzyć moduł zawierający proces, zaimportować go do workflow'a, wywołać go z kanałem wejściowym i opublikować wyniki.

### Co dalej?

Dodaj drugi krok, który pobiera wyjście procesu indeksowania i używa go do uruchomienia wywoływania wariantów.

---

## 2. Dodaj drugi proces, aby uruchomić GATK HaplotypeCaller na zaindeksowanym pliku BAM

Teraz, gdy mamy indeks dla naszego pliku wejściowego, możemy przejść do konfiguracji kroku wywoływania wariantów.

Przypomnij sobie polecenie `gatk HaplotypeCaller` z [Części 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

Polecenie przyjmuje plik BAM (`-I`), genom referencyjny (`-R`) i plik interwałów (`-L`), a następnie tworzy plik VCF (`-O`) wraz z jego indeksem.
Narzędzie oczekuje również, że indeks BAM, indeks referencji i słownik referencji będą zlokalizowane razem z odpowiednimi plikami.
URI kontenera to `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Postępujemy według tych samych trzech etapów co wcześniej:

1. Skonfiguruj wejścia
2. Napisz proces wywoływania wariantów i wywołaj go w workflow
3. Skonfiguruj obsługę wyjścia

### 2.1. Skonfiguruj wejścia

Krok wywoływania wariantów wymaga kilku dodatkowych plików wejściowych.
Musimy zadeklarować dla nich parametry, dodać wartości domyślne do profilu testowego i utworzyć zmienne do ich wczytania.

#### 2.1.1. Dodaj deklaracje parametrów dla plików pomocniczych

Ponieważ nasz nowy proces oczekuje kilku dodatkowych plików, dodaj deklaracje parametrów dla nich w `genomics.nf` w sekcji `Pipeline parameters`:

=== "Po"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Pliki pomocnicze
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

Jak poprzednio, dostarczamy wartości domyślne przez profil testowy, a nie bezpośrednio.

#### 2.1.2. Dodaj wartości domyślne plików pomocniczych do profilu testowego

Tak jak zrobiliśmy dla `reads_bam` w sekcji 1.1.2, dodaj wartości domyślne dla plików pomocniczych do profilu testowego w `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Teraz musimy utworzyć zmienne, które wczytują te ścieżki plików do użycia w workflow.

#### 2.1.3. Utwórz zmienne dla plików pomocniczych

Dodaj zmienne dla ścieżek plików pomocniczych wewnątrz bloku workflow:

=== "Po"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Wczytaj ścieżki plików dla plików pomocniczych (referencja i interwały)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Utwórz plik indeksu dla wejściowego pliku BAM
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Utwórz plik indeksu dla wejściowego pliku BAM
        SAMTOOLS_INDEX(reads_ch)
    ```

Składnia `file()` mówi Nextflow'owi jawnie, aby obsługiwał te wejścia jako ścieżki plików.
Możesz dowiedzieć się więcej na ten temat w Side Quest [Working with files](../../side_quests/working_with_files.md).

### 2.2. Napisz proces wywoływania wariantów i wywołaj go w workflow

Musimy napisać definicję procesu w pliku modułu, zaimportować go do workflow'a za pomocą instrukcji include oraz wywołać go na wejściowych odczytach plus wyjściu kroku indeksowania i plikach pomocniczych.

#### 2.2.1. Wypełnij moduł dla procesu wywoływania wariantów

Otwórz `modules/gatk_haplotypecaller.nf` i przeanalizuj zarys definicji procesu.

Wypełnij definicję procesu samodzielnie, używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Po"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Wywołaj warianty za pomocą GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

Zauważysz, że ten proces ma więcej wejść niż wymaga samo polecenie GATK.
GATK wie, gdzie szukać pliku indeksu BAM i plików pomocniczych genomu referencyjnego na podstawie konwencji nazewnictwa, ale Nextflow jest niezależny od dziedziny i nie zna tych konwencji.
Musimy wymienić je jawnie, aby Nextflow umieścił je w katalogu roboczym w czasie wykonywania; w przeciwnym razie GATK zgłosi błąd o brakujących plikach.

Podobnie wymieniamy plik indeksu wyjściowego VCF (`"${input_bam}.vcf.idx"`) jawnie, aby Nextflow śledził go dla kolejnych kroków.
Używamy składni `emit:`, aby przypisać nazwę do każdego kanału wyjściowego, co stanie się przydatne, gdy podłączymy wyjścia do bloku publish.

Po ukończeniu tego proces jest gotowy.
Aby użyć go w workflow, musisz zaimportować moduł i dodać wywołanie procesu.

#### 2.2.2. Zaimportuj nowy moduł

Zaktualizuj `genomics.nf`, aby zaimportować nowy moduł:

=== "Po"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

Proces jest teraz dostępny w zakresie workflow'a.

#### 2.2.3. Dodaj wywołanie procesu

Dodaj wywołanie procesu w treści workflow, pod `main:`:

=== "Po"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Utwórz plik indeksu dla wejściowego pliku BAM
        SAMTOOLS_INDEX(reads_ch)

        // Wywołaj warianty z zaindeksowanego pliku BAM
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="33"
        // Utwórz plik indeksu dla wejściowego pliku BAM
        SAMTOOLS_INDEX(reads_ch)
    ```

Powinieneś rozpoznać składnię `*.out` ze szkolenia Hello Nextflow; mówimy Nextflow'owi, aby wziął kanał wyjściowy z `SAMTOOLS_INDEX` i podłączył go do wywołania procesu `GATK_HAPLOTYPECALLER`.

!!! note "Uwaga"

    Zauważ, że wejścia są podawane w dokładnie tej samej kolejności w wywołaniu procesu, jak są wymienione w bloku wejściowym procesu.
    W Nextflow'ie wejścia są pozycyjne, co oznacza, że _musisz_ zachować tę samą kolejność; i oczywiście musi być taka sama liczba elementów.

### 2.3. Skonfiguruj obsługę wyjścia

Musimy dodać nowe wyjścia do deklaracji publish i skonfigurować, gdzie mają trafić.

#### 2.3.1. Dodaj cele publikacji dla wyjść wywoływania wariantów

Dodaj wyjścia VCF i indeksu do sekcji `publish:`:

=== "Po"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Teraz musimy powiedzieć Nextflow'owi, gdzie umieścić nowe wyjścia.

#### 2.3.2. Skonfiguruj nowe cele wyjściowe

Dodaj wpisy dla celów `vcf` i `vcf_idx` w bloku `output {}`, publikując oba do podkatalogu `vcf/`:

=== "Po"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

VCF i jego indeks są publikowane jako oddzielne cele, które oba trafiają do podkatalogu `vcf/`.

### 2.4. Uruchom workflow'a

Uruchom rozszerzony workflow, tym razem dodając `-resume`, abyśmy nie musieli ponownie uruchamiać kroku indeksowania.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Teraz, jeśli spojrzymy na wyjście konsoli, widzimy wymienione dwa procesy.

Pierwszy proces został pominięty dzięki buforowaniu, zgodnie z oczekiwaniami, podczas gdy drugi proces został uruchomiony, ponieważ jest zupełnie nowy.

Pliki wyjściowe znajdziesz w katalogu wyników (jako dowiązania symboliczne do katalogu roboczego).

??? abstract "Zawartość katalogu"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Jeśli otworzysz plik VCF, powinieneś zobaczyć tę samą zawartość co w pliku wygenerowanym przez bezpośrednie uruchomienie polecenia GATK w kontenerze.

??? abstract "Zawartość pliku"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

To jest wyjście, które zależy nam na wygenerowaniu dla każdej próbki w naszym badaniu.

### Podsumowanie

Wiesz, jak stworzyć dwukrokowy modularny workflow, który wykonuje rzeczywistą pracę analityczną i jest w stanie radzić sobie z osobliwościami formatów plików genomicznych, takimi jak pliki pomocnicze.

### Co dalej?

Spraw, aby workflow obsługiwał wiele próbek zbiorczo.

---

## 3. Dostosuj workflow'a do uruchamiania na zestawie próbek

To wszystko jest w porządku, aby mieć workflow'a, który może zautomatyzować przetwarzanie pojedynczej próbki, ale co jeśli masz 1000 próbek?
Czy musisz napisać skrypt bash, który przechodzi przez wszystkie Twoje próbki?

Nie, dzięki Bogu! Po prostu wprowadź drobną zmianę w kodzie, a Nextflow również to dla Ciebie obsłuży.

### 3.1. Zaktualizuj wejście, aby wymieniało trzy próbki

Aby uruchomić na wielu próbkach, zaktualizuj profil testowy, aby dostarczał tablicę ścieżek plików zamiast pojedynczej.
To szybki sposób na przetestowanie wykonywania wielopróbkowego; w następnym kroku przejdziemy na bardziej skalowalny sposób przy użyciu pliku wejść.

Najpierw zakomentuj adnotację typu w deklaracji parametru, ponieważ tablice nie mogą używać deklaracji typowanych:

=== "Po"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Wejście główne (tablica trzech próbek)
        reads_bam //: Path
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Następnie zaktualizuj profil testowy, aby wymieniał wszystkie trzy próbki:

=== "Po"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

Fabryka kanałów w treści workflow'a (`.fromPath`) akceptuje wiele ścieżek plików tak samo dobrze jak pojedynczą, więc nie są potrzebne żadne inne zmiany.

### 3.2. Uruchom workflow'a

Spróbuj teraz uruchomić workflow'a, gdy instalacja jest skonfigurowana do uruchamiania na wszystkich trzech próbkach testowych.

```bash
nextflow run genomics.nf -profile test -resume
```

Zabawna rzecz: to _może zadziałać_ LUB _może się nie udać_. Na przykład, oto uruchomienie, które się powiodło:

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Jeśli Twoje uruchomienie workflow'a się powiodło, uruchom je ponownie, aż otrzymasz błąd taki jak ten:

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Jeśli spojrzysz na wyjście błędu polecenia GATK, będzie tam linia taka jak ta:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Cóż, to dziwne, biorąc pod uwagę, że jawnie zaindeksowaliśmy pliki BAM w pierwszym kroku workflow'a. Czy może być coś nie tak z instalacją?

### 3.3. Rozwiąż problem

Sprawdzimy katalogi robocze i użyjemy operatora `view()`, aby dowiedzieć się, co poszło nie tak.

#### 3.3.1. Sprawdź katalogi robocze dla odpowiednich wywołań

Zajrzyj do katalogu roboczego dla nieudanego wywołania procesu `GATK_HAPLOTYPECALLER` wymienionego w wyjściu konsoli.

??? abstract "Zawartość katalogu"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Zwróć szczególną uwagę na nazwy pliku BAM i indeksu BAM, które są wymienione w tym katalogu: `reads_son.bam` i `reads_father.bam.bai`.

Co do diabła? Nextflow umieścił plik indeksu w katalogu roboczym tego wywołania procesu, ale jest to niewłaściwy. Jak to mogło się stać?

#### 3.3.2. Użyj [operatora view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) do sprawdzenia zawartości kanału

Dodaj te dwie linie w treści workflow'a przed wywołaniem procesu `GATK_HAPLOTYPECALLER`, aby wyświetlić zawartość kanału:

=== "Po"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // tymczasowa diagnostyka
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Wywołaj warianty z zaindeksowanego pliku BAM
        GATK_HAPLOTYPECALLER(
    ```

=== "Przed"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Wywołaj warianty z zaindeksowanego pliku BAM
        GATK_HAPLOTYPECALLER(
    ```

Następnie uruchom ponownie polecenie workflow'a.

```bash
nextflow run genomics.nf -profile test
```

Ponownie, to może się udać lub nie. Oto jak wygląda wyjście dwóch wywołań `.view()` dla nieudanego uruchomienia:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Pierwsze trzy linie odpowiadają kanałowi wejściowemu, a drugie - kanałowi wyjściowemu.
Widać, że pliki BAM i pliki indeksów dla trzech próbek nie są wymienione w tej samej kolejności!

!!! note "Uwaga"

    Gdy wywołujesz proces Nextflow'a na kanale zawierającym wiele elementów, Nextflow będzie próbował zrównoleglić wykonywanie tak bardzo, jak to możliwe, i będzie zbierał wyjścia w jakiejkolwiek kolejności staną się dostępne.
    Konsekwencją jest to, że odpowiednie wyjścia mogą być zbierane w innej kolejności niż oryginalne wejścia zostały podane.

Jak obecnie napisany, nasz skrypt workflow'a zakłada, że pliki indeksów wyjdą z kroku indeksowania wymienione w tej samej kolejności matka/ojciec/syn, jak podano wejścia.
Ale nie ma gwarancji, że tak będzie, dlatego czasami (choć nie zawsze) niewłaściwe pliki są sparowane w drugim kroku.

Aby to naprawić, musimy upewnić się, że pliki BAM i ich pliki indeksów podróżują razem przez kanały.

!!! tip "Wskazówka"

    Instrukcje `view()` w kodzie workflow'a nic nie robią, więc nie ma problemu, aby je zostawić.
    Jednak będą zaśmiecać Twoje wyjście konsoli, więc zalecamy ich usunięcie, gdy skończysz rozwiązywać problem.

### 3.4. Zaktualizuj workflow'a, aby poprawnie obsługiwał pliki indeksów

Poprawka polega na spakowaniu każdego pliku BAM z jego indeksem w krotkę, a następnie zaktualizowaniu procesu downstream i instalacji workflow'a, aby pasowały.

#### 3.4.1. Zmień wyjście modułu SAMTOOLS_INDEX na krotkę

Najprostszym sposobem zapewnienia, że plik BAM i jego indeks pozostają ściśle powiązane, jest spakowanie ich razem w krotkę wychodzącą z zadania indeksowania.

!!! note "Uwaga"

    **Krotka** to skończona, uporządkowana lista elementów, która jest powszechnie używana do zwracania wielu wartości z funkcji. Krotki są szczególnie przydatne do przekazywania wielu wejść lub wyjść między procesami przy zachowaniu ich powiązania i kolejności.

Zaktualizuj wyjście w `modules/samtools_index.nf`, aby uwzględnić plik BAM:

=== "Po"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Przed"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

W ten sposób każdy plik indeksu będzie ściśle powiązany z oryginalnym plikiem BAM, a ogólne wyjście kroku indeksowania będzie pojedynczym kanałem zawierającym pary plików.

#### 3.4.2. Zmień wejście modułu GATK_HAPLOTYPECALLER, aby akceptował krotkę

Ponieważ zmieniliśmy "kształt" wyjścia pierwszego procesu, musimy zaktualizować definicję wejścia drugiego procesu, aby pasowała.

Zaktualizuj `modules/gatk_haplotypecaller.nf`:

=== "Po"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Przed"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Teraz musimy zaktualizować workflow, aby odzwierciedlał nową strukturę krotki w wywołaniu procesu i celach publikacji.

#### 3.4.3. Zaktualizuj wywołanie GATK_HAPLOTYPECALLER w workflow

Nie musimy już dostarczać oryginalnego `reads_ch` do procesu `GATK_HAPLOTYPECALLER`, ponieważ plik BAM jest teraz spakowany w kanał wyjściowy przez `SAMTOOLS_INDEX`.

Zaktualizuj wywołanie w `genomics.nf`:

=== "Po"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Na koniec musimy zaktualizować cele publikacji, aby odzwierciedlały nową strukturę wyjścia.

#### 3.4.4. Zaktualizuj cel publikacji dla wyjścia zaindeksowanego BAM

Ponieważ wyjście SAMTOOLS_INDEX jest teraz krotką zawierającą zarówno plik BAM, jak i jego indeks, zmień nazwę celu publikacji z `bam_index` na `indexed_bam`, aby lepiej odzwierciedlić jego zawartość:

=== "Po"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Dzięki tym zmianom BAM i jego indeks mają gwarancję podróżowania razem, więc parowanie zawsze będzie poprawne.

### 3.5. Uruchom poprawiony workflow'a

Uruchom workflow'a ponownie, aby upewnić się, że będzie działał niezawodnie w przyszłości.

```bash
nextflow run genomics.nf -profile test
```

Tym razem (i za każdym razem) wszystko powinno działać poprawnie:

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Katalog wyników zawiera teraz zarówno pliki BAM, jak i BAI dla każdej próbki (z krotki), wraz z wyjściami VCF:

??? abstract "Zawartość katalogu wyników"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

Pakując powiązane pliki w krotki, zapewniliśmy, że właściwe pliki zawsze podróżują razem przez workflow'a.
Workflow teraz przetwarza dowolną liczbę próbek niezawodnie, ale wymienianie ich indywidualnie w konfiguracji nie jest zbyt skalowalne.
W następnym kroku przejdziemy na odczytywanie wejść z pliku.

### Podsumowanie

Wiesz, jak sprawić, aby Twój workflow działał na wielu próbkach (niezależnie).

### Co dalej?

Ułatw obsługę próbek zbiorczo.

---

## 4. Spraw, aby workflow akceptował plik tekstowy zawierający zestaw plików wejściowych

Bardzo powszechnym sposobem dostarczania wielu plików danych wejściowych do workflow'a jest zrobienie tego za pomocą pliku tekstowego zawierającego ścieżki plików.
Może to być tak proste, jak plik tekstowy wymieniający jedną ścieżkę pliku na linię i nic więcej, lub plik może zawierać dodatkowe metadane, w którym to przypadku jest często nazywany arkuszem próbek.

Tutaj pokażemy Ci, jak zrobić prosty przypadek.

### 4.1. Przeanalizuj dostarczony plik tekstowy wymieniający ścieżki plików wejściowych

Już stworzyliśmy plik tekstowy wymieniający ścieżki plików wejściowych, nazwany `sample_bams.txt`, który możesz znaleźć w katalogu `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Jak widać, wymieniliśmy jedną ścieżkę pliku na linię i są to ścieżki bezwzględne.

!!! note "Uwaga"

    Pliki, których tutaj używamy, znajdują się po prostu w lokalnym systemie plików Twojego GitHub Codespaces, ale moglibyśmy również wskazywać na pliki w pamięci chmurowej.
    Jeśli nie używasz dostarczonego środowiska Codespaces, możesz potrzebować dostosować ścieżki plików, aby pasowały do Twojej lokalnej konfiguracji.

### 4.2. Zaktualizuj parametr i profil testowy

Przełącz parametr `reads_bam`, aby wskazywał na plik `sample_bams.txt` zamiast wymieniać poszczególne próbki.

Przywróć adnotację typu w bloku params (ponieważ to znowu pojedyncza ścieżka):

=== "Po"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Wejście główne (plik plików wejściowych, jeden na linię)
        reads_bam: Path
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="10"
        // Wejście główne (tablica trzech próbek)
        reads_bam
    ```

Następnie zaktualizuj profil testowy, aby wskazywał na plik tekstowy:

=== "Po"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

Lista plików nie znajduje się już w kodzie w ogóle, co jest dużym krokiem we właściwym kierunku.

### 4.3. Zaktualizuj fabrykę kanałów, aby odczytywała linie z pliku

Obecnie nasza fabryka kanałów wejściowych traktuje wszystkie pliki, które jej dajemy, jako wejścia danych, które chcemy przekazać do procesu indeksowania.
Ponieważ teraz dajemy jej plik, który wymienia ścieżki plików wejściowych, musimy zmienić jej zachowanie, aby parsowała plik i traktowała ścieżki plików, które zawiera, jako wejścia danych.

Możemy to zrobić, używając tego samego wzorca, którego użyliśmy w [Części 2 Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): stosując operator [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) do parsowania pliku, a następnie operację `map`, aby wybrać pierwsze pole każdej linii.

=== "Po"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Utwórz kanał wejściowy z pliku CSV wymieniającego ścieżki plików wejściowych
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="24"
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Technicznie moglibyśmy to zrobić prościej, używając operatora [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext), ponieważ nasz plik wejściowy obecnie zawiera tylko ścieżki plików.
Jednak używając bardziej wszechstronnego operatora `splitCsv` (uzupełnionego przez `map`), możemy przyszłościować nasz workflow na wypadek, gdybyśmy zdecydowali się dodać metadane do pliku zawierającego ścieżki plików.

!!! tip "Wskazówka"

    Jeśli nie jesteś pewien, że rozumiesz, co robią operatory, to kolejna świetna okazja do użycia operatora `.view()`, aby zobaczyć, jak wygląda zawartość kanału przed i po ich zastosowaniu.

### 4.4. Uruchom workflow'a

Uruchom workflow'a jeszcze raz. Powinno to dać taki sam wynik jak poprzednio, prawda?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Tak! W rzeczywistości Nextflow poprawnie wykrywa, że wywołania procesów są dokładnie takie same i nawet nie zadaje sobie trudu ponownego uruchamiania wszystkiego, ponieważ uruchamialiśmy z `-resume`.

I to wszystko! Nasz prosty workflow wywoływania wariantów ma wszystkie podstawowe funkcje, których chcieliśmy.

### Podsumowanie

Wiesz, jak stworzyć wielokrokowy modularny workflow do indeksowania pliku BAM i zastosowania wywoływania wariantów dla poszczególnych próbek przy użyciu GATK.

Bardziej ogólnie, nauczyłeś się, jak używać podstawowych komponentów i logiki Nextflow'a do budowania prostego pipeline'u genomicznego, który wykonuje rzeczywistą pracę, biorąc pod uwagę osobliwości formatów plików genomicznych i wymagań narzędzi.

### Co dalej?

Świętuj swój sukces i weź dodatkową długą przerwę!

W następnej części tego kursu nauczysz się, jak przekształcić ten prosty workflow wywoływania wariantów dla poszczególnych próbek, aby zastosować wspólne wywoływanie wariantów do danych.
