# Podstawowe Wzorce Skryptowe w Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow to język programowania działający na Java Virtual Machine. Chociaż Nextflow jest zbudowany na [Groovy](http://groovy-lang.org/) i dzieli wiele elementów składni, Nextflow to coś więcej niż tylko "Groovy z rozszerzeniami" -- jest to samodzielny język z w pełni określoną [składnią](https://nextflow.io/docs/latest/reference/syntax.html) i [biblioteką standardową](https://nextflow.io/docs/latest/reference/stdlib.html).

Można pisać dużo kodu w Nextflow, nie wykraczając poza podstawową składnię zmiennych, map i list. Większość tutoriali Nextflow koncentruje się na orkiestracji workflow (kanały, procesy i przepływ danych) i można zajść zaskakująco daleko, używając tylko tego.

Jednak gdy trzeba manipulować danymi, parsować złożone nazwy plików, implementować logikę warunkową lub budować solidne workflow produkcyjne, pomaga myślenie o dwóch odrębnych aspektach kodu: **przepływ danych** (kanały, operatory, procesy i workflow) oraz **skryptowanie** (kod wewnątrz closures, funkcji i skryptów procesów). Choć to rozróżnienie jest nieco arbitralne—to wszystko jest kodem Nextflow—zapewnia użyteczny model mentalny do zrozumienia, kiedy orkiestrujesz pipeline, a kiedy manipulujesz danymi. Opanowanie obu dramatycznie poprawia zdolność pisania przejrzystych, łatwych w utrzymaniu workflow.

### Cele szkolenia

Ten side quest zabiera Cię w praktyczną drogę od podstawowych konceptów do wzorców gotowych do produkcji.
Przekształcimy prosty workflow odczytujący CSV w zaawansowany pipeline bioinformatyczny, rozwijając go krok po kroku przez realistyczne wyzwania:

- **Zrozumienie granic:** Rozróżnienie między operacjami przepływu danych a skryptowaniem i zrozumienie, jak współpracują
- **Manipulacja danymi:** Wyodrębnianie, transformowanie i wybieranie podzbiorów map i kolekcji przy użyciu potężnych operatorów
- **Przetwarzanie ciągów znaków:** Parsowanie złożonych schematów nazewnictwa plików za pomocą wzorców regex i opanowanie interpolacji zmiennych
- **Funkcje wielokrotnego użytku:** Wyodrębnianie złożonej logiki do nazwanych funkcji dla czystszych, łatwiejszych w utrzymaniu workflow
- **Dynamiczna logika:** Budowanie procesów, które dostosowują się do różnych typów wejściowych i używanie closures do dynamicznej alokacji zasobów
- **Warunkowe kierowanie:** Inteligentne kierowanie próbek przez różne procesy na podstawie ich cech metadanych
- **Bezpieczne operacje:** Bezpieczna obsługa brakujących danych z operatorami bezpiecznymi względem null i walidacja wejść z czytelnymi komunikatami błędów
- **Obsługa zdarzeń oparta na konfiguracji:** Użycie handlerów zdarzeń workflow do logowania, powiadomień i zarządzania cyklem życia

### Wymagania wstępne

Przed przystąpieniem do tego side questa powinieneś:

- Ukończyć tutorial [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Czuć się komfortowo z podstawowymi konceptami i mechanizmami Nextflow (procesy, kanały, operatory, praca z plikami, metadane)
- Mieć podstawową znajomość popularnych konstrukcji programistycznych (zmienne, mapy, listy)

Ten tutorial wyjaśni koncepty programistyczne w miarę ich pojawiania się, więc nie potrzebujesz rozległego doświadczenia programistycznego.
Zaczniemy od fundamentalnych konceptów i zbudujemy zaawansowane wzorce.

---

## 0. Rozpoczęcie

#### Otwórz codespace szkoleniowy

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzysz środowisko szkoleniowe zgodnie z opisem w [Konfiguracja Środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki dla tego tutoriala.

```bash
cd side-quests/essential_scripting_patterns
```

#### Przejrzyj materiały

Znajdziesz główny plik workflow oraz katalog `data` zawierający przykładowe pliki danych.

```console title="Zawartość katalogu"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Nasz przykładowy CSV zawiera informacje o próbkach biologicznych wymagających różnego przetwarzania w zależności od ich charakterystyki:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Użyjemy tego realistycznego zbioru danych do eksploracji praktycznych technik programistycznych, z którymi spotkasz się w rzeczywistych workflow bioinformatycznych.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista kontrolna gotowości

Uważasz, że jesteś gotowy do zanurzenia się?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace jest uruchomiony
- [ ] Ustawiłem odpowiednio mój katalog roboczy
<!-- - [ ] I understand the assignment -->

Jeśli możesz zaznaczyć wszystkie pola, możesz zacząć.

---

## 1. Przepływ Danych vs Skryptowanie: Zrozumienie Granic

### 1.1. Identyfikowanie Co Jest Czym

Pisząc workflow w Nextflow, ważne jest rozróżnienie między **przepływem danych** (jak dane przemieszczają się przez kanały i procesy) a **skryptowaniem** (kodem, który manipuluje danymi i podejmuje decyzje). Zbudujmy workflow demonstrujący, jak współpracują.

#### 1.1.1. Podstawowy Workflow Nextflow

Zacznij od prostego workflow, który po prostu odczytuje plik CSV (już to zrobiliśmy dla Ciebie w `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Blok `workflow` definiuje naszą strukturę pipeline, podczas gdy `channel.fromPath()` tworzy kanał ze ścieżki pliku. Operator `.splitCsv()` przetwarza plik CSV i konwertuje każdy wiersz w strukturę danych map.

Uruchom ten workflow, aby zobaczyć surowe dane CSV:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Dodawanie Operatora Map

Teraz dodamy skryptowanie do transformacji danych, używając operatora `.map()`, który prawdopodobnie już znasz. Ten operator przyjmuje 'closure', w którym możemy pisać kod do transformacji każdego elementu.

!!! note "Uwaga"

    **Closure** to blok kodu, który może być przekazywany i wykonywany później. Pomyśl o tym jak o funkcji, którą definiujesz inline. Closures są zapisywane za pomocą nawiasów klamrowych `{ }` i mogą przyjmować parametry. Są fundamentalne dla działania operatorów Nextflow i jeśli od jakiegoś czasu piszesz w Nextflow, możliwe, że już ich używałeś, nie zdając sobie z tego sprawy!

Oto jak wygląda ta operacja map:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

To nasze pierwsze **closure** - funkcja anonimowa, którą możesz przekazać jako argument (podobnie jak lambdy w Python lub funkcje strzałkowe w JavaScript). Closures są niezbędne do pracy z operatorami Nextflow.

Closure `{ row -> return row }` przyjmuje parametr `row` (może być dowolna nazwa: `item`, `sample`, itp.).

Gdy operator `.map()` przetwarza każdy element kanału, przekazuje ten element do Twojego closure. Tutaj `row` przechowuje jeden wiersz CSV na raz.

Zastosuj tę zmianę i uruchom workflow:

```bash
nextflow run main.nf
```

Zobaczysz to samo wyjście co poprzednio, ponieważ po prostu zwracamy wejście bez zmian. To potwierdza, że operator map działa poprawnie. Teraz zacznijmy transformować dane.

#### 1.1.3. Tworzenie Struktury Danych Map

Teraz napiszemy logikę **skryptowania** wewnątrz naszego closure, aby transformować każdy wiersz danych. To tutaj przetwarzamy pojedyncze elementy danych zamiast orkiestrować przepływ danych.

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Skryptowanie do transformacji danych
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

Mapa `sample_meta` to struktura danych klucz-wartość (jak słowniki w Python, obiekty w JavaScript lub hasze w Ruby) przechowująca powiązane informacje: ID próbki, organizm, typ tkanki, głębokość sekwencjonowania i wynik jakości.

Używamy metod manipulacji ciągami znaków jak `.toLowerCase()` i `.replaceAll()` do czyszczenia naszych danych oraz metod konwersji typów jak `.toInteger()` i `.toDouble()` do konwersji danych tekstowych z CSV na odpowiednie typy numeryczne.

Zastosuj tę zmianę i uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Dodawanie Logiki Warunkowej

Teraz dodajmy więcej skryptowania - tym razem używając operatora ternarnego do podejmowania decyzji na podstawie wartości danych.

Wprowadź następującą zmianę:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

Operator ternarny to skrót dla instrukcji if/else zgodnie ze wzorcem `warunek ? wartość_jeśli_prawda : wartość_jeśli_fałsz`. Ta linia oznacza: "Jeśli jakość jest większa niż 40, użyj 'high', w przeciwnym razie użyj 'normal'". Jego kuzyn, **operator Elvis** (`?:`), zapewnia wartości domyślne, gdy coś jest null lub puste - zbadamy ten wzorzec później w tym tutorialu.

Operator dodawania map `+` tworzy **nową mapę** zamiast modyfikować istniejącą. Ta linia tworzy nową mapę zawierającą wszystkie pary klucz-wartość z `sample_meta` plus nowy klucz `priority`.

!!! Note "Uwaga"

    Nigdy nie modyfikuj map przekazanych do closures - zawsze twórz nowe używając `+` (na przykład). W Nextflow te same dane często przepływają przez wiele operacji jednocześnie. Modyfikowanie mapy w miejscu może powodować nieprzewidywalne efekty uboczne, gdy inne operacje odwołują się do tego samego obiektu. Tworzenie nowych map zapewnia, że każda operacja ma własną czystą kopię.

Uruchom zmodyfikowany workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Pomyślnie dodaliśmy logikę warunkową do wzbogacenia naszych metadanych o poziom priorytetu oparty na wynikach jakości.

#### 1.1.5. Wybieranie Podzbiorów Map za pomocą `.subMap()`

Podczas gdy operator `+` dodaje klucze do mapy, czasami musisz zrobić coś odwrotnego - wyodrębnić tylko określone klucze. Metoda `.subMap()` jest do tego idealna.

Dodajmy linię tworzącą uproszczoną wersję naszych metadanych zawierającą tylko pola identyfikacyjne:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Skryptowanie do transformacji danych
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "Tylko pola ID: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Skryptowanie do transformacji danych
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Uruchom zmodyfikowany workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    Tylko pola ID: [id:sample_001, organism:human, tissue:liver]
    Tylko pola ID: [id:sample_002, organism:mouse, tissue:brain]
    Tylko pola ID: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

To pokazuje zarówno pełne metadane wyświetlone przez operację `view()`, jak i wyodrębniony podzbiór, który wydrukowaliśmy za pomocą `println`.

Metoda `.subMap()` przyjmuje listę kluczy i zwraca nową mapę zawierającą tylko te klucze. Jeśli klucz nie istnieje w oryginalnej mapie, po prostu nie jest uwzględniany w wyniku.

Jest to szczególnie przydatne, gdy musisz tworzyć różne wersje metadanych dla różnych procesów - niektóre mogą potrzebować pełnych metadanych, podczas gdy inne potrzebują tylko minimalnych pól identyfikacyjnych.

Teraz usuń te instrukcje println, aby przywrócić workflow do poprzedniego stanu, ponieważ nie będą nam potrzebne dalej.

!!! tip "Podsumowanie Operacji na Mapach"

    - **Dodawanie kluczy**: `map1 + [new_key: value]` - Tworzy nową mapę z dodatkowymi kluczami
    - **Wyodrębnianie kluczy**: `map1.subMap(['key1', 'key2'])` - Tworzy nową mapę z tylko określonymi kluczami
    - **Obie operacje tworzą nowe mapy** - Oryginalne mapy pozostają niezmienione

#### 1.1.6. Łączenie Map i Zwracanie Wyników

Do tej pory zwracaliśmy tylko to, co społeczność Nextflow nazywa 'meta map', i ignorowaliśmy pliki, do których te metadane się odnoszą. Ale jeśli piszesz workflow w Nextflow, prawdopodobnie chcesz coś zrobić z tymi plikami.

Wyprowadźmy strukturę kanału składającą się z krotki 2 elementów: wzbogaconej mapy metadanych i odpowiadającej ścieżki pliku. Jest to powszechny wzorzec w Nextflow do przekazywania danych do procesów.

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Zastosuj tę zmianę i uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Ta struktura krotki `[meta, file]` to powszechny wzorzec w Nextflow do przekazywania zarówno metadanych, jak i powiązanych plików do procesów.

!!! note "Uwaga"

    **Mapy i Metadane**: Mapy są fundamentalne dla pracy z metadanymi w Nextflow. Bardziej szczegółowe wyjaśnienie pracy z mapami metadanych znajdziesz w side queście [Praca z metadanymi](./metadata.md).

Nasz workflow demonstruje główny wzorzec: **operacje przepływu danych** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orkiestrują, jak dane przemieszczają się przez pipeline, podczas gdy **skryptowanie** (mapy `[key: value]`, metody na ciągach znaków, konwersje typów, operatory ternarne) wewnątrz closure `.map()` obsługuje transformację poszczególnych elementów danych.

### 1.2. Zrozumienie Różnych Typów: Channel vs List

Do tej pory wszystko dobrze, możemy rozróżnić operacje przepływu danych od skryptowania. Ale co z przypadkiem, gdy ta sama nazwa metody istnieje w obu kontekstach?

Doskonałym przykładem jest metoda `collect`, która istnieje zarówno dla typów kanałów, jak i typów List w bibliotece standardowej Nextflow. Metoda `collect()` na List transformuje każdy element, podczas gdy operator `collect()` na kanale zbiera wszystkie emisje kanału w kanał jednoelementowy.

Zademonstrujmy to na przykładowych danych, zaczynając od odświeżenia wiedzy o tym, co robi operator `collect()` na kanale. Sprawdź `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - grupuje wiele emisji kanału w jedną
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Pojedynczy element kanału: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "wynik channel.collect(): ${list} (${list.size()} elementów zgrupowanych w 1)" }
```

Kroki:

- Definiujemy listę ID próbek
- Tworzymy kanał za pomocą `fromList()`, który emituje każde ID próbki osobno
- Drukujemy każdy element za pomocą `view()` w miarę przepływu
- Zbieramy wszystkie elementy w jedną listę za pomocą operatora `collect()` kanału
- Drukujemy zebrany wynik (pojedynczy element zawierający wszystkie ID próbek) za pomocą drugiego `view()`

Zmieniliśmy strukturę kanału, ale nie zmieniliśmy samych danych.

Uruchom workflow, aby to potwierdzić:

```bash
nextflow run collect.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Pojedynczy element kanału: sample_001
    Pojedynczy element kanału: sample_002
    Pojedynczy element kanału: sample_003
    wynik channel.collect(): [sample_001, sample_002, sample_003] (3 elementów zgrupowanych w 1)
    ```

`view()` zwraca wyjście dla każdej emisji kanału, więc wiemy, że to pojedyncze wyjście zawiera wszystkie 3 oryginalne elementy zgrupowane w jedną listę.

Teraz zobaczmy metodę `collect` na List w akcji. Zmodyfikuj `collect.nf`, aby zastosować metodę `collect` List do oryginalnej listy ID próbek:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Pojedynczy element kanału: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "wynik channel.collect(): ${list} (${list.size()} elementów zgrupowanych w 1)" }

    // List.collect() - transformuje każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "wynik List.collect(): ${formatted_ids} (${sample_ids.size()} elementów przekształconych w ${formatted_ids.size()})"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Pojedynczy element kanału: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "wynik channel.collect(): ${list} (${list.size()} elementów zgrupowanych w 1)" }
    ```

W tym nowym fragmencie:

- Definiujemy nową zmienną `formatted_ids`, która używa metody `collect` List do transformacji każdego ID próbki w oryginalnej liście
- Drukujemy wynik używając `println`

Uruchom zmodyfikowany workflow:

```bash
nextflow run collect.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    wynik List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elementów przekształconych w 3)
    Pojedynczy element kanału: sample_001
    Pojedynczy element kanału: sample_002
    Pojedynczy element kanału: sample_003
    wynik channel.collect(): [sample_001, sample_002, sample_003] (3 elementów zgrupowanych w 1)
    ```

Tym razem NIE zmieniliśmy struktury danych, wciąż mamy 3 elementy na liście, ale przekształciliśmy każdy element używając metody `collect` List, aby uzyskać nową listę ze zmodyfikowanymi wartościami. To podobne do użycia operatora `map` na kanale, ale operuje na strukturze danych List zamiast na kanale.

`collect` to ekstremalny przypadek, którego używamy tutaj, aby podkreślić punkt. Kluczową lekcją jest to, że pisząc workflow, zawsze rozróżniaj między **strukturami danych** (Lists, Maps, itp.) a **kanałami** (konstrukcje przepływu danych). Operacje mogą mieć te same nazwy, ale zachowują się zupełnie inaczej w zależności od typu, na którym są wywoływane.

### 1.3. Operator Rozprzestrzeniania (`*.`) - Skrót do Wyodrębniania Właściwości

Powiązany z metodą `collect` List jest operator rozprzestrzeniania (`*.`), który zapewnia zwięzły sposób wyodrębniania właściwości z kolekcji. Jest to zasadniczo syntaktyczny cukier dla powszechnego wzorca `collect`.

Dodajmy demonstrację do naszego pliku `collect.nf`:

=== "Po"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Pojedynczy element kanału: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "wynik channel.collect(): ${list} (${list.size()} elementów zgrupowanych w 1)" }

    // List.collect() - transformuje każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "wynik List.collect(): ${formatted_ids} (${sample_ids.size()} elementów przekształconych w ${formatted_ids.size()})"

    // Operator rozprzestrzeniania - zwięzły dostęp do właściwości
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Wynik operatora rozprzestrzeniania: ${all_ids}"
    ```

=== "Przed"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Pojedynczy element kanału: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "wynik channel.collect(): ${list} (${list.size()} elementów zgrupowanych w 1)" }

    // List.collect() - transformuje każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "wynik List.collect(): ${formatted_ids} (${sample_ids.size()} elementów przekształconych w ${formatted_ids.size()})"
    ```

Uruchom zaktualizowany workflow:

```bash title="Testuj operator rozprzestrzeniania"
nextflow run collect.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    wynik List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elementów przekształconych w 3)
    Wynik operatora rozprzestrzeniania: [s1, s2, s3]
    Pojedynczy element kanału: sample_001
    Pojedynczy element kanału: sample_002
    Pojedynczy element kanału: sample_003
    wynik channel.collect(): [sample_001, sample_002, sample_003] (3 elementów zgrupowanych w 1)
    ```

Operator rozprzestrzeniania `*.` jest skrótem dla powszechnego wzorca collect:

```groovy
// Są równoważne:
def ids = samples*.id
def ids = samples.collect { it.id }

// Działa również z wywołaniami metod:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Operator rozprzestrzeniania jest szczególnie przydatny, gdy musisz wyodrębnić pojedynczą właściwość z listy obiektów - jest bardziej czytelny niż pisanie pełnego closure `collect`.

!!! tip "Kiedy Używać Spread vs Collect"

    - **Użyj spread (`*.`)** do prostego dostępu do właściwości: `samples*.id`, `files*.name`
    - **Użyj collect** do transformacji lub złożonej logiki: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Wnioski

W tej sekcji nauczyłeś się:

- **Przepływ danych vs skryptowanie**: Operatory kanałów orkiestrują, jak dane przepływają przez Twój pipeline, podczas gdy skryptowanie transformuje poszczególne elementy danych
- **Zrozumienie typów**: Ta sama nazwa metody (jak `collect`) może zachowywać się inaczej w zależności od typu, na którym jest wywoływana (Channel vs List)
- **Kontekst ma znaczenie**: Zawsze bądź świadomy, czy pracujesz z kanałami (przepływ danych) czy strukturami danych (skryptowanie)

Zrozumienie tych granic jest niezbędne do debugowania, dokumentacji i pisania łatwych w utrzymaniu workflow.

Następnie zagłębimy się w możliwości przetwarzania ciągów znaków, które są niezbędne do obsługi rzeczywistych danych.

---

## 2. Przetwarzanie Ciągów Znaków i Dynamiczne Generowanie Skryptów

Opanowanie przetwarzania ciągów znaków oddziela kruche workflow od solidnych pipelineów. Ta sekcja obejmuje parsowanie złożonych nazw plików, dynamiczne generowanie skryptów i interpolację zmiennych.

### 2.1. Dopasowywanie Wzorców i Wyrażenia Regularne

Pliki bioinformatyczne często mają złożone konwencje nazewnictwa kodujące metadane. Wyodrębnijmy to automatycznie używając dopasowywania wzorców z wyrażeniami regularnymi.

Wrócimy do naszego workflow `main.nf` i dodamy logikę dopasowywania wzorców, aby wyodrębnić dodatkowe informacje o próbkach z nazw plików. Pliki FASTQ w naszym zbiorze danych podążają za konwencjami nazewnictwa Illumina z nazwami jak `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Mogą wyglądać enigmatycznie, ale faktycznie kodują użyteczne metadane jak ID próbki, numer pasa i kierunek odczytu. Użyjemy możliwości regex do parsowania tych nazw.

Wprowadź następującą zmianę do istniejącego workflow `main.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Skryptowanie do transformacji danych
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Skryptowanie do transformacji danych
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

To demonstruje kluczowe **koncepty przetwarzania ciągów znaków**:

1. **Literały wyrażeń regularnych** używające składni `~/wzorzec/` - tworzy wzorzec regex bez konieczności escape'owania ukośników odwrotnych
2. **Dopasowywanie wzorców** z operatorem `=~` - próbuje dopasować ciąg znaków do wzorca regex
3. **Obiekty matcher**, które przechwytują grupy z `[0][1]`, `[0][2]`, itd. - `[0]` odnosi się do całego dopasowania, `[1]`, `[2]`, itd. odnoszą się do przechwyconych grup w nawiasach

Przeanalizujmy wzorzec regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Wzorzec             | Dopasowuje                                  | Przechwytuje                          |
| ------------------- | ------------------------------------------- | ------------------------------------- |
| `^(.+)`             | Nazwa próbki od początku                    | Grupa 1: nazwa próbki                 |
| `_S(\d+)`           | Numer próbki `_S1`, `_S2`, itd.             | Grupa 2: numer próbki                 |
| `_L(\d{3})`         | Numer pasa `_L001`                          | Grupa 3: pas (3 cyfry)                |
| `_(R[12])`          | Kierunek odczytu `_R1` lub `_R2`            | Grupa 4: kierunek odczytu             |
| `_(\d{3})`          | Numer fragmentu `_001`                      | Grupa 5: fragment (3 cyfry)           |
| `\.fastq(?:\.gz)?$` | Rozszerzenie pliku `.fastq` lub `.fastq.gz` | Nieprzechwycone (?: to non-capturing) |

To parsuje konwencje nazewnictwa Illumina, aby automatycznie wyodrębnić metadane.

Uruchom zmodyfikowany workflow:

```bash title="Testuj dopasowywanie wzorców"
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

To pokazuje metadane wzbogacone z nazw plików.

### 2.2. Dynamiczne Generowanie Skryptów w Procesach

Bloki skryptowe procesów są zasadniczo wieloliniowymi ciągami znaków, które są przekazywane do powłoki. Możesz używać **logiki warunkowej** (if/else, operatory ternarne) do dynamicznego generowania różnych ciągów skryptowych na podstawie charakterystyki wejścia. Jest to niezbędne do obsługi różnych typów wejściowych—jak odczyty single-end vs paired-end—bez duplikowania definicji procesów.

Dodajmy proces do naszego workflow, który demonstruje ten wzorzec. Otwórz `modules/fastp.nf` i spójrz:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

Proces przyjmuje pliki FASTQ jako wejście i uruchamia narzędzie `fastp` do przycinania adapterów i filtrowania odczytów niskiej jakości. Niestety, osoba, która napisała ten proces, nie uwzględniła odczytów single-end, które mamy w naszym przykładowym zbiorze danych. Dodajmy go do naszego workflow i zobaczmy, co się stanie:

Najpierw dołącz moduł w samej pierwszej linii Twojego workflow `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Następnie zmodyfikuj blok `workflow`, aby połączyć kanał `ch_samples` z procesem `FASTP`:

=== "Po"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Uruchom ten zmodyfikowany workflow:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Widać, że proces próbuje uruchomić `fastp` z wartością `null` dla drugiego pliku wejściowego, co powoduje błąd. Dzieje się tak, ponieważ nasz zbiór danych zawiera odczyty single-end, ale proces jest zakodowany na stałe do oczekiwania odczytów paired-end (dwa pliki wejściowe na raz).

Napraw to, dodając logikę warunkową do bloku `script:` procesu `FASTP`. Instrukcja if/else sprawdza liczbę plików odczytowych i odpowiednio dostosowuje polecenie.

=== "Po"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Proste wykrywanie single-end vs paired-end
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Teraz workflow może płynnie obsługiwać zarówno odczyty single-end, jak i paired-end. Logika warunkowa sprawdza liczbę plików wejściowych i konstruuje odpowiednie polecenie dla `fastp`. Zobaczmy, czy działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

Wygląda dobrze! Jeśli sprawdzimy faktyczne polecenia, które zostały uruchomione (dostosuj dla Swojego hasha zadania):

```console title="Sprawdź wykonane polecenia"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Możemy zobaczyć, że Nextflow poprawnie wybrał właściwe polecenie dla odczytów single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Inne powszechne zastosowanie dynamicznej logiki skryptowej można zobaczyć w [module Genomics kursu Nextflow for Science](../../nf4science/genomics/02_joint_calling). W tym module wywoływany proces GATK może przyjmować wiele plików wejściowych, ale każdy musi być poprzedzony `-V`, aby utworzyć poprawną linię poleceń. Proces używa skryptowania do transformacji kolekcji plików wejściowych (`all_gvcfs`) w poprawne argumenty polecenia:

```groovy title="manipulacja linią poleceń dla GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Te wzorce używania skryptowania w blokach skryptowych procesów są niezwykle potężne i mogą być zastosowane w wielu scenariuszach - od obsługi zmiennych typów wejściowych po budowanie złożonych argumentów linii poleceń z kolekcji plików, czyniąc Twoje procesy naprawdę adaptacyjnymi do różnorodnych wymagań rzeczywistych danych.

### 2.3. Interpolacja Zmiennych: Zmienne Nextflow i Shell

Skrypty procesów mieszają zmienne Nextflow, zmienne powłoki i podstawienia poleceń, każde z inną składnią interpolacji. Użycie niewłaściwej składni powoduje błędy. Przeanalizujmy to na procesie tworzącym raport przetwarzania.

Spójrz na plik modułu `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Przetwarzanie ${reads}" > ${meta.id}_report.txt
    echo "Próbka: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Ten proces zapisuje prosty raport z ID próbki i nazwą pliku. Teraz uruchommy go, aby zobaczyć, co się dzieje, gdy musimy mieszać różne typy zmiennych.

Dołącz proces w Swoim `main.nf` i dodaj go do workflow:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

Teraz uruchom workflow i sprawdź wygenerowane raporty w `results/reports/`. Powinny zawierać podstawowe informacje o każdej próbce.

<!-- TODO: add the run command -->

??? success "Wyjście polecenia"

    ```console
    <!-- TODO: output -->
    ```

Ale co, jeśli chcemy dodać informacje o tym, kiedy i gdzie przetwarzanie miało miejsce? Z
