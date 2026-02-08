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

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_hawking] DSL2 - revision: 8a9f8d4e61

    executor >  local (6)
    [4d/e1c9f2] process > FASTP (2)           [100%] 3 of 3 ✔
    [2c/7b8a43] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Ale co, jeśli chcemy dodać informacje o tym, kiedy i gdzie przetwarzanie miało miejsce? Zmodyfikujmy proces, aby używał zmiennych **shell** i trochę podstawienia poleceń, aby uwzględnić bieżącego użytkownika, nazwę hosta i datę w raporcie:

=== "Po"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Przetwarzanie ${reads}" > ${meta.id}_report.txt
        echo "Próbka: ${meta.id}" >> ${meta.id}_report.txt
        echo "Przetworzony przez: ${USER}" >> ${meta.id}_report.txt
        echo "Nazwa hosta: $(hostname)" >> ${meta.id}_report.txt
        echo "Data: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Przed"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Przetwarzanie ${reads}" > ${meta.id}_report.txt
        echo "Próbka: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Jeśli uruchomisz to, zauważysz błąd - Nextflow próbuje zinterpretować `${USER}` jako zmienną Nextflow, która nie istnieje.

??? failure "Wyjście polecenia"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Przetworzony przez: ${USER}" >> ${meta.id}_report.txt
    ╰     |                               ^^^^

    ERROR ~ Script compilation failed
    ```

Musimy to zescapować, aby Bash mógł to obsłużyć zamiast Nextflow.

Napraw to, escapując zmienne shell i podstawienia poleceń za pomocą ukośnika wstecznego (`\`):

=== "Po"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Przetwarzanie ${reads}" > ${meta.id}_report.txt
        echo "Próbka: ${meta.id}" >> ${meta.id}_report.txt
        echo "Przetworzony przez: \${USER}" >> ${meta.id}_report.txt
        echo "Nazwa hosta: \$(hostname)" >> ${meta.id}_report.txt
        echo "Data: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Przed"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Przetwarzanie ${reads}" > ${meta.id}_report.txt
        echo "Próbka: ${meta.id}" >> ${meta.id}_report.txt
        echo "Przetworzony przez: ${USER}" >> ${meta.id}_report.txt
        echo "Nazwa hosta: $(hostname)" >> ${meta.id}_report.txt
        echo "Data: $(date)" >> ${meta.id}_report.txt
        """
    ```

Teraz działa! Ukośnik wsteczny (`\`) mówi Nextflow "nie interpretuj tego, przekaż to do Bash."

### Wnioski

W tej sekcji nauczyłeś się technik **przetwarzania ciągów znaków**:

- **Wyrażenia regularne do parsowania plików**: Używanie operatora `=~` i wzorców regex (`~/wzorzec/`) do wyodrębniania metadanych ze złożonych konwencji nazewnictwa plików
- **Dynamiczne generowanie skryptów**: Używanie logiki warunkowej (if/else, operatory ternarne) do generowania różnych ciągów skryptowych na podstawie charakterystyki wejścia
- **Interpolacja zmiennych**: Zrozumienie, kiedy Nextflow interpretuje ciągi znaków, a kiedy robi to powłoka
  - `${var}` - zmienne Nextflow (interpolowane przez Nextflow w czasie kompilacji workflow)
  - `\${var}` - zmienne środowiskowe powłoki (escapowane, przekazywane do bash w czasie uruchomienia)
  - `\$(cmd)` - podstawienie polecenia powłoki (escapowane, wykonywane przez bash w czasie uruchomienia)

Te wzorce przetwarzania i generowania ciągów znaków są niezbędne do obsługi różnorodnych formatów plików i konwencji nazewnictwa, które napotkasz w rzeczywistych workflow bioinformatycznych.

---

## 3. Tworzenie Funkcji Wielokrotnego Użytku

Złożona logika workflow inline w operatorach kanałów lub definicjach procesów zmniejsza czytelność i łatwość w utrzymaniu. **Funkcje** pozwalają wyodrębnić tę logikę do nazwanych, wielokrotnego użytku komponentów.

Nasza operacja map stała się długa i złożona. Wyodrębnijmy ją do funkcji wielokrotnego użytku używając słowa kluczowego `def`.

Aby zilustrować, jak to wygląda z naszym istniejącym workflow, wprowadź modyfikację poniżej, używając `def` do zdefiniowania funkcji wielokrotnego użytku o nazwie `separateMetadata`:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
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

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
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

Wyodrębniając tę logikę do funkcji, zredukowaliśmy faktyczną logikę workflow do czegoś znacznie czystszego:

```groovy title="minimalny workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

To sprawia, że logika workflow jest o wiele łatwiejsza do przeczytania i zrozumienia na pierwszy rzut oka. Funkcja `separateMetadata` enkapsuluje całą złożoną logikę parsowania i wzbogacania metadanych, czyniąc ją wielokrotnego użytku i łatwą do przetestowania.

Uruchom workflow, aby upewnić się, że nadal działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

Wyjście powinno pokazywać, że oba procesy zakończyły się pomyślnie. Workflow jest teraz o wiele czystszy i łatwiejszy w utrzymaniu, z całą złożoną logiką przetwarzania metadanych enkapsulowaną w funkcji `separateMetadata`.

### Wnioski

W tej sekcji nauczyłeś się **tworzenia funkcji**:

- **Definiowanie funkcji za pomocą `def`**: Słowo kluczowe do tworzenia nazwanych funkcji (jak `def` w Python lub `function` w JavaScript)
- **Zakres funkcji**: Funkcje zdefiniowane na poziomie skryptu są dostępne w całym workflow Nextflow
- **Wartości zwracane**: Funkcje automatycznie zwracają ostatnie wyrażenie lub używają jawnego `return`
- **Czystszy kod**: Wyodrębnianie złożonej logiki do funkcji to fundamentalna praktyka inżynierii oprogramowania w każdym języku

Następnie zbadamy, jak używać closures w dyrektywach procesów do dynamicznej alokacji zasobów.

---

## 4. Dynamiczne Dyrektywy Zasobów z Closures

Do tej pory używaliśmy skryptowania w bloku `script` procesów. Ale **closures** (wprowadzone w Sekcji 1.1) są również niezwykle przydatne w dyrektywach procesów, szczególnie do dynamicznej alokacji zasobów. Dodajmy dyrektywy zasobów do naszego procesu FASTP, które dostosowują się na podstawie charakterystyki próbki.

### 4.1. Alokacja zasobów specyficzna dla próbki

Obecnie nasz proces FASTP używa domyślnych zasobów. Uczyńmy go mądrzejszym, alokując więcej CPU dla próbek o wysokiej głębokości. Edytuj `modules/fastp.nf`, aby uwzględnić dynamiczną dyrektywę `cpus` i statyczną dyrektywę `memory`:

=== "Po"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Przed"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

Closure `{ meta.depth > 40000000 ? 2 : 1 }` używa **operatora ternarnego** (omówionego w Sekcji 1.1) i jest ewaluowane dla każdego zadania, umożliwiając alokację zasobów per próbka. Próbki o wysokiej głębokości (>40M odczytów) otrzymują 2 CPU, podczas gdy inne otrzymują 1 CPU.

!!! note "Dostęp do Zmiennych Wejściowych w Dyrektywach"

    Closure może uzyskać dostęp do dowolnych zmiennych wejściowych (jak `meta` tutaj), ponieważ Nextflow ewaluuje te closures w kontekście każdego wykonania zadania.

Uruchom workflow ponownie z opcją `-ansi-log false`, aby ułatwić zobaczenie hashów zadań.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Możesz sprawdzić dokładne polecenie `docker`, które zostało uruchomione, aby zobaczyć alokację CPU dla dowolnego zadania:

```console title="Sprawdź polecenie docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Powinieneś zobaczyć coś takiego:

```bash title="polecenie docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

W tym przykładzie wybraliśmy przykład, który zażądał 2 CPU (`--cpu-shares 2048`), ponieważ była to próbka o wysokiej głębokości, ale powinieneś zobaczyć różne alokacje CPU w zależności od głębokości próbki. Spróbuj tego również dla innych zadań.

### 4.2. Strategie ponownych prób

Innym potężnym wzorcem jest użycie `task.attempt` do strategii ponownych prób. Aby pokazać, dlaczego jest to przydatne, zaczniemy od zmniejszenia alokacji pamięci dla FASTP do mniej niż jest mu potrzebne. Zmień dyrektywę `memory` w `modules/fastp.nf` na `1.GB`:

=== "Po"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Przed"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... i uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

To wskazuje, że proces został zabity za przekroczenie limitów pamięci.

To bardzo powszechny scenariusz w rzeczywistych workflow - czasami po prostu nie wiesz, ile pamięci będzie potrzebne zadaniu, dopóki go nie uruchomisz.

Aby uczynić nasz workflow bardziej odpornym, możemy zaimplementować strategię ponownych prób, która zwiększa alokację pamięci przy każdej próbie, ponownie używając closure Groovy. Zmodyfikuj dyrektywę `memory`, aby pomnożyć pamięć bazową przez `task.attempt`, i dodaj dyrektywy `errorStrategy 'retry'` oraz `maxRetries 2`:

=== "Po"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Przed"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Teraz, jeśli proces zawiedzie z powodu niewystarczającej pamięci, Nextflow ponowi próbę z większą ilością pamięci:

- Pierwsza próba: 1 GB (task.attempt = 1)
- Druga próba: 2.GB (task.attempt = 2)

... i tak dalej, aż do limitu `maxRetries`.

### Wnioski

Dynamiczne dyrektywy z closures pozwalają Ci:

- Alokować zasoby na podstawie charakterystyki wejścia
- Implementować automatyczne strategie ponownych prób z rosnącymi zasobami
- Łączyć wiele czynników (metadane, numer próby, priorytety)
- Używać logiki warunkowej do złożonych obliczeń zasobów

To sprawia, że Twoje workflow są zarówno bardziej wydajne (nie nadmierna alokacja), jak i bardziej odporne (automatyczne ponowienie próby z większymi zasobami).

---

## 5. Logika Warunkowa i Kontrola Procesów

Wcześniej używaliśmy `.map()` ze skryptowaniem do transformacji danych kanału. Teraz użyjemy logiki warunkowej, aby kontrolować, które procesy są wykonywane na podstawie danych—niezbędne dla elastycznych workflow dostosowujących się do różnych typów próbek.

[Operatory przepływu danych](https://www.nextflow.io/docs/latest/reference/operator.html) Nextflow przyjmują closures ewaluowane w czasie wykonania, umożliwiając logice warunkowej kierowanie decyzjami workflow na podstawie zawartości kanału.

### 5.1. Kierowanie za pomocą `.branch()`

Na przykład, powiedzmy, że nasze próbki sekwencjonowania muszą być przycięte za pomocą FASTP tylko wtedy, gdy są próbkami ludzkimi z pokryciem powyżej pewnego progu. Próbki myszy lub próbki o niskim pokryciu powinny być uruchomione z Trimgalore zamiast tego (to wymyślony przykład, ale ilustruje punkt).

Dostarczyliśmy prosty proces Trimgalore w `modules/trimgalore.nf`, spójrz, jeśli chcesz, ale szczegóły nie są ważne dla tego ćwiczenia. Kluczowym punktem jest to, że chcemy kierować próbki na podstawie ich metadanych.

Dołącz nowy moduł z `modules/trimgalore.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... a następnie zmodyfikuj swój workflow `main.nf`, aby rozgałęzić próbki na podstawie ich metadanych i kierować je przez odpowiedni proces przycinania, w ten sposób:

=== "Po"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Uruchom ten zmodyfikowany workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Tutaj użyliśmy małych, ale potężnych wyrażeń warunkowych wewnątrz operatora `.branch{}`, aby kierować próbki na podstawie ich metadanych. Próbki ludzkie o wysokim pokryciu przechodzą przez `FASTP`, podczas gdy wszystkie inne próbki przechodzą przez `TRIMGALORE`.

### 5.2. Używanie `.filter()` z Prawdziwością

Innym potężnym wzorcem do kontrolowania wykonania workflow jest operator `.filter()`, który używa closure do określenia, które elementy powinny kontynuować w pipeline. Wewnątrz closure filtra napiszesz **wyrażenia boolowskie**, które decydują, które elementy przejdą.

Nextflow (jak wiele języków dynamicznych) ma koncepcję **"prawdziwości"**, która określa, jakie wartości są ewaluowane jako `true` lub `false` w kontekstach boolowskich:

- **Prawdziwe**: Wartości nie-null, niepuste ciągi znaków, liczby niezerowe, niepuste kolekcje
- **Fałszywe**: `null`, puste ciągi znaków `""`, zero `0`, puste kolekcje `[]` lub `[:]`, `false`

To oznacza, że samo `meta.id` (bez jawnego `!= null`) sprawdza, czy ID istnieje i nie jest puste. Użyjmy tego do odfiltrowania próbek, które nie spełniają naszych wymagań jakościowych.

Dodaj następujące przed operacją branch:

=== "Po"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Odfiltruj nieprawidłowe lub niskiej jakości próbki
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Ponieważ wybraliśmy filtr, który wyklucza niektóre próbki, wykonano mniej zadań.

Wyrażenie filtra `meta.id && meta.organism && meta.depth >= 25000000` łączy prawdziwość z jawnymi porównaniami:

- `meta.id && meta.organism` sprawdza, czy oba pola istnieją i nie są puste (używając prawdziwości)
- `meta.depth >= 25000000` zapewnia wystarczającą głębokość sekwencjonowania za pomocą jawnego porównania

!!! note "Prawdziwość w Praktyce"

    Wyrażenie `meta.id && meta.organism` jest bardziej zwięzłe niż pisanie:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    To sprawia, że logika filtrowania jest o wiele czystsza i łatwiejsza do odczytania.

### Wnioski

W tej sekcji nauczyłeś się używać logiki warunkowej do kontrolowania wykonania workflow, używając interfejsów closure operatorów Nextflow jak `.branch{}` i `.filter{}`, wykorzystując prawdziwość do pisania zwięzłych wyrażeń warunkowych.

Nasz pipeline teraz inteligentnie kieruje próbki przez odpowiednie procesy, ale workflow produkcyjne muszą obsługiwać nieprawidłowe dane z wdziękiem. Uczyńmy nasz workflow odpornym na brakujące lub null wartości.

---

## 6. Bezpieczna Nawigacja i Operatory Elvis

Nasza funkcja `separateMetadata` obecnie zakłada, że wszystkie pola CSV są obecne i prawidłowe. Ale co się dzieje z niepełnymi danymi? Dowiedzmy się.

### 6.1. Problem: Dostęp do Właściwości, Które Nie Istnieją

Powiedzmy, że chcemy dodać wsparcie dla opcjonalnych informacji o sekwencjonowaniu. W niektórych laboratoriach próbki mogą mieć dodatkowe pole dla ID uruchomienia sekwencjonowania lub numeru partii, ale nasz obecny CSV nie ma tej kolumny. Spróbujmy uzyskać do niej dostęp mimo to.

Zmodyfikuj funkcję `separateMetadata`, aby uwzględnić pole run_id:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Teraz uruchom workflow:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

To powoduje awarię z NullPointerException.

Problem polega na tym, że `row.run_id` zwraca `null`, ponieważ kolumna `run_id` nie istnieje w naszym CSV. Gdy próbujemy wywołać `.toUpperCase()` na `null`, dochodzi do awarii. Tu właśnie operator bezpiecznej nawigacji ratuje sytuację.

### 6.2. Operator Bezpiecznej Nawigacji (`?.`)

Operator bezpiecznej nawigacji (`?.`) zwraca `null` zamiast rzucać wyjątek, gdy jest wywoływany na wartości `null`. Jeśli obiekt przed `?.` jest `null`, całe wyrażenie ewaluuje się do `null` bez wykonywania metody.

Zaktualizuj funkcję, aby używała bezpiecznej nawigacji:

=== "Po"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Uruchom ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fervent_lichterman] DSL2 - revision: 7a8b4c9d2e

    executor >  local (6)
    [4e/2f3a91] process > FASTP (2)           [100%] 2 of 2 ✔
    [7c/8d5b42] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [3a/1e7c93] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Nie ma awarii! Workflow teraz obsługuje brakujące pole z wdziękiem. Gdy `row.run_id` jest `null`, operator `?.` zapobiega wywołaniu `.toUpperCase()`, a `run_id` staje się `null` zamiast powodować wyjątek.

### 6.3. Operator Elvis (`?:`) dla Wartości Domyślnych

Operator Elvis (`?:`) zapewnia wartości domyślne, gdy lewa strona jest "fałszywa" (jak wyjaśniono wcześniej). Nazywa się tak od Elvisa Presleya, ponieważ `?:` wygląda jak jego słynne włosy i oczy, gdy jest oglądany bokiem!

Teraz, gdy używamy bezpiecznej nawigacji, `run_id` będzie `null` dla próbek bez tego pola. Użyjmy operatora Elvis, aby zapewnić wartość domyślną i dodać ją do naszej mapy `sample_meta`:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Dodaj również operator `view()` w workflow, aby zobaczyć wyniki:

=== "Po"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

i uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfekcyjnie! Teraz wszystkie próbki mają pole `run` z ich rzeczywistym ID uruchomienia (wielkimi literami) lub wartością domyślną 'UNSPECIFIED'. Kombinacja `?.` i `?:` zapewnia zarówno bezpieczeństwo (bez awarii), jak i sensowne wartości domyślne.

Usuń teraz operator `.view()`, teraz gdy potwierdziliśmy, że działa.

!!! tip "Łączenie Bezpiecznej Nawigacji i Elvis"

    Wzorzec `value?.method() ?: 'default'` jest powszechny w workflow produkcyjnych:

    - `value?.method()` - Bezpiecznie wywołuje metodę, zwraca `null`, jeśli `value` jest `null`
    - `?: 'default'` - Zapewnia fallback, jeśli wynik jest `null`

    Ten wzorzec obsługuje brakujące/niepełne dane z wdziękiem.

Używaj tych operatorów konsekwentnie w funkcjach, closures operatorów (`.map{}`, `.filter{}`), skryptach procesów i plikach konfiguracyjnych. Zapobiegają awariam podczas obsługi rzeczywistych danych.

### Wnioski

- **Bezpieczna nawigacja (`?.`)**: Zapobiega awariam na wartościach null - zwraca null zamiast rzucać wyjątek
- **Operator Elvis (`?:`)**: Zapewnia wartości domyślne - `value ?: 'default'`
- **Łączenie**: `value?.method() ?: 'default'` to powszechny wzorzec

Te operatory sprawiają, że workflow są odporne na niepełne dane - niezbędne dla rzeczywistej pracy.

---

## 7. Walidacja z `error()` i `log.warn`

Czasami musisz natychmiast zatrzymać workflow, jeśli parametry wejściowe są nieprawidłowe. W Nextflow możesz używać wbudowanych funkcji jak `error()` i `log.warn`, a także standardowych konstrukcji programistycznych jak instrukcje `if` i logika boolowska, aby zaimplementować logikę walidacji. Dodajmy walidację do naszego workflow.

Utwórz funkcję walidacji przed blokiem workflow, wywołaj ją z workflow i zmień tworzenie kanału, aby używało parametru dla ścieżki pliku CSV. Jeśli parametr brakuje lub plik nie istnieje, wywołaj `error()`, aby zatrzymać wykonanie z wyraźnym komunikatem.

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Sprawdź, czy parametr wejściowy jest podany
        if (!params.input) {
            error("Ścieżka pliku CSV wejściowego nie jest podana. Proszę określić --input <file.csv>")
        }

        // Sprawdź, czy plik CSV istnieje
        if (!file(params.input).exists()) {
            error("Plik CSV wejściowy nie został znaleziony: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Teraz spróbuj uruchomić bez pliku CSV:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Ścieżka pliku CSV wejściowego nie jest podana. Proszę określić --input <file.csv>
    ```

Workflow zatrzymuje się natychmiast z wyraźnym komunikatem błędu zamiast zawieść tajemniczo później

Teraz uruchom go z nieistniejącym plikiem:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Plik CSV wejściowy nie został znaleziony: ./data/nonexistent.csv
    ```

Na koniec uruchom go z poprawnym plikiem:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_brahmagupta] DSL2 - revision: 4a8c9d2e5f

    executor >  local (7)
    [5e/3f4a92] process > FASTP (2)           [100%] 2 of 2 ✔
    [8c/6d7b43] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [4a/2e8c94] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Tym razem uruchamia się pomyślnie.

Możesz również dodać walidację wewnątrz funkcji `separateMetadata`. Użyjmy niefatalnego `log.warn` do wydawania ostrzeżeń dla próbek o niskiej głębokości sekwencjonowania, ale nadal pozwolić workflow na kontynuację:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Waliduj, czy dane mają sens
        if (sample_meta.depth < 30000000) {
            log.warn "Niska głębokość sekwencjonowania dla ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Uruchom workflow ponownie z oryginalnym CSV:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Niska głębokość sekwencjonowania dla sample_002: 25000000
    ```

Widzimy ostrzeżenie o niskiej głębokości sekwencjonowania dla jednej z próbek.

### Wnioski

- **`error()`**: Zatrzymuje workflow natychmiast z wyraźnym komunikatem
- **`log.warn`**: Wydaje ostrzeżenia bez zatrzymywania workflow
- **Wczesna walidacja**: Sprawdzaj wejścia przed przetwarzaniem, aby szybko zawieść z pomocnymi błędami
- **Funkcje walidacji**: Twórz wielokrotnego użytku logikę walidacji, którą można wywołać na początku workflow

Odpowiednia walidacja sprawia, że workflow są bardziej odporne i przyjazne dla użytkownika, wychwytując problemy wcześnie z wyraźnymi komunikatami błędów.

---

## 8. Handlery Zdarzeń Workflow

Do tej pory pisaliśmy kod w naszych skryptach workflow i definicjach procesów. Ale jest jeszcze jedna ważna funkcja, o której powinieneś wiedzieć: handlery zdarzeń workflow.

Handlery zdarzeń to closures, które są uruchamiane w określonych punktach cyklu życia Twojego workflow. Są idealne do dodawania logowania, powiadomień lub operacji czyszczenia. Te handlery powinny być zdefiniowane w Twoim skrypcie workflow obok definicji workflow.

### 8.1. Handler `onComplete`

Najczęściej używanym handlerem zdarzeń jest `onComplete`, który jest uruchamiany, gdy workflow się kończy (czy zakończy się sukcesem, czy niepowodzeniem). Dodajmy jeden, aby podsumować wyniki naszego pipeline.

Dodaj handler zdarzeń do Twojego pliku `main.nf`, wewnątrz definicji workflow:

=== "Po"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Podsumowanie wykonania pipeline:"
            println "=========================="
            println "Zakończono o: ${workflow.complete}"
            println "Czas trwania    : ${workflow.duration}"
            println "Sukces     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "status wyjścia : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

To closure uruchamia się, gdy workflow się kończy. Wewnątrz masz dostęp do obiektu `workflow`, który zapewnia użyteczne właściwości dotyczące wykonania.

Uruchom swój workflow, a zobaczysz to podsumowanie pojawiające się na końcu!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Niska głębokość sekwencjonowania dla sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Podsumowanie wykonania pipeline:
    ==========================
    Zakończono o: 2025-10-10T12:14:24.885384+01:00
    Czas trwania    : 2.9s
    Sukces     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    status wyjścia : 0
    ```

Uczyńmy to bardziej użytecznym, dodając logikę warunkową:

=== "Po"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Podsumowanie wykonania pipeline:"
            println "=========================="
            println "Zakończono o: ${workflow.complete}"
            println "Czas trwania    : ${workflow.duration}"
            println "Sukces     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "status wyjścia : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline zakończony pomyślnie!"
            } else {
                println "❌ Pipeline zakończony niepowodzeniem!"
                println "Błąd: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Podsumowanie wykonania pipeline:"
            println "=========================="
            println "Zakończono o: ${workflow.complete}"
            println "Czas trwania    : ${workflow.duration}"
            println "Sukces     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "status wyjścia : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Teraz otrzymujemy jeszcze bardziej informacyjne podsumowanie, w tym komunikat sukcesu/porażki i katalog wyjściowy, jeśli jest określony:

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Niska głębokość sekwencjonowania dla sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Podsumowanie wykonania pipeline:
    ==========================
    Zakończono o: 2025-10-10T12:16:00.522569+01:00
    Czas trwania    : 3.6s
    Sukces     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    status wyjścia : 0

    ✅ Pipeline zakończony pomyślnie!
    ```

Możesz również zapisać podsumowanie do pliku za pomocą operacji na plikach:

```groovy title="main.nf - Zapisywanie podsumowania do pliku"
workflow {
    // ... kod Twojego workflow ...

    workflow.onComplete = {
        def summary = """
        Podsumowanie Wykonania Pipeline
        ===========================
        Zakończono: ${workflow.complete}
        Czas trwania : ${workflow.duration}
        Sukces  : ${workflow.success}
        Polecenie  : ${workflow.commandLine}
        """

        println summary

        // Zapisz do pliku logu
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Handler `onError`

Oprócz `onComplete`, istnieje jeszcze jeden handler zdarzeń, którego możesz użyć: `onError`, który uruchamia się tylko wtedy, gdy workflow zawiedzie:

```groovy title="main.nf - handler onError"
workflow {
    // ... kod Twojego workflow ...

    workflow.onError = {
        println "="* 50
        println "Wykonanie pipeline zakończone niepowodzeniem!"
        println "Komunikat błędu: ${workflow.errorMessage}"
        println "="* 50

        // Zapisz szczegółowy log błędu
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Raport Błędu Workflow
        =====================
        Czas: ${new Date()}
        Błąd: ${workflow.errorMessage}
        Raport błędu: ${workflow.errorReport ?: 'Brak szczegółowego raportu dostępnego'}
        """

        println "Szczegóły błędu zapisane do: ${error_file}"
    }
}
```

Możesz używać wielu handlerów razem w swoim skrypcie workflow:

```groovy title="main.nf - Połączone handlery"
workflow {
    // ... kod Twojego workflow ...

    workflow.onError = {
        println "Workflow zakończony niepowodzeniem: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUKCES ✅" : "NIEPOWODZENIE ❌"

        println """
        Pipeline zakończony: ${status}
        Czas trwania: ${duration_mins} minut
        """
    }
}
```

### Wnioski

W tej sekcji nauczyłeś się:

- **Closures handlerów zdarzeń**: Closures w Twoim skrypcie workflow, które uruchamiają się w różnych punktach cyklu życia
- **Handler `onComplete`**: Do podsumowań wykonania i raportowania wyników
- **Handler `onError`**: Do obsługi błędów i logowania niepowodzeń
- **Właściwości obiektu workflow**: Dostęp do `workflow.success`, `workflow.duration`, `workflow.errorMessage`, itp.

Handlery zdarzeń pokazują, jak możesz używać pełnej mocy języka Nextflow w swoich skryptach workflow, aby dodać zaawansowane możliwości logowania i powiadamiania.

---

## Podsumowanie

Gratulacje, udało Ci się!

W trakcie tego side questa zbudowałeś kompleksowy pipeline przetwarzania próbek, który ewoluował od podstawowej obsługi metadanych do wyrafinowanego, gotowego do produkcji workflow.
Każda sekcja opierała się na poprzedniej, demonstrując, jak konstrukcje programistyczne przekształcają proste workflow w potężne systemy przetwarzania danych, z następującymi korzyściami:

- **Czystszy kod**: Zrozumienie przepływu danych vs skryptowanie pomaga pisać bardziej zorganizowane workflow
- **Odporna obsługa**: Bezpieczna nawigacja i operatory Elvis sprawiają, że workflow są odporne na brakujące dane
- **Elastyczne przetwarzanie**: Logika warunkowa pozwala Twoim workflow odpowiednio przetwarzać różne typy próbek
- **Adaptacyjne zasoby**: Dynamiczne dyrektywy optymalizują wykorzystanie zasobów na podstawie charakterystyki wejścia

Ta progresja odzwierciedla rzeczywistą ewolucję pipeline'ów bioinformatycznych, od prototypów badawczych obsługujących kilka próbek do systemów produkcyjnych przetwarzających tysiące próbek w laboratoriach i instytucjach.
Każde wyzwanie, które rozwiązałeś, i wzorzec, którego się nauczyłeś, odzwierciedla rzeczywiste problemy, z którymi borykają się deweloperzy podczas skalowania workflow Nextflow.

Zastosowanie tych wzorców w Twojej własnej pracy umożliwi Ci budowanie solidnych, gotowych do produkcji workflow.

### Kluczowe wzorce

1.  **Przepływ Danych vs Skryptowanie:** Nauczyłeś się rozróżniać między operacjami przepływu danych (orkiestra kanałów) a skryptowaniem (kod, który manipuluje danymi), w tym kluczowe różnice między operacjami na różnych typach, jak `collect` na Channel vs List.

    - Przepływ danych: orkiestra kanałów

    ````groovy
    channel.fromPath('*.fast```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ````

    - Skryptowanie: przetwarzanie danych w kolekcjach

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Zaawansowane Przetwarzanie Ciągów Znaków**: Opanowałeś wyrażenia regularne do parsowania nazw plików, dynamiczne generowanie skryptów w procesach i interpolację zmiennych (Nextflow vs Bash vs Shell).

    - Dopasowywanie wzorców

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Funkcja z warunkowym zwracaniem

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Kolekcja plików do argumentów poleceń (w bloku skryptowym procesu)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Tworzenie Funkcji Wielokrotnego Użytku**: Nauczyłeś się wyodrębniać złożoną logikę do nazwanych funkcji, które mogą być wywoływane z operatorów kanałów, czyniąc workflow bardziej czytelnymi i łatwiejszymi w utrzymaniu.

    - Zdefiniuj nazwaną funkcję

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* kod ukryty dla zwięzłości */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* kod ukryty dla zwięzłości */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Wywołaj nazwaną funkcję w workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Dynamiczne Dyrektywy Zasobów z Closures**: Zbadałeś używanie closures w dyrektywach procesów do adaptacyjnej alokacji zasobów na podstawie charakterystyki wejścia.

    - Nazwane closures i kompozycja

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures z dostępem do zakresu

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logika Warunkowa i Kontrola Procesów**: Dodałeś inteligentne kierowanie za pomocą operatorów `.branch()` i `.filter()`, wykorzystując prawdziwość do zwięzłych wyrażeń warunkowych.

    - Użyj `.branch()` do kierowania danych przez różne gałęzie workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Ewaluacja boolowska z Groovy Truth

    ```groovy
    if (sample.files) println "Ma pliki"
    ```

    - Użyj `filter()` do wyodrębniania podzbiorów danych z 'prawdziwością'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Bezpieczna Nawigacja i Operatory Elvis**: Uczyniłeś pipeline odpornym na brakujące dane, używając `?.` do bezpiecznego dostępu do właściwości względem null i `?:` do zapewniania wartości domyślnych.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Walidacja z error() i log.warn**: Nauczyłeś się walidować wejścia wcześnie i szybko zawodzić z wyraźnymi komunikatami błędów.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Nieprawidłowe: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Błąd: ${e.message}"
    }
    ```

8.  **Handlery Zdarzeń Konfiguracji**: Nauczyłeś się używać handlerów zdarzeń workflow (`onComplete` i `onError`) do logowania, powiadomień i zarządzania cyklem życia.

    - Używanie `onComplete` do logowania i powiadamiania

    ```groovy
    workflow.onComplete = {
        println "Sukces     : ${workflow.success}"
        println "status wyjścia : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline zakończony pomyślnie!"
        } else {
            println "❌ Pipeline zakończony niepowodzeniem!"
            println "Błąd: ${workflow.errorMessage}"
        }
    }
    ```

    - Używanie `onError` do podejmowania działań specyficznie w przypadku niepowodzenia

    ```groovy
    workflow.onError = {
        // Zapisz szczegółowy log błędu
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Czas: ${new Date()}
        Błąd: ${workflow.errorMessage}
        Raport błędu: ${workflow.errorReport ?: 'Brak szczegółowego raportu dostępnego'}
        """

        println "Szczegóły błędu zapisane do: ${error_file}"
    }
    ```

### Dodatkowe zasoby

- [Nextflow Language Reference](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Syntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standard Library](https://nextflow.io/docs/latest/reference/stdlib.html)

Upewnij się, że sprawdzisz te zasoby, gdy musisz zbadać bardziej zaawansowane funkcje.

Skorzystasz na praktykowaniu i rozwijaniu swoich umiejętności, aby:

- Pisać czystsze workflow z odpowiednim rozdzieleniem między przepływem danych a skryptowaniem
- Opanować interpolację zmiennych, aby uniknąć powszechnych pułapek ze zmiennymi Nextflow, Bash i shell
- Używać dynamicznych dyrektyw zasobów dla wydajnych, adaptacyjnych workflow
- Przekształcać kolekcje plików w poprawnie sformatowane argumenty linii poleceń
- Obsługiwać różne konwencje nazewnictwa plików i formaty wejściowe z wdziękiem, używając regex i przetwarzania ciągów znaków
- Budować wielokrotnego użytku, łatwy w utrzymaniu kod, używając zaawansowanych wzorców closures i programowania funkcyjnego
- Przetwarzać i organizować złożone zbiory danych, używając operacji na kolekcjach
- Dodawać walidację, obsługę błędów i logowanie, aby uczynić swoje workflow gotowymi do produkcji
- Implementować zarządzanie cyklem życia workflow za pomocą handlerów zdarzeń

---

## Co dalej?

Wróć do [menu Side Questów](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
