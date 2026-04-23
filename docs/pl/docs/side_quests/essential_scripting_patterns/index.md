# Podstawowe Wzorce Skryptowania w Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow to język programowania działający na wirtualnej maszynie Java (JVM). Choć Nextflow jest zbudowany na [Groovy](http://groovy-lang.org/) i dzieli z nim dużą część składni, nie jest jedynie „Groovy z rozszerzeniami" — to samodzielny język z w pełni zdefiniowaną [składnią](https://nextflow.io/docs/latest/reference/syntax.html) i [biblioteką standardową](https://nextflow.io/docs/latest/reference/stdlib.html).

Wiele workflow'ów można napisać, nie wychodząc poza podstawową składnię dla zmiennych, map i list. Większość samouczków Nextflow skupia się na orkiestracji workflow'u (kanały, procesy i przepływ danych), a samo to pozwala zajść zaskakująco daleko.

Jednak gdy trzeba manipulować danymi, parsować złożone nazwy plików, implementować logikę warunkową lub budować solidne produkcyjne workflow'y, warto myśleć o dwóch odrębnych aspektach kodu: **dataflow** (kanały, operatory, procesy i workflow'y) oraz **skryptowaniu** (kod wewnątrz domknięć, funkcji i skryptów procesów). To rozróżnienie jest nieco arbitralne — to wciąż ten sam kod Nextflow — ale stanowi użyteczny model mentalny pomagający zrozumieć, kiedy orkiestrujesz pipeline, a kiedy manipulujesz danymi. Opanowanie obu aspektów znacząco poprawia zdolność pisania przejrzystych i łatwych w utrzymaniu workflow'ów.

### Cele szkolenia

Ten side quest zabierze Cię w praktyczną podróż od podstawowych konceptów do wzorców gotowych na środowisko produkcyjne.
Przekształcimy prosty workflow odczytujący CSV w zaawansowany pipeline bioinformatyczny, rozwijając go krok po kroku przez realistyczne wyzwania:

- **Rozumienie granic:** Odróżnianie operacji dataflow od skryptowania i rozumienie, jak ze sobą współpracują
- **Manipulacja danymi:** Wyodrębnianie, transformowanie i podzbiory map i kolekcji przy użyciu potężnych operatorów
- **Przetwarzanie stringów:** Parsowanie złożonych schematów nazewnictwa plików za pomocą wzorców regex i opanowanie interpolacji zmiennych
- **Funkcje wielokrotnego użytku:** Wyodrębnianie złożonej logiki do nazwanych funkcji dla czystszych i łatwiejszych w utrzymaniu workflow'ów
- **Logika dynamiczna:** Budowanie procesów adaptujących się do różnych typów wejść i używanie domknięć do dynamicznej alokacji zasobów
- **Warunkowe trasowanie:** Inteligentne kierowanie próbek przez różne procesy na podstawie ich cech metadanych
- **Bezpieczne operacje:** Obsługa brakujących danych z operatorami null-safe i walidacja wejść z czytelnymi komunikatami błędów
- **Handlery oparte na konfiguracji:** Używanie handlerów zdarzeń workflow do logowania, powiadomień i zarządzania cyklem życia

### Wymagania wstępne

Przed podjęciem tego side questu powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow (procesy, kanały, operatory, praca z plikami, metadane)
- Posiadać podstawową znajomość typowych konstrukcji programistycznych (zmienne, mapy, listy)

Samouczek będzie wyjaśniał koncepty programistyczne w miarę ich napotkania, więc nie potrzebujesz rozległego doświadczenia w programowaniu.
Zaczniemy od podstawowych konceptów i będziemy stopniowo przechodzić do zaawansowanych wzorców.

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/essential_scripting_patterns
```

#### Przejrzyj materiały

Znajdziesz tam główny plik workflow oraz katalog `data` zawierający przykładowe pliki danych.

```console title="Directory contents"
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

Nasz przykładowy plik CSV zawiera informacje o próbkach biologicznych, które wymagają różnego przetwarzania w zależności od ich cech:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Użyjemy tego realistycznego zestawu danych do eksploracji praktycznych technik programistycznych, które napotkasz w prawdziwych workflow'ach bioinformatycznych.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiedni katalog roboczy
<!-- - [ ] I understand the assignment -->

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Dataflow a skryptowanie: rozumienie granic

### 1.1. Identyfikowanie, co jest czym

Pisząc workflow'y w Nextflow, ważne jest rozróżnienie między **dataflow** (jak dane przemieszczają się przez kanały i procesy) a **skryptowaniem** (kod manipulujący danymi i podejmujący decyzje). Zbudujmy workflow demonstrujący, jak ze sobą współpracują.

#### 1.1.1. Podstawowy workflow Nextflow

Zacznij od prostego workflow, który tylko odczytuje plik CSV (zrobiliśmy to już za Ciebie w `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Blok `workflow` definiuje strukturę naszego pipeline'u, a `channel.fromPath()` tworzy kanał ze ścieżki do pliku. Operator `.splitCsv()` przetwarza plik CSV i konwertuje każdy wiersz na strukturę danych map.

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

#### 1.1.2. Dodawanie operatora map

Teraz dodamy skryptowanie do transformacji danych, używając operatora `.map()`, który prawdopodobnie już znasz. Operator ten przyjmuje „domknięcie", w którym możemy napisać kod transformujący każdy element.

!!! note "Uwaga"

    **Domknięcie** to blok kodu, który można przekazywać i wykonywać później. Pomyśl o nim jak o funkcji definiowanej w miejscu użycia. Domknięcia zapisuje się w nawiasach klamrowych `{ }` i mogą przyjmować parametry. Są fundamentem działania operatorów Nextflow — jeśli piszesz w Nextflow od jakiegoś czasu, być może już ich używałeś, nie zdając sobie z tego sprawy!

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

To nasze pierwsze **domknięcie** — anonimowa funkcja, którą można przekazać jako argument (podobna do lambd w Pythonie lub funkcji strzałkowych w JavaScript). Domknięcia są niezbędne do pracy z operatorami Nextflow.

Domknięcie `{ row -> return row }` przyjmuje parametr `row` (może mieć dowolną nazwę: `item`, `sample` itp.).

Gdy operator `.map()` przetwarza każdy element kanału, przekazuje go do Twojego domknięcia. Tutaj `row` przechowuje jeden wiersz CSV na raz.

Zastosuj tę zmianę i uruchom workflow:

```bash
nextflow run main.nf
```

Zobaczysz taki sam wynik jak poprzednio, ponieważ po prostu zwracamy wejście bez zmian. Potwierdza to, że operator map działa poprawnie. Teraz zacznijmy transformować dane.

#### 1.1.3. Tworzenie struktury danych map

Teraz napiszemy logikę **skryptowania** wewnątrz naszego domknięcia, aby transformować każdy wiersz danych. Tu przetwarzamy poszczególne elementy danych, a nie orkiestrujemy przepływ danych.

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

Mapa `sample_meta` to struktura danych klucz-wartość (jak słowniki w Pythonie, obiekty w JavaScript lub hashe w Ruby) przechowująca powiązane informacje: identyfikator próbki, organizm, typ tkanki, głębokość sekwencjonowania i wynik jakości.

Używamy metod manipulacji stringami, takich jak `.toLowerCase()` i `.replaceAll()`, do oczyszczania danych, oraz metod konwersji typów, takich jak `.toInteger()` i `.toDouble()`, do konwersji danych tekstowych z CSV na odpowiednie typy numeryczne.

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

#### 1.1.4. Dodawanie logiki warunkowej

Dodajmy teraz więcej skryptowania — tym razem używając operatora trójargumentowego do podejmowania decyzji na podstawie wartości danych.

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

Operator trójargumentowy to skrócony zapis instrukcji if/else, który ma postać `warunek ? wartość_jeśli_prawda : wartość_jeśli_fałsz`. Ten wiersz oznacza: „Jeśli jakość jest większa niż 40, użyj 'high', w przeciwnym razie użyj 'normal'". Jego kuzyn, **operator Elvis** (`?:`), dostarcza wartości domyślnych, gdy coś jest null lub puste — ten wzorzec poznamy później w tym samouczku.

Operator dodawania map `+` tworzy **nową mapę** zamiast modyfikować istniejącą. Ten wiersz tworzy nową mapę zawierającą wszystkie pary klucz-wartość z `sample_meta` plus nowy klucz `priority`.

!!! Note "Uwaga"

    Nigdy nie modyfikuj map przekazywanych do domknięć — zawsze twórz nowe, używając np. `+`. W Nextflow te same dane często przepływają przez wiele operacji jednocześnie. Modyfikacja mapy w miejscu może powodować nieprzewidywalne efekty uboczne, gdy inne operacje odwołują się do tego samego obiektu. Tworzenie nowych map zapewnia, że każda operacja ma własną czystą kopię.

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

Pomyślnie dodaliśmy logikę warunkową, wzbogacając nasze metadane o poziom priorytetu oparty na wynikach jakości.

#### 1.1.5. Podzbiory map za pomocą `.subMap()`

Podczas gdy operator `+` dodaje klucze do mapy, czasem potrzebujemy zrobić coś odwrotnego — wyodrębnić tylko określone klucze. Metoda `.subMap()` jest do tego idealna.

Dodajmy wiersz tworzący uproszczoną wersję naszych metadanych zawierającą tylko pola identyfikacyjne:

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
                println "ID fields only: ${id_only}"

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

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Widać tu zarówno pełne metadane wyświetlane przez operację `view()`, jak i wyodrębniony podzbiór wydrukowany przez `println`.

Metoda `.subMap()` przyjmuje listę kluczy i zwraca nową mapę zawierającą tylko te klucze. Jeśli klucz nie istnieje w oryginalnej mapie, po prostu nie jest uwzględniany w wyniku.

Jest to szczególnie przydatne, gdy trzeba tworzyć różne wersje metadanych dla różnych procesów — niektóre mogą potrzebować pełnych metadanych, podczas gdy inne wymagają tylko minimalnych pól identyfikacyjnych.

Teraz usuń te instrukcje println, aby przywrócić workflow do poprzedniego stanu, ponieważ nie będziemy ich już potrzebować.

!!! tip "Wskazówka: Podsumowanie operacji na mapach"

    - **Dodawanie kluczy**: `map1 + [new_key: value]` — tworzy nową mapę z dodatkowymi kluczami
    - **Wyodrębnianie kluczy**: `map1.subMap(['key1', 'key2'])` — tworzy nową mapę tylko z określonymi kluczami
    - **Obie operacje tworzą nowe mapy** — oryginalne mapy pozostają niezmienione

#### 1.1.6. Łączenie map i zwracanie wyników

Do tej pory zwracaliśmy tylko to, co społeczność Nextflow nazywa „meta mapą", ignorując pliki, których te metadane dotyczą. Ale jeśli piszesz workflow'y w Nextflow, prawdopodobnie chcesz coś z tymi plikami zrobić.

Wyprowadźmy strukturę kanału składającą się z krotki 2 elementów: wzbogaconej mapy metadanych i odpowiadającej jej ścieżki do pliku. To powszechny wzorzec w Nextflow do przekazywania danych do procesów.

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

    **Mapy i metadane**: Mapy są fundamentem pracy z metadanymi w Nextflow. Bardziej szczegółowe wyjaśnienie pracy z mapami metadanych znajdziesz w side queście [Working with metadata](../metadata/).

Nasz workflow demonstruje podstawowy wzorzec: **operacje dataflow** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orkiestrują przepływ danych przez pipeline, podczas gdy **skryptowanie** (mapy `[key: value]`, metody stringów, konwersje typów, operatory trójargumentowe) wewnątrz domknięcia `.map()` obsługuje transformację poszczególnych elementów danych.

### 1.2. Rozumienie różnych typów: kanał a lista

Do tej pory dobrze — potrafimy odróżnić operacje dataflow od skryptowania. Ale co, gdy ta sama nazwa metody istnieje w obu kontekstach?

Doskonałym przykładem jest metoda `collect`, która istnieje zarówno dla typów kanałów, jak i typów List w bibliotece standardowej Nextflow. Metoda `collect()` na liście transformuje każdy element, podczas gdy operator `collect()` na kanale zbiera wszystkie emisje kanału w jeden element.

Zademonstrujmy to na przykładowych danych, zaczynając od przypomnienia sobie, co robi operator `collect()` na kanale. Sprawdź `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - grupuje wiele emisji kanału w jedną
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Kroki:

- Zdefiniuj listę identyfikatorów próbek
- Utwórz kanał za pomocą `fromList()`, który emituje każdy identyfikator próbki osobno
- Wydrukuj każdy element za pomocą `view()` w miarę przepływu
- Zbierz wszystkie elementy w jedną listę za pomocą operatora `collect()` kanału
- Wydrukuj zebrany wynik (pojedynczy element zawierający wszystkie identyfikatory próbek) za pomocą drugiego `view()`

Zmieniliśmy strukturę kanału, ale nie zmieniliśmy samych danych.

Uruchom workflow, aby to potwierdzić:

```bash
nextflow run collect.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` zwraca wynik dla każdej emisji kanału, więc wiemy, że ten pojedynczy wynik zawiera wszystkie 3 oryginalne elementy zgrupowane w jedną listę.

Teraz zobaczmy metodę `collect` na liście w akcji. Zmodyfikuj `collect.nf`, aby zastosować metodę `collect` listy do oryginalnej listy identyfikatorów próbek:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transformuje każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

W tym nowym fragmencie:

- Definiujemy nową zmienną `formatted_ids`, która używa metody `collect` listy do transformacji każdego identyfikatora próbki na oryginalnej liście
- Drukujemy wynik za pomocą `println`

Uruchom zmodyfikowany workflow:

```bash
nextflow run collect.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Tym razem NIE zmieniliśmy struktury danych — nadal mamy 3 elementy na liście — ale TRANSFORMOWALIŚMY każdy element za pomocą metody `collect` listy, tworząc nową listę ze zmodyfikowanymi wartościami. Jest to podobne do używania operatora `map` na kanale, ale operuje na strukturze danych List, a nie na kanale.

`collect` to skrajny przypadek, którego używamy tu, aby zilustrować pewien punkt. Kluczowa lekcja jest taka: pisząc workflow'y, zawsze odróżniaj **struktury danych** (listy, mapy itp.) od **kanałów** (konstrukcji dataflow). Operacje mogą mieć te same nazwy, ale zachowywać się zupełnie inaczej w zależności od typu, na którym są wywoływane.

### 1.3. Operator spread (`*.`) — skrótowy zapis do wyodrębniania właściwości

Powiązany z metodą `collect` listy jest operator spread (`*.`), który zapewnia zwięzły sposób wyodrębniania właściwości z kolekcji. Jest to w zasadzie cukier składniowy dla powszechnego wzorca `collect`.

Dodajmy demonstrację do naszego pliku `collect.nf`:

=== "Po"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transformuje każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Operator spread - zwięzły dostęp do właściwości
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Przed"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transformuje każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Uruchom zaktualizowany workflow:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Operator spread `*.` to skrótowy zapis powszechnego wzorca collect:

```groovy
// Te zapisy są równoważne:
def ids = samples*.id
def ids = samples.collect { it.id }

// Działa też z wywołaniami metod:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Operator spread jest szczególnie przydatny, gdy trzeba wyodrębnić jedną właściwość z listy obiektów — jest bardziej czytelny niż pisanie pełnego domknięcia `collect`.

!!! tip "Wskazówka: Kiedy używać spread zamiast collect"

    - **Używaj spread (`*.`)** do prostego dostępu do właściwości: `samples*.id`, `files*.name`
    - **Używaj collect** do transformacji lub złożonej logiki: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Dataflow a skryptowanie**: Operatory kanałów orkiestrują przepływ danych przez pipeline, podczas gdy skryptowanie transformuje poszczególne elementy danych
- **Rozumienie typów**: Ta sama nazwa metody (jak `collect`) może zachowywać się inaczej w zależności od typu, na którym jest wywoływana (kanał vs lista)
- **Kontekst ma znaczenie**: Zawsze miej świadomość, czy pracujesz z kanałami (dataflow) czy strukturami danych (skryptowanie)

Rozumienie tych granic jest niezbędne do debugowania, dokumentowania i pisania łatwych w utrzymaniu workflow'ów.

Następnie zagłębimy się w możliwości przetwarzania stringów, które są niezbędne do obsługi danych z prawdziwego świata.

---

## 2. Przetwarzanie stringów i dynamiczne generowanie skryptów

Opanowanie przetwarzania stringów odróżnia kruche workflow'y od solidnych pipeline'ów. Ta sekcja omawia parsowanie złożonych nazw plików, dynamiczne generowanie skryptów i interpolację zmiennych.

### 2.1. Dopasowywanie wzorców i wyrażenia regularne

Pliki bioinformatyczne często mają złożone konwencje nazewnictwa kodujące metadane. Wyodrębnijmy je automatycznie za pomocą dopasowywania wzorców z wyrażeniami regularnymi.

Wrócimy do naszego workflow `main.nf` i dodamy logikę dopasowywania wzorców, aby wyodrębnić dodatkowe informacje o próbce z nazw plików. Pliki FASTQ w naszym zestawie danych stosują konwencje nazewnictwa w stylu Illumina, z nazwami takimi jak `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Mogą wyglądać tajemniczo, ale faktycznie kodują przydatne metadane, takie jak identyfikator próbki, numer ścieżki i kierunek odczytu. Użyjemy możliwości regex do parsowania tych nazw.

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

Demonstruje to kluczowe **koncepty przetwarzania stringów**:

1. **Literały wyrażeń regularnych** używające składni `~/wzorzec/` — tworzy wzorzec regex bez potrzeby escapowania ukośników odwrotnych
2. **Dopasowywanie wzorców** za pomocą operatora `=~` — próbuje dopasować string do wzorca regex
3. **Obiekty Matcher** przechwytujące grupy za pomocą `[0][1]`, `[0][2]` itp. — `[0]` odnosi się do całego dopasowania, `[1]`, `[2]` itp. do grup przechwyconych w nawiasach

Rozłóżmy wzorzec regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Wzorzec             | Dopasowuje                                  | Przechwytuje                                       |
| ------------------- | ------------------------------------------- | -------------------------------------------------- |
| `^(.+)`             | Nazwa próbki od początku                    | Grupa 1: nazwa próbki                              |
| `_S(\d+)`           | Numer próbki `_S1`, `_S2` itp.              | Grupa 2: numer próbki                              |
| `_L(\d{3})`         | Numer ścieżki `_L001`                       | Grupa 3: ścieżka (3 cyfry)                         |
| `_(R[12])`          | Kierunek odczytu `_R1` lub `_R2`            | Grupa 4: kierunek odczytu                          |
| `_(\d{3})`          | Numer fragmentu `_001`                      | Grupa 5: fragment (3 cyfry)                        |
| `\.fastq(?:\.gz)?$` | Rozszerzenie pliku `.fastq` lub `.fastq.gz` | Nie przechwycone (`?:` to grupa nieprzechwytująca) |

Parsuje to konwencje nazewnictwa w stylu Illumina, automatycznie wyodrębniając metadane.

Uruchom zmodyfikowany workflow:

```bash title="Test pattern matching"
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

Widać metadane wzbogacone o informacje z nazw plików.

### 2.2. Dynamiczne generowanie skryptów w procesach

Bloki skryptów procesów to w zasadzie wieloliniowe stringi przekazywane do powłoki. Możesz używać **logiki warunkowej** (if/else, operatory trójargumentowe) do dynamicznego generowania różnych stringów skryptów na podstawie cech wejścia. Jest to niezbędne do obsługi różnorodnych typów wejść — takich jak odczyty single-end vs paired-end — bez duplikowania definicji procesów.

Dodajmy do naszego workflow proces demonstrujący ten wzorzec. Otwórz `modules/fastp.nf` i przyjrzyj się:

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

Proces przyjmuje pliki FASTQ jako wejście i uruchamia narzędzie `fastp` do przycinania adapterów i filtrowania odczytów niskiej jakości. Niestety, osoba, która napisała ten proces, nie uwzględniła odczytów single-end, które mamy w naszym przykładowym zestawie danych. Dodajmy go do naszego workflow i zobaczmy, co się stanie:

Najpierw dołącz moduł w pierwszym wierszu pliku `main.nf`:

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

Widać, że proces próbuje uruchomić `fastp` z wartością `null` dla drugiego pliku wejściowego, co powoduje błąd. Dzieje się tak, ponieważ nasz zestaw danych zawiera odczyty single-end, ale proces jest zakodowany na stałe do obsługi odczytów paired-end (dwa pliki wejściowe jednocześnie).

Napraw to, dodając logikę warunkową do bloku `script:` procesu `FASTP`. Instrukcja if/else sprawdza liczbę plików odczytu i odpowiednio dostosowuje polecenie.

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

Teraz workflow może obsługiwać zarówno odczyty single-end, jak i paired-end. Logika warunkowa sprawdza liczbę plików wejściowych i konstruuje odpowiednie polecenie dla `fastp`. Sprawdźmy, czy działa:

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

Wygląda dobrze! Jeśli sprawdzimy faktyczne wykonane polecenia (dostosuj do swojego hasha zadania):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Widzimy, że Nextflow poprawnie wybrał właściwe polecenie dla odczytów single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Inny przykład dynamicznej logiki skryptów można znaleźć w [module Genomics kursu Nextflow for Science](../../nf4science/genomics/02_joint_calling). W tym module wywoływany proces GATK może przyjmować wiele plików wejściowych, ale każdy musi być poprzedzony `-V`, aby utworzyć poprawną linię poleceń. Proces używa skryptowania do transformacji kolekcji plików wejściowych (`all_gvcfs`) w poprawne argumenty polecenia:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Te wzorce używania skryptowania w blokach skryptów procesów są niezwykle potężne i można je stosować w wielu scenariuszach — od obsługi zmiennych typów wejść po budowanie złożonych argumentów wiersza poleceń z kolekcji plików, czyniąc Twoje procesy naprawdę adaptowalnymi do różnorodnych wymagań danych z prawdziwego świata.

### 2.3. Interpolacja zmiennych: Nextflow i zmienne powłoki

Skrypty procesów mieszają zmienne Nextflow, zmienne powłoki i podstawienia poleceń, każde z inną składnią interpolacji. Użycie złej składni powoduje błędy. Zbadajmy to na przykładzie procesu tworzącego raport przetwarzania.

Przyjrzyj się plikowi modułu `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Ten proces zapisuje prosty raport z identyfikatorem próbki i nazwą pliku. Teraz uruchommy go, aby zobaczyć, co się dzieje, gdy musimy mieszać różne typy zmiennych.

Dołącz proces do `main.nf` i dodaj go do workflow:

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

Ale co, jeśli chcemy dodać informacje o tym, kiedy i gdzie odbyło się przetwarzanie? Zmodyfikujmy proces, aby używał zmiennych **powłoki** i podstawień poleceń, aby uwzględnić w raporcie bieżącego użytkownika, nazwę hosta i datę:

=== "Po"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Przed"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Jeśli to uruchomisz, zauważysz błąd — Nextflow próbuje zinterpretować `${USER}` jako zmienną Nextflow, która nie istnieje.

??? failure "Wyjście polecenia"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Musimy to escapować, aby Bash mógł to obsłużyć.

Napraw to, escapując zmienne powłoki i podstawienia poleceń ukośnikiem odwrotnym (`\`):

=== "Po"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Przed"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

Teraz działa! Ukośnik odwrotny (`\`) mówi Nextflow'owi: „nie interpretuj tego, przekaż to do Basha".

### Podsumowanie

W tej sekcji nauczyłeś się technik **przetwarzania stringów**:

- **Wyrażenia regularne do parsowania plików**: Używanie operatora `=~` i wzorców regex (`~/wzorzec/`) do wyodrębniania metadanych ze złożonych konwencji nazewnictwa plików
- **Dynamiczne generowanie skryptów**: Używanie logiki warunkowej (if/else, operatory trójargumentowe) do generowania różnych stringów skryptów na podstawie cech wejścia
- **Interpolacja zmiennych**: Rozumienie, kiedy Nextflow interpretuje stringi, a kiedy robi to powłoka
  - `${var}` — zmienne Nextflow (interpolowane przez Nextflow w czasie kompilacji workflow'u)
  - `\${var}` — zmienne środowiskowe powłoki (escapowane, przekazywane do basha w czasie wykonania)
  - `\$(cmd)` — podstawienie poleceń powłoki (escapowane, wykonywane przez basha w czasie wykonania)

Te wzorce przetwarzania i generowania stringów są niezbędne do obsługi różnorodnych formatów plików i konwencji nazewnictwa, które napotkasz w prawdziwych workflow'ach bioinformatycznych.

---

## 3. Tworzenie funkcji wielokrotnego użytku

Złożona logika workflow wbudowana w operatory kanałów lub definicje procesów zmniejsza czytelność i łatwość utrzymania. **Funkcje** pozwalają wyodrębnić tę logikę do nazwanych, wielokrotnego użytku komponentów.

Nasza operacja map stała się długa i złożona. Wyodrębnijmy ją do funkcji wielokrotnego użytku za pomocą słowa kluczowego `def`.

Aby zilustrować, jak to wygląda w naszym istniejącym workflow, wprowadź poniższą modyfikację, używając `def` do zdefiniowania funkcji wielokrotnego użytku o nazwie `separateMetadata`:

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

Wyodrębniając tę logikę do funkcji, sprowadziliśmy właściwą logikę workflow do czegoś znacznie czystszego:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Dzięki temu logika workflow jest znacznie łatwiejsza do odczytania i zrozumienia na pierwszy rzut oka. Funkcja `separateMetadata` hermetyzuje całą złożoną logikę parsowania i wzbogacania metadanych, czyniąc ją wielokrotnego użytku i testowalną.

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

Wynik powinien pokazywać oba procesy zakończone pomyślnie. Workflow jest teraz znacznie czystszy i łatwiejszy w utrzymaniu, a cała złożona logika przetwarzania metadanych jest hermetyzowana w funkcji `separateMetadata`.

### Podsumowanie

W tej sekcji nauczyłeś się **tworzenia funkcji**:

- **Definiowanie funkcji za pomocą `def`**: Słowo kluczowe do tworzenia nazwanych funkcji (jak `def` w Pythonie lub `function` w JavaScript)
- **Zakres funkcji**: Funkcje zdefiniowane na poziomie skryptu są dostępne w całym workflow Nextflow
- **Wartości zwracane**: Funkcje automatycznie zwracają ostatnie wyrażenie lub używają jawnego `return`
- **Czystszy kod**: Wyodrębnianie złożonej logiki do funkcji to fundamentalna praktyka inżynierii oprogramowania w każdym języku

Następnie zbadamy, jak używać domknięć w dyrektywach procesów do dynamicznej alokacji zasobów.

---

## 4. Dynamiczne dyrektywy zasobów z domknięciami

Do tej pory używaliśmy skryptowania w bloku `script` procesów. Ale **domknięcia** (wprowadzone w sekcji 1.1) są również niezwykle przydatne w dyrektywach procesów, szczególnie do dynamicznej alokacji zasobów. Dodajmy dyrektywy zasobów do naszego procesu FASTP, które adaptują się na podstawie cech próbki.

### 4.1. Alokacja zasobów specyficzna dla próbki

Obecnie nasz proces FASTP używa domyślnych zasobów. Uczyńmy go mądrzejszym, przydzielając więcej procesorów dla próbek o wysokiej głębokości. Edytuj `modules/fastp.nf`, aby uwzględnić dynamiczną dyrektywę `cpus` i statyczną dyrektywę `memory`:

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

Domknięcie `{ meta.depth > 40000000 ? 2 : 1 }` używa **operatora trójargumentowego** (omówionego w sekcji 1.1) i jest ewaluowane dla każdego zadania, umożliwiając alokację zasobów per próbka. Próbki o wysokiej głębokości (>40M odczytów) otrzymują 2 procesory, podczas gdy pozostałe otrzymują 1.

!!! note "Uwaga: Dostęp do zmiennych wejściowych w dyrektywach"

    Domknięcie może uzyskać dostęp do dowolnych zmiennych wejściowych (jak tutaj `meta`), ponieważ Nextflow ewaluuje te domknięcia w kontekście każdego wykonania zadania.

Uruchom workflow ponownie z opcją `-ansi-log false`, aby łatwiej zobaczyć hashe zadań.

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

Możesz sprawdzić dokładne polecenie `docker`, które zostało uruchomione, aby zobaczyć alokację procesorów dla dowolnego zadania:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Powinieneś zobaczyć coś takiego:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

W tym przykładzie wybraliśmy zadanie, które zażądało 2 procesorów (`--cpu-shares 2048`), ponieważ była to próbka o wysokiej głębokości, ale powinieneś zobaczyć różne alokacje procesorów w zależności od głębokości próbki. Wypróbuj to również dla innych zadań.

### 4.2. Strategie ponawiania

Innym potężnym wzorcem jest używanie `task.attempt` do strategii ponawiania. Aby pokazać, dlaczego jest to przydatne, zaczniemy od zmniejszenia alokacji pamięci dla FASTP do mniej niż potrzebuje. Zmień dyrektywę `memory` w `modules/fastp.nf` na `1.GB`:

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

Wskazuje to, że proces został zabity za przekroczenie limitów pamięci.

To bardzo powszechny scenariusz w prawdziwych workflow'ach — czasem po prostu nie wiesz, ile pamięci będzie potrzebować zadanie, dopóki go nie uruchomisz.

Aby uczynić nasz workflow bardziej solidnym, możemy zaimplementować strategię ponawiania, która zwiększa alokację pamięci przy każdej próbie, ponownie używając domknięcia Groovy. Zmodyfikuj dyrektywę `memory`, aby mnożyła bazową pamięć przez `task.attempt`, i dodaj dyrektywy `errorStrategy 'retry'` i `maxRetries 2`:

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

Teraz, jeśli proces zakończy się niepowodzeniem z powodu niewystarczającej pamięci, Nextflow ponowi próbę z większą ilością pamięci:

- Pierwsza próba: 1 GB (task.attempt = 1)
- Druga próba: 2 GB (task.attempt = 2)

... i tak dalej, aż do limitu `maxRetries`.

### Podsumowanie

Dynamiczne dyrektywy z domknięciami pozwalają:

- Alokować zasoby na podstawie cech wejścia
- Implementować automatyczne strategie ponawiania ze zwiększającymi się zasobami
- Łączyć wiele czynników (metadane, numer próby, priorytety)
- Używać logiki warunkowej do złożonych obliczeń zasobów

Dzięki temu Twoje workflow'y są zarówno bardziej wydajne (bez nadmiernej alokacji), jak i bardziej solidne (automatyczne ponawianie z większymi zasobami).

---

## 5. Logika warunkowa i kontrola procesów

Wcześniej używaliśmy `.map()` ze skryptowaniem do transformacji danych kanałów. Teraz użyjemy logiki warunkowej do kontrolowania, które procesy są wykonywane na podstawie danych — co jest niezbędne dla elastycznych workflow'ów adaptujących się do różnych typów próbek.

[Operatory dataflow](https://www.nextflow.io/docs/latest/reference/operator.html) Nextflow przyjmują domknięcia ewaluowane w czasie wykonania, umożliwiając logice warunkowej sterowanie decyzjami workflow na podstawie zawartości kanałów.

### 5.1. Trasowanie za pomocą `.branch()`

Na przykład, udawajmy, że nasze próbki sekwencjonowania muszą być przycinane za pomocą FASTP tylko wtedy, gdy są to próbki ludzkie z pokryciem powyżej określonego progu. Próbki mysie lub o niskim pokryciu powinny być uruchamiane z Trimgalore (to wymyślony przykład, ale ilustruje punkt).

Dostarczyliśmy prosty proces Trimgalore w `modules/trimgalore.nf` — możesz go przejrzeć, ale szczegóły nie są ważne dla tego ćwiczenia. Kluczowe jest to, że chcemy trasować próbki na podstawie ich metadanych.

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

... a następnie zmodyfikuj blok `workflow` w `main.nf`, aby rozgałęziać próbki na podstawie ich metadanych i kierować je przez odpowiedni proces przycinania:

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

Użyliśmy tu małych, ale potężnych wyrażeń warunkowych wewnątrz operatora `.branch{}` do trasowania próbek na podstawie ich metadanych. Próbki ludzkie o wysokim pokryciu przechodzą przez `FASTP`, podczas gdy wszystkie pozostałe przechodzą przez `TRIMGALORE`.

### 5.2. Używanie `.filter()` z prawdziwością

Innym potężnym wzorcem do kontrolowania wykonania workflow jest operator `.filter()`, który używa domknięcia do określenia, które elementy powinny kontynuować przepływ przez pipeline. Wewnątrz domknięcia filter piszesz **wyrażenia boolowskie** decydujące, które elementy przechodzą.

Nextflow (jak wiele dynamicznych języków) ma koncepcję **„prawdziwości"** (truthiness), która określa, jakie wartości ewaluują się do `true` lub `false` w kontekstach boolowskich:

- **Prawdziwe**: Wartości inne niż null, niepuste stringi, liczby różne od zera, niepuste kolekcje
- **Fałszywe**: `null`, puste stringi `""`, zero `0`, puste kolekcje `[]` lub `[:]`, `false`

Oznacza to, że samo `meta.id` (bez jawnego `!= null`) sprawdza, czy ID istnieje i nie jest puste. Użyjmy tego do filtrowania próbek, które nie spełniają naszych wymagań jakościowych.

Dodaj następujące przed operacją branch:

=== "Po"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtrowanie nieprawidłowych lub niskiej jakości próbek
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

Ponieważ wybraliśmy filtr wykluczający niektóre próbki, wykonano mniej zadań.

Wyrażenie filtrujące `meta.id && meta.organism && meta.depth >= 25000000` łączy prawdziwość z jawnymi porównaniami:

- `meta.id && meta.organism` sprawdza, że oba pola istnieją i nie są puste (używając prawdziwości)
- `meta.depth >= 25000000` zapewnia wystarczającą głębokość sekwencjonowania za pomocą jawnego porównania

!!! note "Uwaga: Prawdziwość w praktyce"

    Wyrażenie `meta.id && meta.organism` jest bardziej zwięzłe niż pisanie:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Dzięki temu logika filtrowania jest znacznie czystsza i łatwiejsza do odczytania.

### Podsumowanie

W tej sekcji nauczyłeś się używać logiki warunkowej do kontrolowania wykonania workflow za pomocą interfejsów domknięć operatorów Nextflow, takich jak `.branch{}` i `.filter{}`, wykorzystując prawdziwość do pisania zwięzłych wyrażeń warunkowych.

Nasz pipeline teraz inteligentnie kieruje próbki przez odpowiednie procesy, ale produkcyjne workflow'y muszą obsługiwać nieprawidłowe dane w sposób kontrolowany. Uczyńmy nasz workflow odpornym na brakujące lub null wartości.

---

## 6. Operatory bezpiecznej nawigacji i Elvis

Nasza funkcja `separateMetadata` obecnie zakłada, że wszystkie pola CSV są obecne i prawidłowe. Ale co się stanie z niekompletnymi danymi? Sprawdźmy.

### 6.1. Problem: dostęp do właściwości, które nie istnieją

Powiedzmy, że chcemy dodać obsługę opcjonalnych informacji o przebiegu sekwencjonowania. W niektórych laboratoriach próbki mogą mieć dodatkowe pole dla identyfikatora przebiegu sekwencjonowania lub numeru partii, ale nasz obecny plik CSV nie ma tej kolumny. Spróbujmy mimo to uzyskać do niej dostęp.

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

Crash z NullPointerException.

Problem polega na tym, że `row.run_id` zwraca `null`, ponieważ kolumna `run_id` nie istnieje w naszym pliku CSV. Gdy próbujemy wywołać `.toUpperCase()` na `null`, program się crashuje. Tu z pomocą przychodzi operator bezpiecznej nawigacji.

### 6.2. Operator bezpiecznej nawigacji (`?.`)

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
    <!-- TODO: output -->
    ```

Żadnego crashu! Workflow obsługuje teraz brakujące pole w sposób kontrolowany. Gdy `row.run_id` jest `null`, operator `?.` zapobiega wywołaniu `.toUpperCase()`, a `run_id` staje się `null` zamiast powodować wyjątek.

### 6.3. Operator Elvis (`?:`) dla wartości domyślnych

Operator Elvis (`?:`) dostarcza wartości domyślnych, gdy lewa strona jest „fałszywa" (jak wyjaśniono wcześniej). Nazwa pochodzi od Elvisa Presleya, ponieważ `?:` wygląda jak jego słynna fryzura i oczy oglądane z boku!

Teraz, gdy używamy bezpiecznej nawigacji, `run_id` będzie `null` dla próbek bez tego pola. Użyjmy operatora Elvis, aby dostarczyć wartość domyślną i dodać ją do naszej mapy `sample_meta`:

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

Dodaj też operator `view()` w workflow, aby zobaczyć wyniki:

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

Wszystkie próbki mają teraz pole `run` z rzeczywistym identyfikatorem przebiegu (pisanym wielkimi literami) lub wartością domyślną 'UNSPECIFIED'. Kombinacja `?.` i `?:` zapewnia zarówno bezpieczeństwo (brak crashy), jak i sensowne wartości domyślne.

Usuń teraz operator `.view()`, skoro potwierdziliśmy, że działa.

!!! tip "Wskazówka: Łączenie bezpiecznej nawigacji i Elvis"

    Wzorzec `value?.method() ?: 'default'` jest powszechny w produkcyjnych workflow'ach:

    - `value?.method()` — bezpiecznie wywołuje metodę, zwraca `null` jeśli `value` jest `null`
    - `?: 'default'` — dostarcza wartość zastępczą, jeśli wynik jest `null`

    Ten wzorzec obsługuje brakujące/niekompletne dane w sposób kontrolowany.

Używaj tych operatorów konsekwentnie w funkcjach, domknięciach operatorów (`.map{}`, `.filter{}`), skryptach procesów i plikach konfiguracyjnych. Zapobiegają crashom podczas obsługi danych z prawdziwego świata.

### Podsumowanie

- **Bezpieczna nawigacja (`?.`)**: Zapobiega crashom na wartościach null — zwraca null zamiast rzucać wyjątek
- **Operator Elvis (`?:`)**: Dostarcza wartości domyślne — `value ?: 'default'`
- **Łączenie**: `value?.method() ?: 'default'` to powszechny wzorzec

Te operatory czynią workflow'y odpornymi na niekompletne dane — co jest niezbędne w pracy z danymi z prawdziwego świata.

---

## 7. Walidacja za pomocą `error()` i `log.warn`

Czasem trzeba natychmiast zatrzymać workflow, jeśli parametry wejściowe są nieprawidłowe. W Nextflow możesz używać wbudowanych funkcji, takich jak `error()` i `log.warn`, a także standardowych konstrukcji programistycznych, takich jak instrukcje `if` i logika boolowska, do implementacji logiki walidacji. Dodajmy walidację do naszego workflow.

Utwórz funkcję walidacji przed blokiem workflow, wywołaj ją z workflow i zmień tworzenie kanału, aby używało parametru dla ścieżki do pliku CSV. Jeśli parametr brakuje lub plik nie istnieje, wywołaj `error()`, aby zatrzymać wykonanie z czytelnym komunikatem.

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Sprawdź, czy parametr wejściowy jest podany
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Sprawdź, czy plik CSV istnieje
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
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
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

Workflow zatrzymuje się natychmiast z czytelnym komunikatem błędu zamiast tajemniczo zawodzić później.

Teraz uruchom z nieistniejącym plikiem:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Na koniec uruchom z poprawnym plikiem:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Wyjście polecenia"

    ```console
    <!-- TODO: output -->
    ```

Tym razem uruchamia się pomyślnie.

Możesz też dodać walidację wewnątrz funkcji `separateMetadata`. Użyjmy nieśmiertelnego `log.warn` do wydawania ostrzeżeń dla próbek o niskiej głębokości sekwencjonowania, ale nadal pozwalając workflow kontynuować:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Walidacja danych
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
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

Uruchom workflow ponownie z oryginalnym plikiem CSV:

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
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Widzimy ostrzeżenie o niskiej głębokości sekwencjonowania dla jednej z próbek.

### Podsumowanie

- **`error()`**: Natychmiast zatrzymuje workflow z czytelnym komunikatem
- **`log.warn`**: Wydaje ostrzeżenia bez zatrzymywania workflow
- **Wczesna walidacja**: Sprawdzaj wejścia przed przetwarzaniem, aby szybko zawodzić z pomocnymi błędami
- **Funkcje walidacji**: Twórz wielokrotnego użytku logikę walidacji, którą można wywoływać na początku workflow

Właściwa walidacja czyni workflow'y bardziej solidnymi i przyjaznymi dla użytkownika, wychwytując problemy wcześnie z czytelnymi komunikatami błędów.

---

## 8. Handlery zdarzeń workflow

Do tej pory pisaliśmy kod w skryptach workflow i definicjach procesów. Jest jednak jeszcze jedna ważna funkcja, o której powinieneś wiedzieć: handlery zdarzeń workflow.

Handlery zdarzeń to domknięcia uruchamiane w określonych punktach cyklu życia Twojego workflow. Są idealne do dodawania logowania, powiadomień lub operacji czyszczenia. Powinny być zdefiniowane w skrypcie workflow obok definicji workflow.

### 8.1. Handler `onComplete`

Najczęściej używanym handlerem zdarzeń jest `onComplete`, który uruchamia się po zakończeniu workflow (niezależnie od tego, czy zakończyło się sukcesem, czy niepowodzeniem). Dodajmy jeden, aby podsumować wyniki naszego pipeline'u.

Dodaj handler zdarzeń do pliku `main.nf`, wewnątrz definicji workflow:

=== "Po"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
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

To domknięcie uruchamia się po zakończeniu workflow. Wewnątrz masz dostęp do obiektu `workflow`, który dostarcza przydatnych właściwości dotyczących wykonania.

Uruchom workflow, a na końcu zobaczysz to podsumowanie!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

Uczyńmy to bardziej użytecznym, dodając logikę warunkową:

=== "Po"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
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
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Teraz otrzymujemy jeszcze bardziej informatywne podsumowanie, w tym komunikat o sukcesie/niepowodzeniu:

<!-- TODO: add run command -->

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

Możesz też zapisać podsumowanie do pliku za pomocą operacji na plikach:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... kod workflow ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // Zapisz do pliku logu
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Handler `onError`

Oprócz `onComplete` istnieje jeszcze jeden handler zdarzeń: `onError`, który uruchamia się tylko w przypadku niepowodzenia workflow:

```groovy title="main.nf - onError handler"
workflow {
    // ... kod workflow ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Zapisz szczegółowy log błędów
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

Możesz używać wielu handlerów razem w skrypcie workflow:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... kod workflow ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Domknięcia handlerów zdarzeń**: Domknięcia w skrypcie workflow uruchamiane w różnych punktach cyklu życia
- **Handler `onComplete`**: Do podsumowań wykonania i raportowania wyników
- **Handler `onError`**: Do obsługi błędów i logowania niepowodzeń
- **Właściwości obiektu workflow**: Dostęp do `workflow.success`, `workflow.duration`, `workflow.errorMessage` itp.

Handlery zdarzeń pokazują, jak możesz używać pełnej mocy języka Nextflow w skryptach workflow, aby dodać zaawansowane możliwości logowania i powiadamiania.

---

## Podsumowanie

Gratulacje, udało Ci się!

W trakcie tego side questu zbudowałeś kompleksowy pipeline przetwarzania próbek, który ewoluował od podstawowej obsługi metadanych do zaawansowanego, gotowego na środowisko produkcyjne workflow.
Każda sekcja budowała na poprzedniej, demonstrując, jak konstrukcje programistyczne transformują proste workflow'y w potężne systemy przetwarzania danych, z następującymi korzyściami:

- **Czystszy kod**: Rozumienie dataflow a skryptowania pomaga pisać bardziej zorganizowane workflow'y
- **Solidna obsługa**: Operatory bezpiecznej nawigacji i Elvis czynią workflow'y odpornymi na brakujące dane
- **Elastyczne przetwarzanie**: Logika warunkowa pozwala workflow'om odpowiednio przetwarzać różne typy próbek
- **Adaptacyjne zasoby**: Dynamiczne dyrektywy optymalizują wykorzystanie zasobów na podstawie cech wejścia

Ta progresja odzwierciedla rzeczywistą ewolucję pipeline'ów bioinformatycznych — od prototypów badawczych obsługujących kilka próbek do systemów produkcyjnych przetwarzających tysiące próbek w laboratoriach i instytucjach.
Każde wyzwanie, które rozwiązałeś, i każdy wzorzec, którego się nauczyłeś, odzwierciedla rzeczywiste problemy, z którymi deweloperzy mierzą się przy skalowaniu workflow'ów Nextflow.

Stosowanie tych wzorców we własnej pracy pozwoli Ci budować solidne, gotowe na środowisko produkcyjne workflow'y.

### Kluczowe wzorce

1.  **Dataflow a skryptowanie:** Nauczyłeś się odróżniać operacje dataflow (orkiestracja kanałów) od skryptowania (kod manipulujący danymi), w tym kluczowe różnice między operacjami na różnych typach, jak `collect` na kanale vs liście.

    - Dataflow: orkiestracja kanałów

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Skryptowanie: przetwarzanie danych na kolekcjach

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Zaawansowane przetwarzanie stringów**: Opanowałeś wyrażenia regularne do parsowania nazw plików, dynamiczne generowanie skryptów w procesach i interpolację zmiennych (Nextflow vs Bash vs powłoka).

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

    - Kolekcja plików do argumentów polecenia (w bloku skryptu procesu)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Tworzenie funkcji wielokrotnego użytku**: Nauczyłeś się wyodrębniać złożoną logikę do nazwanych funkcji, które można wywoływać z operatorów kanałów, czyniąc workflow'y bardziej czytelnymi i łatwymi w utrzymaniu.

    - Definiowanie nazwanej funkcji

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Wywoływanie nazwanej funkcji w workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Dynamiczne dyrektywy zasobów z domknięciami**: Zbadałeś używanie domknięć w dyrektywach procesów do adaptacyjnej alokacji zasobów na podstawie cech wejścia.

    - Nazwane domknięcia ikompozycja

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Domknięcia z dostępem do zakresu

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logika warunkowa i kontrola procesów**: Dodałeś inteligentne trasowanie za pomocą operatorów `.branch()` i `.filter()`, wykorzystując prawdziwość do pisania zwięzłych wyrażeń warunkowych.

    - Używanie `.branch()` do kierowania danych przez różne gałęzie workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Ewaluacja boolowska z prawdziwością Groovy

    ```groovy
    if (sample.files) println "Has files"
    ```

    - Używanie `filter()` do podzbioru danych z „prawdziwością"

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operatory bezpiecznej nawigacji i Elvis**: Uczyniłeś pipeline odpornym na brakujące dane, używając `?.` do bezpiecznego dostępu do właściwości i `?:` do dostarczania wartości domyślnych.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Walidacja za pomocą error() i log.warn**: Nauczyłeś się walidować wejścia wcześnie i szybko zawodzić z czytelnymi komunikatami błędów.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Handlery zdarzeń konfiguracji**: Nauczyłeś się używać handlerów zdarzeń workflow (`onComplete` i `onError`) do logowania, powiadomień i zarządzania cyklem życia.

    - Używanie `onComplete` do logowania i powiadamiania

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Używanie `onError` do działania wyłącznie w przypadku niepowodzenia

    ```groovy
    workflow.onError = {
        // Zapisz szczegółowy log błędów
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Dodatkowe zasoby

- [Nextflow Language Reference](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Syntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standard Library](https://nextflow.io/docs/latest/reference/stdlib.html)

Zajrzyj do tych zasobów, gdy będziesz chciał zgłębić bardziej zaawansowane funkcje.

Ćwiczenie i rozwijanie umiejętności pozwoli Ci:

- Pisać czystsze workflow'y z właściwym rozdzieleniem dataflow od skryptowania
- Opanować interpolację zmiennych, aby unikać typowych pułapek ze zmiennymi Nextflow, Bash i powłoki
- Używać dynamicznych dyrektyw zasobów do wydajnych, adaptacyjnych workflow'ów
- Transformować kolekcje plików w poprawnie sformatowane argumenty wiersza poleceń
- Obsługiwać różne konwencje nazewnictwa plików i formaty wejściowe za pomocą regex i przetwarzania stringów
- Budować wielokrotnego użytku, łatwy w utrzymaniu kod przy użyciu zaawansowanych wzorców domknięć i programowania funkcyjnego
- Przetwarzać i organizować złożone zestawy danych za pomocą operacji na kolekcjach
- Dodawać walidację, obsługę błędów i logowanie, aby Twoje workflow'y były gotowe na środowisko produkcyjne
- Implementować zarządzanie cyklem życia workflow za pomocą handlerów zdarzeń

---

## Co dalej?

Wróć do [menu Side Quests](../) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
