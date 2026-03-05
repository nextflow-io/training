# Podstawowe Wzorce Skryptowe w Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow to język programowania działający na maszynie wirtualnej Java. Choć Nextflow jest zbudowany na [Groovy](http://groovy-lang.org/) i dzieli z nim znaczną część składni, Nextflow to coś więcej niż tylko „Groovy z rozszerzeniami" — jest to samodzielny język z w pełni określoną [składnią](https://nextflow.io/docs/latest/reference/syntax.html) i [biblioteką standardową](https://nextflow.io/docs/latest/reference/stdlib.html).

Możesz napisać wiele kodu w Nextflow, nie wykraczając poza podstawową składnię zmiennych, map i list. Większość tutoriali Nextflow koncentruje się na orkiestracji workflow'ów (kanały, procesy i przepływ danych) i możesz zajść zaskakująco daleko, posługując się tylko tym.

Jednak gdy musisz manipulować danymi, parsować złożone nazwy plików, implementować logikę warunkową lub budować solidne workflow'y produkcyjne, warto myśleć o dwóch odrębnych aspektach Twojego kodu: **przepływie danych** (kanały, operatory, procesy i workflow'y) oraz **skryptowaniu** (kod wewnątrz domknięć, funkcji i skryptów procesów). Choć to rozróżnienie jest nieco arbitralne — to wszystko jest kod Nextflow — dostarcza użytecznego modelu myślowego do zrozumienia, kiedy orkiestrujesz swój pipeline, a kiedy manipulujesz danymi. Opanowanie obu aspektów znacząco poprawia Twoją zdolność do pisania przejrzystych, łatwych w utrzymaniu workflow'ów.

### Cele szkolenia

Ten side quest zabierze Cię w praktyczną podróż od podstawowych koncepcji do wzorców gotowych do produkcji.
Przekształcimy prosty workflow czytający CSV w zaawansowany pipeline bioinformatyczny, rozwijając go krok po kroku przez realistyczne wyzwania:

- **Rozumienie granic:** Rozróżnianie między operacjami przepływu danych a skryptowaniem oraz zrozumienie, jak współpracują ze sobą
- **Manipulacja danymi:** Wyodrębnianie, przekształcanie i wybieranie podzbiorów map i kolekcji przy użyciu potężnych operatorów
- **Przetwarzanie napisów:** Parsowanie złożonych schematów nazewnictwa plików za pomocą wzorców regex i opanowanie interpolacji zmiennych
- **Funkcje wielokrotnego użytku:** Wyodrębnianie złożonej logiki do nazwanych funkcji dla czystszych, łatwiejszych w utrzymaniu workflow'ów
- **Dynamiczna logika:** Budowanie procesów, które dostosowują się do różnych typów wejściowych i używanie domknięć do dynamicznej alokacji zasobów
- **Warunkowe kierowanie:** Inteligentne kierowanie próbek przez różne procesy na podstawie ich cech metadanych
- **Bezpieczne operacje:** Eleganckie radzenie sobie z brakującymi danymi za pomocą operatorów bezpiecznych dla null i walidacja wejść z czytelnymi komunikatami błędów
- **Obsługa oparta na konfiguracji:** Używanie procedur obsługi zdarzeń workflow'a do logowania, powiadomień i zarządzania cyklem życia

### Wymagania wstępne

Przed podjęciem tego side questa powinieneś:

- Ukończyć tutorial [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory, praca z plikami, metadane)
- Mieć podstawową znajomość powszechnych konstrukcji programistycznych (zmienne, mapy, listy)

Ten tutorial będzie wyjaśniał koncepcje programistyczne w miarę ich napotkania, więc nie potrzebujesz rozległego doświadczenia programistycznego.
Zaczniemy od fundamentalnych koncepcji i będziemy budować w kierunku zaawansowanych wzorców.

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe w codespace

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzyłeś środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki dla tego tutoriala.

```bash
cd side-quests/essential_scripting_patterns
```

#### Przejrzyj materiały

Znajdziesz główny plik workflow'a oraz katalog `data` zawierający przykładowe pliki danych.

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

Nasz przykładowy CSV zawiera informacje o próbkach biologicznych, które wymagają różnego przetwarzania w zależności od ich charakterystyk:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Użyjemy tego realistycznego zbioru danych do eksploracji praktycznych technik programistycznych, które napotkasz w rzeczywistych workflow'ach bioinformatycznych.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy
<!-- - [ ] I understand the assignment -->

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Przepływ danych vs Skryptowanie: Rozumienie Granic

### 1.1. Identyfikowanie Czym Jest Co

Pisząc workflow'y w Nextflow, ważne jest rozróżnienie między **przepływem danych** (jak dane przemieszczają się przez kanały i procesy) a **skryptowaniem** (kod, który manipuluje danymi i podejmuje decyzje). Zbudujmy workflow demonstrujący, jak współpracują ze sobą.

#### 1.1.1. Podstawowy Workflow w Nextflow

Zacznij od prostego workflow'a, który po prostu odczytuje plik CSV (zrobiliśmy to już dla Ciebie w `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Blok `workflow` definiuje strukturę naszego pipeline'u, podczas gdy `channel.fromPath()` tworzy kanał ze ścieżki pliku. Operator `.splitCsv()` przetwarza plik CSV i konwertuje każdy wiersz na strukturę danych mapy.

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

Teraz dodamy skryptowanie do przekształcenia danych, używając operatora `.map()`, który prawdopodobnie już znasz. Ten operator przyjmuje 'domknięcie', w którym możemy napisać kod do przekształcenia każdego elementu.

!!! note "Uwaga"

    **Domknięcie** to blok kodu, który może być przekazywany i wykonywany później. Pomyśl o nim jak o funkcji, którą definiujesz inline. Domknięcia są zapisywane za pomocą nawiasów klamrowych `{ }` i mogą przyjmować parametry. Są fundamentalne dla działania operatorów Nextflow i jeśli piszesz w Nextflow od jakiegoś czasu, możliwe, że już ich używałeś, nie zdając sobie z tego sprawy!

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

To nasze pierwsze **domknięcie** — funkcja anonimowa, którą możesz przekazać jako argument (podobnie do lambd w Pythonie lub funkcji strzałkowych w JavaScript). Domknięcia są niezbędne do pracy z operatorami Nextflow.

Domknięcie `{ row -> return row }` przyjmuje parametr `row` (może to być dowolna nazwa: `item`, `sample`, itp.).

Gdy operator `.map()` przetwarza każdy element kanału, przekazuje ten element do Twojego domknięcia. Tutaj `row` przechowuje jeden wiersz CSV na raz.

Zastosuj tę zmianę i uruchom workflow:

```bash
nextflow run main.nf
```

Zobaczysz to samo wyjście co wcześniej, ponieważ po prostu zwracamy wejście bez zmian. To potwierdza, że operator map działa poprawnie. Teraz zacznijmy przekształcać dane.

#### 1.1.3. Tworzenie Struktury Danych Mapy

Teraz napiszemy logikę **skryptowania** wewnątrz naszego domknięcia, aby przekształcić każdy wiersz danych. To tutaj przetwarzamy poszczególne elementy danych, zamiast orkiestrować przepływ danych.

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Skryptowanie do przekształcenia danych
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

Mapa `sample_meta` to struktura danych klucz-wartość (jak słowniki w Pythonie, obiekty w JavaScript lub hasze w Ruby) przechowująca powiązane informacje: ID próbki, organizm, typ tkanki, głębokość sekwencjonowania i wynik jakości.

Używamy metod manipulacji napisami takich jak `.toLowerCase()` i `.replaceAll()` do czyszczenia naszych danych oraz metod konwersji typów takich jak `.toInteger()` i `.toDouble()` do konwersji danych napisowych z CSV na odpowiednie typy numeryczne.

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

Teraz dodajmy więcej skryptowania — tym razem używając operatora trójargumentowego do podejmowania decyzji na podstawie wartości danych.

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

Operator trójargumentowy to skrót dla instrukcji if/else, który podąża za wzorcem `warunek ? wartość_jeśli_prawda : wartość_jeśli_fałsz`. Ta linia oznacza: „Jeśli jakość jest większa niż 40, użyj 'high', w przeciwnym razie użyj 'normal'". Jego kuzyn, **operator Elvisa** (`?:`), dostarcza wartości domyślnych, gdy coś jest null lub puste — zbadamy ten wzorzec później w tym tutorialu.

Operator dodawania map `+` tworzy **nową mapę** zamiast modyfikować istniejącą. Ta linia tworzy nową mapę, która zawiera wszystkie pary klucz-wartość z `sample_meta` plus nowy klucz `priority`.

!!! Note "Uwaga"

    Nigdy nie modyfikuj map przekazanych do domknięć — zawsze twórz nowe używając `+` (na przykład). W Nextflow te same dane często przepływają przez wiele operacji jednocześnie. Modyfikowanie mapy w miejscu może powodować nieprzewidywalne efekty uboczne, gdy inne operacje odwołują się do tego samego obiektu. Tworzenie nowych map zapewnia, że każda operacja ma swoją własną czystą kopię.

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

#### 1.1.5. Wybieranie Podzbiorów Map za Pomocą `.subMap()`

Podczas gdy operator `+` dodaje klucze do mapy, czasami musisz zrobić coś przeciwnego — wyodrębnić tylko określone klucze. Metoda `.subMap()` jest do tego idealna.

Dodajmy linię, aby utworzyć uproszczoną wersję naszych metadanych, która zawiera tylko pola identyfikacyjne:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Skryptowanie do przekształcenia danych
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
                // Skryptowanie do przekształcenia danych
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

To pokazuje zarówno pełne metadane wyświetlone przez operację `view()`, jak i wyodrębniony podzbiór, który wydrukowaliśmy za pomocą `println`.

Metoda `.subMap()` przyjmuje listę kluczy i zwraca nową mapę zawierającą tylko te klucze. Jeśli klucz nie istnieje w oryginalnej mapie, po prostu nie jest uwzględniany w wyniku.

Jest to szczególnie użyteczne, gdy musisz utworzyć różne wersje metadanych dla różnych procesów — niektóre mogą potrzebować pełnych metadanych, podczas gdy inne potrzebują tylko minimalnych pól identyfikacyjnych.

Teraz usuń te instrukcje println, aby przywrócić workflow do poprzedniego stanu, ponieważ nie będą nam potrzebne w dalszej części.

!!! tip "Podsumowanie Operacji na Mapach"

    - **Dodawanie kluczy**: `map1 + [new_key: value]` - Tworzy nową mapę z dodatkowymi kluczami
    - **Wyodrębnianie kluczy**: `map1.subMap(['key1', 'key2'])` - Tworzy nową mapę tylko z określonymi kluczami
    - **Obie operacje tworzą nowe mapy** - Oryginalne mapy pozostają niezmienione

#### 1.1.6. Łączenie Map i Zwracanie Wyników

Do tej pory zwracaliśmy tylko to, co społeczność Nextflow nazywa 'mapą meta', i ignorowaliśmy pliki, do których te metadane się odnoszą. Ale jeśli piszesz workflow'y w Nextflow, prawdopodobnie chcesz coś zrobić z tymi plikami.

Wyprodukujmy strukturę kanału składającą się z krotki 2 elementów: wzbogaconej mapy metadanych i odpowiadającej ścieżki pliku. To powszechny wzorzec w Nextflow do przekazywania danych do procesów.

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

    **Mapy i Metadane**: Mapy są fundamentalne dla pracy z metadanymi w Nextflow. Aby uzyskać bardziej szczegółowe wyjaśnienie pracy z mapami metadanych, zobacz side quest [Working with metadata](./metadata.md).

Nasz workflow demonstruje podstawowy wzorzec: **operacje przepływu danych** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orkiestrują, jak dane przemieszczają się przez pipeline, podczas gdy **skryptowanie** (mapy `[key: value]`, metody napisów, konwersje typów, operatory trójargumentowe) wewnątrz domknięcia `.map()` obsługuje przekształcenie poszczególnych elementów danych.

### 1.2. Rozumienie Różnych Typów: Channel vs List

Do tej pory dobrze nam idzie, możemy rozróżnić między operacjami przepływu danych a skryptowaniem. Ale co z sytuacją, gdy ta sama nazwa metody istnieje w obu kontekstach?

Doskonałym przykładem jest metoda `collect`, która istnieje zarówno dla typów kanałów, jak i typów List w bibliotece standardowej Nextflow. Metoda `collect()` na Liście przekształca każdy element, podczas gdy operator `collect()` na kanale zbiera wszystkie emisje kanału w kanał jednoelementowy.

Zademonstrujmy to na przykładowych danych, zaczynając od odświeżenia sobie, co robi operator `collect()` kanału. Sprawdź `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - grupuje wiele emisji kanału w jedną
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Kroki:

- Zdefiniuj Listę ID próbek
- Utwórz kanał za pomocą `fromList()`, który emituje każde ID próbki osobno
- Wydrukuj każdy element za pomocą `view()` w miarę przepływu
- Zbierz wszystkie elementy w jedną listę za pomocą operatora `collect()` kanału
- Wydrukuj zebrany wynik (pojedynczy element zawierający wszystkie ID próbek) za pomocą drugiego `view()`

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

`view()` zwraca wyjście dla każdej emisji kanału, więc wiemy, że to pojedyncze wyjście zawiera wszystkie 3 oryginalne elementy zgrupowane w jedną listę.

Teraz zobaczmy metodę `collect` na Liście w akcji. Zmodyfikuj `collect.nf`, aby zastosować metodę `collect` Listy do oryginalnej listy ID próbek:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - przekształca każdy element, zachowuje strukturę
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

- Definiujemy nową zmienną `formatted_ids`, która używa metody `collect` Listy do przekształcenia każdego ID próbki w oryginalnej liście
- Drukujemy wynik używając `println`

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

Tym razem NIE zmieniliśmy struktury danych, nadal mamy 3 elementy na liście, ale przekształciliśmy każdy element używając metody `collect` Listy, aby wyprodukować nową listę ze zmodyfikowanymi wartościami. Jest to podobne do użycia operatora `map` na kanale, ale operuje na strukturze danych List zamiast na kanale.

`collect` to ekstremalny przypadek, którego używamy tutaj, aby zilustrować punkt. Kluczową lekcją jest to, że pisząc workflow'y, zawsze rozróżniaj między **strukturami danych** (Listy, Mapy, itp.) a **kanałami** (konstrukcje przepływu danych). Operacje mogą dzielić nazwy, ale zachowują się zupełnie inaczej w zależności od typu, na którym są wywoływane.

### 1.3. Operator Rozproszenia (`*.`) - Skrót do Wyodrębniania Właściwości

Powiązany z metodą `collect` Listy jest operator rozproszenia (`*.`), który dostarcza zwięzłego sposobu wyodrębniania właściwości z kolekcji. Jest to zasadniczo cukier składniowy dla powszechnego wzorca `collect`.

Dodajmy demonstrację do naszego pliku `collect.nf`:

=== "Po"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - grupuje wiele emisji kanału w jedną
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - przekształca każdy element, zachowuje strukturę
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Operator rozproszenia - zwięzły dostęp do właściwości
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

    // List.collect() - przekształca każdy element, zachowuje strukturę
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

Operator rozproszenia `*.` to skrót dla powszechnego wzorca collect:

```groovy
// Te są równoważne:
def ids = samples*.id
def ids = samples.collect { it.id }

// Działa również z wywołaniami metod:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Operator rozproszenia jest szczególnie użyteczny, gdy musisz wyodrębnić pojedynczą właściwość z listy obiektów — jest bardziej czytelny niż wypisywanie pełnego domknięcia `collect`.

!!! tip "Kiedy Używać Rozproszenia vs Collect"

    - **Użyj rozproszenia (`*.`)** dla prostego dostępu do właściwości: `samples*.id`, `files*.name`
    - **Użyj collect** dla przekształceń lub złożonej logiki: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Przepływ danych vs skryptowanie**: Operatory kanałów orkiestrują, jak dane przepływają przez Twój pipeline, podczas gdy skryptowanie przekształca poszczególne elementy danych
- **Rozumienie typów**: Ta sama nazwa metody (jak `collect`) może zachowywać się inaczej w zależności od typu, na którym jest wywoływana (Channel vs List)
- **Kontekst ma znaczenie**: Zawsze bądź świadomy, czy pracujesz z kanałami (przepływ danych) czy strukturami danych (skryptowanie)

Rozumienie tych granic jest niezbędne do debugowania, dokumentacji i pisania łatwych w utrzymaniu workflow'ów.

Następnie zagłębimy się w możliwości przetwarzania napisów, które są niezbędne do obsługi rzeczywistych danych.

---

## 2. Przetwarzanie Napisów i Dynamiczne Generowanie Skryptów

Opanowanie przetwarzania napisów oddziela kruche workflow'y od solidnych pipeline'ów. Ta sekcja obejmuje parsowanie złożonych nazw plików, dynamiczne generowanie skryptów i interpolację zmiennych.

### 2.1. Dopasowywanie Wzorców i Wyrażenia Regularne

Pliki bioinformatyczne często mają złożone konwencje nazewnictwa kodujące metadane. Wyodrębnijmy to automatycznie używając dopasowywania wzorców z wyrażeniami regularnymi.

Wrócimy do naszego workflow'a `main.nf` i dodamy logikę dopasowywania wzorców, aby wyodrębnić dodatkowe informacje o próbkach z nazw plików. Pliki FASTQ w naszym zbiorze danych podążają za konwencjami nazewnictwa w stylu Illumina z nazwami takimi jak `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Mogą wyglądać tajemniczo, ale faktycznie kodują użyteczne metadane, takie jak ID próbki, numer ścieżki i kierunek odczytu. Użyjemy możliwości regex do parsowania tych nazw.

Wprowadź następującą zmianę do istniejącego workflow'a `main.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Skryptowanie do przekształcenia danych
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
                // Skryptowanie do przekształcenia danych
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

To demonstruje kluczowe **koncepcje przetwarzania napisów**:

1. **Literały wyrażeń regularnych** używające składni `~/wzorzec/` - tworzy wzorzec regex bez potrzeby escapowania ukośników wstecznych
2. **Dopasowywanie wzorców** za pomocą operatora `=~` - próbuje dopasować napis do wzorca regex
3. **Obiekty dopasowujące**, które przechwytują grupy za pomocą `[0][1]`, `[0][2]`, itp. - `[0]` odnosi się do całego dopasowania, `[1]`, `[2]`, itp. odnoszą się do przechwyconych grup w nawiasach

Rozbijmy wzorzec regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Wzorzec             | Dopasowuje                                  | Przechwytuje                              |
| ------------------- | ------------------------------------------- | ----------------------------------------- |
| `^(.+)`             | Nazwa próbki od początku                    | Grupa 1: nazwa próbki                     |
| `_S(\d+)`           | Numer próbki `_S1`, `_S2`, itp.             | Grupa 2: numer próbki                     |
| `_L(\d{3})`         | Numer ścieżki `_L001`                       | Grupa 3: ścieżka (3 cyfry)                |
| `_(R[12])`          | Kierunek odczytu `_R1` lub `_R2`            | Grupa 4: kierunek odczytu                 |
| `_(\d{3})`          | Numer fragmentu `_001`                      | Grupa 5: fragment (3 cyfry)               |
| `\.fastq(?:\.gz)?$` | Rozszerzenie pliku `.fastq` lub `.fastq.gz` | Nieprzechwycone (?: to nieprzechwytujące) |

To parsuje konwencje nazewnictwa w stylu Illumina, aby automatycznie wyodrębnić metadane.

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

To pokazuje metadane wzbogacone z nazw plików.

### 2.2. Dynamiczne Generowanie Skryptów w Procesach

Bloki skryptów procesów są zasadniczo wieloliniowymi napisami, które są przekazywane do powłoki. Możesz użyć **logiki warunkowej** (if/else, operatory trójargumentowe) do dynamicznego generowania różnych napisów skryptowych na podstawie charakterystyk wejściowych. Jest to niezbędne do obsługi różnorodnych typów wejściowych — takich jak odczyty single-end vs paired-end — bez duplikowania definicji procesów.

Dodajmy proces do naszego workflow'a, który demonstruje ten wzorzec. Otwórz `modules/fastp.nf` i spójrz:

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

Proces przyjmuje pliki FASTQ jako wejście i uruchamia narzędzie `fastp` do przycinania adapterów i filtrowania odczytów niskiej jakości. Niestety osoba, która napisała ten proces, nie uwzględniła odczytów single-end, które mamy w naszym przykładowym zbiorze danych. Dodajmy go do naszego workflow'a i zobaczmy, co się stanie:

Najpierw dołącz moduł w pierwszej linii Twojego workflow'a `main.nf`:

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

Widać, że proces próbuje uruchomić `fastp` z wartością `null` dla drugiego pliku wejściowego, co powoduje jego niepowodzenie. Dzieje się tak, ponieważ nasz zbiór danych zawiera odczyty single-end, ale proces jest zakodowany na stałe, aby oczekiwać odczytów paired-end (dwa pliki wejściowe na raz).

Napraw to, dodając logikę warunkową do bloku `script:` procesu `FASTP`. Instrukcja if/else sprawdza liczbę plików odczytów i odpowiednio dostosowuje polecenie.

=== "Po"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Prosta detekcja single-end vs paired-end
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

Teraz workflow może elegancko obsługiwać zarówno odczyty single-end, jak i paired-end. Logika warunkowa sprawdza liczbę plików wejściowych i konstruuje odpowiednie polecenie dla `fastp`. Zobaczmy, czy działa:

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

Wygląda dobrze! Jeśli sprawdzimy faktyczne polecenia, które zostały uruchomione (dostosuj dla Twojego hasha zadania):

```console title="Check commands executed"
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

Inne powszechne użycie dynamicznej logiki skryptowej można zobaczyć w [module Nextflow for Science Genomics](../../nf4science/genomics/02_joint_calling). W tym module wywoływany proces GATK może przyjmować wiele plików wejściowych, ale każdy musi być poprzedzony `-V`, aby utworzyć poprawną linię poleceń. Proces używa skryptowania do przekształcenia kolekcji plików wejściowych (`all_gvcfs`) w poprawne argumenty polecenia:

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

Te wzorce używania skryptowania w blokach skryptów procesów są niezwykle potężne i mogą być stosowane w wielu scenariuszach — od obsługi zmiennych typów wejściowych po budowanie złożonych argumentów linii poleceń z kolekcji plików, czyniąc Twoje procesy naprawdę adaptacyjnymi do różnorodnych wymagań rzeczywistych danych.

### 2.3. Interpolacja Zmiennych: Zmienne Nextflow i Powłoki

Skrypty procesów mieszają zmienne Nextflow, zmienne powłoki i podstawienia poleceń, każde z inną składnią interpolacji. Użycie niewłaściwej składni powoduje błędy. Zbadajmy to za pomocą procesu, który tworzy raport przetwarzania.

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
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Ten proces zapisuje prosty raport z ID próbki i nazwą pliku. Teraz uruchommy go, aby zobaczyć, co się dzieje, gdy musimy mieszać różne typy zmiennych.

Dołącz proces w swoim `main.nf` i dodaj go do workflow'a:

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

Ale co, jeśli chcemy dodać informacje o tym, kiedy i gdzie nastąpiło przetwarzanie? Zmodyfikujmy proces, aby użyć zmiennych **powłoki** i trochę podstawienia poleceń, aby uwzględnić bieżącego użytkownika, nazwę hosta i datę w raporcie:

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

Musimy to escapować, aby Bash mógł to obsłużyć zamiast tego.

Napraw to, escapując zmienne powłoki i podstawienia poleceń ukośnikiem wstecznym (`\`):

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

Teraz działa! Ukośnik wsteczny (`\`) mówi Nextflow „nie interpretuj tego, przekaż to do Bash".

### Podsumowanie

W tej sekcji nauczyłeś się technik **przetwarzania napisów**:

- **Wyrażenia regularne do parsowania plików**: Używanie operatora `=~` i wzorców regex (`~/wzorzec/`) do wyodrębniania metadanych ze złożonych konwencji nazewnictwa plików
- **Dynamiczne generowanie skryptów**: Używanie logiki warunkowej (if/else, operatory trójargumentowe) do generowania różnych napisów skryptowych na podstawie charakterystyk wejściowych
- **Interpolacja zmiennych**: Rozumienie, kiedy Nextflow interpretuje napisy vs kiedy robi to powłoka
  - `${var}` - zmienne Nextflow (interpolowane przez Nextflow w czasie kompilacji workflow'a)
  - `\${var}` - zmienne środowiskowe powłoki (escapowane, przekazywane do bash w czasie wykonania)
  - `\$(cmd)` - podstawienie polecenia powłoki (escapowane, wykonywane przez bash w czasie wykonania)

Te wzorce przetwarzania i generowania napisów są niezbędne do obsługi różnorodnych formatów plików i konwencji nazewnictwa, które napotkasz w rzeczywistych workflow'ach bioinformatycznych.

---

## 3. Tworzenie Funkcji Wielokrotnego Użytku

Złożona logika workflow'a inline w operatorach kanałów lub definicjach procesów zmniejsza czytelność i łatwość utrzymania. **Funkcje** pozwalają wyodrębnić tę logikę do nazwanych, wielokrotnego użytku komponentów.

Nasza operacja map urosła i stała się złożona. Wyodrębnijmy ją do funkcji wielokrotnego użytku używając słowa kluczowego `def`.

Aby zilustrować, jak to wygląda z naszym istniejącym workflow'em, wprowadź poniższą modyfikację, używając `def` do zdefiniowania funkcji wielokrotnego użytku o nazwie `separateMetadata`:

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

Wyodrębniając tę logikę do funkcji, zredukowaliśmy faktyczną logikę workflow'a do czegoś znacznie czystszego:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

To sprawia, że logika workflow'a jest znacznie łatwiejsza do odczytania i zrozumienia na pierwszy rzut oka. Funkcja `separateMetadata` enkapsuluje całą złożoną logikę parsowania i wzbogacania metadanych, czyniąc ją wielokrotnego użytku i testowalną.

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

Wyjście powinno pokazywać oba procesy kończące się pomyślnie. Workflow jest teraz znacznie czystszy i łatwiejszy w utrzymaniu, z całą złożoną logiką przetwarzania metadanych enkapsulowaną w funkcji `separateMetadata`.

### Podsumowanie

W tej sekcji nauczyłeś się **tworzenia funkcji**:

- **Definiowanie funkcji za pomocą `def`**: Słowo kluczowe do tworzenia nazwanych funkcji (jak `def` w Pythonie lub `function` w JavaScript)
- **Zakres funkcji**: Funkcje zdefiniowane na poziomie skryptu są dostępne w całym Twoim workflow'ie Nextflow
- **Wartości zwracane**: Funkcje automatycznie zwracają ostatnie wyrażenie lub używają jawnego `return`
- **Czystszy kod**: Wyodrębnianie złożonej logiki do funkcji to fundamentalna praktyka inżynierii oprogramowania w każdym języku

Następnie zbadamy, jak używać domknięć w dyrektywach procesów do dynamicznej alokacji zasobów.

---

## 4. Dynamiczne Dyrektywy Zasobów z Domknięciami

Do tej pory używaliśmy skryptowania w bloku `script` procesów. Ale **domknięcia** (wprowadzone w Sekcji 1.1) są również niezwykle użyteczne w dyrektywach procesów, szczególnie do dynamicznej alokacji zasobów. Dodajmy dyrektywy zasobów do naszego procesu FASTP, które dostosowują się na podstawie charakterystyk próbki.

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

Domknięcie `{ meta.depth > 40000000 ? 2 : 1 }` używa **operatora trójargumentowego** (omówionego w Sekcji 1.1) i jest ewaluowane dla każdego zadania, umożliwiając alokację zasobów per próbka. Próbki o wysokiej głębokości (>40M odczytów) otrzymują 2 CPU, podczas gdy inne otrzymują 1 CPU.

!!! note "Dostęp do Zmiennych Wejściowych w Dyrektywach"

    Domknięcie może uzyskać dostęp do dowolnych zmiennych wejściowych (jak `meta` tutaj), ponieważ Nextflow ewaluuje te domknięcia w kontekście każdego wykonania zadania.

Uruchom workflow ponownie z opcją `-ansi-log false`, aby łatwiej było zobaczyć hasze zadań.

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

Możesz sprawdzić dokładne polecenie `docker`, które zostało uruchomione, aby zobaczyć alokację CPU dla danego zadania:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Powinieneś zobaczyć coś takiego:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

W tym przykładzie wybraliśmy przykład, który zażądał 2 CPU (`--cpu-shares 2048`), ponieważ była to próbka o wysokiej głębokości, ale powinieneś zobaczyć różne alokacje CPU w zależności od głębokości próbki. Spróbuj tego również dla innych zadań.

### 4.2. Strategie ponawiania

Innym potężnym wzorcem jest używanie `task.attempt` dla strategii ponawiania. Aby pokazać, dlaczego jest to użyteczne, zaczniemy od zmniejszenia alokacji pamięci dla FASTP do mniej niż potrzebuje. Zmień dyrektywę `memory` w `modules/fastp.nf` na `1.GB`:

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

To bardzo powszechny scenariusz w rzeczywistych workflow'ach — czasami po prostu nie wiesz, ile pamięci będzie potrzebować zadanie, dopóki go nie uruchomisz.

Aby uczynić nasz workflow bardziej solidnym, możemy zaimplementować strategię ponawiania, która zwiększa alokację pamięci przy każdej próbie, ponownie używając domknięcia Groovy. Zmodyfikuj dyrektywę `memory`, aby mnożyć bazową pamięć przez `task.attempt`, i dodaj dyrektywy `errorStrategy 'retry'` oraz `maxRetries 2`:

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

Teraz jeśli proces zawiedzie z powodu niewystarczającej pamięci, Nextflow ponowi próbę z większą pamięcią:

- Pierwsza próba: 1 GB (task.attempt = 1)
- Druga próba: 2.GB (task.attempt = 2)

... i tak dalej, aż do limitu `maxRetries`.

### Podsumowanie

Dynamiczne dyrektywy z domknięciami pozwalają:

- Alokować zasoby na podstawie charakterystyk wejściowych
- Implementować automatyczne strategie ponawiania ze zwiększającymi się zasobami
- Łączyć wiele czynników (metadane, numer próby, priorytety)
- Używać logiki warunkowej do złożonych obliczeń zasobów

To sprawia, że Twoje workflow'y są zarówno bardziej wydajne (nie nadmiernie alokujące), jak i bardziej solidne (automatyczne ponawianie z większymi zasobami).

---

## 5. Logika Warunkowa i Kontrola Procesów

Wcześniej używaliśmy `.map()` ze skryptowaniem do przekształcania danych kanału. Teraz użyjemy logiki warunkowej do kontrolowania, które procesy są wykonywane na podstawie danych — niezbędne dla elastycznych workflow'ów dostosowujących się do różnych typów próbek.

[Operatory przepływu danych](https://www.nextflow.io/docs/latest/reference/operator.html) Nextflow przyjmują domknięcia ewaluowane w czasie wykonania, umożliwiając logice warunkowej kierowanie decyzjami workflow'a na podstawie zawartości kanału.

### 5.1. Kierowanie za pomocą `.branch()`

Na przykład, udajmy, że nasze próbki sekwencjonowania muszą być przycinane za pomocą FASTP tylko jeśli są próbkami ludzkimi z pokryciem powyżej pewnego progu. Próbki mysie lub próbki o niskim pokryciu powinny być uruchamiane z Trimgalore zamiast tego (to jest wymyślony przykład, ale ilustruje punkt).

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

... a następnie zmodyfikuj swój workflow `main.nf`, aby rozgałęziać próbki na podstawie ich metadanych i kierować je przez odpowiedni proces przycinania, w ten sposób:

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

Innym potężnym wzorcem do kontrolowania wykonania workflow'a jest operator `.filter()`, który używa domknięcia do określenia, które elementy powinny kontynuować w pipeline. Wewnątrz domknięcia filtra napiszesz **wyrażenia boolowskie**, które decydują, które elementy przechodzą.

Nextflow (jak wiele języków dynamicznych) ma koncepcję **"prawdziwości"**, która określa, jakie wartości ewaluują się do `true` lub `false` w kontekstach boolowskich:

- **Prawdziwe**: Wartości nie-null, niepuste napisy, niezerowe liczby, niepuste kolekcje
- **Fałszywe**: `null`, puste napisy `""`, zero `0`, puste kolekcje `[]` lub `[:]`, `false`

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

- `meta.id && meta.organism` sprawdza, że oba pola istnieją i nie są puste (używając prawdziwości)
- `meta.depth >= 25000000` zapewnia wystarczającą głębokość sekwencjonowania za pomocą jawnego porównania

!!! note "Prawdziwość w Praktyce"

    Wyrażenie `meta.id && meta.organism` jest bardziej zwięzłe niż pisanie:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    To sprawia, że logika filtrowania jest znacznie czystsza i łatwiejsza do odczytania.

### Podsumowanie

W tej sekcji nauczyłeś się używać logiki warunkowej do kontrolowania wykonania workflow'a za pomocą interfejsów domknięć operatorów Nextflow takich jak `.branch{}` i `.filter{}`, wykorzystując prawdziwość do pisania zwięzłych wyrażeń warunkowych.

Nasz pipeline teraz inteligentnie kieruje próbki przez odpowiednie procesy, ale workflow'y produkcyjne muszą elegancko obsługiwać nieprawidłowe dane. Uczyńmy nasz workflow odpornym na brakujące lub null wartości.

---

## 6. Bezpieczna Nawigacja i Operatory Elvisa

Nasza funkcja `separateMetadata` obecnie zakłada, że wszystkie pola CSV są obecne i prawidłowe. Ale co się dzieje z niekompletnymi danymi? Dowiedzmy się.

### 6.1. Problem: Dostęp do Właściwości, Które Nie Istnieją

Powiedzmy, że chcemy dodać wsparcie dla opcjonalnych informacji o przebiegu sekwencjonowania. W niektórych laboratoriach próbki mogą mieć dodatkowe pole dla ID przebiegu sekwencjonowania lub numeru partii, ale nasz obecny CSV nie ma tej kolumny. Spróbujmy uzyskać do niej dostęp mimo to.

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

Problem polega na tym, że `row.run_id` zwraca `null`, ponieważ kolumna `run_id` nie istnieje w naszym CSV. Gdy próbujemy wywołać `.toUpperCase()` na `null`, następuje awaria. To tutaj operator bezpiecznej nawigacji ratuje sytuację.

### 6.2. Operator Bezpiecznej Nawigacji (`?.`)

Operator bezpiecznej nawigacji (`?.`) zwraca `null` zamiast rzucać wyjątek, gdy jest wywoływany na wartości `null`. Jeśli obiekt przed `?.` jest `null`, całe wyrażenie ewaluuje się do `null` bez wykonywania metody.

Zaktualizuj funkcję, aby używać bezpiecznej nawigacji:

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

Brak awarii! Workflow teraz elegancko obsługuje brakujące pole. Gdy `row.run_id` jest `null`, operator `?.` zapobiega wywołaniu `.toUpperCase()`, a `run_id` staje się `null` zamiast powodować wyjątek.

### 6.3. Operator Elvisa (`?:`) dla Wartości Domyślnych

Operator Elvisa (`?:`) dostarcza wartości domyślnych, gdy lewa strona jest "fałszywa" (jak wyjaśniono wcześniej). Jest nazwany na cześć Elvisa Presleya, ponieważ `?:` wygląda jak jego słynne włosy i oczy, gdy jest oglądany z boku!

Teraz, gdy używamy bezpiecznej nawigacji, `run_id` będzie `null` dla próbek bez tego pola. Użyjmy operatora Elvisa, aby dostarczyć wartość domyślną i dodać ją do naszej mapy `sample_meta`:

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

Idealnie! Teraz wszystkie próbki mają pole `run` z ich faktycznym ID przebiegu (wielkimi literami) lub wartością domyślną 'UNSPECIFIED'. Kombinacja `?.` i `?:` zapewnia zarówno bezpieczeństwo (brak awarii), jak i sensowne wartości domyślne.

Usuń teraz operator `.view()`, teraz gdy potwierdziliśmy, że działa.

!!! tip "Łączenie Bezpiecznej Nawigacji i Elvisa"

    Wzorzec `value?.method() ?: 'default'` jest powszechny w workflow'ach produkcyjnych:

    - `value?.method()` - Bezpiecznie wywołuje metodę, zwraca `null`, jeśli `value` jest `null`
    - `?: 'default'` - Dostarcza wartość zastępczą, jeśli wynik jest `null`

    Ten wzorzec elegancko obsługuje brakujące/niekompletne dane.

Używaj tych operatorów konsekwentnie w funkcjach, domknięciach operatorów (`.map{}`, `.filter{}`), skryptach procesów i plikach konfiguracyjnych. Zapobiegają awariom podczas obsługi rzeczywistych danych.

### Podsumowanie

- **Bezpieczna nawigacja (`?.`)**: Zapobiega awariom na wartościach null - zwraca null zamiast rzucać wyjątek
- **Operator Elvisa (`?:`)**: Dostarcza wartości domyślnych - `value ?: 'default'`
- **Łączenie**: `value?.method() ?: 'default'` to powszechny wzorzec

Te operatory sprawiają, że workflow'y są odporne na niekompletne dane — niezbędne do pracy w rzeczywistych warunkach.

---

## 7. Walidacja za pomocą `error()` i `log.warn`

Czasami musisz natychmiast zatrzymać workflow, jeśli parametry wejściowe są nieprawidłowe. W Nextflow możesz użyć wbudowanych funkcji takich jak `error()` i `log.warn`, a także standardowych konstrukcji programistycznych takich jak instrukcje `if` i logika boolowska, aby zaimplementować logikę walidacji. Dodajmy walidację do naszego workflow'a.

Utwórz funkcję walidacji przed blokiem workflow, wywołaj ją z workflow'a i zmień tworzenie kanału, aby używać parametru dla ścieżki pliku CSV. Jeśli parametr brakuje lub plik nie istnieje, wywołaj `error()`, aby zatrzymać wykonanie z jasnym komunikatem.

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Sprawdź, czy parametr wejściowy jest dostarczony
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

Workflow zatrzymuje się natychmiast z jasnym komunikatem błędu zamiast zawodzić tajemniczo później

Teraz uruchom go z nieistniejącym plikiem:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Na koniec uruchom go z poprawnym plikiem:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Wyjście polecenia"

    ```console
    <!-- TODO: output -->
    ```

Tym razem działa pomyślnie.

Możesz również dodać walidację wewnątrz funkcji `separateMetadata`. Użyjmy niefatalnego `log.warn` do wydawania ostrzeżeń dla próbek o niskiej głębokości sekwencjonowania, ale nadal pozwalając workflow'owi kontynuować:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Waliduj, czy dane mają sens
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
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Widzimy ostrzeżenie o niskiej głębokości sekwencjonowania dla jednej z próbek.

### Podsumowanie

- **`error()`**: Natychmiast zatrzymuje workflow z jasnym komunikatem
- **`log.warn`**: Wydaje ostrzeżenia bez zatrzymywania workflow'a
- **Wczesna walidacja**: Sprawdź wejścia przed przetwarzaniem, aby szybko zawieść z pomocnymi błędami
- **Funkcje walidacji**: Twórz logikę walidacji wielokrotnego użytku, którą można wywołać na początku workflow'a

Właściwa walidacja sprawia, że workflow'y są bardziej solidne i przyjazne dla użytkownika, wychwytując problemy wcześnie z jąsnymi komunikatami błędów.

---

## 8. Procedury Obsługi Zdarzeń Workflow'a

Do tej pory pisaliśmy kod w naszych skryptach workflow'a i definicjach procesów. Ale jest jeszcze jedna ważna funkcja, o której powinieneś wiedzieć: procedury obsługi zdarzeń workflow'a.

Procedury obsługi zdarzeń to domknięcia, które uruchamiają się w określonych punktach cyklu życia Twojego workflow'a. Są idealne do dodawania logowania, powiadomień lub operacji czyszczenia. Te procedury powinny być zdefiniowane w Twoim skrypcie workflow'a obok definicji workflow'a.

### 8.1. Procedura Obsługi `onComplete`

Najczęściej używaną procedurą obsługi zdarzeń jest `onComplete`, która uruchamia się, gdy Twój workflow się kończy (niezależnie od tego, czy zakończył się sukcesem, czy niepowodzeniem). Dodajmy jedną, aby podsumować wyniki naszego pipeline'u.

Dodaj procedurę obsługi zdarzeń do pliku `main.nf`, wewnątrz definicji workflow'a:

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

To domknięcie uruchamia się, gdy workflow się kończy. Wewnątrz masz dostęp do obiektu `workflow`, który dostarcza użytecznych właściwości o wykonaniu.

Uruchom swój workflow, a zobaczysz to podsumowanie pojawiające się na końcu!

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

Teraz otrzymujemy jeszcze bardziej informacyjne podsumowanie, w tym komunikat o sukcesie/niepowodzeniu i katalog wyjściowy, jeśli został określony:

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

Możesz również zapisać podsumowanie do pliku używając operacji na plikach:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... Twój kod workflow'a ...

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

### 8.2. Procedura Obsługi `onError`

Oprócz `onComplete`, jest jeszcze jedna procedura obsługi zdarzeń, której możesz użyć: `onError`, która uruchamia się tylko jeśli workflow zawiedzie:

```groovy title="main.nf - onError handler"
workflow {
    // ... Twój kod workflow'a ...

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

Możesz używać wielu procedur obsługi razem w swoim skrypcie workflow'a:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... Twój kod workflow'a ...

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

- **Domknięcia procedur obsługi zdarzeń**: Domknięcia w Twoim skrypcie workflow'a, które uruchamiają się w różnych punktach cyklu życia
- **Procedura obsługi `onComplete`**: Do podsumowań wykonania i raportowania wyników
- **Procedura obsługi `onError`**: Do obsługi błędów i logowania niepowodzeń
- **Właściwości obiektu workflow**: Dostęp do `workflow.success`, `workflow.duration`, `workflow.errorMessage`, itp.

Procedury obsługi zdarzeń pokazują, jak możesz użyć pełnej mocy języka Nextflow w swoich skryptach workflow'a, aby dodać zaawansowane możliwości logowania i powiadamiania.

---

## Podsumowanie

Gratulacje, udało Ci się!

W trakcie tego side questa zbudowałeś kompleksowy pipeline przetwarzania próbek, który ewoluował od podstawowej obsługi metadanych do zaawansowanego, gotowego do produkcji workflow'a.
Każda sekcja budowała na poprzedniej, demonstrując, jak konstrukcje programistyczne przekształcają proste workflow'y w potężne systemy przetwarzania danych, z następującymi korzyściami:

- **Czystszy kod**: Rozumienie przepływu danych vs skryptowania pomaga pisać bardziej zorganizowane workflow'y
- **Solidna obsługa**: Bezpieczna nawigacja i operatory Elvisa sprawiają, że workflow'y są odporne na brakujące dane
- **Elastyczne przetwarzanie**: Logika warunkowa pozwala Twoim workflow'om odpowiednio przetwarzać różne typy próbek
- **Adaptacyjne zasoby**: Dynamiczne dyrektywy optymalizują użycie zasobów na podstawie charakterystyk wejściowych

Ta progresja odzwierciedla rzeczywistą ewolucję pipeline'ów bioinformatycznych, od prototypów badawczych obsługujących kilka próbek do systemów produkcyjnych przetwarzających tysiące próbek w laboratoriach i instytucjach.
Każde wyzwanie, które rozwiązałeś, i wzorzec, którego się nauczyłeś, odzwierciedla rzeczywiste problemy, z którymi borykają się deweloperzy podczas skalowania workflow'ów Nextflow.

Stosowanie tych wzorców w Twojej własnej pracy umożliwi Ci budowanie solidnych, gotowych do produkcji workflow'ów.

### Kluczowe wzorce

1.  **Przepływ danych vs Skryptowanie:** Nauczyłeś się rozróżniać między operacjami przepływu danych (orkiestracja kanałów) a skryptowaniem (kod, który manipuluje danymi), w tym kluczowe różnice między operacjami na różnych typach, takich jak `collect` na Channel vs List.

    - Przepływ danych: orkiestracja kanałów

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Skryptowanie: przetwarzanie danych na kolekcjach

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Zaawansowane Przetwarzanie Napisów**: Opanowałeś wyrażenia regularne do parsowania nazw plików, dynamiczne generowanie skryptów w procesach oraz interpolację zmiennych (Nextflow vs Bash vs Powłoka).

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

3.  **Tworzenie Funkcji Wielokrotnego Użytku**: Nauczyłeś się wyodrębniać złożoną logikę do nazwanych funkcji, które mogą być wywoływane z operatorów kanałów, czyniąc workflow'y bardziej czytelnymi i łatwiejszymi w utrzymaniu.

    - Zdefiniuj nazwaną funkcję

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

    - Wywołaj nazwaną funkcję w workflow'ie

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Dynamiczne Dyrektywy Zasobów z Domknięciami**: Zbadałeś używanie domknięć w dyrektywach procesów do adaptacyjnej alokacji zasobów na podstawie charakterystyk wejściowych.

    - Nazwane domknięcia i kompozycja

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Domknięcia z dostępem do zakresu

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logika Warunkowa i Kontrola Procesów**: Dodałeś inteligentne kierowanie używając operatorów `.branch()` i `.filter()`, wykorzystując prawdziwość do zwięzłych wyrażeń warunkowych.

    - Użyj `.branch()` do kierowania danych przez różne gałęzie workflow'a

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
    if (sample.files) println "Has files"
    ```

    - Użyj `filter()` do wybierania podzbiorów danych z 'prawdziwością'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Bezpieczna Nawigacja i Operatory Elvisa**: Uczyniłeś pipeline odpornym na brakujące dane używając `?.` do bezpiecznego dostępu do właściwości null i `?:` do dostarczania wartości domyślnych.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Walidacja za pomocą error() i log.warn**: Nauczyłeś się walidować wejścia wcześnie i szybko zawodzić z jąsnymi komunikatami błędów.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Procedury Obsługi Zdarzeń Konfiguracji**: Nauczyłeś się używać procedur obsługi zdarzeń workflow'a (`onComplete` i `onError`) do logowania, powiadomień i zarządzania cyklem życia.

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

    - Używanie `onError` do podejmowania działań specyficznie w przypadku niepowodzenia

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

- [Dokumentacja Języka Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operatory Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Składnia Skryptów Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteka Standardowa Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Upewnij się, że sprawdzisz te zasoby, gdy musisz zbadać bardziej zaawansowane funkcje.

Skorzystasz z praktykowania i rozwijania swoich umiejętności, aby:

- Pisać czystsze workflow'y z właściwym rozdzieleniem między przepływem danych a skryptowaniem
- Opanować interpolację zmiennych, aby uniknąć powszechnych pułapek ze zmiennymi Nextflow, Bash i powłoki
- Używać dynamicznych dyrektyw zasobów dla wydajnych, adaptacyjnych workflow'ów
- Przekształcać kolekcje plików w poprawnie sformatowane argumenty linii poleceń
- Elegancko obsługiwać różne konwencje nazewnictwa plików i formaty wejściowe używając regex i przetwarzania napisów
- Budować kod wielokrotnego użytku i łatwy w utrzymaniu używając zaawansowanych wzorców domknięć i programowania funkcyjnego
- Przetwarzać i organizować złożone zbiory danych używając operacji na kolekcjach
- Dodawać walidację, obsługę błędów i logowanie, aby uczynić Twoje workflow'y gotowymi do produkcji
- Implementować zarządzanie cyklem życia workflow'a za pomocą procedur obsługi zdarzeń

---

## Co dalej?

Wróć do [menu Side Questów](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
